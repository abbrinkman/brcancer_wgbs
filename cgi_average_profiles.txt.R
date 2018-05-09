library(GenomicRanges)
library(ggplot2)
library(reshape2)
library(parallel)
library(data.table)
library(rtracklayer)

meth <- get(load("~/BiSeq_BASIS/methcounts/meth.gr.RData"))
pmds <- get(load("~/BiSeq_BASIS/PMDs_all_samples/pmdmeth.basis.RData"))


pat <- intersect(names(pmds), gsub("\\.[MT]","", colnames(mcols(meth))))
pat <- pat[!pat %in% c("PD9590a")]


# get CGIs, promoters, genes
session <- browserSession()
genome(session) <- "hg19"
query <- ucscTableQuery(session, "CpG Islands")
cgi <- GRanges(track(query))
cgi$name <- NULL
cgi$elem <- "CGI"
cgi <- cgi[,"elem"]

# make CGI bins
makeBinsUndirectional <- function(gr, up.down, inside, win.size) {
  # up.down:  number of bp upstream and downstream from the used elements
  # inside:   number of bp into (inside) the used elements (overshooting is trimmed)
  # win.size: window size (bin size) 
  d1 <- as.data.frame(gr[,0])[,1:3]
  
  # upstream bins
  l1 <- list()
  for (i in seq(-up.down, inside, by=win.size)) {
     l1[[paste0("up.", i)]] <- d1
     l1[[paste0("up.", i)]]$start <- d1$start + i - win.size
     l1[[paste0("up.", i)]]$end <- d1$start + i
  }
  # downstream bins
  for (i in -seq(-inside, up.down, by=win.size)) {
     l1[[paste0("down.", i)]] <- d1
     l1[[paste0("down.", i)]]$start <- d1$end - i 
     l1[[paste0("down.", i)]]$end <- d1$end - i + win.size
  }
  for (i in names(l1)) {
     l1[[i]]$pos <- sapply(strsplit(i, "\\."), function(x) {x[1]})
     l1[[i]]$dist <- as.numeric(sapply(strsplit(i, "\\."), function(x) {x[2]}))
  }
  gr1 <- unlist(GRangesList(lapply(l1, makeGRangesFromDataFrame, keep.extra.columns=T)))
  
  # take out all 'inside' bins (positive distance) that actually overshoot the CGI
  gr1 <- gr1[!c(gr1$dist > 0 & !overlapsAny(gr1, gr))]
}
cgi.bins <- makeBinsUndirectional(cgi, up.down=10000, inside=5000, win.size=100)


# function for calculation of weightedmean in regions
weightedMethWholeSet <- function(REGIONS, WGBS, NCPG, NCORES) {
  # REGIONS: GRanges object
  # WGBS: GRanges object with columns (PDxxxx.T , PDxxxx.M) or (xxxx_cov, xxxx_meth) as patients
  # get unique sample names
  samples <- as.list(unique(gsub("_cov","",gsub("_meth","",gsub("\\.[TM]$","", 
      colnames(mcols(WGBS)))))))
  names(samples) <- unlist(samples)
  wgbs.tmp <- WGBS
  ov <- findOverlaps(REGIONS, wgbs.tmp)
    getMeans <- function(x) {
    message(x)
    cols.t <- which(grepl(x, colnames(mcols(wgbs.tmp))) & c(grepl("\\.T$", 
        colnames(mcols(wgbs.tmp))) | grepl("_cov", colnames(mcols(wgbs.tmp)))))
    cols.m <- which(grepl(x, colnames(mcols(wgbs.tmp))) & c(grepl("\\.M$", 
        colnames(mcols(wgbs.tmp))) | grepl("_meth", colnames(mcols(wgbs.tmp)))))

    dt1 <- data.table("regions"=queryHits(ov),
        "T"=mcols(wgbs.tmp[subjectHits(ov)])[,cols.t],
        "M"=mcols(wgbs.tmp[subjectHits(ov)])[,cols.m])
    if (identical(range(dt1$M, na.rm=T), c(0,1))) {
      dt1$M  <- round(dt1$M * dt1$T, 0)
    }
    setkey(dt1, regions)

    # calculate weighted mean
    dtT <- dt1[,sum(T, na.rm=T),by=regions]
    dtM <- dt1[,sum(M, na.rm=T),by=regions]
    dtC <- dt1[,length(M),by=regions]

    df.mean <- data.frame(
        "regions"=as.numeric(1:length(REGIONS)),
        "meth"=as.numeric(rep(NaN, length(REGIONS))),
        "ncpg"=as.numeric(rep(0, length(REGIONS))))
    df.mean$meth[df.mean$regions %in% dtM$regions] <- round(dtM$V1/dtT$V1, 2)
    df.mean$ncpg[df.mean$regions %in% dtM$regions] <- dtC$V1
    df.mean$meth[df.mean$ncpg < NCPG] <- NaN
    colnames(df.mean) <- gsub("meth", x, colnames(df.mean))

    # return data
    df.mean[,x]
  }
  means.list <- mclapply(samples, getMeans, mc.cores=NCORES)
  gr4 <- REGIONS[,0]
  gr4@elementMetadata@listData <- means.list
  gr4
}

# calculate mean methylation per bin
cgi.bins.meth <- weightedMethWholeSet(REGIONS=cgi.bins, WGBS=meth, NCPG=2, NCORES=2)
save(cgi.bins.meth, file="/tmp/cgi.bins.meth.RData")
cgi.bins.meth <- cgi.bins.meth[,pat]

# determine for each bin whether it's within a PMD (per patient)
cgi.bins.pmd <- cgi.bins[,0]
mcols(cgi.bins.pmd) <- as.data.frame(sapply(pat, function(x) {
    ifelse(overlapsAny(cgi.bins, pmds[[x]], type="within"), "in", "out")}))


meanMethBinPMD <- function(pd) {
  message(pd)
  meth.in <- tapply(
      mcols(cgi.bins.meth)[[pd]][mcols(cgi.bins.pmd)[[pd]]=="in"], 
      names(cgi.bins.meth)[mcols(cgi.bins.pmd)[[pd]]=="in"], function(x) {mean(x, na.rm=T)})
  meth.out <- tapply(
      mcols(cgi.bins.meth)[[pd]][mcols(cgi.bins.pmd)[[pd]]=="out"],
      names(cgi.bins.meth)[mcols(cgi.bins.pmd)[[pd]]=="out"], function(x) {mean(x, na.rm=T)})
  list("in"=meth.in, "out"=meth.out)
}
cgi.bins.meth.mean <- sapply(pat, meanMethBinPMD, simplify=F, USE.NAMES=T)

# shape the data
d1 <- melt(cgi.bins.meth.mean)
colnames(d1) <- c("pos.dist","meth","PMD","patient")
d1$pos <- sapply(strsplit(as.character(d1$pos.dist), "\\."), function(x) {x[1]})
d1$dist <- as.numeric(sapply(strsplit(as.character(d1$pos.dist), "\\."), function(x) {x[2]}))
d1$pos.PMD <- factor(paste(d1$pos, d1$PMD, sep="_"), levels=c("up_out","down_out","up_in","down_in"))

# prepare data for ribbonplot
d2 <- d1[!d1$patient %in% c("HMEC","MCF7"),]
d2 <- do.call(rbind,lapply(split(d2, d2$pos.dist), function(x) {melt(tapply(x$meth, x$PMD, summary))}))
d2 <- data.frame(d2, t(sapply(strsplit(as.character(d2$value), ","), function(x) {as.numeric(gsub("[ )(c]","",x))})))
d2$value <- NULL
colnames(d2) <- c("PMD", "min", "q1", "median", "mean", "q3", "max")
d2$pos <- gsub("\\..+$","", rownames(d2))
d2$dist <- as.numeric(gsub("\\..+$","",gsub("down\\.","", gsub("up\\.","", rownames(d2)))))
d2$pos.PMD <- factor(paste(d2$pos, d2$PMD, sep="_"), levels=c("up_out","down_out","up_in","down_in"))

p1 <- ggplot(d2[d2$dist >= -4000 & d2$dist <= 800,], aes(dist, median)) + 
      geom_ribbon(aes(ymax=max, ymin=min, fill=PMD), alpha=0.3) + 
      geom_line(aes(color=PMD)) + facet_wrap(~pos) + 
      scale_fill_manual(values=c("in"="red", "out"="black")) + 
      scale_color_manual(values=c("in"="red","out"="black")) +
      theme_classic() +
      theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1))
ggsave(p1, file="ribbonplot_CGImeth_InOutPMD.pdf", height=3.3, width=7, scale=0.55) # not used for publication

# repeat the ribbonplot over all patients, but without first creating means per patient
l4 <- list()
l4[["in"]] <- data.frame()
l4[["out"]] <- data.frame()
for (pd in pat) {
  message(pd)
  l4[["in"]] <- rbind(l4[["in"]],
           data.frame("pos"=names(cgi.bins.meth)[mcols(cgi.bins.pmd)[[pd]]=="in"],
                     "meth"=mcols(cgi.bins.meth)[[pd]][mcols(cgi.bins.pmd)[[pd]]=="in"]))
  l4[["out"]] <- rbind(l4[["out"]],
          data.frame("pos"=names(cgi.bins.meth)[mcols(cgi.bins.pmd)[[pd]]=="out"],
                     "meth"=mcols(cgi.bins.meth)[[pd]][mcols(cgi.bins.pmd)[[pd]]=="out"]))
}
d4.in <- as.data.frame(do.call(rbind, tapply(l4$`in`$meth, l4$`in`$pos, summary, na.rm=T)))
d4.out <- as.data.frame(do.call(rbind, tapply(l4$out$meth, l4$out$pos, summary, na.rm=T)))
d4.in$PMD <- "in"
d4.out$PMD <- "out"

d4.in$`NA's` <- NULL
colnames(d4.in) <- c("min","q25","median","mean","q75", "max", "PMD")
d4.in$dist <- as.numeric(gsub("\\..+$","",gsub("down\\.","", gsub("up\\.","", rownames(d4.in)))))
d4.in$pos <- gsub("\\..+$","", rownames(d4.in))

d4.out$`NA's` <- NULL
colnames(d4.out) <- c("min","q25","median","mean","q75", "max", "PMD")
d4.out$dist <- as.numeric(gsub("\\..+$","",gsub("down\\.","", gsub("up\\.","", rownames(d4.out)))))
d4.out$pos <- gsub("\\..+$","", rownames(d4.out))

d4 <- rbind(d4.in, d4.out)

p3 <- ggplot(d4[d4$dist >= -4000 & d4$dist <= 800,], aes(dist, median)) +
      geom_ribbon(aes(ymax=q75, ymin=q25, fill=PMD), alpha=0.3) +
      geom_line(aes(color=PMD)) + facet_wrap(~pos) +
      scale_fill_manual(values=c("in"="red", "out"="black")) +
      scale_color_manual(values=c("in"="red","out"="black")) +
      theme_classic() +
      theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1))
ggsave(p3, file="ribbonplot_CGImeth_InOutPMD_1.pdf", height=3.3, width=7, scale=0.55) #### Figure 3C ####


