library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(data.table)
library(matrixStats)
library(mixtools)
library(ggplot2)
library(reshape2)

meth <- get(load("~/BiSeq_BASIS/methcounts/meth.gr.RData"))
meth <- meth[, !colnames(mcols(meth)) %in% c("PD9590a")]

# get the sequence context (tetranucleotides) for each CpG
meth$context <- as.data.frame(getSeq(Hsapiens, 
    promoters(GenomicRanges::shift(meth, shift=1),upstream=1, downstream=3)))$x

# make a restricted dataset with only the solo-WCGW CpGs
meth.solo <- meth[grep("[AT]CG[AT]", meth$context)]
meth.solo$context <- NULL

#   # are the solo-CpGs the same CpGs as those within CGIs?
#   session <- browserSession()
#   genome(session) <- "hg19"
#   query <- ucscTableQuery(session, "CpG Islands")
#   cgi <- GRanges(track(query))
#   meth$cgi <- ifelse(overlapsAny(meth, cgi), "CGI", "noCGI")
#   meth$solo <- ifelse(grepl("[AT]CG[AT]", meth$context), "Solo", "noSolo")
#   cgi.solo <- table(meth$cgi, meth$solo)
#   
#   #  cgi.solo
#   #
#   #          noSolo     Solo
#   #  CGI    1946190   143398
#   #  noCGI 19413601  6714259


# make bins of 100 kb
# get the hg19 genome for constructing the tiles
hg19 <- BSgenome.Hsapiens.UCSC.hg19
hg19 <- as.data.frame(seqlengths(hg19))
colnames(hg19) <- "end"
hg19$start <- 1
hg19$chr <- rownames(hg19)
hg19 <- makeGRangesFromDataFrame(hg19)
hg19 <- hg19[!grepl("chrUn", seqnames(hg19))]
hg19 <- hg19[!grepl("random", seqnames(hg19))]
hg19 <- hg19[!grepl("hap", seqnames(hg19))]
hg19 <- hg19[!grepl("chrM", seqnames(hg19))]
#hg19 <- hg19[!grepl("chrX", seqnames(hg19))] # leave in chrX here
hg19 <- hg19[!grepl("chrY", seqnames(hg19))]

# remove centromers from genome
centro <- read.delim(pipe("cat ~/Datarepository/hg19_cytoBand.bed |grep -e acen -e gvar -e stalk"),
    header=F, sep="\t")
colnames(centro) <- c("chr","start","end","arm","type")
centro <- makeGRangesFromDataFrame(centro, keep.extra.columns=T)
hg19 <- setdiff(hg19, centro)

# make the tiles
tile.width <- 100000
tiles <- unlist(tile(hg19, width=tile.width))

# determine mean meth (solo-CpGs) per tile 
weightedMethWholeSet <- function(REGIONS, WGBS, NCPG, NCORES=1) {
  # REGIONS: GRanges object
  # WGBS: GRanges object with columns (xx.T , xx.M) as patients
  # get unique sample names
  samples <- as.list(unique(gsub("\\.[TM]$","", colnames(mcols(WGBS)))))
  names(samples) <- unlist(samples)
  ov <- findOverlaps(REGIONS, WGBS)
    getMeans <- function(x) {
    message(x)
    cols.t <- which(grepl(x, colnames(mcols(WGBS))) & grepl("\\.T$", colnames(mcols(WGBS))))
    cols.m <- which(grepl(x, colnames(mcols(WGBS))) & grepl("\\.M$", colnames(mcols(WGBS))))
    dt1 <- data.table("regions"=queryHits(ov),
        "T"=mcols(WGBS[subjectHits(ov)])[,cols.t],
        "M"=mcols(WGBS[subjectHits(ov)])[,cols.m])
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
tiles.meth.solo <- weightedMethWholeSet(REGIONS=tiles, WGBS=meth.solo[, grep("^PD", colnames(mcols(meth.solo)))], 
    NCPG=5, NCORES=2)

tiles.meth.all <- weightedMethWholeSet(REGIONS=tiles, WGBS=meth[, grep("^PD", colnames(mcols(meth)))],
    NCPG=5, NCORES=1)

save(tiles.meth.solo, file="tiles.meth.solo.RData")
save(tiles.meth.all, file="tiles.meth.all.RData")

# calculate the StDev for each tile
tiles$sd.solo <- rowSds(as.matrix(mcols(tiles.meth.solo)))
tiles$mean.solo <- rowMeans(as.matrix(mcols(tiles.meth.solo)))

tiles$sd.all <- rowSds(as.matrix(mcols(tiles.meth.all)))
tiles$mean.all <- rowMeans(as.matrix(mcols(tiles.meth.all)))


# mixed Gaussion model with default parameters
mdl <- normalmixEM(tiles$sd.solo[!is.na(tiles$sd.solo)], )

pdf("densityplots_mixedGaussian.pdf") ##### Supplemental Figure 4E #####
par(mfrow=c(2,2))
plot(mdl, density=T, xlab2="StDev (solo-CpGs in 100 kb tiles)")

# overplot with original density line
lines(density(tiles$sd.solo, na.rm=T), col="blue", lty=2, lwd=2)

# our cutoff
pmd.cutoff <- 0.062

# display our dutoff
lines(c(pmd.cutoff, pmd.cutoff), c(par("usr")[3], par("usr")[4]), col="brown", lwd=2)

dev.off()

# determine PMDs according Berman's method
tiles$pmd <- ifelse(tiles$sd.solo > pmd.cutoff, "PMD", "noPMD")

# total PMD bases for this method
 sum(width(tiles[which(tiles$pmd=="PMD")]))
# [1] 1875682239

# fraction of the genome PMDs for this method
sum(width(tiles[which(tiles$pmd=="PMD")]))/3.2e09
# [1] 0.5861507

pmds.solo <- reduce(tiles[which(tiles$pmd=="PMD")])

# overlap with the MethylSeekR PMDs called on all CpGs
pmds <- get(load("~/BiSeq_BASIS/PMDs_all_samples/pmdmeth.basis.RData"))

# remove chrY (makes no sense in female tumors)
pmds <- lapply(pmds, function(x) {x[seqnames(x) != "chrY"]})

# fraction overlap that the solo-PMDs have with each individual PMD-track
fractionOverlap <- function(x) { # x = GRanges with PMDs (PDxxxx)
  sum(width(intersect(x, pmds.solo)))/sum(width(x))
}
pmd.solo.overlap <- sapply(pmds[!names(pmds) %in% c("HMEC","MCF7")], fractionOverlap)

d1 <- data.frame("fraction"=pmd.solo.overlap, "sample"=names(pmd.solo.overlap))
p1 <- ggplot(d1, aes(reorder(sample,fraction), fraction)) + 
  geom_bar(stat="identity", fill="black") + 
  theme_classic() + 
  theme(axis.text.x=element_text(color="black", angle=45, hjust=1, vjust=1), 
      axis.text.y=element_text(color="black")) + 
  xlab("sample") + 
  ylab("fraction overlap with solo-WCGW common PMDs")

ggsave(p1, file="barplot_overlap_PMDs_soloPMDs.pdf", height=4, width=7) ##### Supplemental Figure 4F #####

# fraction overlap that the solo-PMDs have with the merged (union) of all PMD-tracks
pmds.union <- reduce(unlist(GRangesList(lapply(pmds[!names(pmds) %in% c("HMEC","MCF7")], 
    function(x) {colnames(mcols(x)) <- "meth" ; x}))))

sum(width(intersect(pmds.union, pmds.solo)))/sum(width(pmds.union))
# [1] 0.8232483

