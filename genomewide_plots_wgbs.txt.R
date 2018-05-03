library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(rtracklayer)
library(reshape2)
library(ggplot2)
library(data.table)
library(RColorBrewer)

# create genomic tiles of 10 kb
tiles.pmd <- unlist(tileGenome(seqlengths=seqlengths(Hsapiens)[1:23], ntile=10000))

# get the WGBS data (GRanges object, with CpG positions as rows, and
# colnames PDxxxx.T = total reads, PDxxxx.M = methylated reads)
wgbs <- get(load("~/BiSeq_BASIS/methcounts/meth.gr.RData"))
wgbs <- wgbs[, -grep("PD9590a", colnames(mcols(wgbs)))]

# for calculation of mean meth per tile, remove CGIs, CGI shores, and promoters
session <- browserSession()
genome(session) <- "hg19"
query <- ucscTableQuery(session, "CpG Islands")
cgi <- GRanges(track(query))
cgi <- cgi[seqnames(cgi) %in% seqnames(Hsapiens)[1:25]]
seqlevels(cgi) <- seqnames(Hsapiens)[1:25]
shore.size <- 2000
cgi.shores <- cgi
start(cgi.shores) <- start(cgi.shores)-shore.size
end(cgi.shores) <- end(cgi.shores)+shore.size

query <- ucscTableQuery(session, "GENCODE Genes V19")
genes <- GRanges(track(query))
genes <- genes[seqnames(genes) %in% seqnames(Hsapiens)[1:25]]
seqlevels(genes) <- seqnames(Hsapiens)[1:25]
prom.size <- 1000
prom <- promoters(genes, upstream=prom.size/2, downstream=prom.size/2)

cgi.shores.prom <- c(cgi.shores[,0], prom[,0])

# take out the promoter and CpG island CpGs
wgbs1 <- wgbs[-queryHits(findOverlaps(wgbs, cgi.shores.prom))]


## calculate mean methylation in each tile
weightedMethWholeSet <- function(REGIONS, WGBS, NCPG, NCORES=1) {
  # REGIONS: GRanges object
  # WGBS: GRanges object with columns (PDxxxx.T , PDxxxx.M) as patients
  # NCPG: minimal number of CpGs per tile
  # NCORES: number of parallel processes

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
tiles.mean <- weightedMethWholeSet(REGIONS=tiles.pmd, WGBS=wgbs1, NCPG=5)


mcols(tiles.mean)$pos <- round(start(tiles.mean)+(width(tiles.mean)/2),0)
mcols(tiles.mean)$chr <- seqnames(tiles.mean)

d6 <- melt(as.data.frame(mcols(tiles.mean[,!colnames(mcols(tiles.mean)) %in% c("HMEC","MCF7")])), 
    id.vars=c("chr", "pos"),variable.name="patient", value.name="meth")
d6$pos <- as.factor(d6$pos)

# order the samples by clustering on the tiled means, use only the autosomes
m1 <- as.matrix(mcols(tiles.mean[!seqnames(tiles.mean) %in% c("chrX","chrY"),
    -grep("[pc][oh][sr]", colnames(mcols(tiles.mean)))]))
m1 <- m1[complete.cases(m1),]
hc <- hclust(as.dist(1-cor(m1)), method="ward.D")
d6$patient.hc <- factor(d6$patient, levels=hc$labels[hc$order])

# make the genome-wide tiled plot, clustered based on the clustering above
p <- ggplot(d6, aes(pos, patient.hc)) +
  geom_tile(aes(colour=meth, fill=meth)) +
  scale_fill_gradientn(colours=c("white", "white", "yellow","blue"), values=c(0, 0.25, 0.5, 1), na.value="white") +
  scale_colour_gradientn(colours=c("white", "white", "yellow","blue"), values=c(0, 0.25, 0.5, 1), na.value="white") +
  facet_grid(~chr, scales="free_x", space="free_x") +
  theme(strip.text = element_text(angle=90), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
      axis.text.y=element_text(colour="black"), axis.ticks.y=element_line(colour="black"),
      panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
      panel.border=element_blank()) +
  xlab("") +
  ylab("patient")
ggsave(p, file="tilemaps_PMDmeth_clustered.png", height=4, width=10) #### Figure 1A / Supplemental Figure 2A ####
ggsave(p, file="tilemaps_PMDmeth_clustered.pdf", height=4, width=10)

# function for plotting only one chromosome (clustered, as above)
plotPMDmeth_chromosome <- function(x) {
  p <- ggplot(d6[d6$chr==x,], aes(pos, patient.hc)) +
    geom_tile(aes(colour=meth, fill=meth)) +
    scale_fill_gradientn(colours=c("white","white","yellow","blue"), values=c(0, 0.25, 0.5, 1), na.value="white") +
    scale_colour_gradientn(colours=c("white","white","yellow","blue"), values=c(0, 0.25, 0.5, 1), na.value="white") +
    facet_grid(~chr, scales="free_x", space="free_x") +
    theme(strip.text = element_text(angle=90), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_text(colour="black"), axis.ticks.y=element_line(colour="black"),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.border=element_blank()) +
    xlab("") +
    ylab("patient")
  ggsave(p, file=paste0("tilemaps_PMDmeth_per_chromosome/tilemaps_PMDmeth_clustered_",x,".png"), height=4, width=10)
  ggsave(p, file=paste0("tilemaps_PMDmeth_per_chromosome/tilemaps_PMDmeth_clustered_", x, ".pdf"), height=4, width=10)
}
system("mkdir -p tilemaps_PMDmeth_per_chromosome")
sapply(unique(d6$chr), plotPMDmeth_chromosome) #### Figure 1B / Supplemental Figure 2B ####


# check significance of association between ER and overall PMDmeth within the tiles
source("~/tools/BASIS_common_functions.txt.R")
d6.er <- as.data.frame(tapply(d6$meth, d6$patient, mean, na.rm=T))
names(d6.er) <- "meth.mean"
d6.er$patient <- rownames(d6.er)
d6.er$ER <- clin$ER.FPKM[match(d6.er$patient, paste0(clin$sample_name, "a"))]

# check if data is normally distributed
shapiro.test(d6.er$meth.mean)
#  
#          Shapiro-Wilk normality test
#  
#  data:  d6.er$meth.mean
#  W = 0.97017, p-value = 0.5438
qqnorm(d6.er$meth.mean, col=d6.er$ER, pch=20, cex=3)
qqline(d6.er$meth.mean)

# OK to do a t-test
t.test(meth.mean ~ ER, data=d6.er)
#  
#          Welch Two Sample t-test
#  
#  data:  meth.mean by ER
#  t = -1.6971, df = 4.9051, p-value = 0.1516
#  alternative hypothesis: true difference in means is not equal to 0
#  95 percent confidence interval:
#   -0.13384471  0.02776713
#  sample estimates:
#  mean in group negative mean in group positive
#               0.6825972              0.7356360
#  


