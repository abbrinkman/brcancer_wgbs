library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(reshape2)
library(ggplot2)

source("~/tools/BASIS_common_functions.txt.R")

pmdmeth.basis <- get(load("~/BiSeq_BASIS/PMDs_all_samples/pmdmeth.basis.RData"))
pmdmeth.basis <- pmdmeth.basis[!names(pmdmeth.basis) %in% c("HMEC","MCF7")]

# get the CNVs from BASIS
cnv <- read.table("~/Infinium/BASIS_data_all_TN_QC_analysis/suppl_nik-zainal_landscape_nature17676-s3/Supplementary Table 4.Copy.Number.Segments.txt", header=T)
cnv$chr <- paste0("chr", cnv$chr)
cnv$chr <- gsub("chr23","chrX", cnv$chr)
cnv <- makeGRangesFromDataFrame(cnv, keep.extra.columns=T)
mcols(cnv)$pd <- gsub("a2","a",gsub("a_2","a", mcols(cnv)$SampleID)) # fix names only for the "a" samples (=tumors)
sum(names(pmdmeth.basis) %in% mcols(cnv)$pd) # should be 25

# using genomic tiles, calculate the mean CN per tile, over 560 patients
tiles.cn <- unlist(tileGenome(seqlengths=seqlengths(Hsapiens)[1:23], ntile=10000))
cnv.l <- split(cnv, mcols(cnv)$pd)
countCN <- function(PD) {
  gr1 <- cnv.l[[PD]]
  ov <- findOverlaps(tiles.cn, gr1, type="within")
  gr2 <- tiles.cn
  mcols(gr2)[[PD]] <- rep(NA, length(gr2))
  mcols(gr2)[[PD]][queryHits(ov)] <- mcols(gr1)$seg.mean[subjectHits(ov)]
  mcols(gr2)[[PD]]
}

# count, for only the 25 overlapping patients with WGBS
l3 <- sapply(names(pmdmeth.basis)[names(pmdmeth.basis) %in% mcols(cnv)$pd], function(PD) { message(PD)
    countCN(PD)}, USE.NAMES=T, simplify=F) # use only the 'PDxxxxa' patients ('PDxxxxb' = normal)

cn.count <- tiles.cn
mcols(cn.count)$cn.mean <- apply((do.call(cbind, l3)), 1, mean, na.rm=T)
mcols(cn.count)$pos <- round(start(cn.count)+(width(cn.count)/2),0)
mcols(cn.count)$chr <- seqnames(cn.count)
d3 <- melt(as.data.frame(mcols(cn.count)), id.vars=c("chr", "pos"),value.name="CN" )

chr.cols <- rep(c("#d95f02", "#7570b3"), length(levels(d3$chr)))[1:length(levels(d3$chr))]
names(chr.cols) <- levels(d3$chr)

p <- ggplot(d3, aes(pos, CN)) + 
  geom_line(aes(colour=chr)) + 
  facet_grid(~chr, scales="free_x", space="free_x") + 
  scale_colour_manual(values=chr.cols) + 
  theme(strip.text = element_text(angle=90), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
      panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  xlab("") +
  ylab("mean CN (25 patients)")
ggsave(p, file="scatterplots_CNcounts.pdf", height=3, width=10) #### Supplemental Figure 1C ####

# same for the PMDs, calculate PMD frequency per tile, all 30 patients
tiles.pmd <- unlist(tileGenome(seqlengths=seqlengths(Hsapiens)[1:23], ntile=10000))
countPMD <- function(PD) {
  countOverlaps(tiles.pmd, pmdmeth.basis[[PD]], type="within")
}
l4 <- sapply(names(pmdmeth.basis), function(PD) { message(PD)
    countPMD(PD)}, USE.NAMES=T, simplify=F)

# plot PMD frequency
pmd.count <- tiles.pmd
mcols(pmd.count)$pmd.freq <- apply((do.call(cbind, l4)), 1, sum, na.rm=T)
mcols(pmd.count)$pos <- round(start(pmd.count)+(width(pmd.count)/2),0)
mcols(pmd.count)$chr <- seqnames(pmd.count)
runmean <- function(x, N) {
  l <- split(x, unlist(lapply(seq(1,length(x), N), function(y) {rep(y, N)}))[1:length(x)])
  unlist(lapply(l, function(z) {rep(mean(z, na.rm=T), length(z))}))
}
pmd.count.1 <- pmd.count
mcols(pmd.count.1)$pmd.freq <- runmean(mcols(pmd.count)$pmd.freq, 15)
d4 <- melt(as.data.frame(mcols(pmd.count.1)), id.vars=c("chr", "pos"),value.name="PMD.freq" )



# is PMDfreq dependent on CNVs?
d5 <- data.frame(as.data.frame(mcols(cn.count)), as.data.frame(mcols(pmd.count)) )
d5 <- d5[complete.cases(d5),]

# Pearson correlation
cor(d5$pmd.freq, d5$cn.mean)
  # [1] 0.1745826

