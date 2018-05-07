library(reshape2)
library(GenomicRanges)
library(rtracklayer)

pmds <- get(load("~/BiSeq_BASIS/PMDs_all_samples/pmdmeth.basis.RData"))
pmds <- pmds[!names(pmds) %in% c("HMEC","MCF7")]


# use multiIntersecBed
tmp <- tempdir()
sapply(names(pmds), function(x) {export(pmds[[x]], con=paste0(tmp, "/", x), format="BED")})
pmds.m <- read.delim(pipe(paste0("multiIntersectBed -header -i ", 
    paste0(list.files(path=tmp, full.names=T), collapse=' '))), header=T)
colnames(pmds.m) <- gsub("^X.+\\.PD","PD", colnames(pmds.m))
pmds.m <- makeGRangesFromDataFrame(pmds.m, keep.extra.columns=T)

# take out chrY, we have only females
pmds.m <- pmds.m[seqnames(pmds.m) != "chrY"]

# calculate number of bases inside PMDs, per patient
pmds.bases <- sapply(colnames(mcols(pmds.m))[grep("^PD", colnames(mcols(pmds.m)))], 
    function(x) {sum(width(pmds.m[mcols(pmds.m)[[x]]==1]))})

# same, but as fraction of total genome size
pmds.frac <- pmds.bases/3.2e09

# boxplot/stripchart PMD bases as fraction of the genome
pdf("boxplot_PMD_bases_per_patient.pdf") #### Figure 2A ####
par(mfrow=c(3,4), mar=c(5,5,4,0))
boxplot(pmds.frac, boxwex=0.5, names="patients", ylab="fraction of genome\n covered by PMDs")
stripchart(pmds.frac, vertical=T, pch=1, col="red", method="jitter",add=T)
dev.off()

# determine PMD frequency: the number of patients in which a genomic regions is (inside) a PMD
patient.cols <- grep("^PD", colnames(mcols(pmds.m)))
pmds.m$pmd.freq <- apply(as.matrix(mcols(pmds.m[,patient.cols])), 1, sum)

# per PMD frequency, determine the number of bases
pmds.freq.bases <- sapply(as.character(unique(pmds.m$pmd.freq)[order(unique(pmds.m$pmd.freq))]), function(x) {
    sum(as.numeric(width(pmds.m[pmds.m$pmd.freq >= as.numeric(x)])))}, USE.NAMES=T)

# same, but as a fraction of total genome size
pmds.freq.frac <- pmds.freq.bases/3.2e09

# dotplot PMD bases (fraction of genome) as a function of PMD frequency
pdf("dotplot_PMD_bases_vs_PMD_frequency.pdf") #### Figure 2B ####
par(mfrow=c(3,3), mar=c(5,5,4,0))
plot(pmds.freq.frac, pch=20, col="red", xlab= "PMD frequency\n(number of tumors)", 
    ylab="fraction of genome\ncovered by common PMDs")
dev.off()
