library(GenomicRanges)
library(ggplot2)
library(reshape2)
library(stringr)

pmdmeth <- get(load("~/BiSeq_BASIS/PMDs_all_samples/pmdmeth.basis.RData"))
pmdmeth <- pmdmeth[!names(pmdmeth) %in% c("HMEC", "MCF7")]

# calculate the mean for each sample
getMean <- function(x) {
  sapply(x, function(y) {mean(y@elementMetadata@listData[[1]], na.rm=T)})
}
pmdmeth <- getMean(pmdmeth)

source("~/tools/BASIS_common_functions.txt.R")
clinj  <- clinj[paste0(clinj$sample_name, "a") %in% names(pmdmeth),]

d1  <- data.frame(
  "aberrant_cells"=clinj$AberrantCells[match(names(pmdmeth), paste0(clinj$sample_name, "a"))],
  "PMDmeth"=pmdmeth)

pdf("pmd_mean_meth_vs_aberrant_cell_fraction.pdf") ##### Supplemental Figure 1E #####
par(mfrow=c(3,3))
plot(PMDmeth ~ aberrant_cells, data=d1, ylim=c(0,1), xlim=c(0,1),
     ylab="mean PMD meth", xlab="aberrant cell fraction",
     pch=20, col="blue")
dev.off()

