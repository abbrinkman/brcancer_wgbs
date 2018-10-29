library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ConsensusClusterPlus)
library(rtracklayer)
library(cluster)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ggthemes)
library(ggdendro)
library(grid)
library(gridExtra)


# get hg19 tiles
tile.size <- 5000 # basepairs
chr.lengths <- seqlengths(Hsapiens)
chr.lengths <- chr.lengths[1:22]
tiles <- tileGenome(chr.lengths, tilewidth=tile.size,
        cut.last.tile.in.chrom=T)

# get BASIS PMDs
basis.files <- list.files(pattern="PMDs_PD[0-9]+a_noCent.bed$", path="~/BiSeq_BASIS/PMDs_all_samples/", full.names=T)
basis.files <- basis.files[!grepl("PD9590", basis.files)]
basis.files <- as.list(basis.files)
names(basis.files) <- sapply(strsplit(basename(unlist(basis.files)), "_"), function(x) {x[2]})

pmds.basis <- lapply(basis.files, function(x) {read.table(x, col.names=c("chr","start","end"))})
pmds.basis <- lapply(pmds.basis, makeGRangesFromDataFrame)

# get TCGA PMDs
pmds.tcga <- get(load("~/BiSeq_BASIS/PMDs_other_celltypes/tcga/tcga.pmds.RData"))

# take out the one sample for which PMD calling was faulty
pmds.tcga <- pmds.tcga[!names(pmds.tcga) %in% "BLCA_NIC1254A93_tumor"]

# take out the normals
pmds.tcga <- pmds.tcga[!grepl("normal", names(pmds.tcga))]

# get lymphoma PMDs
pmds.lymph <- get(load("~/BiSeq_BASIS/PMDs_other_celltypes/lymphoma/lymph.pmds.RData"))

# take out the normal germinal centre B-cells from the lymphoma PMDs
pmds.lymph <- pmds.lymph[-grep("GCB",names(pmds.lymph))]

# remove the two lymphoma samples for which PMD calling was bad
pmds.lymph <- pmds.lymph[!names(pmds.lymph) %in% c("FL_4158726", "FL_4159170")]


# score whether a PMD completely overlaps a tile (1) or not (0)
scorePMDs <- function(PMD_GR_LIST, LABEL) { # LABEL can be used for cluster labelling later
  for (i in names(PMD_GR_LIST)) {
    message(i)
    pmds <- PMD_GR_LIST[[i]]
    scores <- ifelse(countOverlaps(tiles, pmds, type="within") >= 1, 1 ,0)
    pmd.scores@elementMetadata@listData[[paste0(LABEL, "_", i)]] <- scores
  }
  pmd.scores
}

pmd.scores <- tiles

pmd.scores <- scorePMDs(pmds.basis, "basis")
pmd.scores <- scorePMDs(pmds.tcga, "tcga")
pmd.scores <- scorePMDs(pmds.lymph, "lymph")


# create distance table from Jaccard index
pmd.scores.jacc <- matrix(nc=ncol(mcols(pmd.scores)), nr=ncol(mcols(pmd.scores)))
colnames(pmd.scores.jacc) <- colnames(mcols(pmd.scores))
rownames(pmd.scores.jacc) <- colnames(mcols(pmd.scores))
for (r in 1:nrow(pmd.scores.jacc)) {
  for (c in 1:ncol(pmd.scores.jacc)) {
    message(rownames(pmd.scores.jacc)[r], " vs. ", colnames(pmd.scores.jacc)[c])
    jacc <- sum(mcols(pmd.scores[,r])==mcols(pmd.scores[,c]))/length(pmd.scores)
    pmd.scores.jacc[r,c] <- jacc
  }
}

# hierarchical clustering
hc1 <- hclust(as.dist(1-pmd.scores.jacc))

# get the silhouette widths for this clustering
sapply(2:5, function(X) {summary(silhouette(cutree(hc1, k=X), 1-pmd.scores.jacc))$avg.width})
# [1] 0.2467066 0.1385606 0.1069243 0.1123475 -> 2 rather than 3 clusters

# plot hierarchical clustering
ddata <- dendro_data(hc1, type="rectangle")
dd <- ggplot(segment(ddata))
dd <- dd + geom_segment(aes(x = x, y = y, xend = xend, yend = yend))
dd <- dd + theme(axis.text.x = element_blank(),
                 axis.text.y = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.ticks.y = element_blank(),
                 axis.title.x = element_blank(),
                 axis.title.y = element_blank(),
                 legend.title = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 legend.position = "none",
                 plot.margin=unit(c(0.2,0.26,-0.35,0.34), "cm"))

# plot clinical parameters
source("~/tools/BASIS_common_functions.txt.R")
hc1.d <- data.frame("name"=hc1$labels[hc1$order])
hc1.d$value <- rep(0, nrow(hc1.d))
hc1.d$type <- gsub("^PD.+","BRCA",sapply(strsplit(as.character(hc1.d$name), "_"), function(x) {x[2]}))
hc1.d$name <- factor(hc1.d$name, levels=as.character(hc1.d$name))
hc1.d$ER <- clinj$ER.FPKM[match(gsub("a$","", gsub("^basis_","",hc1.d$name)), clinj$sample_name)]
p <- ggplot(data=hc1.d, aes(name, value))
p <- p + geom_tile(aes(fill=type, colour=type))
p <- p + scale_fill_manual(values=c(brewer.pal(9,"Set1"), c("darkblue", "black") ))
p <- p + scale_colour_manual(values=c(brewer.pal(9,"Set1"), c("darkblue", "black") ))
p <- p + theme_classic()

p1 <- p + theme(axis.text.x = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.x = element_blank(),
               axis.ticks.y = element_blank(),
               axis.title.x = element_blank(),
               axis.title.y = element_blank(),
               legend.title = element_blank(),
               legend.position = "none",
               plot.margin=unit(c(-0.3,0.9,-0.3,1), "cm"))
 
p2 <- p + theme(axis.text.x = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.x = element_blank(),
               axis.ticks.y = element_blank(),
               axis.title.x = element_blank(),
               axis.title.y = element_blank(),
               plot.margin=unit(c(-0.3,0.9,-0.3,1), "cm"))

p3 <- ggplot(data=hc1.d, aes(name, value))
p3 <- p3 + geom_tile(aes(fill=ER, colour=ER))
p3 <- p3 + scale_fill_brewer(palette="Set1")
p3 <- p3 + scale_colour_brewer(palette="Set1")
p3 <- p3 + theme_classic()
p3 <- p3 + theme(axis.text.x = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.x = element_blank(),
               axis.ticks.y = element_blank(),
               axis.title.x = element_blank(),
               axis.title.y = element_blank(),
               legend.title = element_blank(),
               legend.position = "none",
               plot.margin=unit(c(-0.3,0.9,-0.3,1), "cm"))

pdf("hclust_PMDs_other_celltypes_jaccard_withColors.pdf", onefile=T) #### Figure 4C ####
n <- 30
grid.arrange(dd,
             p1,
             p3,
             nrow=n,
             heights=c(0.2, rep(0.8/(n-1), n-1)))
p2
dev.off()


