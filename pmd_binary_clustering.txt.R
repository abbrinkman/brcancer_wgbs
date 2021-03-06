library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(rtracklayer)
library(cluster)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ggthemes)
library(ggdendro)
library(grid)
library(gridExtra)


# get the hg19 genome
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
hg19 <- hg19[!grepl("chrX", seqnames(hg19))]
hg19 <- hg19[!grepl("chrY", seqnames(hg19))]

# remove centromers from genome
centro <- read.delim(pipe("cat ~/Datarepository/hg19_cytoBand.bed |grep -e acen -e gvar -e stalk"),
    header=F, sep="\t")
colnames(centro) <- c("chr","start","end","arm","type")
centro <- makeGRangesFromDataFrame(centro, keep.extra.columns=T)
hg19 <- setdiff(hg19, centro)

# make the hg19 tiles
tile.size <- 5000 # basepairs
tiles <- unlist(tile(hg19, width=tile.size))


# get PMDs
pmdmeth.select <- get(load("~/BiSeq_BASIS/PMDs_other_celltypes/pmdmeth.select.RData"))

# select only the tumors for this
pmdmeth.tumors <- pmdmeth.select[grep("^tumor_", names(pmdmeth.select))]

# score whether a PMD completely overlaps a tile (1) or not (0)
scorePMDs <- function(PMD_GR_LIST) { # LABEL can be used for cluster labelling later
  for (i in names(PMD_GR_LIST)) {
    message(i)
    pmds <- PMD_GR_LIST[[i]]
    scores <- ifelse(countOverlaps(tiles, pmds, type="within") >= 1, 1 ,0)
    pmd.scores@elementMetadata@listData[[i]] <- scores
  }
  pmd.scores
}

pmd.scores <- tiles
pmd.scores <- scorePMDs(pmdmeth.tumors)

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
load("~/BiSeq_BASIS/PMDs_other_celltypes/PMD_clustering/pmd.scores.jacc_1.RData")

# hierarchical clustering
hc1 <- hclust(as.dist(1-pmd.scores.jacc))

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
#hc1.d$type <- gsub("^PD.+","BRCA",sapply(strsplit(as.character(hc1.d$name), "_"), function(x) {x[2]}))
hc1.d$type <- sapply(strsplit(sapply(strsplit(as.character(hc1.d$name), "_"), 
    function(x) {x[3]}), "\\."), function(x) {x[1]})
hc1.d$name <- factor(hc1.d$name, levels=as.character(hc1.d$name))
#hc1.d$ER <- clinj$ER.FPKM[match(gsub("a$","", gsub("^basis_","",hc1.d$name)), clinj$sample_name)]
hc1.d$ER <- clinj$ER.FPKM[match(gsub("a$","", gsub("^.+\\.","",hc1.d$name)), clinj$sample_name)]
cluster.colors <- c(brewer.pal(12,"Paired"),  
    c("darkblue", "black", "pink", "yellow", "grey", "magenta", "darkslategrey") )

p <- ggplot(data=hc1.d, aes(name, value))
p <- p + geom_tile(aes(fill=type, colour=type))
p <- p + scale_fill_manual(values=cluster.colors)
p <- p + scale_colour_manual(values=cluster.colors)
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

pdf("hclust_PMDs_other_celltypes_jaccard_withColors.pdf", onefile=T) ##### Figure 4B #####
n <- 30
grid.arrange(dd,
             p1,
             p3,
             nrow=n,
             heights=c(0.2, rep(0.8/(n-1), n-1)))
p2
dev.off()

