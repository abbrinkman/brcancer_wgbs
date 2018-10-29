library(GenomicRanges)
library(ggplot2)
library(reshape2)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19)
library(gplots)

# get the PMDs
pmds <- get(load("~/BiSeq_BASIS/PMDs_all_samples/pmdmeth.basis.RData"))
pmds <- pmds[!c(names(pmds) %in% c("HMEC","MCF7"))]



## use two different annotations: 1: ENSEMBL (from the 'fpkm' object); 2: RefSeq

# get the refseq gene annotation
session <- browserSession()
genome(session) <- "hg19"
query <- ucscTableQuery(session, "RefSeq Genes")
refseq <- getTable(query)
colnames(refseq) <- gsub("chrom","chr",gsub("^tx","",colnames(refseq)))
refseq <- makeGRangesFromDataFrame(refseq, keep.extra.columns=T)

# get the used gene annotation for the RNA-seq cohort
fpkm <- get(load("~/Infinium/BASIS_data_all_TN_QC_analysis/fpkm.gr.RData"))

# clean up the ENSEMBL annotation
out.1 <- grepl("^[ABF][JECPLX]\\d{5,}\\.\\d{1,2}$", fpkm$Name)   # e.g. AC008175.6, AE000662.5, AF250324.1
out.2 <- grepl("^RP\\d{1,2}-.+\\.\\d{1,2}$", fpkm$Name)          # e.g. RP11-501J20.2
out.3 <- grepl("^Z\\d{5}\\.\\d{1,2}$", fpkm$Name)                # e.g. Z84723.1
out.4 <- grepl("^CT[ABCD]-\\d+[A-Z]\\d+", fpkm$Name)             # e.g. CTD-2135D7.1
out.5 <- grepl("^AJ\\d{6}\\.\\d{1,2}$", fpkm$Name)               # e.g. AJ003147.2
out.6 <- grepl("^RNU\\d-\\d{1,2}", fpkm$Name)                    # e.g. RNU7-54P
out.7 <- grepl("^LA16c-\\d{1,2}", fpkm$Name)                     # e.g. LA16c-3G11.5
out.8 <- grepl("^XXbac-", fpkm$Name)                             # e.g. XXbac-B444P24.8
out.9 <- grepl("^LL22.+\\-.+\\.", fpkm$Name)                     # e.g. LL22NC03-75H12.2
out.10 <- grepl("^.+-.+\\.", fpkm$Name)                          # e.g. KB-1592A4.13, WI2-85898F10.1
fpkm <- fpkm[!c(out.1 | out.2 | out.3 | out.4 | out.5 | out.6 | out.7 | out.8 | out.9 | out.10)]


# convert into one nonoverlapping set
fpkm.1 <- reduce(fpkm, ignore.strand=T)
refseq.1 <- reduce(refseq, ignore.strand=T)

# convert into one nonoverlapping set, but keep directionality (allows overlapping genes)
fpkm.2 <- reduce(fpkm, ignore.strand=F)
refseq.2 <- reduce(refseq, ignore.strand=F)

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
hg19 <- hg19[!grepl("chrX", seqnames(hg19))]
hg19 <- hg19[!grepl("chrY", seqnames(hg19))]

# remove centromers from genome
centro <- read.delim(pipe("cat ~/Datarepository/hg19_cytoBand.bed |grep -e acen -e gvar -e stalk"),
    header=F, sep="\t")
colnames(centro) <- c("chr","start","end","arm","type")
centro <- makeGRangesFromDataFrame(centro, keep.extra.columns=T)
hg19 <- setdiff(hg19, centro)

# make the tiles
tile.width <- 5000 # try various tile widths: 5 kb gives good resolution
tiles <- unlist(tile(hg19, width=tile.width))

# count the pmds per tile
# tiles:  |=======1========||=======2========||=======3========|
# PMDs:      <------------------------------------->
# 2 is a hit ("in" PMD), 1 and 3 not
tiles$pmd.freq <- apply(sapply(pmds, function(x) {ifelse(overlapsAny(tiles, x, type="within"), 1, 0)}), 1, sum)
tiles$pmd.freq.bin <- cut(tiles$pmd.freq, breaks=c(-0.01, seq(0,30, by=3)))

# split the tiles into PMDfreq classes
tiles.freq <- split(tiles, tiles$pmd.freq)
tiles.freq.bin <- split(tiles, tiles$pmd.freq.bin)

# for each PMD frequency (bin) calculate the 
# - total tile bases
# - total gene bases
getTileBases_geneBases <- function(gr.tiles, gr.genes) { 
  gene.bp <- sum(width(intersect(gr.tiles, gr.genes)))
  tile.bp <- sum(width(gr.tiles))
  data.frame("gene.bp"=gene.bp, "tile.bp"=tile.bp, "frac"=gene.bp/tile.bp)
}

## RNA-seq (ENSEMBL) annotation
# get the basecount for PMDfreq
d1 <- do.call(rbind, lapply(tiles.freq, function(x) {getTileBases_geneBases(x, fpkm.1)}))
d1$pmd.freq <- factor(rownames(d1), levels=as.character(rownames(d1)))

# get the basecount for PMDfreq bins
d2 <- do.call(rbind, lapply(tiles.freq.bin, function(x) {getTileBases_geneBases(x, fpkm.1)}))
d2$pmd.freq.bin <- rownames(d2)
d2$pmd.freq.bin <- factor(d2$pmd.freq.bin, levels=d2$pmd.freq.bin)

## RefSeq annotation
# get the basecount for PMDfreq
d3 <- do.call(rbind, lapply(tiles.freq, function(x) {getTileBases_geneBases(x, refseq.1)}))
d3$pmd.freq <- factor(rownames(d3), levels=as.character(rownames(d3)))

# get the basecount for PMDfreq bins
d4 <- do.call(rbind, lapply(tiles.freq.bin, function(x) {getTileBases_geneBases(x, refseq.1)}))
d4$pmd.freq.bin <- rownames(d4)
d4$pmd.freq.bin <- factor(d4$pmd.freq.bin, levels=d4$pmd.freq.bin)

# for each PMD frequency (bin) calculate the
# - number of genes (regard the tss start sites to avoid over-counting)
getTileGenesNr <- function(gr.tiles, gr.genes) {
  ov <- findOverlaps(promoters(gr.genes, upstream=0, downstream=1), gr.tiles, type="within")
  gene.nr <- length(subjectHits(ov))
  gene.nr
}
d1$gene.nr <- sapply(tiles.freq, function(x) {getTileGenesNr(x, fpkm.2)})
d2$gene.nr <- sapply(tiles.freq.bin, function(x) {getTileGenesNr(x, fpkm.2)})
d3$gene.nr <- sapply(tiles.freq, function(x) {getTileGenesNr(x, refseq.2)})
d4$gene.nr <- sapply(tiles.freq.bin, function(x) {getTileGenesNr(x, refseq.2)})


# for each PMD frequency (bin) calculate the
# - gene length(regard the tss start sites to avoid over-counting)
getTileGeneLength <- function(gr.tiles, gr.genes) {
  ov <- findOverlaps(promoters(gr.genes, upstream=0, downstream=1), gr.tiles, type="within")
  gene.ln <- width(gr.genes[queryHits(ov)])
  gene.ln
}
d5 <- melt(lapply(tiles.freq, function(x) {getTileGeneLength(x, fpkm.2)}), value.name="gene.length")
d5$pmd.freq <- factor(d5$L1, levels=unique(d5$L1))
d6 <- melt(lapply(tiles.freq.bin, function(x) {getTileGeneLength(x, fpkm.2)}), value.name="gene.length")
d6$pmd.freq.bin <- factor(d6$L1, levels=unique(d6$L1))
d7 <- melt(lapply(tiles.freq, function(x) {getTileGeneLength(x, refseq.2)}), value.name="gene.length")
d7$pmd.freq <- factor(d7$L1, levels=unique(d7$L1))
d8 <- melt(lapply(tiles.freq.bin, function(x) {getTileGeneLength(x, refseq.2)}), value.name="gene.length")
d8$pmd.freq.bin <- factor(d8$L1, levels=unique(d8$L1))

# also determine the number of genes per PMD frequency class, independent of the tiles
fpkm$pmd.freq <- apply(sapply(pmds, function(x) {overlapsAny(fpkm, x, type="within")}), 1, sum)
fpkm$pmd.freq.bin <- cut(fpkm$pmd.freq, breaks=c(-0.01, seq(0,30, by=3)))
d9 <- as.data.frame(table(fpkm$pmd.freq))
colnames(d9) <- c("pmd.freq","count")
d9$frac <- d9$count/sum(d9$count)

d10 <- as.data.frame(table(fpkm$pmd.freq.bin))
colnames(d10) <- c("pmd.freq.bin","count")
d10$frac <- d10$count/sum(d10$count)
d10$pmd.freq.bin <- factor(d10$pmd.freq.bin, levels=d10$pmd.freq.bin)

refseq$pmd.freq <- apply(sapply(pmds, function(x) {overlapsAny(refseq, x, type="within")}), 1, sum)
refseq$pmd.freq.bin <- cut(refseq$pmd.freq, breaks=c(-0.01, seq(0,30, by=3)))
d11 <- as.data.frame(table(refseq$pmd.freq))
colnames(d11) <- c("pmd.freq","count")
d11$frac <- d11$count/sum(d11$count)

d12 <- as.data.frame(table(refseq$pmd.freq.bin))
colnames(d12) <- c("pmd.freq.bin","count")
d12$frac <- d12$count/sum(d12$count)
d12$pmd.freq.bin <- factor(d12$pmd.freq.bin, levels=d12$pmd.freq.bin)

## plotting
# gene density (bp/bp)
p1 <- ggplot(d1, aes(pmd.freq, frac)) +
      geom_bar(stat="identity", aes(fill=pmd.freq), color="black") +
      scale_fill_manual(values=colorpanel(nrow(d1), "white", "red")) +
      theme_classic() + 
      theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"), 
          axis.text.x=element_text(angle=45, hjust=1, vjust=1), legend.position="none") + 
      ylab("relative gene fraction (bp/bp)") + 
      xlab("PMD frequency") + 
      ggtitle("RNA-seq annotation (ENSEMBL)")
ggsave(p1, file="relative_gene_density_ENSEMBL_vs_PMDfrequency.pdf", width=4, height=2.8, scale=1.5)

p2 <- ggplot(d2, aes(pmd.freq.bin, frac)) +
      geom_bar(stat="identity", aes(fill=pmd.freq.bin), color="black") +
      scale_fill_manual(values=colorpanel(nrow(d2), "white", "red")) +
      theme_classic() +
      theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1), legend.position="none") +
      ylab("relative gene fraction (bp/bp)") +
      xlab("PMD frequency") +
      ggtitle("RNA-seq annotation (ENSEMBL)")
ggsave(p2, file="relative_gene_density_ENSEMBL_vs_PMDfrequencyBin.pdf", width=4, height=2.8, scale=1) #### Supplemental Figure 3A ####

p3 <- ggplot(d3, aes(pmd.freq, frac)) +
      geom_bar(stat="identity", aes(fill=pmd.freq), color="black") +
      scale_fill_manual(values=colorpanel(nrow(d3), "white", "red")) +
      theme_classic() +
      theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1), legend.position="none") +
      ylab("relative gene fraction (bp/bp)") +
      xlab("PMD frequency") +
      ggtitle("RefSeq annotation")
ggsave(p3, file="relative_gene_density_RefSeq_vs_PMDfrequency.pdf", width=4, height=2.8, scale=1.5)

p4 <- ggplot(d4, aes(pmd.freq.bin, frac)) +
      geom_bar(stat="identity", aes(fill=pmd.freq.bin), color="black") +
      scale_fill_manual(values=colorpanel(nrow(d4), "white", "red")) +
      theme_classic() +
      theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1), legend.position="none") +
      ylab("relative gene fraction (bp/bp)") +
      xlab("PMD frequency") +
      ggtitle("RefSeq annotation")
ggsave(p4, file="relative_gene_density_RefSeq_vs_PMDfrequencyBin.pdf", width=4, height=2.8, scale=1)

# number of genes
p5 <- ggplot(d1, aes(pmd.freq, gene.nr/1000)) +
      geom_bar(stat="identity", aes(fill=pmd.freq), color="black") +
      scale_fill_manual(values=colorpanel(nrow(d1), "white", "red")) +
      theme_classic() +
      theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1), legend.position="none") +
      ylab("#genes (x1000)") +
      xlab("PMD frequency") +
      ggtitle("RNA-seq annotation (ENSEMBL)")
ggsave(p5, file="numberOfGenes_ENSEMBL_vs_PMDfrequency.pdf", width=4, height=2.8, scale=1.5)


p6 <- ggplot(d2, aes(pmd.freq.bin, gene.nr/1000)) +
      geom_bar(stat="identity", aes(fill=pmd.freq.bin), color="black") +
      scale_fill_manual(values=colorpanel(nrow(d2), "white", "red")) +
      theme_classic() +
      theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1), legend.position="none") +
      ylab("#genes (x1000)") +
      xlab("PMD frequency") +
      ggtitle("RNA-seq annotation (ENSEMBL)")
ggsave(p6, file="numberOfGenes_ENSEMBL_vs_PMDfrequencyBin.pdf", width=4, height=2.8, scale=1)

p7 <- ggplot(d3, aes(pmd.freq, gene.nr/1000)) +
      geom_bar(stat="identity", aes(fill=pmd.freq), color="black") +
      scale_fill_manual(values=colorpanel(nrow(d3), "white", "red")) +
      theme_classic() +
      theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1), legend.position="none") +
      ylab("#genes (x1000)") +
      xlab("PMD frequency") +
      ggtitle("RefSeq annotation")
ggsave(p7, file="numberOfGenes_RefSeq_vs_PMDfrequency.pdf", width=4, height=2.8, scale=1.5)


p8 <- ggplot(d4, aes(pmd.freq.bin, gene.nr/1000)) +
      geom_bar(stat="identity", aes(fill=pmd.freq.bin), color="black") +
      scale_fill_manual(values=colorpanel(nrow(d4), "white", "red")) +
      theme_classic() +
      theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1), legend.position="none") +
      ylab("#genes (x1000)") +
      xlab("PMD frequency") +
      ggtitle("RefSeq annotation")
ggsave(p8, file="numberOfGenes_RefSeq_vs_PMDfrequencyBin.pdf", width=4, height=2.8, scale=1)

# gene lengths
p9 <- ggplot(d5, aes(pmd.freq, gene.length/1000)) +
      geom_boxplot(aes(fill=pmd.freq), outlier.shape=NA) +
      scale_fill_manual(values=colorpanel(length(levels(d5$pmd.freq)), "white", "red")) +
      theme_classic() +
      coord_cartesian(ylim=c(0, 150)) +
      theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1), legend.position="none") +
      ylab("gene length (kb)") +
      xlab("PMD frequency") +
      ggtitle("RNA-seq annotation (ENSEMBL)")
ggsave(p9, file="gene_length_ENSEMBL_vs_PMDfrequency.pdf", width=4, height=2.8, scale=1.5)

p10 <- ggplot(d6, aes(pmd.freq.bin, gene.length/1000)) +
      geom_boxplot(aes(fill=pmd.freq.bin), outlier.shape=NA) +
      scale_fill_manual(values=colorpanel(length(levels(d6$pmd.freq.bin)), "white", "red")) +
      theme_classic() +
      coord_cartesian(ylim=c(0, 150)) +
      theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1), legend.position="none") +
      ylab("gene length (kb)") +
      xlab("PMD frequency") +
      ggtitle("RNA-seq annotation (ENSEMBL)")
ggsave(p10, file="gene_length_ENSEMBL_vs_PMDfrequencyBin.pdf", width=4, height=2.8, scale=1)

p11 <- ggplot(d7, aes(pmd.freq, gene.length/1000)) +
      geom_boxplot(aes(fill=pmd.freq), outlier.shape=NA) +
      scale_fill_manual(values=colorpanel(length(levels(d7$pmd.freq)), "white", "red")) +
      theme_classic() +
      coord_cartesian(ylim=c(0, 150)) +
      theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1), legend.position="none") +
      ylab("gene length (kb)") +
      xlab("PMD frequency") +
      ggtitle("RefSeq annotation")
ggsave(p11, file="gene_length_RefSeq_vs_PMDfrequency.pdf", width=4, height=2.8, scale=1.5)

p12 <- ggplot(d8, aes(pmd.freq.bin, gene.length/1000)) +
      geom_boxplot(aes(fill=pmd.freq.bin), outlier.shape=NA) +
      scale_fill_manual(values=colorpanel(length(levels(d8$pmd.freq.bin)), "white", "red")) +
      theme_classic() +
      coord_cartesian(ylim=c(0, 150)) +
      theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1), legend.position="none") +
      ylab("gene length (kb)") +
      xlab("PMD frequency") +
      ggtitle("RefSeq annotation")
ggsave(p12, file="gene_length_RefSeq_vs_PMDfrequencyBin.pdf", width=4, height=2.8, scale=1)

# plot the gene fraction, but gene-based (not count-based)
p13 <- ggplot(d9, aes(pmd.freq, frac)) +
      geom_bar(stat="identity", aes(fill=pmd.freq), color="black") + 
      scale_fill_manual(values=colorpanel(length(levels(d9$pmd.freq)), "white", "red")) +
      theme_classic() +
      theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1), legend.position="none") +
      ylab("gene fraction") +
      xlab("PMD frequency") +
      ggtitle("RNA-seq annotation (ENSEMBL)")
ggsave(p13, file="gene_fraction_ENSEMBL_vs_PMDfrequency.pdf", width=4, height=2.8, scale=1.5)

p14 <- ggplot(d10, aes(pmd.freq.bin, frac)) +
      geom_bar(stat="identity", aes(fill=pmd.freq.bin), color="black") +
      scale_fill_manual(values=colorpanel(length(levels(d10$pmd.freq.bin)), "white", "red")) +
      theme_classic() +
      theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1), legend.position="none") +
      ylab("gene fraction") +
      xlab("PMD frequency") +
      ggtitle("RNA-seq annotation (ENSEMBL)")
ggsave(p14, file="gene_fraction_ENSEMBL_vs_PMDfrequencyBin.pdf", width=4, height=2.8, scale=0.9) #### Figure 2E ####

p15 <- ggplot(d11, aes(pmd.freq, frac)) +
      geom_bar(stat="identity", aes(fill=pmd.freq), color="black") +
      scale_fill_manual(values=colorpanel(length(levels(d11$pmd.freq)), "white", "red")) +
      theme_classic() +
      theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1), legend.position="none") +
      ylab("gene fraction") +
      xlab("PMD frequency") +
      ggtitle("RefSeq annotation")
ggsave(p15, file="gene_fraction_RefSeq_vs_PMDfrequency.pdf", width=4, height=2.8, scale=1.5)

p16 <- ggplot(d12, aes(pmd.freq.bin, frac)) +
      geom_bar(stat="identity", aes(fill=pmd.freq.bin), color="black") +
      scale_fill_manual(values=colorpanel(length(levels(d12$pmd.freq.bin)), "white", "red")) +
      theme_classic() +
      theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1), legend.position="none") +
      ylab("gene fraction") +
      xlab("PMD frequency") +
      ggtitle("RefSeq annotation")
ggsave(p16, file="gene_fraction_RefSeq_vs_PMDfrequencyBin.pdf", width=4, height=2.8, scale=0.9)


