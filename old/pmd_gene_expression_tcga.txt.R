library(GenomicRanges)
library(rtracklayer)
library(reshape2)
library(ggplot2)

# get PMDs
pmds <- get(load("~/BiSeq_BASIS/PMDs_all_samples/pmdmeth.basis.RData"))
pmds$HMEC <- NULL
pmds$MCF7 <- NULL

# get CGIs for determining whether genes have CGI promoters
cgi <- get(load("~/BiSeq_BASIS/PMDs_all_samples/cgimeth.basis.RData"))[[1]][,0]

# get the TCGA 'BRCA' RNAseq expression data,
# downloaded from the TCGA portal (level 3, =readcounts)
# Example file fragment from 1 patient:
# gene    raw_counts      median_length_normalized        RPKM
# ?|100130426     0       0       0
# ?|100133144     87      3.59354304635762        0.79460879941902
# ?|100134869     32      1.00376411543287        0.221494380870373
# ?|10357 294     23.1132075471698        5.10024767460536
# ?|10431 7819    349.627906976744        77.1632411408924
# ?|136542        0       0       0
# ?|155060        781     13.4048931771192        2.96930409703286

# get the files
md <- read.table(
    "~/BRCA_RNAseq/METADATA/UNC__IlluminaHiSeq_RNASeq/unc.edu_BRCA.IlluminaHiSeq_RNASeq.2.sdrf.txt", 
    sep="\t", header=T)
colnames(md) <- gsub('\\.$','', gsub('\\.{1,}','\\.', colnames(md)))
dir <-  "~/BRCA_RNAseq/RNASeq/UNC__IlluminaHiSeq_RNASeq/Level_3"
files <- list.files(path=dir, pattern="gene.quantification.txt", full.names=T)
files <- as.list(files)
names(files) <- sapply(strsplit(sapply(strsplit(unlist(files), "/"), 
  function(X) X[11]), "\\."), function(X) X[2])

# read the data
readData <- function(FILE) {
  message(FILE)
  d <- read.table(FILE, header=T)
  d$genename <- sapply(strsplit(as.character(d$gene), "\\|"), function(X) X[1])
  d$geneid <- sapply(strsplit(as.character(d$gene), "\\|"), function(X) X[2])
  d
}
l1 <- mclapply(files, readData, mc.cores=8)

# get gene position data (RefSeq table was downloaded from UCSC Table browser)
refseq <- read.delim(pipe("cat ~/Datarepository/RefSeq_hg19_v20140606.full.table |sed 's/^#//'"), 
    sep="\t", header=T)
refseq <- refseq[,colnames(refseq) %in% c("chrom","txStart","txEnd","name2", "strand")]
colnames(refseq) <- c("chr","strand","start","end", "name2")
refseq$start <- as.numeric(refseq$start)
refseq$end <- as.numeric(refseq$end)
refseq$strand <- as.character(refseq$strand)
refseq <- refseq[!refseq$name2=="",] 
refseq <- makeGRangesFromDataFrame(refseq, keep.extra.columns=T)
refseq$cgi.prom <- overlapsAny(promoters(refseq, upstream=1, downstream=1), cgi)
cgi.prom.genes <- unique(as.character(refseq$name2[refseq$cgi.prom]))
refseq.l <- split(refseq[, "name2"], refseq$name2)
refseq.l <- reduce(refseq.l)
refseq <- unlist(refseq.l)
refseq$name2 <- names(refseq)
refseq$cgi.prom <- ifelse(refseq$name2 %in% cgi.prom.genes, TRUE, FALSE)
refseq$pmdfreq <- apply(sapply(pmds, function(x) {overlapsAny(refseq, x, type="within")}), 1, sum)

expr.l <- lapply(l1, function(x) {x$RPKM[match(refseq$name2,x$genename)]})
refseq@elementMetadata@listData <- c(refseq@elementMetadata@listData, expr.l)
refseq <- refseq[complete.cases(as.matrix(mcols(refseq)))]

# log2 transformation 
refseq.log <- refseq
for (i in colnames(mcols(refseq.log))[grep("^TCGA", colnames(mcols(refseq.log)))]) {
  message(i)
  mcols(refseq.log)[[i]] <- log2(mcols(refseq.log)[[i]]+1)
}

# define tumors and normals
# Tumor types range from 01 - 09, normal types from 10 - 19 and control samples from 20 - 29
# take only one (the first) vial (=technical replicate): "A" at the "vial" position
tumors <- colnames(mcols(refseq))[grep(".+-.+-.+-0[1-9]A-.+-.+", colnames(mcols(refseq)))]
normals <- colnames(mcols(refseq))[grep(".+-.+-.+-1[0-9]A-.+-.+", colnames(mcols(refseq)))]

# calculate gene-means 
l1.tumor <- split(as.data.frame(mcols(refseq.log[,c("cgi.prom", tumors)])),refseq.log$pmdfreq)
l1.normal <- split(as.data.frame(mcols(refseq.log[,c("cgi.prom",normals)])),refseq.log$pmdfreq)
l2.tumor <- lapply(l1.tumor, function(x) {data.frame("cgi.prom"=x$cgi.prom, 
    "mean"=apply(x[,!colnames(x) %in% "cgi.prom"], 1, mean))})
l2.normal <- lapply(l1.normal, function(x) {data.frame("cgi.prom"=x$cgi.prom, 
    "mean"=apply(x[,!colnames(x) %in% "cgi.prom"], 1, mean))})

# plot all gene-means, tumor vs. normal, and per PMD frequency class
l3 <- list(tumor=l2.tumor, normal=l2.normal)
d3 <- melt(l3)
d3$variable <- NULL
colnames(d3) <- c("cgi.prom", "expression.log2","PMD.frequency","normal.tumor")
d3$PMD.frequency <- factor(d3$PMD.frequency, levels=unique(d3$PMD.frequency))
d3$PMD.frequency.bin <- cut(as.numeric(as.character(d3$PMD.frequency)), breaks=c(-0.01, seq(0,30, by=3)))


p1 <- ggplot(d3, aes(PMD.frequency.bin, expression.log2)) +
      geom_boxplot(aes(fill=normal.tumor), outlier.shape=NA) +
      scale_fill_manual(values=c("normal"="green","tumor"="red")) +
      coord_cartesian(ylim=c(0, 10)) +
      theme_classic() +
      theme(panel.background=element_blank(), axis.text=element_text(color="black"), 
          axis.text.x=element_text(angle=45, hjust=1, vjust=1))
ggsave(p1, file="boxplot_expression_gene-means_vs_PMDfreqBin_allPromoters.pdf", width=4, 
    height=2.8, scale=1.1) #### Figure 3H ####

# is gene expression associated with PMD frequency?
summary(lm(expression.log2 ~ as.numeric(PMD.frequency), data=d3))
# 
# Call:
# lm(formula = expression.log2 ~ as.numeric(PMD.frequency), data = d3)
# 
# Residuals:
#     Min      1Q  Median      3Q     Max
# -2.7205 -1.5284 -0.1218  1.1589 10.1559
# 
# Coefficients:
#                            Estimate Std. Error t value Pr(>|t|)
# (Intercept)                2.815944   0.010920   257.9   <2e-16 ***
# as.numeric(PMD.frequency) -0.095494   0.001218   -78.4   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.877 on 39974 degrees of freedom
# Multiple R-squared:  0.1333,    Adjusted R-squared:  0.1333


# is gene expression associated with tumor vs. healthy?
summary(lm(expression.log2 ~ normal.tumor, data=d3))
# 
# Call:
# lm(formula = expression.log2 ~ normal.tumor, data = d3)
# 
# Residuals:
#     Min      1Q  Median      3Q     Max
# -2.3811 -1.9639 -0.1428  1.3959  9.5644
# 
# Coefficients:
#                    Estimate Std. Error t value Pr(>|t|)
# (Intercept)        2.381054   0.014263 166.944   <2e-16 ***
# normal.tumortumor -0.004588   0.020170  -0.227     0.82
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 2.016 on 39974 degrees of freedom
# Multiple R-squared:  1.294e-06, Adjusted R-squared:  -2.372e-05
# F-statistic: 0.05174 on 1 and 39974 DF,  p-value: 0.8201 



# repeat the analysis using matched normal/tumor samples (i.e. from the same patient)
patients <- unique(sapply(strsplit(c(normals, tumors), '-',), function(x) {x[3]}))
patients <- as.list(patients)
names(patients) <- unlist(patients)
tumor.normal.pairs <- lapply(patients, function(x) {
    c("tumor"=tumors[grep(paste0(".+-.+-", x, "-01A-.+-.+.+"), tumors)], 
      "normal"=normals[grep(paste0(".+-.+-", x, "-11A-.+-.+.+"), normals)])
})
tumor.normal.pairs <- tumor.normal.pairs[lapply(tumor.normal.pairs, length)==2]
tumor.normal.ratios <- lapply(tumor.normal.pairs, function(x) {
   mcols(refseq.log)[[x["tumor"]]] - mcols(refseq.log)[[x["normal"]]]
})

# from here, split into genes with and without CGI-promoter

### CGI-promoters
l5.cgi.prom <- split( as.data.frame(mcols(refseq.log[refseq.log$cgi.prom,-c(1,3)])), 
    refseq.log$pmdfreq[refseq.log$cgi.prom])
l5.cgi.prom <- lapply(l5.cgi.prom, function(x) {
    tmp <- x ; colnames(tmp) <- gsub("\\.","-", colnames(tmp)) ; tmp})

# for each PMD frequence class, calculate the ratio tumor/normal per gene ('minus' since in log2)
l6.cgi.prom <- lapply(l5.cgi.prom, function(y) {
    lapply(tumor.normal.pairs, function(x) {
        y[,x["tumor"]] - y[,x["normal"]]})})

# calculate the mean ratio tumor/normal for each gene
l6.cgi.prom <- lapply(l6.cgi.prom, function(x) {apply(do.call(cbind, x), 1, mean)})
d6.cgi.prom <- melt(l6.cgi.prom)
colnames(d6.cgi.prom) <- c("log2.ratio.tumor.normal","PMD.frequency")
d6.cgi.prom$PMD.frequency <- factor(d6.cgi.prom$PMD.frequency, levels=unique(d6.cgi.prom$PMD.frequency))
d6.cgi.prom$PMD.frequency.bin <- cut(as.numeric(as.character(d6.cgi.prom$PMD.frequency)), 
    breaks=c(-0.01, seq(0,30, by=3)))

### non-CGI promoters
l5.noncgi.prom <- split( as.data.frame(mcols(refseq.log[!refseq.log$cgi.prom,-c(1,3)])),
    refseq.log$pmdfreq[!refseq.log$cgi.prom])
l5.noncgi.prom <- lapply(l5.noncgi.prom, function(x) {
    tmp <- x ; colnames(tmp) <- gsub("\\.","-", colnames(tmp)) ; tmp})

# for each PMD frequence class, calculate the ratio tumor/normal per gene ('minus' since in log2)
l6.noncgi.prom <- lapply(l5.noncgi.prom, function(y) {
    lapply(tumor.normal.pairs, function(x) {
        y[,x["tumor"]] - y[,x["normal"]]})})

# calculate the mean ratio tumor/normal for each gene
l6.noncgi.prom <- lapply(l6.noncgi.prom, function(x) {apply(do.call(cbind, x), 1, mean)})
d6.noncgi.prom <- melt(l6.noncgi.prom)
colnames(d6.noncgi.prom) <- c("log2.ratio.tumor.normal","PMD.frequency")
d6.noncgi.prom$PMD.frequency <- factor(d6.noncgi.prom$PMD.frequency, levels=unique(d6.noncgi.prom$PMD.frequency))
d6.noncgi.prom$PMD.frequency.bin <- cut(as.numeric(as.character(d6.noncgi.prom$PMD.frequency)),
     breaks=c(-0.01, seq(0,30, by=3)))

# plot these ratios (tumor/normal MATCHED)
p5.1 <- ggplot(d6.cgi.prom, aes(PMD.frequency.bin, log2.ratio.tumor.normal)) +
      geom_boxplot(outlier.shape=NA) +
      coord_cartesian(ylim=c(-1.75, 1.75)) +
      theme_classic() +
      theme(panel.background=element_blank(), axis.text=element_text(color="black"),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1))
ggsave(p5.1, file="boxplot_expression_mean-ratios_matched_vs_PMDfreqBin_CGIpromoters.pdf",
     width=4, height=2.8, scale=1.1) #### Figure 3I ####

p5.2 <- ggplot(d6.noncgi.prom, aes(PMD.frequency.bin, log2.ratio.tumor.normal)) +
      geom_boxplot(outlier.shape=NA) +
      coord_cartesian(ylim=c(-1.75, 1.75)) +
      theme_classic() +
      theme(panel.background=element_blank(), axis.text=element_text(color="black"),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1))
ggsave(p5.2, file="boxplot_expression_mean-ratios_matched_vs_PMDfreqBin_nonCGIpromoters.pdf",
     width=4, height=2.8, scale=1.1) #### Supplemental Figure 4B ####



##########################################################################
########## validation of individual genes from other analyses ############
##########################################################################

d9 <- as.data.frame(mcols(refseq.log[,-c(1:2)]))
rownames(d9) <- make.unique(refseq.log$name2)
tumor.status <- ifelse(colnames(d9) %in% gsub("-","\\.", tumors), "tumor", "normal")

system("mkdir -p validation_genes")

boxplotNormalTumor <- function(genename) {
  par(mar=c(5,4,5,2))
  pval <- signif(t.test(t(d9[genename,]) ~ tumor.status)$p.value,2)
  boxplot(t(d9[genename,]) ~ tumor.status, boxwex=0.5, las=1, ylab="FPKM (log2)") 
  p <- par("usr")
  lines(c(1, 2), rep(p[4]+(p[4]-p[3])*0.07, 2), xpd=T)
  lines(c(1,1), c(p[4]+(p[4]-p[3])*0.03, p[4]+(p[4]-p[3])*0.07), xpd=T)
  lines(c(2,2), c(p[4]+(p[4]-p[3])*0.03, p[4]+(p[4]-p[3])*0.07), xpd=T)
  text(1.5, p[4]+(p[4]-p[3])*0.15, labels=paste0("p = ", pval), xpd=T)
  text(p[1], p[4]+(p[4]-p[3])*0.27, labels=genename, xpd=T, cex=1.2, adj=c(0,0))
}

boxplotNormalTumorMatched <- function(genename) {
  par(mar=c(5,5,5,2))
  d.tmp <- t(sapply(tumor.normal.pairs, function(x) {d9[genename, gsub("-",".",x)]}))
  fold.diff <- unlist(d.tmp[,1])-unlist(d.tmp[,2]) # because log
  pval <- signif(t.test(fold.diff)$p.value,2)
  boxplot(fold.diff, boxwex=0.5, las=1, ylab="fold change\ntumor/normal (log2)")
  abline(h=0, col="darkgrey")
  p <- par("usr")
  text(1, p[4]+(p[4]-p[3])*0.15, labels=paste0("p = ", pval), xpd=T)
  text(p[1], p[4]+(p[4]-p[3])*0.27, labels=genename, xpd=T, cex=1.2, adj=c(0,0))

} 

val.genes <- c( # these genes are depicted in Suppl. Fig. 5A
  "FAT1",
  "GRIN2A",
  "PTPRC",
  "WT1",
  "EBF1",
  "IKZF1",
  "ATP2B3",
  "PRF1",
  "SFRP4",
  "ELF4",
  "BCORL1",
  "FOXL2",
  "GPC3",
  "GATA1",
  "MED12",
  "FOXO4",
  "RBM10",
  "ATRX",
  "STAG2",
  "IRF4",
  "RPL10",
  "PHF6",
  "RUNXT1",
  "BTK",
  "RHOH",
  "BCL9L",
  "SLC34A2",
  "APOBEC3B",
  "EGFR",
  "PREX2",
  "PDGFRA",
  "CLDN8",
  "PSTPIP1",
  "FGG",
  "AZGP1",
  "CD3G",
  "VGLL1",
  "TSPAN8",
  "NEUROD1",
  "THPO",
  "TNFRSF17",
  "CA6",
  "RBP4",
  "PIP",
  "CPA3",
  "TRPA1",
  "CEACAM6",
  "SERPINB7",
  "TMPRSS11D",
  "UGT8",
  "CXCL9",
  "INSM1",
  "CTAG2",
  "CXCL13",
  "AFP")

pdf("validation_genes/boxplots_TSG_validation_genes_normal_vs_tumor_unmatched_matched.pdf") #### Supplemental Figure 5D ####
par(mfrow=c(3,4))
for (i in val.genes) { message(i) ; try(boxplotNormalTumor(i)) ; try(boxplotNormalTumorMatched(i))}
dev.off()

