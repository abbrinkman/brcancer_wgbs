library(FactoMineR)
library(GenomicRanges)
library(reshape2)
library(ggplot2)
library(matrixStats)

# get the WGBS data (GRanges object, with CpG positions as rows, and
# colnames PDxxxx.T = total reads, PDxxxx.M = methylated reads)
meth <- get(load("~/BiSeq_BASIS/methcounts/meth.gr.RData"))
meth <- meth[!seqnames(meth) %in% c("chrX", "chrY")]
meth <- meth[,-grep("^[HM][MC][EF][C7]",colnames(mcols(meth)))]
meth <- meth[,-grep("PD9590a",colnames(mcols(meth)))]

# set a coverage cutoff, select only CpGs covered >= 10 reads in EACH sample
cov.thr <- 10
cols.t <- grep("^PD\\d+a\\.T", colnames(mcols(meth)))
cols.m <- grep("^PD\\d+a\\.M", colnames(mcols(meth)))
cov.rows <- apply(as.matrix(mcols(meth[,cols.t])), 1, function(x) {sum(x >= cov.thr)==length(cols.t)})
meth <- meth[cov.rows]
gc()
for (i in unique(gsub("\\.[TM]","", colnames(mcols(meth))))) {
  message(i)
  mcols(meth)[[i]] <- round(mcols(meth)[[paste0(i, ".M")]]/mcols(meth)[[paste0(i, ".T")]],3)
}
meth <- meth[,-c(cols.t,cols.m)]
gc()

# set a stdev cutoff: select top 5% most variable CpGs
meth.sd <- rowSds(as.matrix(mcols(meth)))
gc()
d1 <- as.data.frame(mcols(meth[meth.sd > quantile(meth.sd, p=0.95)]))
rownames(d1) <- paste(seqnames(meth[meth.sd > quantile(meth.sd, p=0.95)]),
    start(meth[meth.sd > quantile(meth.sd, p=0.95)]), sep="_")
d1 <- t(d1)

# add clinical data
source("~/tools/BASIS_common_functions.txt.R")
clinj <- clinj[match(rownames(d1), paste0(clinj$sample_name, "a")),]

clinj$paper_inclusion <- NULL
clinj$For_analyses <- NULL
clinj$Assay <- NULL
clinj$reps <- NULL
clinj$Nprobes <- NULL
clinj$wgbs <- NULL

# remove some of the unrelevant features from the clinical data
for (i in colnames(clinj)[grep("Signature.+", colnames(clinj))]) {
  clinj[,i] <- NULL
}
for (i in colnames(clinj)[grep("mutMatColor\\..+", colnames(clinj))]) {
  clinj[,i] <- NULL
}
for (i in colnames(clinj)[grep("ConsensusCluster[0-9]+", colnames(clinj))]) {
  clinj[,i] <- NULL
}
for (i in colnames(clinj)[grep("mutMatConsequence\\..+", colnames(clinj))]) {
  clinj[,i] <- NULL
}

# add rearrangement data
colnames(rearr) <- gsub("Cluster.ID","rearr.cluster", gsub("Signature","rearr.sig", colnames(rearr)))
clinj <- data.frame(clinj, rearr[match(clinj$sample_name, rearr$Sample.Name),
    !colnames(rearr) %in% c("Sample.Name","Total")])

# fix staging nomenclature inconsitencies
clinj$N_stage <- as.factor(gsub("X","x", clinj$N_stage))
clinj$M_stage <- as.factor(gsub("X","x", clinj$M_stage))


# add signatures data
colnames(sig) <- gsub("Cluster.ID","mut.cluster", gsub("Signature","mut.sig", colnames(sig)))
sig$mut.cluster[grepl("Singleton", sig$mut.cluster)] <- NA
sig$mut.cluster <- factor(
    sig$mut.cluster, levels=as.numeric(unique(sig$mut.cluster))[order(as.numeric(unique(sig$mut.cluster)))])
clinj <- data.frame(clinj, sig[match(clinj$sample_name, sig$Sample.Name),
    !colnames(sig) %in% c("Sample.Name","Accuracy")])
for (i in colnames(clinj)[grepl("^mut.sig.[0-9]+", colnames(clinj)) | grepl("^rearr.sig.[0-9]+", colnames(clinj))]){
  clinj[,i] <- as.numeric(clinj[,i])
}
colnames(clinj) <- gsub("^mutMat\\.","", colnames(clinj))

# combine methylation and clinical data
d2 <- cbind(d1, clinj)
d2$sample_name <- NULL

# remove clinical features (factor data) with too few or too many levels (e.g. gender)
for (i in colnames(d2)[grep("^chr", colnames(d2), invert=T)]) {
  if ( length(unique(d2[,i]))==1 | sum(is.na(d2[,i]))==nrow(d2)  ) {
    message(i)
    d2[,i] <- NULL
  }
}

# perform PCA on a small test set first (d4) 
# d4 <- d2[,c(sample(max(grep("^chr", colnames(d2))), 1000, replace=F), grep("^chr", colnames(d2), invert=T))]

# take the full set (top 5%)
d4 <- d2

# take out all features that have only one level or are only NA
cols.cg <- grep("^chr", colnames(d4))
cols.sup.quali <- unlist(sapply((max(cols.cg)+1):ncol(d4),
    function(x) {if (class(d4[,x])=="factor") {x}}))
cols.sup.quanti <- unlist(sapply((max(cols.cg)+1):ncol(d4),
     function(x) {if (class(d4[,x])=="numeric") {x}}))

# save some space
rm(meth, meth.gr) 
gc()

## start the PCA
# do the PCA (for the complete dataset (top 5%) + clinical data this takes ~16 hours)
# this full set contains 642858 CpGs and 192 clinical/genomic features
system("touch start.log")
d4.pca <- PCA(d4, scale.unit=TRUE, ncp=20,
      quali.sup=cols.sup.quali, quanti.sup=cols.sup.quanti, graph=F)
save(d4.pca, file="d4.pca.RData")

# determine the associations with the clinical data
d4.dimdesc <- dimdesc(d4.pca, axes=c(1:20), proba=1)
save(d4.dimdesc, file="d4.dimdesc.RData")
system("touch end.log")

# get pvals for the association of traits with the dimenstions (PCs)
getAssocTraits <- function(DIM) { 
  l <- d4.dimdesc[[DIM]]
  d1 <- as.data.frame(l$quanti[grep("^chr", rownames(l$quanti), invert=T),])
  d2 <- as.data.frame(l$quali[grep("^chr", rownames(l$quali), invert=T),])
  colnames(d1) <- c("cor-R2", "p.value")
  colnames(d2) <- c("cor-R2", "p.value")
  d3 <- rbind(d1, d2)
  out <- d3$p.value
  names(out) <- rownames(d3)
  out
}
d4.pvals <- sapply(names(d4.dimdesc), getAssocTraits, simplify=F, USE.NAMES=T)
d4.pvals <- lapply(d4.pvals, function(x) {x[unique(unlist(lapply(d4.pvals, names)))]}) #fix same order 
d4.pvals <- do.call(cbind, d4.pvals)
d4.pvals <- d4.pvals[order(apply(d4.pvals, 1, min), decreasing=T),]

# adjust pvalues
d4.padj <- apply(d4.pvals, 2, p.adjust)
d4.padj <- d4.padj[complete.cases(d4.padj),]

# get the explained variance
var.expl <- signif(d4.pca$eig[1:length(d4.dimdesc),"percentage of variance"],2)


# scatterplot of PC1 vs PC2
d5 <- as.data.frame(d4.pca$ind$coord[,1:2])
colnames(d5) <- gsub("Dim\\.","PC",colnames(d5))
d5$patient <- rownames(d5)
d5$ER <- d4.pca$call$X$ER.FPKM
d5$PMDmeth <- d4.pca$call$X$PMD_meth_wgbs
d5$AberrantCells <- d4.pca$call$X$AberrantCells

p <- ggplot(d5[!is.na(d5$ER),], aes(PC1, PC2)) +
     geom_point(aes(color=PMDmeth, shape=ER, size=2)) +
     scale_color_gradient(low="yellow", high="blue") +
     theme_classic() +
     theme(axis.text=element_text(color="black"))
ggsave(p, file="PCA_WGBS_PC1_vs_PC2.pdf", width=3.5, height=2.5) #### Figure 1D ####




## plot heatmap of PCA scores vs. clinical features (association)
# function for textual formatting of clinical parameters
formatParams <- function(M) {
  rownames(M) <- gsub('ER.FPKM','ER',rownames(M))
  rownames(M) <- gsub('Signatures\\.','',rownames(M))
  rownames(M) <- gsub('PMD_meth_450k','mean_PMD_meth',
      rownames(M), ignore.case=T)
  rownames(M) <- gsub('nbr','#',rownames(M), ignore.case=T)
  rownames(M) <- gsub('([TNM])_stage','\\1-stage ',rownames(M))
  rownames(M) <- gsub('Signatures\\.Signature\\.(\\d{1,}.{0,})','Signature \\1 (#Subs)', 
      rownames(M))
  rownames(M) <- gsub('AberrantCells','aberrant_cells',rownames(M))
  rownames(M) <- gsub('AIMS','AIMS expr type',rownames(M))
  rownames(M) <- gsub('^Sig','sig',rownames(M))
  rownames(M) <- gsub('^Central','central',rownames(M))
  rownames(M) <- gsub('^Sig','sig',rownames(M))
  rownames(M) <- gsub('Rearr','rearr',rownames(M))
  rownames(M) <- gsub('^Lymph','lymph',rownames(M))
  rownames(M) <- gsub('^Sentrix','sentrix',rownames(M))
  rownames(M) <- gsub('^Invasive','invasive',rownames(M))
  rownames(M) <- gsub('Sub','sub',rownames(M))
  rownames(M) <- gsub('New','new',rownames(M))
  rownames(M) <- gsub('Normal','normal',rownames(M))
  rownames(M) <- gsub('Stroma','stroma',rownames(M))
  rownames(M) <- gsub('Adipose','adipose',rownames(M))
  rownames(M) <- gsub('Sample','sample',rownames(M))
  rownames(M) <- gsub('Ploidy','ploidy',rownames(M))
  rownames(M) <- gsub('Necrosis','necrosis',rownames(M))
  rownames(M) <- gsub('_',' ',rownames(M))
  rownames(M) <- gsub('\\.',' ',rownames(M))
  M
}

# provide function for filtering clinical parameters, 
# these provide no relevant/usable information for the heatmap, or are repeats of other features
selectParams <- function(M) {
   exclude.parameters <- c(
   "Project",
   "paper_inclusion",
   "study_ID",
   "donor_age_at_last_follow.up",
   "donor_vital_status",
   "disease_status_last_follow.up",
   "donor_interval_of_last_follow.up_in_DAYS",
   "For_analyses",
   "grade", #also present as "tumour grade"
   "Assay",
   "OSbin",
   "source_of_normal",
   "SignatureCluster",
   "ER", # use ER.FPKM instead
   )
  M <- M[!rownames(M) %in% exclude.parameters,]
  M
}

# function for pval cutoff
cutPval <- function(M, PVALUE=1.5e-01) {
  M <- M[apply(d4.padj, 1, min, na.rm=T) <= PVALUE,]
  M
}

plotHeatmap <- function(PVALS, NAME, BIN=FALSE, HEIGHT=10, WIDTH=5) { 
  #PVALS: matrix with pvals ; BIN: make pval bins ; NAME: pdf name
  br <- c(0, -3, -6, -10, -20, -30, -50, -110)
  m1 <- PVALS
  if (BIN) {
    m1 <- apply(log10(m1), 2, function(x) {out <- cut(x, breaks=br)
      names(out) <- names(x) ; out})
  }
  m1.m <- melt(m1)
  colnames(m1.m) <- c("param","PC","padj")
  m1.m$PC <- gsub("Dim\\.","", m1.m$PC)
  m1.m$PC <- factor(m1.m$PC, levels=c(1:length(d4.dimdesc)))
  m1.m$PC.var <- paste0("(", var.expl, ")  ", levels(m1.m$PC))[as.numeric(m1.m$PC)]
  pc.ord <- unique(m1.m$PC.var)[order(as.numeric(sapply(strsplit(unique(m1.m$PC.var), " "), 
      function(x) {x[length(x)]})))]
  m1.m$PC <- factor(m1.m$PC.var, levels=pc.ord)
  if (BIN) {
    ord <- order(as.numeric(sapply(strsplit(gsub('\\(','',levels(m1.m$padj)), ","),
        function(x) {x[1]})))
    m1.m$padj <- factor(m1.m$padj, levels=levels(m1.m$padj)[rev(ord)])
  }
  p <- ggplot(m1.m, aes(PC, param))
  if (BIN) {
    p <- p + geom_tile(aes(fill=padj, colour=padj))
    p <- p + scale_fill_brewer(palette="YlGnBu")
    p <- p + scale_colour_brewer(palette="YlGnBu")
  } else {
    p <- p + geom_tile(aes(fill=-log10(padj), colour=-log10(padj)))
    p <- p + scale_fill_gradientn(colours=brewer.pal(5, "YlGnBu"))
    p <- p + scale_colour_gradientn(colours=brewer.pal(5, "YlGnBu"))
  }
  p <- p + theme_classic()
  p <- p + theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))
  p <- p + xlab("PC (%variance explained)")
  p <- p + ylab("")
  ggsave(p, file=NAME, width=WIDTH, height=HEIGHT)
}
plotHeatmap(formatParams(d4.padj), "PCA_factominer_vs_clinj_padj.pdf", WIDTH=5.75, HEIGHT=18)
plotHeatmap(formatParams((cutPval(d4.padj))), "PCA_factominer_vs_clinj_padj_cutPval.pdf", WIDTH=5.75, HEIGHT=2.75)



