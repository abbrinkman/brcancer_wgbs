library(GenomicRanges)
library(ggplot2)
library(reshape2)
library(rtracklayer)
library(parallel)
library(data.table)

# get the genome positions of all CGIs and CGI shores
session <- browserSession()
genome(session) <- "hg19"
query <- ucscTableQuery(session, "CpG Islands")
cgi <- GRanges(track(query))
cgi$name <- NULL
shore.size <- 2000
cgi.shores <- cgi
start(cgi.shores) <- start(cgi.shores)-shore.size
end(cgi.shores) <- end(cgi.shores)+shore.size

# get the genome positions of all promoters
query <- ucscTableQuery(session, "GENCODE Genes V19")
genes <- GRanges(track(query))
prom.size <- 2000
prom <- promoters(genes, upstream=prom.size/2, downstream=prom.size/2)

# GRanges containing the CGIs plus shores, and promoters
cgi.shores.prom <- c(cgi.shores[,0], prom[,0])

# function for calculation of mean methylation (weighted mean) in regions
weightedMeth <- function(REGIONS, WGBS, COLNAME, PMDS=TRUE) {
  # REGIONS: GRanges object
  # WGBS: GRanges object with 2 columns (xx.T , xx.M) as patients
  # COLNAME: name the the mcols column in output GRanges object
  # PMDS: setting that determines whether CGIs/promoters should be taken out
  message(COLNAME)

  if (PMDS) { # delete CpGs in promoters, CGIs+shores
    wgbs1 <- WGBS[-queryHits(findOverlaps(WGBS, cgi.shores.prom))]
  } else {
    wgbs1 <- WGBS
  }

  # get unique sample names
  ov <- findOverlaps(REGIONS, wgbs1)
  cols.t <- which(grepl("\\T$", colnames(mcols(wgbs1))) | grepl("_cov$", colnames(mcols(wgbs1))))
  cols.m <- which(grepl("\\M$", colnames(mcols(wgbs1))) | grepl("_meth$", colnames(mcols(wgbs1))))
  if (identical(range(mcols(WGBS)[[cols.m]], na.rm=T), c(0,1))) {
    mcols(wgbs1)[[cols.m]]  <- round(mcols(wgbs1)[[cols.m]] * mcols(wgbs1)[[cols.t]], 0) 
  }

  dt1 <- data.table("pmds"=queryHits(ov),
      "T"=mcols(wgbs1[subjectHits(ov)])[,cols.t],
      "M"=mcols(wgbs1[subjectHits(ov)])[,cols.m])
  setkey(dt1, pmds)

  # calculate weighted mean
  dtT <- dt1[,sum(T, na.rm=T),by=pmds]
  dtM <- dt1[,sum(M, na.rm=T),by=pmds]
  gr1 <- REGIONS[dtT$pmds]
  gr1@elementMetadata@listData[[COLNAME]] <- round(dtM$V1/dtT$V1, 2)

  # include also all tiles that had no overlap with CpGs
  gr2 <- REGIONS[!c(1:length(REGIONS)) %in% queryHits(ov)]
  gr2@elementMetadata@listData[[COLNAME]] <- as.numeric(rep(NA, length(gr2)))
  gr3 <- c(gr1, gr2)
  gr3 <- sort(gr3)

  # just to make absolutely sure that the order of each gr3 is the same as REGIONS
  ov1 <- findOverlaps(REGIONS, gr3)
  gr3 <- gr3[subjectHits(ov1)]

  # return data
  gr3
}

# basis
basis.files <- list.files(pattern="PMDs_.+_noCent.bed$", 
    path="~/BiSeq_BASIS/PMDs_all_samples/", full.names=T)
basis.files <- basis.files[!grepl("PD9590", basis.files)]
basis.files <- as.list(basis.files)
names(basis.files) <- sapply(strsplit(basename(unlist(basis.files)), "_"), function(x) {x[2]})
pmds.basis <- lapply(basis.files, function(x) {read.table(x, col.names=c("chr","start","end"))})
pmds.basis <- lapply(pmds.basis, makeGRangesFromDataFrame)
wgbs.basis <- get(load("~/BiSeq_BASIS/methcounts/meth.gr.RData"))
samples.basis <- as.list(names(pmds.basis))
names(samples.basis) <- unlist(samples.basis)

pmdmeth.basis <- mclapply(samples.basis, function(x) {
    weightedMeth(pmds.basis[[x]], wgbs.basis[,grep(x, colnames(mcols(wgbs.basis)))], x) }, mc.cores=2)
save(pmdmeth.basis, file="pmdmeth.basis.RData")

cgimeth.basis <- mclapply(samples.basis, function(x) {
    weightedMeth(cgi, wgbs.basis[,grep(x, colnames(mcols(wgbs.basis)))], x, PMDS=FALSE) }, mc.cores=2)
save(cgimeth.basis, file="cgimeth.basis.RData")

# lymphoma
pmds.lymph <- get(load("/home/arjen/BiSeq_BASIS/PMDs_other_celltypes/lymphoma/lymph.pmds.RData"))
ispmd.lymph <- read.table(
    "/home/arjen/BiSeq_BASIS/PMDs_other_celltypes/lymphoma/lymphoma_samples_with_PMDs.txt")
wgbs.lymph <- get(load("/home/arjen/BiSeq_BASIS/PMDs_other_celltypes/lymphoma/lymph.wgbs.RData"))
samples.lymph <- as.list(as.character(ispmd.lymph$V1[ispmd.lymph$V3=="OK" & ispmd.lymph$V2=="PMDs"]))
names(samples.lymph) <- unlist(samples.lymph)

pmdmeth.lymph <- mclapply(samples.lymph, function(x) {
    weightedMeth(pmds.lymph[[x]], wgbs.lymph[,grep(x, colnames(mcols(wgbs.lymph)))], x) }, mc.cores=2)
save(pmdmeth.lymph, file="pmdmeth.lymph.RData")

cgimeth.lymph <- mclapply(samples.lymph, function(x) {
    weightedMeth(cgi, wgbs.lymph[,grep(x, colnames(mcols(wgbs.lymph)))], x, PMDS=FALSE) }, mc.cores=2)
save(cgimeth.lymph, file="cgimeth.lymph.RData")

# tcga
tcga.files <- as.list(list.files(path="~/BiSeq_BASIS/PMDs_other_celltypes/tcga", 
    pattern="gr.RData", full.names=T))
names(tcga.files) <- gsub("\\.gr\\.RData$","",basename(unlist(tcga.files)))
pmds.tcga <- get(load("/home/arjen/BiSeq_BASIS/PMDs_other_celltypes/tcga/tcga.pmds.RData"))
samples.tcga <- names(pmds.tcga)
samples.tcga <- as.list(samples.tcga[!grepl("BLCA_NIC1254A93", samples.tcga) & !grepl("normal", 
    samples.tcga)])
names(samples.tcga) <- unlist(samples.tcga)

pmdmeth.tcga <- mclapply(samples.tcga, function(x) {
    weightedMeth(pmds.tcga[[x]], get(load(tcga.files[[x]])), x) }, mc.cores=2)
save(pmdmeth.tcga, file="pmdmeth.tcga.RData")

cgimeth.tcga <- mclapply(samples.tcga, function(x) {
    weightedMeth(cgi, get(load(tcga.files[[x]])), x, PMDS=FALSE) }, mc.cores=2)
save(cgimeth.tcga, file="cgimeth.tcga.RData")

# roadmap
wgbs.rdmp <- get(load("/home/arjen/BiSeq_BASIS/enhancers/roadmap.wgbs.gr.RData"))
pmds.rdmp <- get(load("/home/arjen/BiSeq_BASIS/enhancers/roadmap/LMRs_UMRs/roadmap.pmds.RData")) 
ispmd.rdmp <- read.table("/home/arjen/BiSeq_BASIS/enhancers/roadmap/LMRs_UMRs/samples_with_PMDs.txt")
samples.rdmp <- as.list(as.character(ispmd.rdmp$V1[ispmd.rdmp$V2=="PMDs"]))
names(samples.rdmp) <- unlist(samples.rdmp)

pmdmeth.rdmp <- mclapply(samples.rdmp, function(x) {
    weightedMeth(pmds.rdmp[[x]], wgbs.rdmp[,grep(x, colnames(mcols(wgbs.rdmp)))], x) }, mc.cores=2)
save(pmdmeth.rdmp, file="pmdmeth.roadmap.RData")

cgimeth.rdmp <- mclapply(samples.rdmp, function(x) {
    weightedMeth(cgi, wgbs.rdmp[,grep(x, colnames(mcols(wgbs.rdmp)))], x, PMDS=FALSE) }, mc.cores=1)
save(cgimeth.rdmp, file="cgimeth.roadmap.RData")

# calculate the mean for each sample
pmdmeth <- c(
  sapply(pmdmeth.basis, function(x) {mean(x@elementMetadata@listData[[1]], na.rm=T)}),
  sapply(pmdmeth.lymph, function(x) {mean(x@elementMetadata@listData[[1]], na.rm=T)}),
  sapply(pmdmeth.tcga, function(x) {mean(x@elementMetadata@listData[[1]], na.rm=T)}),
  sapply(pmdmeth.rdmp, function(x) {mean(x@elementMetadata@listData[[1]], na.rm=T)})
)
pmdmeth <- data.frame("sample"=names(pmdmeth), "pmd.meth"=pmdmeth)

cgimeth <- c(
  sapply(cgimeth.basis, function(x) {mean(x@elementMetadata@listData[[1]], na.rm=T)}),
  sapply(cgimeth.lymph, function(x) {mean(x@elementMetadata@listData[[1]], na.rm=T)}),
  sapply(cgimeth.tcga, function(x) {mean(x@elementMetadata@listData[[1]], na.rm=T)}),
  sapply(cgimeth.rdmp, function(x) {mean(x@elementMetadata@listData[[1]], na.rm=T)})
)
cgimeth <- data.frame("sample"=names(cgimeth), "cgi.meth"=cgimeth)

pmd.cgi.meth <- merge(pmdmeth, cgimeth, by.x="sample", by.y="sample")



# clean up names of tissues
pmd.cgi.meth$tissue <- pmd.cgi.meth$sample
pmd.cgi.meth$tissue <- gsub("NIC.+_","", pmd.cgi.meth$tissue)
pmd.cgi.meth$tissue <- gsub("^GCB_\\d+","GCB", pmd.cgi.meth$tissue)
pmd.cgi.meth$tissue <- gsub("^FL_\\d+","FL_tumor", pmd.cgi.meth$tissue)
pmd.cgi.meth$tissue <- gsub("^BL_\\d+","BL_tumor", pmd.cgi.meth$tissue)
pmd.cgi.meth$tissue <- gsub("^PD\\d+a","BRCA_tumor", pmd.cgi.meth$tissue)

# get the tissues for roadmap data
tissues <- read.csv("/home/arjen/BiSeq_BASIS/enhancers/roadmap/jul2013.roadmapData.qc.csv")
pmd.cgi.meth$tissue[grepl("^E\\d\\d\\d$", pmd.cgi.meth$tissue)] <- paste(tissues$ANATOMY, 
    tissues$TYPE, sep="_")[match(pmd.cgi.meth$tissue[grepl("^E\\d\\d\\d$", pmd.cgi.meth$tissue)], 
    tissues$Epigenome.ID..EID.)]
pmd.cgi.meth <- pmd.cgi.meth[!grepl("ESC_DERIVED", pmd.cgi.meth$tissue),]
pmd.cgi.meth$name <- make.unique(pmd.cgi.meth$tissue)
pmd.cgi.meth$name <- factor(pmd.cgi.meth$name, levels=as.character(pmd.cgi.meth$name)[order(pmd.cgi.meth$pmd.meth)] )

# clean up names of tissues
pmd.cgi.meth$tissue <- gsub(".+_PrimaryCulture","Primary_culture", pmd.cgi.meth$tissue)
pmd.cgi.meth$tissue <- gsub(".+_PrimaryTissue","Primary_tissue", pmd.cgi.meth$tissue)
pmd.cgi.meth$tissue <- gsub("GCB","Primary_tissue", pmd.cgi.meth$tissue)
pmd.cgi.meth$tissue <- gsub("HMEC","Primary_culture", pmd.cgi.meth$tissue)
pmd.cgi.meth$tissue <- gsub(".+_CellLine","Celline", pmd.cgi.meth$tissue)
pmd.cgi.meth$tissue <- gsub("MCF7","Celline", pmd.cgi.meth$tissue)
pmd.cgi.meth$tissue <- gsub(".+_tumor","Tumor", pmd.cgi.meth$tissue)
pmd.cgi.meth$tissue <- factor(pmd.cgi.meth$tissue, levels=c("Primary_tissue","Primary_culture",
    "Celline","Tumor"))
pmd.cgi.meth$sample <- NULL
pmd.cgi.meth.m <- melt(pmd.cgi.meth, variable.name="element", value.name="meth") 
pmd.cgi.meth.m$element <- toupper(gsub("\\.meth","", pmd.cgi.meth.m$element))
pmd.cgi.meth.m$name <- factor(pmd.cgi.meth.m$name, 
 levels=pmd.cgi.meth.m$name[pmd.cgi.meth.m$element=="PMD"][order(pmd.cgi.meth.m$meth[pmd.cgi.meth.m$element=="PMD"])])
pmd.cgi.meth.m$name1 <- pmd.cgi.meth.m$name
levels(pmd.cgi.meth.m$name1) <- rev(levels(pmd.cgi.meth.m$name1))

pmd.cgi.meth.m$name2 <- gsub("^(.+)_Primary.+$","\\L\\1",gsub("LUNG_CellLine","lung", 
    gsub("GI_","", gsub("_tumor","", gsub(".\\d+","",pmd.cgi.meth.m$name)))), perl=T)
pmd.cgi.meth.m$BASIS <- ifelse(pmd.cgi.meth.m$name2=="BRCA", "\u2190","")

# boxplot of PMDmeth, tumor vs normal tissues -> #### Supplemental Figure 7A ####
pmd.cgi.meth$is.tumor <- ifelse(pmd.cgi.meth$tissue %in% c("Tumor","Celline"), "Tumor", "Normal")
pdf("boxplot_PMDmeth_tumors_vs_normal_tissues.pdf")
par(mfrow=c(3,3))
boxplot(pmd.meth ~ is.tumor, data=pmd.cgi.meth, ylab="mean PMD meth")
legend("bottomleft", bty="n", legend=paste0("p = ", signif(t.test(pmd.meth ~ is.tumor, data=pmd.cgi.meth)$p.value, 3)))
dev.off()


p1 <- ggplot(pmd.cgi.meth.m[pmd.cgi.meth.m$element=="PMD",], aes(name, meth)) + 
geom_point(aes(color=tissue), size=2) + 
scale_color_manual(values=c("Tumor"="red","Primary_tissue"="darkgreen","Primary_culture"="green","Celline"="orange")) +
geom_text(aes(label=name2, color=tissue), size=2.7, hjust=1, vjust=0.5, nudge_y=-0.012, angle=90, fontface="bold") + 
geom_text(aes(label=BASIS), hjust=0, vjust=0, angle=90, nudge_y=0.012, nudge_x=0.5, fontface="bold") + 
theme_classic() + 
theme(axis.text=element_text(color="black")) +
xlab("mean methylation in PMDs")
pdf("dotplot_PMD_methylation_BASIS_TCGA_lymph_roadmap.pdf", width=12, height=3) #### Figure 4B ####
p1
dev.off()


# in how many samples can we detect PMDs?
#  basis 
#    tumor: 30/30
#  lymph
#    normal (GCB): 4/4
#    tumor (BL/FL): 20/21
#  roadmap
#    normal: 21/37
#  tcga
#    tumor: 38/38
#  blueprint
#    ?
#  schulz
#    normal: 0/36

pmd.pres.basis <- data.frame(
  "pmd"=rep("detected", 30), "tumor"=rep("tumor", 30)
)
pmd.pres.lymph <- data.frame(
  "pmd"=c(rep("detected", 24), rep("not_detected", 1)), "tumor"=c(rep("normal", 4),rep("tumor", 21))
)
pmd.pres.roadmap <- data.frame(
  "pmd"=c(rep("detected", 21), rep("not_detected", 16)), "tumor"=rep("normal", 37)
)
pmd.pres.tcga <- data.frame(
  "pmd"=rep("detected", 38), "tumor"=rep("tumor", 38)
)
pmd.pres.schultz <- data.frame(
  "pmd"=rep("not_detected", 36), "tumor"=rep("normal", 36)
)

pmd.pres <- rbind(pmd.pres.basis, pmd.pres.lymph, pmd.pres.roadmap, pmd.pres.tcga, pmd.pres.schultz)

pdf("PMD_presence_normal_tumor.pdf") #### Figure 4A ####
par(mfrow=c(2,3))
barplot(table(pmd.pres), col=c("red","white"), ylab="#samples")
barplot(apply(table(pmd.pres), 2, function(x) {x/sum(x)}), col=c("red","white"), ylab="fraction of samples", 
    ylim=c(0,1))
dev.off()

