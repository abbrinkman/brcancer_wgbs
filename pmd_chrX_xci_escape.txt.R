library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(reshape2)
library(ggplot2)
library(openxlsx)

# get the PMD data, as a list per patient
pmdmeth.basis <- get(load("~/BiSeq_BASIS/PMDs_all_samples/pmdmeth.basis.RData"))
pmdmeth.basis <- pmdmeth.basis[!names(pmdmeth.basis) %in% c("HMEC","MCF7")]

# get the ploidy data
ploidy <- read.xlsx("~/Infinium/BASIS_data_all_TN_QC_analysis/suppl_nik-zainal_landscape_nature17676-s3/Supplementary Table 5.Ploidy.AberrantCellFraction.090402015.v1.xlsx")
ploidy$Sample <- paste0(ploidy$Sample, "a")
ploidy <- split(ploidy, ploidy$Sample)
 
# get the CNV data as a list per patient
cnv <- read.table("~/Infinium/BASIS_data_all_TN_QC_analysis/suppl_nik-zainal_landscape_nature17676-s3/Supplementary Table 4.Copy.Number.Segments.txt", header=T)
cnv$chr <- paste0("chr", cnv$chr)
cnv$chr <- gsub("chr23","chrX", cnv$chr)
cnv <- makeGRangesFromDataFrame(cnv, keep.extra.columns=T)
mcols(cnv)$pd <- gsub("a2","a",gsub("a_2","a", mcols(cnv)$SampleID)) # fix names only for the "a" samples (=tumors)
cnv <- cnv[cnv$pd %in% names(pmdmeth.basis)]
sum(names(pmdmeth.basis) %in% mcols(cnv)$pd) # should be 25
cnv <- split(cnv, cnv$pd)

# make the CNV segments 'relative': for expression, you will not see the difference between a 2N or a 4N genome,
# but if SEGMENTS in the genome deviate from this N, it will matter
# for this, take the ploidy into account, and use Johan Staaf's guidelines:
# gains=ploidy+0.6, amplications=total copy number >2*ploidy, deletion=total copy number < tumor.ploidy-0.6
relativeCN <- function(PD) {
  gr1 <- cnv[[PD]]
  pl <- ploidy[[PD]]$Ploidy
  gain <- which(gr1$seg.mean > (pl+0.6))
  loss <- which(gr1$seg.mean < (pl-0.6))
  ampl <- which(gr1$seg.mean > (2*pl))
  #del <- which((gr1$seg.mean - (pl-0.6)) < 0)
  del <- which(gr1$seg.mean==0)
  gr1$cn.rel <- rep("unchanged", length(gr1))
  gr1$cn.rel[gain] <- "gain"
  gr1$cn.rel[loss] <- "loss"
  gr1$cn.rel[ampl] <- "amplified"
  gr1$cn.rel[del] <- "deleted"
  gr1
}
cnv <- sapply(names(cnv), relativeCN, USE.NAMES=T, simplify=F)

# get the escape/XCI data from Balaton et al.
# checked randomly 5 genes: this is hg19 (not clear from the paper)
esc <- read.xlsx("~/BiSeq_BASIS/PMDs_all_samples/PMDs_chrX/13293_2015_53_MOESM1_ESM.xlsx") 


# get the expression
fpkm <- get(load("~/Infinium/BASIS_data_all_TN_QC_analysis/fpkm.gr.RData"))
esc1 <- esc[!esc$Gene.Name %in% esc$Gene.Name[duplicated(esc$Gene.Name)], ]
fpkm1 <- fpkm[fpkm$Name %in% intersect(esc1$Gene.Name, fpkm$Name)]
esc1 <- esc1[esc1$Gene.Name %in% intersect(esc1$Gene.Name, fpkm$Name),]
fpkm1 <- fpkm1[match(esc1$Gene.Name, fpkm1$Name)]
fpkm1 <- fpkm1[,c(1:4, which(colnames(mcols(fpkm1)) %in% names(pmdmeth.basis)))]

# compile all data into a single list; one data.frame per patient
pds <- intersect(intersect(names(pmdmeth.basis), names(cnv)), colnames(mcols(fpkm1))) # 22 patients

getFPKM_PMD_CNV <- function(PD) {
  gr <- fpkm1[,c(4, which(colnames(mcols(fpkm1))==PD))]
  colnames(mcols(gr)) <- c("gene", "FPKM")
  mcols(gr)$XCI.status <- esc1$Balaton.consensus.calls[match(gr$gene, esc1$Gene.Name)]
  gr$PMD <- ifelse(overlapsAny(gr, pmdmeth.basis[[PD]], type="within"), "in", "out")
  #gr$CN <- rep(2, length(gr)) 
  gr$CN <- rep(NA, length(gr)) 
  gr$CN.rel <- rep(NA, length(gr)) 
  ov <- findOverlaps(gr, cnv[[PD]], type="within")
  gr$CN[queryHits(ov)] <- cnv[[PD]]$seg.mean[subjectHits(ov)]
  gr$CN.rel[queryHits(ov)] <- cnv[[PD]]$cn.rel[subjectHits(ov)]
  as.data.frame(mcols(gr))
}

l1 <- sapply(pds, function(x) {message(x) ; getFPKM_PMD_CNV(x)}, USE.NAMES=T, simplify=F)
d1 <- do.call(rbind, l1)
d1$patient <- gsub("\\..+$","", rownames(d1))
d1 <- d1[!d1$XCI.status %in% c("Discordant", "No call"),]
d1$XCI.status <- factor(d1$XCI.status, levels=c("S","Mostly S","E","Mostly E","VE","Mostly VE", "PAR"))
d1$XCI.status.1 <- d1$XCI.status
d1$XCI.status <- as.factor(as.character(gsub("Mostly ", "", d1$XCI.status))) # aggregate "S" / "Mostly S", etc...
d1$XCI.status <- gsub("^S$","(mostly) S", d1$XCI.status)
d1$XCI.status <- gsub("^E$","(mostly) E", d1$XCI.status)
d1$XCI.status <- gsub("^VE$","(mostly) VE", d1$XCI.status)
d1$XCI.status <- factor(d1$XCI.status, levels=c("(mostly) S","(mostly) E","(mostly) VE","PAR"))
d1$CN.bin <- as.character(d1$CN)
d1$CN.bin[d1$CN > 4] <- "5+"
d1$CN.bin <- as.factor(d1$CN.bin)
d1 <- d1[!is.na(d1$CN.bin),]
rownames(d1) <- NULL

p1 <- ggplot(d1, aes(PMD, FPKM)) +
  geom_boxplot(aes(fill=PMD), outlier.size=0.05) +
  scale_fill_manual(values=c("in"="red", "out"="white")) +
  facet_grid(~XCI.status.1) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(),
      axis.text=element_text(colour="black"), axis.text.x=element_text(angle=45, vjust=1, hjust=1)) +
  ylab("FPKM (log2)")
ggsave(p1, file="boxplots_FPKM_vs_PMD_XCIgenes.pdf", width=6, height=2.5) ##### Figure 3L #####

p2 <- ggplot(d1, aes(CN.rel, FPKM)) +
  geom_boxplot(aes(fill=PMD), outlier.size=0.05) +
  facet_grid(~XCI.status) +
  scale_fill_manual(values=c("in"="red", "out"="white")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(),
      axis.text=element_text(colour="black"), axis.text.x=element_text(angle=45, vjust=1, hjust=1)) +
  xlab("rel. copy number") +
  ylab("FPKM (log2)")
ggsave(p2, file="boxplots_FPKM_vs_CNrelative_vs_PMD_XCIgenes.pdf", width=8, height=2.5) ##### Supplemental Figure 6E #####


# clearly expression of Escape genes in PMDs are more severely affected by PMD formation
# significance of this:
tmp <- d1[d1$PMD=="in",]
tmp$esc <- ifelse(tmp$XCI.status.1=="E", "esc", "not.esc")
t.test(tmp$FPKM ~ tmp$esc)$p.value
  # [1] [1] 2.551283e-40


