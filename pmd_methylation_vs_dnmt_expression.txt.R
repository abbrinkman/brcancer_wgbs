library(GenomicRanges)

# get RNAseq data (FPKM) from Smid et al. Nature Communications 7, 12910 (2016).
# GRanges opject with genes as rows and patients as columns
# Example:
#       seqnames           ranges strand |    UNIQID         Ensembl
#          <Rle>        <IRanges>  <Rle> | <integer>        <factor>
#   [1]    chr16 [ 61554,  64090]      + |     38384 ENSG00000233614
#   [2]    chr16 [103009, 107669]      + |      6475 ENSG00000161981
#   [3]    chr16 [127005, 135852]      + |     29072 ENSG00000103152
#   [4]    chr16 [202685, 204502]      + |     31303 ENSG00000130656
#   [5]    chr16 [203890, 216767]      + |     17157 ENSG00000206177
#                     Source     Name  PD10010a  PD10011a  PD10014a
#                   <factor> <factor> <numeric> <numeric> <numeric>
#   [1] processed_transcript DDX11L10   -0.3454   -0.3780   -0.1590
#   [2]       protein_coding  SNRNP25    2.4905    2.8318    3.2312
#   [3]       protein_coding      MPG    0.7188    2.9728    2.7242
#   [4]       protein_coding      HBZ      <NA>      <NA>      <NA>
#   [5]       protein_coding      HBM      <NA>      <NA>      <NA>
expr <- get(load("fpkm.gr.RData"))
source("~/tools/BASIS_common_functions.txt.R")


# genes of interest for DNAme
genes <- c("DNMT3A","DNMT3L","DNMT1","DNMT3B", "TET1", "TET2", "TET3",
           "UHRF1BP1", "UHRF2", "UHRF1", "UHRF1BP1L", "MBD1", "MBD2", "MBD3", "MBD4", "MECP2",
           "IDH1", "IDH2")

d1 <- as.data.frame(mcols(expr[expr$Name %in% genes]))
rownames(d1) <- paste(d1$Name, d1$UNIQID, d1$Ensembl, sep="_")
d1 <- d1[, grep("^PD", colnames(d1))]
d1 <- t(d1)
d2 <- data.frame(
  d1, "PMD_meth_wgbs"=clinj$PMD_meth_wgbs[match(rownames(d1), paste0(clinj$sample_name, "a"))])

d2 <- d2[!is.na(d2$PMD_meth_wgbs),]

# plot mean PMD methylation as a function of gene expression
plotFunc <- function(GENE, PMD) {
    na <- is.na(d2[,GENE]) | is.na(d2[,PMD])
    r <- round(cor(d2[!na,GENE], d2[!na,PMD]),2)
    mdl <- lm(d2[!na,GENE] ~ d2[!na,PMD])
    p <- signif(summary(mdl)$coefficients[,"Pr(>|t|)"][2], 2)
    r2 <- round(summary(mdl)$r.squared, 2)
    plot(d2[,GENE], d2[,PMD], xlab="FPKM (log2)", ylab=PMD, 
        cex=0.5, cex.main=0.8,col="blue")
    legend("topleft", bty="n", legend=c(paste("R =", r), paste("R2 = ", r2), paste("p = ", p)), cex=0.8)
    legend(par("usr")[1]-0.23*(par("usr")[2]-par("usr")[1]), par("usr")[4]+0.23*(par("usr")[4]-par("usr")[3]), 
        bty="n", legend=gsub("_.+$","",colnames(d2)[GENE]), xpd=T, xjust=0, cex=1.1)
    
}

pdf("PMDmeth_vs_DNMTs.pdf") #### Supplemental Figure 1D ####
par(mfrow=c(3,4), mar=c(5,4,5,2))
sapply(1:length(genes), function(x) {plotFunc(x, "PMD_meth_wgbs")})
dev.off()

