library(GenomicRanges)
library(parallel)
library(ggplot2)
library(reshape2)
library(openxlsx)
library(ggrepel)

# use the TSG genes from Census (COSMIC), which also specifies the type of cancer 
tsg <- read.table("~/BiSeq_BASIS/PMDs_all_samples/TSG_overlap/Census_allMon Oct 23 13-54-39 2017.tsv", 
    sep="\t", header=T)

# also do the analysis for the 93 drivergenes published by Nik-Zainal
driv <- read.xlsx(
"suppl_nik-zainal_landscape_nature17676-s3/Supplementary Table 14.Driver.Events.By.Mutation.Type.01052015.v2.xlsx", sheet=4)

driv.genes <- unlist(strsplit(gsub("\\)","",gsub("Chr8:\\(","", unique(driv$Gene))), "/"))


# get the PMDs
pmd.files <- list.files(pattern="PMDs_PD[0-9]+a_noCent.bed$",
    path="~/BiSeq_BASIS/PMDs_all_samples/", full.names=T)
pmd.files <- pmd.files[!grepl("PD9590", pmd.files)]
pmd.files <- as.list(pmd.files)
names(pmd.files) <- sapply(strsplit(basename(unlist(pmd.files)), "_"), function(x) {x[2]})

pmds <- lapply(pmd.files, function(x) {read.table(x, col.names=c("chr","start","end"))})
pmds <- lapply(pmds, makeGRangesFromDataFrame)

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
fpkm <- get(load("fpkm.gr.RData"))

# clean up the RNA-seq annotation
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

# determine for each gene the PMD frequency, and PMD frequency bin
inPMD <- function(PD) {
  ov <- findOverlaps(fpkm, pmds[[PD]], type="within", select="all")
  in.pmd <- ifelse(1:length(fpkm) %in% queryHits(ov), 1, 0)
  in.pmd
}
pmd.freq <- sapply(names(pmds), inPMD)
pmd.freq <- apply(pmd.freq, 1, sum)
fpkm$pmd.freq <- pmd.freq
fpkm$pmd.freq.bin <- cut(fpkm$pmd.freq, breaks=c(-0.01, seq(0,30, by=3)))

# for each gene, mark whether it is a TSG, breast-specific TSG, or driver (Nik-Zainal)
census.tsg <- as.character(tsg$Gene.Symbol[grepl("TSG", tsg$Role.in.Cancer)])
census.tsg.breast <- as.character(tsg$Gene.Symbol[grepl("TSG", tsg$Role.in.Cancer) & 
  (grepl("breast", tsg$Tumour.Types.Somatic., ignore.case=T) | grepl("breast", tsg$Tumour.Types.Germline., 
   ignore.case=T))])

fpkm$tsg <- ifelse(fpkm$Name %in% census.tsg, "TSG", "noTSG") 
fpkm$tsg.breast <- ifelse(fpkm$Name %in% census.tsg.breast, "TSG", "noTSG") 
fpkm$driv <- ifelse(fpkm$Name %in% driv.genes, "Driv", "noDriv")
fpkm$other <- rep("other", length(fpkm))
fpkm$other[fpkm$tsg=="TSG" | fpkm$tsg.breast=="TSG" | fpkm$driv=="Driv"] <- NA
fpkm <- fpkm[!duplicated(fpkm$Name)]


# per PMD frequency bin, count the number of TSGs (of all classes), in absolute counts or as a fraction
tsg.counts <- table(fpkm$tsg, fpkm$pmd.freq.bin)["TSG",]
tsg.breast.counts <- table(fpkm$tsg.breast, fpkm$pmd.freq.bin)["TSG",]
driv.counts <- table(fpkm$driv, fpkm$pmd.freq.bin)["Driv",]
other.counts <- table(fpkm$other, fpkm$pmd.freq.bin)["other",]

d.counts <- t(data.frame(
    "other"=other.counts,
    "TSG_all_cancers"=tsg.counts, 
    "TSG_breast_cancer"=tsg.breast.counts, 
    "driver_mutated_genes_breast_cancer"=driv.counts))

d.frac <- t(apply(d.counts, 1, function(x) {x/sum(x)}))


# plot these TSG counts
pdf("histograms_TSG_driver_fractions_vs_PMDfreq_Census.pdf", height=6, width=6) ##### Figure 3K #####
par(mfrow=c(3,2), mar=c(5,4,3,0))
barplot(d.frac, beside=T, col=c("darkgreen","red","orange", "brown"), ylim=c(0,1),
    names.arg=rep("", ncol(d.frac)), ylab="gene fraction", xlab="PMD frequency", border=NA, 
    las=1) -> bp
text(apply(bp, 2, mean), par("usr")[3]-(0.05*par("usr")[4]), srt = 45, adj = 1, xpd = T, labels =
    gsub("0.99-","", unlist(lapply(strsplit(gsub("\\(","",gsub("]","", colnames(d.counts))), ","), 
        function(x) {x <- as.numeric(x) ; x[1] <- x[1]+1 ; paste(x, collapse="-")}))))
legend("topright", bty="n", legend=rownames(d.frac),
    fill=c("darkgreen","red","orange", "brown"), border=NA)

# and for each category separately, to distinguish y-scales (inset)
par(mfrow=c(3,2), mar=c(5,4,3,0))

# other
barplot(d.counts, beside=T, col=c("darkgreen","red","orange", "brown"),
    names.arg=rep("", ncol(d.frac)), ylab="gene fraction", xlab="PMD frequency", 
    border=NA, las=1, ylim=c(0,max(d.counts["other",]))) -> bp
text(apply(bp, 2, mean), par("usr")[3]-(0.05*par("usr")[4]), srt = 45, adj = 1, xpd = T, labels =
    gsub("0.99-","", unlist(lapply(strsplit(gsub("\\(","",gsub("]","", colnames(d.counts))), ","),
        function(x) {x <- as.numeric(x) ; x[1] <- x[1]+1 ; paste(x, collapse="-")}))))
legend("topright", bty="n", legend=rownames(d.counts),
    fill=c("darkgreen","red","orange", "brown"), border=NA)

# tsg
barplot(d.counts, beside=T, col=c("darkgreen","red","orange", "brown"),
    names.arg=rep("", ncol(d.frac)), ylab="gene fraction", xlab="PMD frequency",     border=NA, las=1, ylim=c(0,max(d.counts["TSG_all_cancers",]))) -> bp
text(apply(bp, 2, mean), par("usr")[3]-(0.05*par("usr")[4]), srt = 45, adj = 1, xpd = T, labels =
    gsub("0.99-","", unlist(lapply(strsplit(gsub("\\(","",gsub("]","", colnames(d.counts))), ","),
        function(x) {x <- as.numeric(x) ; x[1] <- x[1]+1 ; paste(x, collapse="-")}))))
legend("topright", bty="n", legend=rownames(d.counts),
    fill=c("darkgreen","red","orange", "brown"), border=NA)

# tsg breast
barplot(d.counts, beside=T, col=c("darkgreen","red","orange", "brown"),
    names.arg=rep("", ncol(d.frac)), ylab="gene fraction", xlab="PMD frequency",     border=NA, las=1, ylim=c(0,max(d.counts["TSG_breast_cancer",]))) -> bp
text(apply(bp, 2, mean), par("usr")[3]-(0.05*par("usr")[4]), srt = 45, adj = 1, xpd = T, labels =
    gsub("0.99-","", unlist(lapply(strsplit(gsub("\\(","",gsub("]","", colnames(d.counts))), ","),
        function(x) {x <- as.numeric(x) ; x[1] <- x[1]+1 ; paste(x, collapse="-")}))))
legend("topright", bty="n", legend=rownames(d.counts),
    fill=c("darkgreen","red","orange", "brown"), border=NA)

# Nik-Zainal breast drivers
barplot(d.counts, beside=T, col=c("darkgreen","red","orange", "brown"),
    names.arg=rep("", ncol(d.frac)), ylab="gene fraction", xlab="PMD frequency",     border=NA, las=1, ylim=c(0,max(d.counts["driver_mutated_genes_breast_cancer",]))) -> bp
text(apply(bp, 2, mean), par("usr")[3]-(0.05*par("usr")[4]), srt = 45, adj = 1, xpd = T, labels =
    gsub("0.99-","", unlist(lapply(strsplit(gsub("\\(","",gsub("]","", colnames(d.counts))), ","),
        function(x) {x <- as.numeric(x) ; x[1] <- x[1]+1 ; paste(x, collapse="-")}))))
legend("topright", bty="n", legend=rownames(d.counts),
    fill=c("darkgreen","red","orange", "brown"), border=NA)

dev.off()


## hypergeometric test for whether TSG genes (or driver genes) are outside of PMDs
d1 <- as.data.frame(mcols(fpkm[, c("pmd.freq", "pmd.freq.bin", "tsg", "tsg.breast", "driv", "other")]))
d1$in.pmd <- ifelse(d1$pmd.freq==0, "out", "in")

doHyperGeom <- function(feat, label) { # feat: colname (e.g. 'tsg'), label: e.g. 'TSG'
  phyper(
    q=sum(d1[,feat]==label & d1$in.pmd=="out", na.rm=T), # white balls drawn
    m=sum(d1[,"in.pmd"]=="out", na.rm=T), # number of white balls present
    n=sum(d1[,"in.pmd"]=="in", na.rm=T), # number of black balls present
    k=sum(d1[,feat]==label, na.rm=T),  # number of balls drawn
    lower.tail=F)
}
doHyperGeom('tsg', 'TSG')
# [1] 8.777996e-16
doHyperGeom('tsg.breast', 'TSG')
# [1] 3.452291e-06
doHyperGeom('driv', 'Driv')
# [1] 1.977179e-11
doHyperGeom('other', 'other')
# [1] 1


###################################################################################
############### consequences of PMD over a (tumor suppressor) gene ################
###################################################################################

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

fpkm <- get(load("fpkm.gr.RData"))
fpkm <- fpkm[seqnames(fpkm) != "chrY"]

# clean up the RNA-seq annotation
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



pat <- intersect(names(pmds), colnames(mcols(fpkm)))
fpkm$tsg <- ifelse(fpkm$Name %in% census.tsg, "yes", "no")
fpkm$tsg.breast <- ifelse(fpkm$Name %in% census.tsg.breast, "yes", "no")
fpkm$driv.breast <- ifelse(fpkm$Name %in% driv.genes, "yes", "no")
fpkm <- fpkm[!duplicated(fpkm$Name)]
pmd.freq <- sapply(names(pmds), inPMD)
pmd.freq <- apply(pmd.freq, 1, sum)
fpkm$pmd.freq <- pmd.freq
fpkm$pmd.freq.bin <- cut(fpkm$pmd.freq, breaks=c(-0.01, seq(0,30, by=3)))


l1 <- list(
  tsg = as.data.frame(mcols(fpkm[, c("Name", "tsg","tsg.breast","driv.breast")])),
  fpkm = as.data.frame(mcols(fpkm[,pat])),
  pmd = sapply(pat, function(x) {ifelse(overlapsAny(fpkm, pmds[[x]], type="within"), "in", "out")})
)

getMeanInOut <- function(x) {
  c("in"=mean(as.numeric(l1$fpkm[x,l1$pmd[x,]=="in"]), na.rm=T), 
    "out"=mean(as.numeric(l1$fpkm[x,l1$pmd[x,]=="out"]), na.rm=T))
}

l1$mean.in.out <- as.data.frame(t(sapply(1:nrow(l1$fpkm), getMeanInOut)))
l1$mean.in.out$in.out <- l1$mean.in.out$`in`-l1$mean.in.out$out

d3 <- do.call(cbind, l1[c("tsg","mean.in.out")])
colnames(d3) <- gsub("^tsg\\.","", gsub("^mean\\.in\\.out\\.","",colnames(d3)))
d3$label <- gsub(" ","\\/", gsub(" $","", gsub("^ ","", gsub(" +"," ", 
    paste(ifelse(d3$tsg=="yes", "TSG", ""), ifelse(d3$tsg.breast=="yes", "TSG_breast", ""), 
    ifelse(d3$driv.breast=="yes", "driver_breast", ""))))))

# plot expression fold change (in/out PMD), for all of them together
p1 <- ggplot(d3[c(d3$tsg=="yes" | d3$tsg.breast=="yes" | d3$driv.breast=="yes") & !is.na(d3$in.out),], 
    aes(reorder(Name, in.out), in.out)) +
      geom_bar(stat="identity", color=NA, aes(fill=label)) +
      scale_fill_manual(values=c("TSG"="red","driver_breast"="brown", "TSG/driver_breast"="brown",
          "TSG/TSG_breast"="orange")) +
      theme_classic() +
      theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1), axis.text=element_text(color="black")) +
      ylab("fold change (log2)") +
      xlab("gene")
ggsave(p1, file="barplot_expression_TSGs_in_out_PMDs_all.pdf", width=18, height=8, scale=0.4) ##### Supplemental Figure 7A ##### 


# add other features to the data frame with mean gene expression changes
d3$type <- ifelse(d3$tsg=="yes", "tsg", "none")
d3$type[d3$tsg.breast=="yes"] <- "tsg.breast"
d3$type[d3$driv.breast=="yes"] <- "driv.breast"
d3 <- d3[!is.na(d3$in.out),]
d3 <- d3[order(d3$in.out),]
d3$rank <- 1:nrow(d3)
d3$pmd.freq <- fpkm$pmd.freq[match(d3$Name, fpkm$Name)]
d3$pmd.freq.bin <- fpkm$pmd.freq.bin[match(d3$Name, fpkm$Name)]
cgi <- get(load("~/BiSeq_BASIS/PMDs_all_samples/cgimeth.basis.RData"))[[1]][,0]
d3$prom.in.cgi <- d3$Name %in% as.character(fpkm$Name[overlapsAny(promoters(fpkm, upstream=1, downstream=1), cgi)])

## t-test: is there generally downregulation inside PMDs?
# for CGI-promoter genes
t.test(d3$in.out[d3$prom.in.cgi])
#  
#          One Sample t-test
#  
#  data:  d3$in.out[d3$prom.in.cgi]
#  t = -6.2798, df = 1204, p-value = 4.725e-10
#  alternative hypothesis: true mean is not equal to 0
#  95 percent confidence interval:
#   -0.9962376 -0.5219313
#  sample estimates:
#   mean of x
#  -0.7590845
#  
t.test(d3$in.out[!d3$prom.in.cgi])
#  
#          One Sample t-test
#  
#  data:  d3$in.out[!d3$prom.in.cgi]
#  t = -2.5182, df = 4194, p-value = 0.01183
#  alternative hypothesis: true mean is not equal to 0
#  95 percent confidence interval:
#   -0.49288419 -0.06137045
#  sample estimates:
#   mean of x
#  -0.2771273

# plot expression change (in/out PMDs), by PMDfreq bin, promoter-CGI genes
p2.1 <- ggplot(d3[d3$prom.in.cgi,], aes(pmd.freq.bin, in.out)) +
        geom_boxplot(outlier.shape=NA, aes(fill=pmd.freq.bin)) +
        geom_hline(yintercept=0, color="darkgrey") +
        scale_fill_manual(values=colorpanel(11, "white","red")[-1]) +
        coord_cartesian(ylim=c(-4,4)) +
        theme_classic() +
        theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"),
             axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
        ylab("fold change inside/outside PMDs (log2 FPKM)") +
        xlab("PMD frequency")
ggsave(p2.1, file="fold_change_inOutPMDs_CGIpromGenes.pdf", width=4, height=2.5, scale=1) ##### Figure 3H #####

# plot expression change (in/out PMDs), by PMDfreq bin, non-CGI-promoter genes
p2.2 <- ggplot(d3[!d3$prom.in.cgi,], aes(pmd.freq.bin, in.out)) +
        geom_boxplot(outlier.shape=NA, aes(fill=pmd.freq.bin)) +
        geom_hline(yintercept=0, color="darkgrey") +
        scale_fill_manual(values=colorpanel(11, "white","red")[-1]) +
        coord_cartesian(ylim=c(-4,4)) +
        theme_classic() +
        theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"),
             axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
        ylab("fold change inside/outside PMDs (log2 FPKM)") +
        xlab("PMD frequency")
ggsave(p2.2, file="fold_change_inOutPMDs_nonCGIpromGenes.pdf", width=4, height=2.5, scale=1) ##### Supplemental Figure 6A #####


# for the affected genes, what is then the correlation DNAme vs. expression?
tmp <- get(load("~/BiSeq_BASIS/PMDs_all_samples/cgimeth.basis.RData"))
cgi.meth <- tmp[[1]][,0]
mcols(cgi.meth) <- as.data.frame((sapply(tmp, function(x) {x@elementMetadata@listData[[1]]})))
ov <- findOverlaps(cgi.meth, promoters(fpkm[,0], upstream=500, downstream=500))

d4 <- data.frame(
   "Name"=fpkm[subjectHits(ov)]$Name,
   as.data.frame(mcols(fpkm[subjectHits(ov), pat])),
   as.data.frame(mcols(cgi.meth[queryHits(ov), pat]))
)
d4 <- d4[complete.cases(d4),]
rownames(d4) <- make.unique(as.character(d4$Name), sep="__")
d4$Name <- NULL
meth.cols <- grep("\\.1$", colnames(d4))
fpkm.cols <- grep("a$", colnames(d4))
m4 <- as.matrix(d4)
corr <- sapply(1:nrow(m4), function(x) {cor(m4[x, meth.cols], m4[x, fpkm.cols])})
names(corr) <- rownames(m4)

d5 <- data.frame(
  "R"=corr,
  "pmd.freq"=fpkm$pmd.freq[match(gsub("__.+$","",names(corr)), fpkm$Name)],
  "pmd.freq.bin"=fpkm$pmd.freq.bin[match(gsub("__.+$","",names(corr)), fpkm$Name)],
  "tsg"=fpkm$tsg[match(gsub("__.+$","",names(corr)), fpkm$Name)],
  "tsg.breast"=fpkm$tsg.breast[match(gsub("__.+$","",names(corr)), fpkm$Name)],
  "driv.breast"=fpkm$driv.breast[match(gsub("__.+$","",names(corr)), fpkm$Name)])
d5$pmd <- ifelse(d5$pmd.freq > 0 , "in", "out")

d6 <- d5[!is.na(d5$R) & d5$pmd.freq > 0 & d5$pmd.freq < 30,]
d6 <- d6[order(d6$R),]
d6$rank <- 1:nrow(d6)
d6$type <- ifelse(d6$tsg=="yes", "tsg", "none")
d6$type[d6$tsg.breast=="yes"] <- "tsg.breast"
d6$type[d6$driv.breast=="yes"] <- "driv.breast"
d6$gene <- rownames(d6)
d6$gene <- gsub("__.+$","", d6$gene)


p3 <- ggplot(d6, aes(rank, R)) + 
      geom_point(color="darkgreen", size=1.75, shape=20) + 
      geom_point(data=d6[!d6$type=="none",] , aes(color=type, size=type, shape=type)) + 
      geom_text_repel(data=d6[!d6$type=="none",], aes(label=gene)) +
      scale_shape_manual(values=c("driv.breast"=17,"tsg"=15)) +
      scale_color_manual(values=c("driv.breast"="brown","tsg"="red")) +
      scale_size_manual(values=c("driv.breast"=6,"tsg"=5)) +
      geom_hline(yintercept = 0) +
      theme_classic() +
      ylab("Pearson R CGImeth vs. expr.")
ggsave(p3, file="correlation_inOutPMDs_allGenes_vs_TSGs.pdf", height=3, width=4, scale=1.5) ##### Supplemental Figure 7C #####


##################################################################################
############### enrichment of GO terms of PMD-downregulated genes ################
##################################################################################

# for which genes does PMD overlap has the strongest effect? -> GO categories?
most.down <- -2.5 # cutoff for most downregulated genes (= log2 scale)
in.out.down <- d3$Name[d3$in.out <= most.down]

# remove pseudogenes and processed transcripts
in.out.down <- as.data.frame(mcols(fpkm[match(in.out.down, fpkm$Name), 3:4]))
in.out.down <- as.character(in.out.down$Name[!in.out.down$Source %in% c("pseudogene","processed_transcript")])
in.out.down <- in.out.down[!grepl("^[ABF][CPLX]\\d{5,}\\.\\d{1,2}$", in.out.down)]
in.out.down <- in.out.down[!grepl("^RP\\d{1,2}-.+\\.\\d$", in.out.down)]
in.out.down.pvals <- sapply(in.out.down, function(x) {
    wilcox.test(as.numeric(l1$fpkm[l1$tsg$Name==x,]) ~ l1$pmd[l1$tsg$Name==x])$p.value})

write.table(in.out.down, quote=F, col.names=F, row.names=F, file=paste0("in.out.down_", most.down, ".txt"))
write.xlsx(in.out.down, file=paste0("in.out.down_", most.down, ".xlsx")) ##### Supplemental Table 4 #####

pdf(paste0("boxplots_in.out.down.", most.down, "fold.p0.05.pdf"))
par(mfrow=c(3,3))
for (i in in.out.down[in.out.down.pvals < 0.05]) {
  boxplot(as.numeric(l1$fpkm[l1$tsg$Name==i,]) ~ l1$pmd[l1$tsg$Name==i], main=i)
}
dev.off()

# gs.files: Obtained from GUI GSEA analysis (http://software.broadinstitute.org/gsea/msigdb/annotate.jsp),
# using the xlsx file above (Suppl Table) as input. Results were exported to text file and imported here.
gs.files <- as.list(list.files(pattern="GSEA_.+\\.txt", 
    path="~/BiSeq_BASIS/PMDs_all_samples/TSG_overlap", full.names=T))
names(gs.files) <- gsub("GSEA_in.out.down_","",gsub("\\.txt","",basename(unlist(gs.files))))

readGSEA <- function(x) {
  read.delim(pipe(paste0("cat ", x, " |grep -A 1000 'Gene Set Name' ",
         " | grep -B 1000 'Gene/Gene Set Overlap Matrix' ",
         " |grep -v 'Gene/Gene Set Overlap Matrix'")), sep="\t", header=T)
}
gs <- lapply(gs.files, readGSEA)
go.terms <- unique(as.character(unlist(lapply(gs, function(x) {x$Gene.Set.Name}))))

# for all GO terms and categories, get the fraction of genes, pval, fdr, for each PMDfreq category
getStats <- function(x) {
  sapply(go.terms, function(y) {
    if (!c(y %in% x$Gene.Set.Name)) {
       frac <- NA
       p <- NA
       q <- NA
    } else {
       s <- x$Gene.Set.Name==y
       frac <- x$k.K[s]
       p <- x$p.value[s]
       q <- x$FDR.q.value[s]
    }
    data.frame("frac"=frac,"pval"=p,"fdr"=q)
  })
}
gs.stats <- lapply(gs, getStats)
for (i in names(gs.stats)) {
  gs.stats[[i]] <- as.data.frame(t(gs.stats[[i]]))
  gs.stats[[i]]$frac <- unlist(gs.stats[[i]]$frac)
  gs.stats[[i]]$pval <- unlist(gs.stats[[i]]$pval)
  gs.stats[[i]]$fdr <- unlist(gs.stats[[i]]$fdr)
  gs.stats[[i]]$term <- rownames(gs.stats[[i]])
}
gs.stats <- do.call(rbind, gs.stats)
gs.stats$cutoff <- sapply(strsplit(rownames(gs.stats), "\\."), function(x) {paste(x[-length(x)], collapse=".")})
rownames(gs.stats) <- NULL

pdf("barplots_GSEA_in.out.down.pdf", height=9, width=3) ##### Supplemental Figure 8A #####
for ( i in unique(gs.stats$cutoff)) {
   tmp <- gs.stats[gs.stats$cutoff==i, c("fdr","term")]
   tmp <- tmp[complete.cases(tmp),]
   tmp <- tmp[grepl("^GO_", tmp$term) | grepl("SMID", tmp$term) | grepl("MODULE", tmp$term),]
   tmp$term <- gsub("_"," ", gsub("^go_","",tolower(tmp$term)))
   p <- ggplot(tmp, aes(reorder(term, fdr), -log10(fdr))) +
        geom_bar(stat="identity", fill="black", color="white") +
        theme_classic() +
        theme(axis.text.x=element_text(color="black", angle=90, hjust=1, vjust=0.5, size=18),
            axis.text.y=element_text(color="black", angle=90, hjust=0.5, vjust=1, size=15),
            axis.ticks=element_line(color="black")) +
        ggtitle(i)
        coord_flip()

   print(p)
}

dev.off()
##################################################################################
################# survival analysis of PMD-downregulated genes ###################
##################################################################################

## clustering of expression data with in.out.down genes --> survival analysis on the obtained clusters
in.out.down.fpkm <- fpkm[fpkm$Name %in% in.out.down, c(2,4, grep("PD", colnames(mcols(fpkm))))]
names(in.out.down.fpkm) <- paste0(in.out.down.fpkm$Name, "_", in.out.down.fpkm$Ensembl)
in.out.down.fpkm$Ensembl <- NULL
in.out.down.fpkm$Name <- NULL
m1 <- as.matrix(mcols(in.out.down.fpkm))
rownames(m1) <- names(in.out.down.fpkm)
m1 <- m1[complete.cases(m1),]

# rank patients based on gene means (in.out.down genes) and make  groups --> survival
in.out.down.fpkm.medians <- apply(as.matrix(mcols(in.out.down.fpkm)), 2, median, na.rm=T)
summr <- summary(in.out.down.fpkm.medians)

# 2-group split
in.out.down.fpkm.ranks <- data.frame(
    "patient"=names(in.out.down.fpkm.medians),
    "cluster"=cut(in.out.down.fpkm.medians, breaks=c(summr[1], summr[3], summr[6]), labels=F),
     "median"=in.out.down.fpkm.medians)
write.xlsx(in.out.down.fpkm.ranks,
    file=paste0("ranks_in.out.down.fpkm_fold", most.down, "_2groups.xlsx"))
pdf(paste0("boxplot_ranks_in.out.down.fpkm_fold", most.down, "_2groups.pdf"))
par(mfrow=c(2,2))
boxplot(in.out.down.fpkm.ranks$median ~ in.out.down.fpkm.ranks$cluster, ylab="median expression")
dev.off()
system(paste0("./survival.R -i ranks_in.out.down.fpkm_fold", most.down, "_2groups.xlsx")) ##### Supplemental Figure 8C #####

