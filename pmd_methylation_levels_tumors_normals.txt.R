library(GenomicRanges)
library(ggplot2)
library(reshape2)
library(stringr)

pmdmeth.select <- get(load("~/BiSeq_BASIS/PMDs_other_celltypes/pmdmeth.select.RData"))
cgimeth.select <- get(load("~/BiSeq_BASIS/PMDs_other_celltypes/cgimeth.select.RData"))

# calculate the mean for each sample
getMean <- function(x) {
  sapply(x, function(y) {mean(y@elementMetadata@listData[[1]], na.rm=T)})
}
pmdmeth <- getMean(pmdmeth.select)
cgimeth <- getMean(cgimeth.select)

pmdmeth <- data.frame("sample"=names(pmdmeth), "pmd.meth"=pmdmeth)
cgimeth <- data.frame("sample"=names(cgimeth), "cgi.meth"=cgimeth)
d1 <- merge(pmdmeth, cgimeth, by.x="sample", by.y="sample")

d2 <- melt(d1, variable.name="element", value.name="meth") 
d2$element <- toupper(gsub("\\.meth","", d2$element))
d2$status <- sapply(str_split(d2$sample, "_"), function(x) {x[1]})
d2$studyID <- sapply(str_split(d2$sample, "_"), function(x) {x[2]})
d2$sampleID <- sapply(str_split(d2$sample, "_"), function(x) {x[3]})
d2$tissue <- sapply(str_split(d2$sampleID, "\\."), function(x) {x[1]})

# organize by tissue type
d2$tissue2 <- d2$tissue
d2$tissue2[d2$studyID=="Blueprint" & d2$status=="normal"] <- "Heamatopoietic"
d2$tissue2[grepl("ventricle", d2$tissue, ignore.case=T)] <- "Heart"
d2$tissue2[grepl("atrium", d2$tissue, ignore.case=T)] <- "Heart"
d2$tissue2[grepl("germinalcenter", d2$tissue, ignore.case=T)] <- "Heamatopoietic"
d2$tissue2[grepl("colon", d2$tissue, ignore.case=T)] <- "Colon"
d2$tissue2[grepl("intestine", d2$tissue, ignore.case=T)] <- "Intestine"
d2$tissue2[grepl("bowel", d2$tissue, ignore.case=T)] <- "Intestine"
d2$tissue2[grepl("stomach", d2$tissue, ignore.case=T)] <- "Stomach"
d2$tissue2[grepl("esophagus", d2$tissue, ignore.case=T)] <- "Esophagus"
d2$tissue2[grepl("psoas", d2$tissue, ignore.case=T)] <- "Muscle"
d2$tissue2[grepl("lung", d2$tissue, ignore.case=T)] <- "Lung"
d2$tissue2 <- gsub("liver","Liver", d2$tissue2)
d2$tissue2.meth <- sapply(1:nrow(d2), function(x) {mean(d2$meth[d2$tissue2==d2$tissue2[x]])})

# for the record, how much do we have for each dataset, tumor/normal?
table(sapply(strsplit(as.character(unique(d2$sample)), "_"), function(x) {
    paste(x[1:2], collapse="_")}))
# normal_Blueprint  normal_Kretzmer        normal_Li   normal_Roadmap
#              120                4                7               20
#   normal_Schultz  tumor_Blueprint   tumor_Kretzmer         tumor_Li
#               36               41               19                5
#       tumor_TCGA  tumor_ThisStudy
#               38               30
table(sapply(strsplit(as.character(unique(d2$sample)), "_"), function(x) {x[1]}))
# normal  tumor
#    187    133

table(d2$tissue2[d2$status=="tumor"])
#    ALL    AML     BL   BLCA breast    CLL   COAD     FL    GBM  Liver   LUAD
#     14     38     26     12     60     12      4     12     12      4     12
#   Lung   LUSC    MCL     MM   READ   STAD    TPL   UCEC
#      6     10     10      2      6     10      6     10


p1 <- ggplot(d2[d2$element=="PMD",], aes(reorder(tissue2, tissue2.meth), meth)) + 
  geom_jitter(aes(color=status), width=0.2, size=1.5) +
  scale_color_manual(values=c("tumor"="red","normal"="darkgreen")) +
  facet_wrap(~status, scales="free_x") +
  theme_linedraw() +
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1, color="black"), 
      axis.text.y=element_text(color="black")) +
  xlab("")
ggsave(p1, file="jitterplot_PMD_methylation_BASIS_TCGA_lymph_roadmap_li_Blueprint_byTissue.pdf", 
     height=2.5, width=8, scale=1.2) ##### Figure 4A #####


# do the statistics
wilcox.test(d2[d2$element=="PMD","meth"] ~ d2[d2$element=="PMD","status"])
#  
#          Wilcoxon rank sum test with continuity correction
#  
#  data:  d2[d2$element == "PMD", "meth"] by d2[d2$element == "PMD", "status"]
#  W = 20848, p-value < 2.2e-16
#  alternative hypothesis: true location shift is not equal to 0

t.test(d2[d2$element=="PMD","meth"] ~ d2[d2$element=="PMD","status"])
#         Welch Two Sample t-test
# 
# data:  d2[d2$element == "PMD", "meth"] by d2[d2$element == "PMD", "status"]
# t = 12.574, df = 209.4, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  0.09033241 0.12392510
# sample estimates:
# mean in group normal  mean in group tumor
#            0.8080934            0.7009647


# boxplot of PMDmeth, tumor vs normal tissues
pdf("boxplot_PMDmeth_tumors_vs_normal_tissues.pdf") ##### Supplemental Figure 9A #####
par(mfrow=c(3,3))
boxplot(meth ~ status, data=d2[d2$element=="PMD",], ylab="mean PMD meth")
legend("bottomleft", bty="n", legend=paste0("p = ", signif(t.test(meth ~ status, data=d2[d2$element=="PMD",])$p.value, 3)))
dev.off()


