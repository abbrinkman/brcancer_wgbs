library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(reshape2)
library(ggplot2)
library(caret)

pmdmeth.basis <- get(load("~/BiSeq_BASIS/PMDs_all_samples/pmdmeth.basis.RData"))
pmdmeth.basis <- pmdmeth.basis[!names(pmdmeth.basis) %in% c("HMEC","MCF7")]

# fraction of chromosome in PMDs
getPMDlengthFracChromosome <- function(GR) {
  pmd.frac <- sapply(seqlevels(GR)[-grep("chrY", seqlevels(GR))], 
      function(x) {sum(width(GR[seqnames(GR)==x])) / seqlengths(Hsapiens)[x]})
  names(pmd.frac) <- gsub("\\..+$","",names(pmd.frac))
  pmd.frac
}
d1 <- sapply(pmdmeth.basis, getPMDlengthFracChromosome)
d1.1 <- as.data.frame(t(apply(d1, 2, function(x) {
	s <- rownames(d1) == "chrX" ; c("autosomes"=mean(x[!s]), x[s])})))

# is there a link between PMDs on chrX and expression of some gene(s)?
fpkm <- get(load("~/Infinium/BASIS_data_all_TN_QC_analysis/fpkm.z.gr.RData"))
pds <- intersect(names(pmdmeth.basis), colnames(mcols(fpkm)))
fpkm <- fpkm[,c("UNIQID","Ensembl","Source","Name", pds)]
pmd.frac.chrX <- d1["chrX",]

# manually define which genes to check: all related to XCI related to XCI
xci.genes <- c("XIST", "TSIX", "XACT", "JARID2", "EED", "SUZ12", "EZH1", "EZH2", "RBBP7", "RBBP4")
# note: TSIX and XACT are not present in the fpkm file

d.xci <- as.data.frame(fpkm[mcols(fpkm)$Name %in% xci.genes])[,-c(1:8)]
rownames(d.xci) <- d.xci$Name
d.xci$Name <- NULL
d.xci <- data.frame("pmd.frac.chrX"=pmd.frac.chrX[rownames(t(d.xci))], t(d.xci))
d.xci.ctrl <- trainControl(method = "repeatedcv", repeats = 5)
d.xci.lmFit <- train(pmd.frac.chrX ~ ., data=d.xci, method="lm", tuneLength=15, trControl=d.xci.ctrl)
d.xci.lmFit
  # Linear Regression
  # 
  # 25 samples
  #  8 predictors
  # 
  # No pre-processing
  # Resampling: Cross-Validated (10 fold, repeated 5 times)
  # Summary of sample sizes: 23, 22, 22, 22, 23, 22, ...
  # Resampling results:
  # 
  #   RMSE        Rsquared
  #   0.06523443  0.8463505
  # 
  # Tuning parameter 'intercept' was held constant at a value of TRUE

summary(lm(pmd.frac.chrX ~ ., data=d.xci))
  # Call:
  # lm(formula = pmd.frac.chrX ~ ., data = d.xci)
  # 
  # Residuals:
  #       Min        1Q    Median        3Q       Max
  # -0.174017 -0.018186  0.001051  0.025108  0.091667
  # 
  # Coefficients:
  #             Estimate Std. Error t value Pr(>|t|)
  # (Intercept)  0.81602    0.02507  32.552 4.74e-16 ***
  # JARID2       0.02440    0.02318   1.053  0.30819
  # RBBP7       -0.03064    0.02256  -1.358  0.19321
  # XIST         0.08569    0.02844   3.013  0.00825 **
  # EED          0.05391    0.02184   2.469  0.02520 *
  # EZH2        -0.05171    0.02311  -2.238  0.03983 *
  # SUZ12        0.03206    0.04174   0.768  0.45370
  # EZH1        -0.05004    0.01776  -2.818  0.01238 *
  # RBBP4        0.01975    0.02427   0.814  0.42782
  # ---
  # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  # 
  # Residual standard error: 0.06486 on 16 degrees of freedom
  # Multiple R-squared:  0.7409,    Adjusted R-squared:  0.6113
  # F-statistic: 5.719 on 8 and 16 DF,  p-value: 0.001534

# same, but selecting only the genes with highest variable importance (see below)
summary(lm(pmd.frac.chrX ~ ., data=d.xci[,c(1,4,5,6,8)]))
  # 
  # Call:
  # lm(formula = pmd.frac.chrX ~ ., data = d.xci[, c(1, 4, 5, 6,
  #     8)])
  # 
  # Residuals:
  #      Min       1Q   Median       3Q      Max
  # -0.17630 -0.03963  0.01189  0.03710  0.07696
  # 
  # Coefficients:
  #             Estimate Std. Error t value Pr(>|t|)
  # (Intercept)  0.82564    0.01764  46.807  < 2e-16 ***
  # XIST         0.09313    0.02194   4.245 0.000397 ***
  # EED          0.03828    0.01713   2.234 0.037041 *
  # EZH2        -0.04252    0.01793  -2.372 0.027845 *
  # EZH1        -0.05025    0.01652  -3.042 0.006431 **
  # ---
  # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  # 
  # Residual standard error: 0.06249 on 20 degrees of freedom
  # Multiple R-squared:  0.6993,    Adjusted R-squared:  0.6392
  # F-statistic: 11.63 on 4 and 20 DF,  p-value: 4.824e-05

# check how good the model is if you randomly select 10 genes (bootstrapping)
lmRandomGenes <- function() {
  s <- complete.cases(as.data.frame(mcols(fpkm)))
  d.ran <- as.data.frame(mcols(fpkm[s][sample(1:length(fpkm[s]), 10, replace=T), -c(1:3)]))
  rownames(d.ran) <- make.unique(as.character(d.ran$Name))
  d.ran$Name <- NULL
  d.ran <- data.frame("pmd.frac.chrX"=pmd.frac.chrX[rownames(t(d.ran))], t(d.ran))
  d.ran.ctrl <- trainControl(method = "repeatedcv", repeats = 5)
  d.ran.lmFit <- train(pmd.frac.chrX ~ ., data=d.ran, method="lm", tuneLength=15, trControl=d.ran.ctrl)
  d.ran.lmFit$results
}
lm.random <- do.call(rbind, mclapply(1:1000, function(x) {lmRandomGenes()}, mc.cores=10))
mean(lm.random$RMSE <= d.xci.lmFit$results$RMSE)
  # [1] 0
mean(lm.random$Rsquared >= d.xci.lmFit$results$Rsquared)
  # [1] 0.007

# check whether the same model can predict PMD formation on the autosomes
d.xci.auto <- data.frame("pmd.frac.auto"=d1.1[rownames(d.xci),]$autosomes, d.xci[,!colnames(d.xci)=="pmd.frac.chrX"])
d.xci.auto.ctrl <- trainControl(method = "repeatedcv", repeats = 5)
d.xci.auto.lmFit <- train(pmd.frac.auto ~ ., data=d.xci.auto, method="lm", tuneLength=15, trControl=d.xci.auto.ctrl)
d.xci.auto.lmFit
  # Linear Regression
  # 
  # 25 samples
  #  8 predictors
  # 
  # No pre-processing
  # Resampling: Cross-Validated (10 fold, repeated 5 times)
  # Summary of sample sizes: 23, 23, 23, 23, 22, 22, ...
  # Resampling results:
  # 
  #   RMSE       Rsquared
  #   0.1214093  0.7542935
  # 
  # Tuning parameter 'intercept' was held constant at a value of TRUE

summary(lm(pmd.frac.auto ~ ., data=d.xci.auto))
  # Call:
  # lm(formula = pmd.frac.auto ~ ., data = d.xci.auto)
  # 
  # Residuals:
  #       Min        1Q    Median        3Q       Max
  # -0.179946 -0.064111  0.003981  0.056971  0.153834
  # 
  # Coefficients:
  #             Estimate Std. Error t value Pr(>|t|)
  # (Intercept)  0.26589    0.03926   6.772 4.48e-06 ***
  # JARID2       0.06201    0.03631   1.708   0.1070
  # RBBP7       -0.02294    0.03533  -0.649   0.5254
  # XIST        -0.04041    0.04455  -0.907   0.3778
  # EED          0.07042    0.03420   2.059   0.0561 .
  # EZH2        -0.02831    0.03620  -0.782   0.4456
  # SUZ12        0.13542    0.06538   2.071   0.0549 .
  # EZH1        -0.04982    0.02781  -1.791   0.0922 .
  # RBBP4        0.01164    0.03802   0.306   0.7634
  # ---
  # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  # 
  # Residual standard error: 0.1016 on 16 degrees of freedom
  # Multiple R-squared:  0.4627,    Adjusted R-squared:  0.1941
  # F-statistic: 1.722 on 8 and 16 DF,  p-value: 0.1687

## so, this is significant; 
## XIST and PRC2 gene expression is linked to PMD formation on chrX, not on the autosomes

# variable importance
pdf("variableImportance_PMDfracChrX_vs_XCIgeneExpr.pdf", width=2.5, height=3)
plot(varImp(d.xci.lmFit, scale=F)) ##### Supplemental Figure 6C #####
dev.off()


pmd.frac.threshold <- 0.8
d5 <- d.xci[, c("pmd.frac.chrX", "XIST","EED","EZH1","EZH2")]
d5 <- melt(d5, id.vars="pmd.frac.chrX", variable.name="gene", value.name="expression")
d5$gene <- factor(d5$gene, levels=c("XIST","EED","EZH2","EZH1"))
d5$high.low <- ifelse(d5$pmd.frac.chrX > pmd.frac.threshold, "high", "low")

p <- ggplot(d5, aes(gene, expression)) +
   geom_boxplot(aes(fill=high.low)) + 
   scale_fill_manual(values=c("high"="orange","low"="blue"), 
       name=paste0("PMD fraction, (thr=", pmd.frac.threshold, ")")) +
   facet_wrap(~high.low)+
   theme_classic() +
   theme(axis.text.x=element_text(color="black", angle=45, vjust=1, hjust=1), 
       axis.text.y=element_text(color="black"))
ggsave(p, file="boxplots_PMDfracChrX_vs_XCIgeneExpr.pdf", height=5, width=7.7, scale=0.5) ##### Supplemental Figure 6D #####


