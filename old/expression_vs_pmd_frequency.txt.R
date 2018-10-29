library(GenomicRanges)
library(ggplot2)
library(reshape2)
library(gplots)

# get the PMDs
pmds <- get(load("~/BiSeq_BASIS/PMDs_all_samples/pmdmeth.basis.RData"))
pmds <- pmds[!names(pmds) %in% c("HMEC","MCF7")]

# get RNAseq data
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



# determine for each gene the PMD frequency
inPMD <- function(PD) {
  ov <- findOverlaps(fpkm, pmds[[PD]], type="within", select="all")
  in.pmd <- ifelse(1:length(fpkm) %in% queryHits(ov), 1, 0)
  in.pmd
}
pmd.freq <- sapply(names(pmds), inPMD)
pmd.freq <- apply(pmd.freq, 1, sum)
pmd.freq.wgbs <- sapply(names(pmds)[names(pmds) %in% colnames(mcols(fpkm))], inPMD)
pmd.freq.wgbs <- apply(pmd.freq.wgbs, 1, sum)

fpkm$pmd.freq <- pmd.freq
fpkm$pmd.freq.bin <- cut(fpkm$pmd.freq, breaks=c(-0.01, seq(0,30, by=3)))
fpkm$pmd.freq.wgbs <- pmd.freq.wgbs

# match breaks to the number of bins in the meth-mut analysis
wgbs.breaks <- c(-0.01, seq(0, 30, by=3))
wgbs.breaks[wgbs.breaks==24] <- 25
wgbs.breaks <- wgbs.breaks[wgbs.breaks <= 25]
fpkm$pmd.freq.bin.wgbs <- cut(fpkm$pmd.freq.wgbs, breaks=wgbs.breaks)

# determine the mean expression for each gene, in samples from the WGBS cohort only
fpkm$mean.expr.wgbs <- apply(as.matrix(mcols(fpkm[,colnames(mcols(fpkm)) %in% names(pmds)])), 1, mean, na.rm=T)

# determine the standard deviation for each gene, in samples from the WGBS cohort only
fpkm$stdev.wgbs <- apply(as.matrix(mcols(fpkm[,colnames(mcols(fpkm)) %in% names(pmds)])), 1, sd, na.rm=T)

# determine the mean expression for each gene, in all samples from the RNA-seq cohort
fpkm$mean.expr.all <- apply(as.matrix(mcols(fpkm[,grep("^PD",colnames(mcols(fpkm)))])), 1, mean, na.rm=T)

# determine the standard deviation for each gene, in all samples from the RNA-seq cohort
fpkm$stdev.all <- apply(as.matrix(mcols(fpkm[,grep("^PD",colnames(mcols(fpkm)))])), 1, sd, na.rm=T)

 #### Supplemental Figure 3B ####
p1 <- ggplot( as.data.frame(mcols(fpkm)), aes(as.factor(pmd.freq.bin.wgbs), mean.expr.wgbs)) +
      geom_boxplot(outlier.shape=NA, aes(fill=pmd.freq.bin.wgbs)) +
      scale_fill_manual(values=colorpanel(length(unique(fpkm$pmd.freq.bin.wgbs)), "white", "red")) +
      coord_cartesian(ylim=c(-12,12)) +
      theme_classic() + 
      theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"), 
          axis.text.x=element_text(angle=45, hjust=1, vjust=1)) + 
      ylab("mean expression (log2 FPKM)") + 
      xlab("PMD frequency") + 
      ggtitle("25 WGBS samples")
ggsave(p1, file="mean_expression_vs_PMDfrequency_WGBS.pdf", width=4, height=2.8, scale=1)

 #### Figure 2F ####
p2 <- ggplot( as.data.frame(mcols(fpkm)), aes(as.factor(pmd.freq.bin), mean.expr.all)) +
      geom_boxplot(outlier.shape=NA, aes(fill=pmd.freq.bin)) +
      scale_fill_manual(values=colorpanel(length(unique(fpkm$pmd.freq.bin)), "white", "red")) +
      coord_cartesian(ylim=c(-12,12)) +
      theme_classic() +
      theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"),
           axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
      ylab("mean expression (log2 FPKM)") +
      xlab("PMD frequency") +
      ggtitle("all RNAseq_samples")
ggsave(p2, file="mean_expression_vs_PMDfrequency_all.pdf", width=4, height=2.8, scale=1)

 #### Supplemental Figure 3B ####
p3 <- ggplot( as.data.frame(mcols(fpkm)), aes(as.factor(pmd.freq.bin.wgbs), stdev.wgbs)) +
      geom_boxplot(outlier.shape=NA, aes(fill=pmd.freq.bin.wgbs)) +
      scale_fill_manual(values=colorpanel(length(unique(fpkm$pmd.freq.bin.wgbs)), "white", "red")) +
      coord_cartesian(ylim=c(0,5)) +
      theme_classic() +
      theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"), 
            axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
      ylab("mean stdev (log2 FPKM)") +
      xlab("PMD frequency") +
      ggtitle("25 WGBS samples")
ggsave(p3, file="stdev_expression_vs_PMDfrequency_WGBS.pdf", width=4, height=2.8, scale=1)

 #### Figure 2F ####
p4 <- ggplot( as.data.frame(mcols(fpkm)), aes(as.factor(pmd.freq.bin), stdev.all)) +
      geom_boxplot(outlier.shape=NA, aes(fill=pmd.freq.bin)) +
      scale_fill_manual(values=colorpanel(length(unique(fpkm$pmd.freq.bin)), "white", "red")) +
      coord_cartesian(ylim=c(0,5)) +
      theme_classic() +
      theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"),
            axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
      ylab("mean stdev (log2 FPKM)") +
      xlab("PMD frequency") +
      ggtitle("all RNAseq_samples")
ggsave(p4, file="stdev_expression_vs_PMDfrequency_all.pdf", width=4, height=3, scale=1)

## do the statistics
summary(lm(mean.expr.wgbs ~ pmd.freq, data=as.data.frame(mcols(fpkm))))
#  
#  Call:
#  lm(formula = mean.expr.wgbs ~ pmd.freq, data = as.data.frame(mcols(fpkm)))
#  
#  Residuals:
#       Min       1Q   Median       3Q      Max
#  -273.890   -2.133    0.546    2.308   16.846
#  
#  Coefficients:
#              Estimate Std. Error t value Pr(>|t|)
#  (Intercept)  1.16136    0.03347   34.70   <2e-16 ***
#  pmd.freq    -0.18472    0.00385  -47.98   <2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  
#  Residual standard error: 4.347 on 20898 degrees of freedom
#    (1624 observations deleted due to missingness)
#  Multiple R-squared:  0.09922,   Adjusted R-squared:  0.09918
#  F-statistic:  2302 on 1 and 20898 DF,  p-value: < 2.2e-16
#  
t.test(mean.expr.wgbs ~ pmd.freq > 0, data=as.data.frame(mcols(fpkm)))
#  
#          Welch Two Sample t-test
#  
#  data:  mean.expr.wgbs by pmd.freq > 0
#  t = 51.688, df = 8765.6, p-value < 2.2e-16
#  alternative hypothesis: true difference in means is not equal to 0
#  95 percent confidence interval:
#   3.746074 4.041411
#  sample estimates:
#  mean in group FALSE  mean in group TRUE
#             1.744915           -2.148827
#  
summary(lm(stdev.wgbs ~ pmd.freq, data=as.data.frame(mcols(fpkm))))
#  
#  Call:
#  lm(formula = stdev.wgbs ~ pmd.freq, data = as.data.frame(mcols(fpkm)))
#  
#  Residuals:
#     Min     1Q Median     3Q    Max
#   -2.20  -0.61  -0.32   0.15 337.49
#  
#  Coefficients:
#              Estimate Std. Error t value Pr(>|t|)
#  (Intercept) 1.224989   0.037328  32.817  < 2e-16 ***
#  pmd.freq    0.032600   0.004471   7.291 3.18e-13 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  
#  Residual standard error: 4.825 on 20455 degrees of freedom
#    (2067 observations deleted due to missingness)
#  Multiple R-squared:  0.002592,  Adjusted R-squared:  0.002544
#  F-statistic: 53.16 on 1 and 20455 DF,  p-value: 3.179e-13
#  
summary(lm(stdev.wgbs ~ pmd.freq > 0, data=as.data.frame(mcols(fpkm))))
#  
#  Call:
#  lm(formula = stdev.wgbs ~ pmd.freq > 0, data = as.data.frame(mcols(fpkm)))
#  
#  Residuals:
#     Min     1Q Median     3Q    Max
#   -1.87  -0.56  -0.29   0.13 337.31
#  
#  Coefficients:
#                   Estimate Std. Error t value Pr(>|t|)
#  (Intercept)       1.09506    0.04084   26.82   <2e-16 ***
#  pmd.freq > 0TRUE  0.77107    0.07223   10.68   <2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  
#  Residual standard error: 4.818 on 20455 degrees of freedom
#    (2067 observations deleted due to missingness)
#  Multiple R-squared:  0.005541,  Adjusted R-squared:  0.005492
#  F-statistic:   114 on 1 and 20455 DF,  p-value: < 2.2e-16



# is the increased expression variability inside PMDs due to low expression or due to being inside a PMD?
d1 <- as.data.frame(mcols(fpkm[,c(
    "pmd.freq", "pmd.freq.bin", "mean.expr.wgbs", "stdev.wgbs", "mean.expr.all", "stdev.all")]))
d1$pmd <- ifelse(d1$pmd.freq > 0, "in", "out")

mdl1 <- lm(stdev.wgbs ~ mean.expr.wgbs + pmd.freq + pmd, data=d1)

summary(mdl1)
#  
#  Call:
#  lm(formula = stdev.wgbs ~ mean.expr.wgbs + pmd.freq + pmd, data = d1)
#  
#  Residuals:
#      Min      1Q  Median      3Q     Max
#   -7.879  -1.381   0.110   1.249 172.747
#  
#  Coefficients:
#                  Estimate Std. Error  t value Pr(>|t|)
#  (Intercept)     0.571583   0.068936    8.292  < 2e-16 ***
#  mean.expr.wgbs -0.768135   0.006038 -127.225  < 2e-16 ***
#  pmd.freq       -0.026935   0.004621   -5.829 5.67e-09 ***
#  pmdout          1.857829   0.077510   23.969  < 2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  
#  Residual standard error: 3.6 on 20453 degrees of freedom
#    (2067 observations deleted due to missingness)
#  Multiple R-squared:  0.4449,    Adjusted R-squared:  0.4448
#  F-statistic:  5463 on 3 and 20453 DF,  p-value: < 2.2e-16
  

mdl2 <- lm(stdev.wgbs ~ mean.expr.wgbs, data=d1)

anova(mdl1, mdl2)
#  Analysis of Variance Table
#  
#  Model 1: stdev.wgbs ~ mean.expr.wgbs + pmd.freq + pmd
#  Model 2: stdev.wgbs ~ mean.expr.wgbs
#    Res.Df    RSS Df Sum of Sq      F    Pr(>F)
#  1  20453 265017
#  2  20455 282909 -2    -17892 690.41 < 2.2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# So: being inside a PMD is heavily associated with higher StDev

# simple MWU-test on stDevs, in/out PMDs
wilcox.test(d1$stdev.wgbs ~ d1$pmd)
#  
#          Wilcoxon rank sum test with continuity correction
#  
#  data:  d1$stdev.wgbs by d1$pmd
#  W = 65712000, p-value < 2.2e-16
#  alternative hypothesis: true location shift is not equal to 0



