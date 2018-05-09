# CIMP definitions:
# there is no universal standard or consensus with respect to defining CIMP ... (Hughes)
# 1. Increased prevalence of CpG island promoter methylation (Hughes)
# 2. CIMP-tumors represents a subtype of cancers that exhibit concurrent hypermethylation of multiple CGIs (Suzuki)

# So, in our dataset we count the number of hypermethylated CGIs for each sample --> CIMP-low/high 
# What is hypermethylated? Try different thresholds, ranging between 0.2 - 0.4

library(GenomicRanges)
library(ggplot2)
library(reshape2)
library(betareg)
library(gplots)
library(openxlsx)

pmds <- get(load("~/BiSeq_BASIS/PMDs_all_samples/pmdmeth.basis.RData"))
cgi <- get(load("~/BiSeq_BASIS/PMDs_all_samples/cgimeth.basis.RData"))

#also provide cgi as an unlisted GRanges object
cgi.1 <- cgi[[1]] 
for (i in names(cgi)) {mcols(cgi.1)[[i]] <- mcols(cgi[[i]])[[1]]}

pat <- intersect(names(pmds), names(cgi))
pat <- pat[!pat %in% c("HMEC","MCF7","PD9590a")]
pat <- pat[!pat %in% c("PD9590a")]
pat <- as.list(pat)
names(pat) <- unlist(pat)

# dataframe with mean CGI meth, mean PMD meth, PMD fraction, CGI methylation in bins
d1 <- data.frame(
  "cgi.mean"=sapply(pat, function(x) {mean(mcols(cgi[[x]])[[1]], na.rm=T)}),
  as.data.frame(t(sapply(pat, function(x) {table(cut(mcols(cgi[[x]])[[1]], 
      breaks=c(-0.01, 0, 0.2, 0.4, 0.6, 0.8, 1)))}))),
  "cgi.in.pmd"=sapply(pat, function(x) {sum(overlapsAny(cgi[[x]], pmds[[x]]))}),
  "pmd.mean"=sapply(pat, function(x) {mean(mcols(pmds[[x]])[[1]], na.rm=T)}),
  "pmd.frac"=sapply(pat, function(x) {sum(as.numeric(width(pmds[[x]])))/3.2e09})
)
colnames(d1) <- gsub("\\.$","", gsub(".0.01","",gsub("^X\\.+","cgi.", colnames(d1))))

# What is the association between the number of CGIs inside PMDs and the number of hypermethylated CGIs?  
getCIMPvalue_CGIinPMDs <- function(pd, hyper.min) { # hyper.min: minimal CGI methylation to call it 'hypermethylated'
  cgi.hyper <- mcols(cgi[[pd]])[[1]] >= hyper.min
  cgi.pmd <- overlapsAny(cgi[[pd]], pmds[[pd]], type="within")
  ncgi.hyper <- mean(cgi.hyper, na.rm=T)
  ncgi.in.pmd <- mean(cgi.pmd, na.rm=T)
  c("cimp"=ncgi.hyper, "cgi.in.pmd"=ncgi.in.pmd)
}
pat.cimp <- sapply(as.character(c(0.2, 0.25, 0.3, 0.4)), function(x) {
    as.data.frame(t(sapply(pat, getCIMPvalue_CGIinPMDs, as.numeric(x))))}, USE.NAMES=T, simplify=F)

pdf("scatterplots_CIMP_vs_CGI-PMDs.pdf") #### Figure 3D ####
par(mfrow=c(3,3), mar=c(5,5,4,2))
lapply(names(pat.cimp), function(x) {
  plot(pat.cimp[[x]]$cimp ~ pat.cimp[[x]]$cgi.in.pmd, xlab=" fraction CGIs in PMDs", 
      ylab="CIMP\n(fraction hypermethylated CGIs)", main=paste0("CGI hyper meth = ",x), pch=20, col="blue")
  mdl <- betareg(pat.cimp[[x]]$cimp ~ pat.cimp[[x]]$cgi.in.pmd)
  pval <- signif(summary(mdl)$coefficients$mean[2,"Pr(>|z|)"],2)
  r2 <- round(summary(mdl)$pseudo.r.squared, 2)
  legend("topleft", bty="n", legend=paste0("R2 = ", r2, "\np = ", pval))
  mdl.int <- mdl$coefficients$mean[1]
  mdl.slp <- mdl$coefficients$mean[2]
  abline(a=exp(mdl.int)/(1+exp(mdl.int)), b=exp(mdl.slp)/(1+exp(mdl.slp)), col="darkgrey")
  summary(mdl)
})
dev.off()
# Conclusion: this clearly shows that the more CGIs end up inside PMDs, the higher the CIMP value.
# This association is strong, for hypermethylation from low values onwards (0.2-0.4)
# However, you can explain only about half of the hypermethylated CGIs in this way. What about the rest?


## Determine the mean CGI methylation and the StDev vs PMD frequency, 
#  to show that CGI methylation outside PMDs is rather stable
cgi.1$pmd.freq <- apply(sapply(pat, function(x) {overlapsAny(cgi[[x]], pmds[[x]], type="within")}), 1, sum)
cgi.1$pmd.freq.bin <- cut(cgi.1$pmd.freq, breaks=c(-0.01, seq(0, 30, by=3)))
cgi.1$cgi.mean <- apply(as.matrix(mcols(cgi.1[,-grep("^pmd", colnames(mcols(cgi.1)))])), 1, mean, na.rm=T)
cgi.1$cgi.sd <- apply(as.matrix(mcols(cgi.1[,-grep("^pmd", colnames(mcols(cgi.1)))])), 1, sd, na.rm=T)
cgi.1 <- cgi.1[!seqnames(cgi.1) %in% "chrY"]

# export CGIs with their PMD frequencies to a file as Suppl Table
cgi.pmdfreq.suppl <- data.frame(cgi.1[,"pmd.freq"])
cgi.pmdfreq.suppl$strand <- NULL
cgi.pmdfreq.suppl$width <- NULL
colnames(cgi.pmdfreq.suppl) <- c("chr","start","end","PMD_frequency")
wb <- createWorkbook()
addWorksheet(wb, "CpG island PMD frequency")
writeData(wb, sheet="CpG island PMD frequency", x=cgi.pmdfreq.suppl)
saveWorkbook(wb, file="cgi.pmdfreq.suppl.xlsx", overwrite=T) #### Supplemental Table 2 ####


# mean CGI methylation
p1.1 <- ggplot(as.data.frame(mcols(cgi.1)), aes(pmd.freq.bin, cgi.mean)) +
        geom_boxplot(outlier.shape=NA, aes(fill=pmd.freq.bin)) +
        scale_fill_manual(values=colorpanel(length(unique(cgi.1$pmd.freq.bin)), "white","red")) +
        theme_classic() +
        theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"),
             axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
        ylab("CGI meth") +
        xlab("PMD frequency")
ggsave(p1.1, file="boxplot_mean_CGImeth_inOutPMDs_noOutliers.pdf", width=4, height=2.5, scale=1)

# StDev CGI methylation
p2.1 <- ggplot(as.data.frame(mcols(cgi.1)), aes(pmd.freq.bin, cgi.sd)) +
        geom_boxplot(outlier.shape=NA, aes(fill=pmd.freq.bin)) +
        scale_fill_manual(values=colorpanel(length(unique(cgi.1$pmd.freq.bin)), "white","red")) +
        theme_classic() +
        theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"),
             axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
        ylab("StDev CGI meth") +
        xlab("PMD frequency")
ggsave(p2.1, file="boxplot_StDev_CGImeth_inOutPMDs_noOutliers.pdf", width=4, height=2.5, scale=1) #### Figure 3E ####

