library(GenomicRanges)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(betareg)
library(gplots)

pmds <- c(get(load("~/BiSeq_BASIS/PMDs_all_samples/pmdmeth.tcga.RData")),
    get(load("~/BiSeq_BASIS/PMDs_all_samples/pmdmeth.lymph.RData")),
    get(load("~/BiSeq_BASIS/PMDs_all_samples/pmdmeth.basis.RData")))
cgi <- c(get(load("~/BiSeq_BASIS/PMDs_all_samples/cgimeth.tcga.RData")),
    get(load("~/BiSeq_BASIS/PMDs_all_samples/cgimeth.lymph.RData")),
    get(load("~/BiSeq_BASIS/PMDs_all_samples/cgimeth.basis.RData")))

pat <- intersect(names(pmds), names(cgi))
pmds <- pmds[pat]
cgi <- cgi[pat]
names(pmds) <- make.unique(gsub("^PD\\d+a$","BASIS", gsub("_.+$","", names(pmds))), sep="_")
names(cgi) <- make.unique(gsub("^PD\\d+a$","BASIS", gsub("_.+$","", names(cgi))), sep="_")

pat <- intersect(names(pmds), names(cgi))
pat <- as.list(pat)
names(pat) <- unlist(pat)

#also provide cgi as an unlisted GRanges object
cgi.1 <- cgi[[1]] 
for (i in names(cgi)) {mcols(cgi.1)[[i]] <- mcols(cgi[[i]])[[1]]}

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

# Conclusion: when taking all tumors (different tissues) together, 
# the trends are not as clear as with breast cancer (=single tissue), although associations are still
# significant (R2 only up to 0.28)
# Likely, with different tumors it's like comparing apples with oranges
# So, split into groups of the same tumor:

pat.cimp.l <- lapply(pat.cimp, function(x) {split(x, gsub("_\\d+$","", rownames(x)))})

plotCIMP_BetaReg <- function(df, hyper.min.cutoff, tumor.type) {
  plot(df$cimp ~ df$cgi.in.pmd, xlab=" fraction CGIs in PMDs",
        ylab="CIMP\n(fraction hypermethylated CGIs)", main=paste0(tumor.type, 
        "\nCGI hyper meth = ",hyper.min.cutoff), pch=20, col="blue")
  mdl <- betareg(df$cimp ~ df$cgi.in.pmd)
  pval <- signif(summary(mdl)$coefficients$mean[2,"Pr(>|z|)"],2)
  r2 <- round(summary(mdl)$pseudo.r.squared, 2)
  legend("topleft", bty="n", legend=paste0("R2 = ", r2, "\np = ", pval))
  mdl.int <- mdl$coefficients$mean[1]
  mdl.slp <- mdl$coefficients$mean[2]
  abline(a=exp(mdl.int)/(1+exp(mdl.int)), b=exp(mdl.slp)/(1+exp(mdl.slp)), col="darkgrey")
  c("n"=nrow(df), "pval"=pval, "r2"=r2, "hyper.min.cutoff"=hyper.min.cutoff, "tumor.type"=tumor.type)
}

pat.cimp.l.mdl <- list()
n <- 1
pdf("scatterplots_CIMP_vs_CGI-PMDs_perTumorType.pdf")
par(mfrow=c(3,3), mar=c(5,5,4,2))
for (i in names(pat.cimp.l)) {
  message(i)
  for (j in names(pat.cimp.l[[i]])) {
    message("  ",j)
    if (nrow(pat.cimp.l[[i]][[j]]) > 2) {
      pat.cimp.l.mdl[[n]] <- plotCIMP_BetaReg(pat.cimp.l[[i]][[j]], i, j)
      n <- n+1
    }
  }  
}
dev.off()

# combine the regressions into one plot
pat.cimp.l.mdl <- as.data.frame(do.call(rbind, pat.cimp.l.mdl))
pat.cimp.l.mdl$pval <- as.numeric(as.character(pat.cimp.l.mdl$pval))
pat.cimp.l.mdl$hyper.min.cutoff <- as.numeric(as.character(pat.cimp.l.mdl$hyper.min.cutoff))
pat.cimp.l.mdl$r2 <- as.numeric(as.character(pat.cimp.l.mdl$r2))
pat.cimp.l.mdl$n <- as.numeric(as.character(pat.cimp.l.mdl$n))
pat.cimp.l.mdl$tumor.type.n <- paste0(pat.cimp.l.mdl$tumor.type, " (n=", pat.cimp.l.mdl$n, ")")

p0 <- ggplot(pat.cimp.l.mdl[pat.cimp.l.mdl$hyper.min.cutoff==0.3 & pat.cimp.l.mdl$tumor.type != "GCB",], 
          aes(jitter(r2), -log10(pval))) + 
      geom_point(color="black", fill="blue", shape=21, size=3) + 
      geom_text_repel(aes(label=tumor.type.n), size=1) +
      theme_classic() +
      theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black")) +
      xlab("R2 (variation explained)") 
ggsave(file="regressions_CIMP_vs_CGI-PMDs_perTumorType_statsPlot.pdf", width=7, height=7, scale=0.4) #### Figure 3F ####


