library(GenomicRanges)
library(reshape2)
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg19)
library(gplots)

# get the PMDs
pmds <- get(load("~/BiSeq_BASIS/PMDs_all_samples/pmdmeth.basis.RData"))
pmds <- pmds[grep("^PD", names(pmds))]
pmds <- lapply(pmds, function(x) {x[!seqnames(x) %in% c("chrX","chrY")]})

# get the hg19 genome
hg19 <- BSgenome.Hsapiens.UCSC.hg19
hg19 <- as.data.frame(seqlengths(hg19))
colnames(hg19) <- "end"
hg19$start <- 1
hg19$chr <- rownames(hg19)
hg19 <- makeGRangesFromDataFrame(hg19)
hg19 <- hg19[!grepl("chrUn", seqnames(hg19))]
hg19 <- hg19[!grepl("random", seqnames(hg19))]
hg19 <- hg19[!grepl("hap", seqnames(hg19))]
hg19 <- hg19[!grepl("chrM", seqnames(hg19))]
hg19 <- hg19[!grepl("chrX", seqnames(hg19))]
hg19 <- hg19[!grepl("chrY", seqnames(hg19))]

# remove centromers from genome
centro <- read.delim(pipe("cat ~/Datarepository/hg19_cytoBand.bed |grep -e acen -e gvar -e stalk"),
    header=F, sep="\t")
colnames(centro) <- c("chr","start","end","arm","type")
centro <- makeGRangesFromDataFrame(centro, keep.extra.columns=T)
hg19 <- setdiff(hg19, centro)

# get the mutation data
read.vcf <- function(FILE) {
  str <- paste0("grep -v ^## ", FILE, " | sed 's/^#//'")
  v <- read.delim(pipe(str), header=T, sep="\t")
  colnames(v) <- gsub("-",".", colnames(v))
  v
}
vcf.sub <- read.vcf(
  "~/Infinium/BASIS_data_all_TN_QC_analysis/Caveman_560_20Nov14_clean.txt")
vcf.indel <- read.vcf(
  "~/Infinium/BASIS_data_all_TN_QC_analysis/Pindel_560_20Nov14_noFemaleY_SUMMS_REP.txt.1000gFiltered")
merge.cols <- intersect(colnames(vcf.indel), colnames(vcf.sub))
vcf.all <- rbind(vcf.sub[,merge.cols], vcf.indel[,merge.cols])

#filter for PASS mutations
vcf.all <- vcf.all[vcf.all$Filter=="PASS",]



# split into one df per patient
vcf.all.pat <- split(vcf.all, vcf.all$Sample)
unique.names <- gsub("a[12]$","a", names(vcf.all.pat))
if (sum(duplicated(unique.names)) == 0) {
  names(vcf.all.pat) <- unique.names
}
rm(unique.names)

# get all overlapping patients
patients <- intersect(names(vcf.all.pat), names(pmds))
patients <- as.list(patients)
names(patients) <- unlist(patients)

# take data only for relevant patients (=WGBS cohort)
mut <- vcf.all.pat[names(patients)]
pmds <- pmds[names(patients)]

# convert mutation data to GRanges
mutToGRanges <- function(x) {
  x$chr <- paste0('chr', x$Chrom)
  x$Chrom <- NULL
  x$start <- x$Pos
  x$end <- x$Pos+1
  x$Pos <- NULL
  makeGRangesFromDataFrame(x, keep.extra.columns=T)
}
mut <- lapply(mut, mutToGRanges)

# remove unnecessary stuff
rm(vcf.all.pat, vcf.all, vcf.sub, vcf.indel)
gc()

# tile the genome
tile.size <- 30000
tiles <- unlist(tile(hg19, width=tile.size))

# determine the PMD frequency for each tile
tiles$pmd.freq <- apply(sapply(pmds, function(x) {ifelse(overlapsAny(tiles, x, type="within"), 1, 0)}), 1, sum)
tiles$pmd.freq.bin <- cut(tiles$pmd.freq, breaks=c(-0.01, seq(0, length(patients), by=3)))
tiles$length <- width(tiles)


# count for each patient the mutations in each tile
countMut <- function(pd, type) {
  countOverlaps(tiles, mut[[pd]][mut[[pd]]$Type==type])
}
l1 <- lapply(patients, function(x) {
    as.data.frame(sapply(c("Sub", "Complex", "Del", "Ins"), function(y) {countMut(x, y)}, USE.NAMES=T, simplify=T))})

# sum the mutations per tile
sumMut <- function(pd) {
  d.out <- as.data.frame(sapply(c("Sub", "Ins", "Del", "Complex"), function(x) {
      tapply(l1[[pd]][,x], tiles$pmd.freq.bin, sum)}))
  d.out$bases <- as.numeric(tapply(tiles$length, tiles$pmd.freq.bin, sum))
  d.out.1 <- as.data.frame(sapply(c("Sub", "Ins", "Del", "Complex"), function(x) {d.out[,x]/(d.out$bases/1e06)}))
  d.out.1$pmd.freq.bin <- rownames(d.out)
  d.out.1
}

l2 <- lapply(patients, sumMut)
d2 <- melt(l2, variable.name="Type", value.name="mut.Mb")
d2$pmd.freq.bin <- factor(d2$pmd.freq, levels=unique(d2$pmd.freq.bin))

p.sub <- ggplot(d2[d2$Type=="Sub",], aes(pmd.freq.bin, mut.Mb)) + 
         geom_boxplot(aes(fill=pmd.freq.bin), outlier.shape=NA) + 
         scale_fill_manual(values=colorpanel(length(unique(d2$pmd.freq.bin)), "white", "red")) + 
         coord_cartesian(ylim=c(0,8)) +
         theme_classic() +
         theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"),
              axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
         ylab("mean mut/Mb") +
         xlab("PMD frequency") +
         ggtitle("Subst, 25 WGBS samples")
ggsave(p.sub, file="boxplots_Sub_perMb_vs_PMDfreq_WGBS.pdf", width=4, height=2.8, scale=1) #### Supplemental Figure 3C ####

p.ins <- ggplot(d2[d2$Type=="Ins",], aes(pmd.freq.bin, mut.Mb)) +
         geom_boxplot(aes(fill=pmd.freq.bin), outlier.shape=NA) +
         scale_fill_manual(values=colorpanel(length(unique(d2$pmd.freq.bin)), "white", "red")) +
         coord_cartesian(ylim=c(0,0.1)) +
         theme_classic() +
         theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"),
              axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
         ylab("mean mut/Mb") +
         xlab("PMD frequency") +
         ggtitle("Ins, 25 WGBS samples")
ggsave(p.ins, file="boxplots_Ins_perMb_vs_PMDfreq_WGBS.pdf", width=4, height=2.8, scale=1) #### Supplemental Figure 3C ####

p.del <- ggplot(d2[d2$Type=="Del",], aes(pmd.freq.bin, mut.Mb)) +
         geom_boxplot(aes(fill=pmd.freq.bin), outlier.shape=NA) +
         scale_fill_manual(values=colorpanel(length(unique(d2$pmd.freq.bin)), "white", "red")) +
         coord_cartesian(ylim=c(0,0.15)) +
         theme_classic() +
         theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"),
              axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
         ylab("mean mut/Mb") +
         xlab("PMD frequency") +
         ggtitle("Del, 25 WGBS samples")
ggsave(p.del, file="boxplots_Del_perMb_vs_PMDfreq_WGBS.pdf", width=4, height=2.8, scale=1) #### Supplemental Figure 3C ####

p.complex <- ggplot(d2[d2$Type=="Complex",], aes(pmd.freq.bin, mut.Mb)) +
         geom_boxplot(aes(fill=pmd.freq.bin), outlier.shape=NA) +
         scale_fill_manual(values=colorpanel(length(unique(d2$pmd.freq.bin)), "white", "red")) +
         coord_cartesian(ylim=c(0,0.015)) +
         theme_classic() +
         theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"),
              axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
         ylab("mean mut/Mb") +
         xlab("PMD frequency") +
         ggtitle("Complex, 25 WGBS samples")
ggsave(p.complex, file="boxplots_Complex_perMb_vs_PMDfreq_WGBS.pdf", width=4, height=2.8, scale=1)

### statistics on the tile counts
d3 <- melt(lapply(l1, function(x) {
    x[,"pmd.freq"] <- tiles$pmd.freq ; 
    x[,"pmd.freq.bin"] <- tiles$pmd.freq.bin ; 
    x[,"width"] <- tile.size ; x}), 
        id.vars=c("Sub","Complex","Ins","Del", "pmd.freq", "pmd.freq.bin", "width"))

# Kruskal Wallis test
kruskal.test(d3$Sub, d3$pmd.freq)
#  
#          Kruskal-Wallis rank sum test
#  
#  data:  d3$Sub and d3$pmd.freq
#  Kruskal-Wallis chi-squared = 1338.8, df = 25, p-value < 2.2e-16
#  
kruskal.test(d3$Sub, d3$pmd.freq.bin)
#  
#          Kruskal-Wallis rank sum test
#  
#  data:  d3$Sub and d3$pmd.freq.bin
#  Kruskal-Wallis chi-squared = 1301.9, df = 8, p-value < 2.2e-16


## logistic regression
# Subst
summary(glm(formula=cbind(Sub, width-Sub) ~ pmd.freq, family = binomial(link = "probit"), data=d3))
#  
#  Call:
#  glm(formula = cbind(Sub, width - Sub) ~ pmd.freq, family = binomial(link = "probit"),
#      data = d3)
#  
#  Deviance Residuals:
#      Min       1Q   Median       3Q      Max
#  -0.4726  -0.4353  -0.4199  -0.4155  18.4003
#  
#  Coefficients:
#                Estimate Std. Error  z value Pr(>|z|)
#  (Intercept) -4.535e+00  6.332e-04 -7162.29   <2e-16 ***
#  pmd.freq     2.184e-03  6.809e-05    32.08   <2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  
#  (Dispersion parameter for binomial family taken to be 1)
#  
#      Null deviance: 1087435  on 2218399  degrees of freedom
#  Residual deviance: 1086425  on 2218398  degrees of freedom
#  AIC: 1441253
#  
#  Number of Fisher Scoring iterations: 6

summary(glm(formula=cbind(Sub, width-Sub) ~ pmd.freq==0, family = binomial(link = "probit"), data=d3))
#  
#  Call:
#  glm(formula = cbind(Sub, width - Sub) ~ pmd.freq == 0, family = binomial(link = "probit"),
#      data = d3)
#  
#  Deviance Residuals:
#      Min       1Q   Median       3Q      Max
#  -0.4355  -0.4355  -0.4355  -0.4155  18.4002
#  
#  Coefficients:
#                      Estimate Std. Error  z value Pr(>|z|)
#  (Intercept)       -4.5153027  0.0005712 -7904.54   <2e-16 ***
#  pmd.freq == 0TRUE -0.0198567  0.0009993   -19.87   <2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  
#  (Dispersion parameter for binomial family taken to be 1)
#  
#      Null deviance: 1087435  on 2218399  degrees of freedom
#  Residual deviance: 1087036  on 2218398  degrees of freedom
#  AIC: 1441864
#  
#  Number of Fisher Scoring iterations: 6


# Ins
summary(glm(formula=cbind(Ins, width-Ins) ~ pmd.freq, family = binomial(link = "probit"), data=d3))
#  
#  Call:
#  glm(formula = cbind(Ins, width - Ins) ~ pmd.freq, family = binomial(link = "probit"),
#      data = d3)
#  
#  Deviance Residuals:
#      Min       1Q   Median       3Q      Max
#  -0.0879  -0.0879  -0.0877  -0.0873   5.8262
#  
#  Coefficients:
#                Estimate Std. Error   z value Pr(>|z|)
#  (Intercept) -5.1520642  0.0026897 -1915.469   <2e-16 ***
#  pmd.freq    -0.0002772  0.0003059    -0.906    0.365
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  
#  (Dispersion parameter for binomial family taken to be 1)
#  
#      Null deviance: 95082  on 2218399  degrees of freedom
#  Residual deviance: 95081  on 2218398  degrees of freedom
#  AIC: 111843
#  
#  Number of Fisher Scoring iterations: 9

summary(glm(formula=cbind(Ins, width-Ins) ~ pmd.freq==0, family = binomial(link = "probit"), data=d3))
#  
#  Call:
#  glm(formula = cbind(Ins, width - Ins) ~ pmd.freq == 0, family = binomial(link = "probit"),
#      data = d3)
#  
#  Deviance Residuals:
#      Min       1Q   Median       3Q      Max
#  -0.0898  -0.0898  -0.0863  -0.0863   5.8439
#  
#  Coefficients:
#                     Estimate Std. Error   z value Pr(>|z|)
#  (Intercept)       -5.158821   0.002546 -2026.153  < 2e-16 ***
#  pmd.freq == 0TRUE  0.014537   0.004229     3.437 0.000588 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  
#  (Dispersion parameter for binomial family taken to be 1)
#  
#      Null deviance: 95082  on 2218399  degrees of freedom
#  Residual deviance: 95070  on 2218398  degrees of freedom
#  AIC: 111832
#  
#  Number of Fisher Scoring iterations: 9

# Del
summary(glm(formula=cbind(Del, width-Del) ~ pmd.freq, family = binomial(link = "probit"), data=d3))
#  
#  Call:
#  glm(formula = cbind(Del, width - Del) ~ pmd.freq, family = binomial(link = "probit"),
#      data = d3)
#  
#  Deviance Residuals:
#      Min       1Q   Median       3Q      Max
#  -0.1270  -0.1270  -0.1261  -0.1239   6.4613
#  
#  Coefficients:
#                Estimate Std. Error  z value Pr(>|z|)
#  (Intercept) -5.0123762  0.0019168 -2614.97  < 2e-16 ***
#  pmd.freq    -0.0009556  0.0002217    -4.31 1.63e-05 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  
#  (Dispersion parameter for binomial family taken to be 1)
#  
#      Null deviance: 171870  on 2218399  degrees of freedom
#  Residual deviance: 171851  on 2218398  degrees of freedom
#  AIC: 205034
#  
#  Number of Fisher Scoring iterations: 8
#  
summary(glm(formula=cbind(Del, width-Del) ~ pmd.freq==0, family = binomial(link = "probit"), data=d3))
#  
#  Call:
#  glm(formula = cbind(Del, width - Del) ~ pmd.freq == 0, family = binomial(link = "probit"),
#      data = d3)
#  
#  Deviance Residuals:
#      Min       1Q   Median       3Q      Max
#  -0.1292  -0.1292  -0.1230  -0.1230   6.4942
#  
#  Coefficients:
#                     Estimate Std. Error   z value Pr(>|z|)
#  (Intercept)       -5.024520   0.001831 -2743.834  < 2e-16 ***
#  pmd.freq == 0TRUE  0.018744   0.003024     6.198 5.71e-10 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  
#  (Dispersion parameter for binomial family taken to be 1)
#  
#      Null deviance: 171870  on 2218399  degrees of freedom
#  Residual deviance: 171832  on 2218398  degrees of freedom
#  AIC: 205015
#  
#  Number of Fisher Scoring iterations: 8
#  

