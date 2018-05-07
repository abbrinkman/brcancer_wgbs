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

# get all patients
patients <- names(vcf.all.pat)
patients <- as.list(patients)
names(patients) <- unlist(patients)

# take data only for relevant patients
mut <- vcf.all.pat[names(patients)]

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
tiles$pmd.freq.bin <- factor(tiles$pmd.freq.bin, 
  levels=unique(tiles$pmd.freq.bin)[order(as.numeric(gsub(",.+$","",gsub("\\(","", unique(tiles$pmd.freq.bin)))))])
tiles$length <- width(tiles)


# count for each patient the mutations in each tile
countMut <- function(pd, type) {
  message(pd)
  countOverlaps(tiles, mut[[pd]][mut[[pd]]$Type==type])
}
l1 <- lapply(patients, function(x) {
    as.data.frame(sapply(c("Sub", "Complex", "Del", "Ins"), function(y) {
    message("    ", y)
    countMut(x, y)}, USE.NAMES=T, simplify=T))})

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
d2$pmd.freq.bin <- factor(d2$pmd.freq, levels=levels(tiles$pmd.freq.bin))

p.sub <- ggplot(d2[d2$Type=="Sub",], aes(pmd.freq.bin, mut.Mb)) + 
         geom_boxplot(aes(fill=pmd.freq.bin), outlier.shape=NA) + 
         scale_fill_manual(values=colorpanel(length(unique(d2$pmd.freq.bin)), "white", "red")) + 
         coord_cartesian(ylim=c(0,8)) +
         theme_classic() +
         theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"),
              axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
         ylab("mean mut/Mb") +
         xlab("PMD frequency") +
         ggtitle("Subst, all samples")
ggsave(p.sub, file="boxplots_Sub_perMb_vs_PMDfreq_all.pdf", width=4, height=2.8, scale=1) #### Figure 2G ####

p.ins <- ggplot(d2[d2$Type=="Ins",], aes(pmd.freq.bin, mut.Mb)) +
         geom_boxplot(aes(fill=pmd.freq.bin), outlier.shape=NA) +
         scale_fill_manual(values=colorpanel(length(unique(d2$pmd.freq.bin)), "white", "red")) +
         coord_cartesian(ylim=c(0,0.1)) +
         theme_classic() +
         theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"),
              axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
         ylab("mean mut/Mb") +
         xlab("PMD frequency") +
         ggtitle("Ins, all samples")
ggsave(p.ins, file="boxplots_Ins_perMb_vs_PMDfreq_all.pdf", width=4, height=2.8, scale=1) #### Figure 2G ####

p.del <- ggplot(d2[d2$Type=="Del",], aes(pmd.freq.bin, mut.Mb)) +
         geom_boxplot(aes(fill=pmd.freq.bin), outlier.shape=NA) +
         scale_fill_manual(values=colorpanel(length(unique(d2$pmd.freq.bin)), "white", "red")) +
         coord_cartesian(ylim=c(0,0.15)) +
         theme_classic() +
         theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"),
              axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
         ylab("mean mut/Mb") +
         xlab("PMD frequency") +
         ggtitle("Del, all samples")
ggsave(p.del, file="boxplots_Del_perMb_vs_PMDfreq_all.pdf", width=4, height=2.8, scale=1) #### Figure 2G ####

p.complex <- ggplot(d2[d2$Type=="Complex",], aes(pmd.freq.bin, mut.Mb)) +
         geom_boxplot(aes(fill=pmd.freq.bin), outlier.shape=NA) +
         scale_fill_manual(values=colorpanel(length(unique(d2$pmd.freq.bin)), "white", "red")) +
         coord_cartesian(ylim=c(0,0.015)) +
         theme_classic() +
         theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"),
              axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
         ylab("mean mut/Mb") +
         xlab("PMD frequency") +
         ggtitle("Complex, all samples")
ggsave(p.complex, file="boxplots_Complex_perMb_vs_PMDfreq_all.pdf", width=4, height=2.8, scale=1)

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
#  Kruskal-Wallis chi-squared = 81575, df = 30, p-value < 2.2e-16
#  
kruskal.test(d3$Sub, d3$pmd.freq.bin)
#  
#          Kruskal-Wallis rank sum test
#  
#  data:  d3$Sub and d3$pmd.freq.bin
#  Kruskal-Wallis chi-squared = 79569, df = 10, p-value < 2.2e-16


#################################################################################
############################# logistic regression ###############################
#################################################################################

################################### whole set ###################################

s <- sample(nrow(d3), 10000000) # to prevent a too lengthy computation time, pick a random sample

summary(glm(formula=cbind(Sub, width-Sub) ~ pmd.freq, family = binomial(link = "probit"),
     data=d3[s,]))
#  
#  Call:
#  glm(formula = cbind(Sub, width - Sub) ~ pmd.freq, family = binomial(link = "probit"),
#      data = d3[s, ])
#  
#  Deviance Residuals:
#      Min       1Q   Median       3Q      Max
#  -0.4364  -0.3698  -0.3446  -0.3356  23.7439
#  
#  Coefficients:
#                Estimate Std. Error z value Pr(>|z|)
#  (Intercept) -4.625e+00  3.552e-04  -13019   <2e-16 ***
#  pmd.freq     3.668e-03  2.957e-05     124   <2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  
#  (Dispersion parameter for binomial family taken to be 1)
#  
#      Null deviance: 3791751  on 9999999  degrees of freedom
#  Residual deviance: 3776950  on 9999998  degrees of freedom
#  AIC: 4947100
#  
#  Number of Fisher Scoring iterations: 6
#  
summary(glm(formula=cbind(Sub, width-Sub) ~ pmd.freq==0, family = binomial(link = "probit"),
     data=d3[s,]))
#  
#  Call:
#  glm(formula = cbind(Sub, width - Sub) ~ pmd.freq == 0, family = binomial(link = "probit"),
#      data = d3[s, ])
#  
#  Deviance Residuals:
#      Min       1Q   Median       3Q      Max
#  -0.3746  -0.3746  -0.3746  -0.3271  23.5994
#  
#  Coefficients:
#                      Estimate Std. Error   z value Pr(>|z|)
#  (Intercept)       -4.5787607  0.0003082 -14854.05   <2e-16 ***
#  pmd.freq == 0TRUE -0.0564387  0.0005726    -98.57   <2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  
#  (Dispersion parameter for binomial family taken to be 1)
#  
#      Null deviance: 3791751  on 9999999  degrees of freedom
#  Residual deviance: 3781722  on 9999998  degrees of freedom
#  AIC: 4951873
#  
#  Number of Fisher Scoring iterations: 6
#  
#  
summary(glm(formula=cbind(Ins, width-Ins) ~ pmd.freq, family = binomial(link = "probit"),
     data=d3[s,]))
#  
#  Call:
#  glm(formula = cbind(Ins, width - Ins) ~ pmd.freq, family = binomial(link = "probit"),
#      data = d3[s, ])
#  
#  Deviance Residuals:
#      Min       1Q   Median       3Q      Max
#  -0.0669  -0.0606  -0.0576  -0.0566   7.3875
#  
#  Coefficients:
#                Estimate Std. Error  z value Pr(>|z|)
#  (Intercept) -5.3145199  0.0018634 -2852.09   <2e-16 ***
#  pmd.freq     0.0020373  0.0001607    12.67   <2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  
#  (Dispersion parameter for binomial family taken to be 1)
#  
#      Null deviance: 222713  on 9999999  degrees of freedom
#  Residual deviance: 222557  on 9999998  degrees of freedom
#  AIC: 256746
#  
#  Number of Fisher Scoring iterations: 10
#  
summary(glm(formula=cbind(Ins, width-Ins) ~ pmd.freq==0, family = binomial(link = "probit"),
     data=d3[s,]))
#  
#  Call:
#  glm(formula = cbind(Ins, width - Ins) ~ pmd.freq == 0, family = binomial(link = "probit"),
#      data = d3[s, ])
#  
#  Deviance Residuals:
#      Min       1Q   Median       3Q      Max
#  -0.0607  -0.0607  -0.0607  -0.0557   7.4050
#  
#  Coefficients:
#                     Estimate Std. Error  z value Pr(>|z|)
#  (Intercept)       -5.289330   0.001664 -3179.29   <2e-16 ***
#  pmd.freq == 0TRUE -0.031075   0.002995   -10.38   <2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  
#  (Dispersion parameter for binomial family taken to be 1)
#  
#      Null deviance: 222713  on 9999999  degrees of freedom
#  Residual deviance: 222604  on 9999998  degrees of freedom
#  AIC: 256793
#  
#  Number of Fisher Scoring iterations: 10
#  
#  
summary(glm(formula=cbind(Del, width-Del) ~ pmd.freq, family = binomial(link = "probit"),
     data=d3[s,]))
#  
#  Call:
#  glm(formula = cbind(Del, width - Del) ~ pmd.freq, family = binomial(link = "probit"),
#      data = d3[s, ])
#  
#  Deviance Residuals:
#      Min       1Q   Median       3Q      Max
#  -0.1101  -0.1028  -0.0993  -0.0982   8.5753
#  
#  Coefficients:
#                Estimate Std. Error  z value Pr(>|z|)
#  (Intercept) -5.111e+00  1.120e-03 -4562.49   <2e-16 ***
#  pmd.freq     1.456e-03  9.868e-05    14.75   <2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  
#  (Dispersion parameter for binomial family taken to be 1)
#  
#      Null deviance: 557352  on 9999999  degrees of freedom
#  Residual deviance: 557139  on 9999998  degrees of freedom
#  AIC: 649803
#  
#  Number of Fisher Scoring iterations: 9
#  
summary(glm(formula=cbind(Del, width-Del) ~ pmd.freq==0, family = binomial(link = "probit"),
     data=d3[s,]))
#  
#  Call:
#  glm(formula = cbind(Del, width - Del) ~ pmd.freq == 0, family = binomial(link = "probit"),
#      data = d3[s, ])
#  
#  Deviance Residuals:
#      Min       1Q   Median       3Q      Max
#  -0.1032  -0.1032  -0.1032  -0.0965   8.5991
#  
#  Coefficients:
#                     Estimate Std. Error  z value Pr(>|z|)
#  (Intercept)       -5.091611   0.001014 -5022.07   <2e-16 ***
#  pmd.freq == 0TRUE -0.025364   0.001804   -14.06   <2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  
#  (Dispersion parameter for binomial family taken to be 1)
#  
#      Null deviance: 557352  on 9999999  degrees of freedom
#  Residual deviance: 557151  on 9999998  degrees of freedom
#  AIC: 649816
#  
#  Number of Fisher Scoring iterations: 9


################################ WGBS set ##################################
d4 <- d3[d3$L1 %in% names(pmds),]

summary(glm(formula=cbind(Sub, width-Sub) ~ pmd.freq, family = binomial(link = "probit"),
     data=d4))
#  
#  Call:
#  glm(formula = cbind(Sub, width - Sub) ~ pmd.freq, family = binomial(link = "probit"),
#      data = d4)
#  
#  Deviance Residuals:
#      Min       1Q   Median       3Q      Max
#  -0.4716  -0.4354  -0.4192  -0.4157  18.3991
#  
#  Coefficients:
#                Estimate Std. Error  z value Pr(>|z|)
#  (Intercept) -4.535e+00  6.271e-04 -7231.73   <2e-16 ***
#  pmd.freq     1.787e-03  5.503e-05    32.48   <2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  
#  (Dispersion parameter for binomial family taken to be 1)
#  
#      Null deviance: 1087435  on 2218399  degrees of freedom
#  Residual deviance: 1086401  on 2218398  degrees of freedom
#  AIC: 1441228
#  
#  Number of Fisher Scoring iterations: 6
  
summary(glm(formula=cbind(Sub, width-Sub) ~ pmd.freq==0, family = binomial(link = "probit"),
     data=d4))
#  
#  Call:
#  glm(formula = cbind(Sub, width - Sub) ~ pmd.freq == 0, family = binomial(link = "probit"),
#      data = d4)
#  
#  Deviance Residuals:
#      Min       1Q   Median       3Q      Max
#  -0.4354  -0.4354  -0.4354  -0.4156  18.4000
#  
#  Coefficients:
#                      Estimate Std. Error  z value Pr(>|z|)
#  (Intercept)       -4.5153878  0.0005702 -7918.55   <2e-16 ***
#  pmd.freq == 0TRUE -0.0197470  0.0010012   -19.72   <2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  
#  (Dispersion parameter for binomial family taken to be 1)
#  
#      Null deviance: 1087435  on 2218399  degrees of freedom
#  Residual deviance: 1087042  on 2218398  degrees of freedom
#  AIC: 1441870
#  
#  Number of Fisher Scoring iterations: 6
  
summary(glm(formula=cbind(Ins, width-Ins) ~ pmd.freq, family = binomial(link = "probit"),
     data=d4))
#  
#  Call:
#  glm(formula = cbind(Ins, width - Ins) ~ pmd.freq, family = binomial(link = "probit"),
#      data = d4)
#  
#  Deviance Residuals:
#      Min       1Q   Median       3Q      Max
#  -0.0879  -0.0879  -0.0877  -0.0873   5.8261
#  
#  Coefficients:
#                Estimate Std. Error   z value Pr(>|z|)
#  (Intercept) -5.1520914  0.0026637 -1934.209   <2e-16 ***
#  pmd.freq    -0.0002257  0.0002475    -0.912    0.362
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
  
summary(glm(formula=cbind(Ins, width-Ins) ~ pmd.freq==0, family = binomial(link = "probit"),
     data=d4))
#  
#  Call:
#  glm(formula = cbind(Ins, width - Ins) ~ pmd.freq == 0, family = binomial(link = "probit"),
#      data = d4)
#  
#  Deviance Residuals:
#      Min       1Q   Median       3Q      Max
#  -0.0899  -0.0899  -0.0863  -0.0863   5.8446
#  
#  Coefficients:
#                     Estimate Std. Error   z value Pr(>|z|)
#  (Intercept)       -5.159061   0.002543 -2028.957   <2e-16 ***
#  pmd.freq == 0TRUE  0.015307   0.004234     3.615    3e-04 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  
#  (Dispersion parameter for binomial family taken to be 1)
#  
#      Null deviance: 95082  on 2218399  degrees of freedom
#  Residual deviance: 95069  on 2218398  degrees of freedom
#  AIC: 111831
#  
#  Number of Fisher Scoring iterations: 9
  
summary(glm(formula=cbind(Del, width-Del) ~ pmd.freq, family = binomial(link = "probit"),
     data=d4))
#  
#  Call:
#  glm(formula = cbind(Del, width - Del) ~ pmd.freq, family = binomial(link = "probit"),
#      data = d4)
#  
#  Deviance Residuals:
#      Min       1Q   Median       3Q      Max
#  -0.1270  -0.1270  -0.1259  -0.1239   6.4606
#  
#  Coefficients:
#                Estimate Std. Error   z value Pr(>|z|)
#  (Intercept) -5.0125125  0.0018984 -2640.352  < 2e-16 ***
#  pmd.freq    -0.0007717  0.0001795    -4.299 1.71e-05 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  
#  (Dispersion parameter for binomial family taken to be 1)
#  
#      Null deviance: 171870  on 2218399  degrees of freedom
#  Residual deviance: 171851  on 2218398  degrees of freedom
#  AIC: 205035
#  
#  Number of Fisher Scoring iterations: 8
#  
summary(glm(formula=cbind(Del, width-Del) ~ pmd.freq==0, family = binomial(link = "probit"),
     data=d4))
#  
#  Call:
#  glm(formula = cbind(Del, width - Del) ~ pmd.freq == 0, family = binomial(link = "probit"),
#      data = d4)
#  
#  Deviance Residuals:
#      Min       1Q   Median       3Q      Max
#  -0.1293  -0.1293  -0.1230  -0.1230   6.4944
#  
#  Coefficients:
#                     Estimate Std. Error  z value Pr(>|z|)
#  (Intercept)       -5.024583   0.001828 -2748.77  < 2e-16 ***
#  pmd.freq == 0TRUE  0.019047   0.003028     6.29 3.17e-10 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  
#  (Dispersion parameter for binomial family taken to be 1)
#  
#      Null deviance: 171870  on 2218399  degrees of freedom
#  Residual deviance: 171831  on 2218398  degrees of freedom
#  AIC: 205014
#  
#  Number of Fisher Scoring iterations: 8
