library(GenomicRanges)
library(reshape2)
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg19)
library(gplots)
library(openxlsx)

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

# get the rearrangment data
rearr <- read.xlsx(
    "/home/arjen/Infinium/BASIS_data_all_TN_QC_analysis/SV_probabilities_corrected.xlsx")

# reorganize into individual breakpoints so that they can be counted
rearr.1 <- rearr[,c(1,2,5:ncol(rearr))]
colnames(rearr.1) <- gsub("^pos\\.\\d$","start",gsub("^Chromosome\\.\\d$","chr", colnames(rearr.1)))
rearr.1$end <- rearr.1$start+1
rearr.2 <- rearr[,c(3,4,5:ncol(rearr))]
colnames(rearr.2) <- gsub("^pos\\.\\d$","start",gsub("^Chromosome\\.\\d$","chr", colnames(rearr.2)))
rearr.2$end <- rearr.2$start+1
rearr.1$chr <- paste0("chr", rearr.1$chr) 
rearr.2$chr <- paste0("chr", rearr.2$chr) 
rearr.all <- makeGRangesFromDataFrame(rbind(rearr.1, rearr.2), keep.extra.columns=T)
rearr.wgbs <- rearr.all[rearr$sample %in% names(pmds)] 

# split into one gr per patient
rearr.all.pat <- split(rearr.all, rearr.all$sample) # all WGS patients

patients <- intersect(names(pmds), names(rearr.all.pat)) # only the WGBS patients

# tile the genome
tile.size <- 30000
tiles <- unlist(tile(hg19, width=tile.size))

# determine the PMD frequency for each tile
tiles$pmd.freq <- apply(sapply(pmds, function(x) {ifelse(overlapsAny(tiles, x, type="within"), 1, 0)}), 1, sum)
tiles$pmd.freq.bin <- cut(tiles$pmd.freq, breaks=c(-0.01, seq(0, length(pmds), by=3)))
tiles$pmd.freq.wgbs <- apply(sapply(pmds[patients], function(x) {
    ifelse(overlapsAny(tiles, x, type="within"), 1, 0)} ), 1, sum)
tiles$pmd.freq.bin.wgbs <- cut(tiles$pmd.freq, breaks=c(-0.01, seq(0, length(patients), by=3)))

tiles$length <- width(tiles)



# count for each patient the breakpoints in each tile
l1 <- sapply(names(rearr.all.pat), function(x) {message(x) ;
    countOverlaps(tiles, rearr.all.pat[[x]])}, 
    USE.NAMES=T, simplify=F)

l2 <- sapply(patients, function(x) {message(x) ;
    countOverlaps(tiles, rearr.all.pat[[x]])},
    USE.NAMES=T, simplify=F)


# sum the breakpoints per tile
d1 <- as.data.frame(sapply(names(rearr.all.pat), function(x) {
    tapply(l1[[x]], tiles$pmd.freq.bin, sum)}))
d1$pmd.freq.bin <- rownames(d1)
d1$Mb <- sapply(d1$pmd.freq.bin, function(x) {
    sum(width(tiles[which(tiles$pmd.freq.bin==x)]))})/1e06

d2 <- as.data.frame(sapply(patients, function(x) {
    tapply(l2[[x]], tiles$pmd.freq.bin.wgbs, sum)}))
d2$pmd.freq.bin <- rownames(d2)
d2$Mb <- sapply(d2$pmd.freq.bin, function(x) {
    sum(width(tiles[which(tiles$pmd.freq.bin.wgbs==x)]))})/1e06

# convert to breakpoints/Mb
for (i in colnames(d1)[grep("^PD", colnames(d1))]) { d1[,i] <- d1[,i]/d1$Mb}
for (i in colnames(d2)[grep("^PD", colnames(d2))]) { d2[,i] <- d2[,i]/d2$Mb}

d3 <- melt(d1, variable.name="patient", value.name="breakpoints.Mb")
d3$pmd.freq.bin <- factor(d3$pmd.freq.bin, levels=unique(d3$pmd.freq.bin))

d4 <- melt(d2, variable.name="patient", value.name="breakpoints.Mb")
d4$pmd.freq.bin <- factor(d4$pmd.freq, levels=unique(d4$pmd.freq.bin))



p1 <- ggplot(d3, aes(pmd.freq.bin, breakpoints.Mb)) + 
         geom_boxplot(aes(fill=pmd.freq.bin), outlier.shape=NA) + 
         scale_fill_manual(values=colorpanel(length(unique(d3$pmd.freq.bin)), "white", "red")) + 
         coord_cartesian(ylim=c(0,0.45)) +
         theme_classic() +
         theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"),
              axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
         ylab("mean breakpoints/Mb") +
         xlab("PMD frequency") +
         ggtitle("all (544) samples")
ggsave(p1, file="boxplots_breakpoints_perMb_vs_PMDfreq_all.pdf", width=4, height=2.8, scale=1) #### Figure 2G ####


p2 <- ggplot(d4, aes(pmd.freq.bin, breakpoints.Mb)) +
         geom_boxplot(aes(fill=pmd.freq.bin), outlier.shape=NA) +
         scale_fill_manual(values=colorpanel(length(unique(d4$pmd.freq.bin)), "white", "red")) +
         coord_cartesian(ylim=c(0,0.45)) +
         theme_classic() +
         theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"),
              axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
         ylab("mean breakpoints/Mb") +
         xlab("PMD frequency") +
         ggtitle("24 WGBS samples")
ggsave(p2, file="boxplots_breakpoints_perMb_vs_PMDfreq_WGBS.pdf", width=4, height=2.8, scale=1) #### Supplemental Figure 3C ####


### statistics on the tile counts
# full cohort
d5 <- as.data.frame(sapply(names(rearr.all.pat), function(x) {
    tapply(l1[[x]], tiles$pmd.freq, sum)}))
d5$pmd.freq <- rownames(d5)
d5$width <- sapply(d5$pmd.freq, function(x) {
    sum(width(tiles[which(tiles$pmd.freq==x)]))})
d5 <- melt(d5, id.vars=c("pmd.freq","width"))
d5$pmd.freq <- as.numeric(d5$pmd.freq)

# WGBS cohort
d6 <- as.data.frame(sapply(patients, function(x) {
    tapply(l2[[x]], tiles$pmd.freq.wgbs, sum)}))
d6$pmd.freq <- rownames(d6)
d6$width <- sapply(d6$pmd.freq, function(x) {
    sum(width(tiles[which(tiles$pmd.freq.wgbs==x)]))})
d6 <- melt(d6, id.vars=c("pmd.freq","width"))
d6$pmd.freq <- as.numeric(d6$pmd.freq)

## logistic regression
# full cohort
summary(glm(formula=cbind(value, width-value) ~ pmd.freq, family = binomial(link = "probit"), 
    data=d5))
#  Call:
#  glm(formula = cbind(value, width - value) ~ pmd.freq, family = binomial(link = "probit"),
#      data = d5)
#  
#  Deviance Residuals:
#      Min       1Q   Median       3Q      Max
#  -14.448   -2.428   -1.037    0.865   47.141
#  
#  Coefficients:
#                Estimate Std. Error  z value Pr(>|z|)
#  (Intercept) -5.175e+00  6.139e-04 -8428.71   <2e-16 ***
#  pmd.freq    -2.972e-03  6.255e-05   -47.51   <2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  
#  (Dispersion parameter for binomial family taken to be 1)
#  
#      Null deviance: 181709  on 16863  degrees of freedom
#  Residual deviance: 179340  on 16862  degrees of freedom
#  AIC: 220819
#  
#  Number of Fisher Scoring iterations: 5

#WGBS cohort
summary(glm(formula=cbind(value, width-value) ~ pmd.freq, family = binomial(link = "probit"),
    data=d6))
#  Call:
#  glm(formula = cbind(value, width - value) ~ pmd.freq, family = binomial(link = "probit"),
#      data = d6)
#  
#  Deviance Residuals:
#       Min        1Q    Median        3Q       Max
#  -12.2696   -2.9546   -1.7282    0.4168   29.9287
#  
#  Coefficients:
#                Estimate Std. Error   z value Pr(>|z|)
#  (Intercept) -5.2239450  0.0033362 -1565.842  < 2e-16 ***
#  pmd.freq    -0.0025237  0.0004138    -6.099 1.07e-09 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  
#  (Dispersion parameter for binomial family taken to be 1)
#  
#      Null deviance: 8832.1  on 599  degrees of freedom
#  Residual deviance: 8793.8  on 598  degrees of freedom
#  AIC: 10109
#  
#  Number of Fisher Scoring iterations: 5


