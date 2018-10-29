library(reshape2)
library(ggplot2)
library(gplots)
library(parallel)
library(RColorBrewer)

# get all files with feature counts in genomic windows
files <- c(
  list.files(path="~/BiSeq_BASIS/PMDs/PMD_characteristics/CTCF", pattern="readcounts", full.names=T),
  list.files(path="~/BiSeq_BASIS/PMDs/PMD_characteristics/LADs", pattern="counts", full.names=T),
  list.files(path="~/BiSeq_BASIS/PMDs/PMD_characteristics/TDs", pattern="counts", full.names=T), 
  list.files(path="~/BiSeq_BASIS/PMDs/PMD_characteristics/repliseq", pattern="counts", full.names=T) 
)
files <- as.list(files)
names(files) <- gsub("Uw","", gsub("\\.counts","",gsub("list","s", gsub("HiCCUPS_","", 
    gsub("Arrowhead_","",gsub("GSE63525_","", gsub("UwTfbs","",gsub("HaibTfbs","", 
    gsub("wgEncode","", gsub("\\.readcounts","", 
    gsub("hg19.Win1000.bed__","",basename(unlist(files)), ignore.case=T)))))))))))

readFiles <- function(FILE) {
  read.delim(pipe(paste0("cat ", FILE, " | grep -v chrM")), header=F)
}
l1 <- mclapply(files, readFiles, mc.cores=8)

# get the genomic window positions
# distance of genomic windows to PMDs was calculated using:
# bedtools closest -D "b" -a <windows bedfile> -b <PMDs bedfile> 
pos <- readFiles("~/BiSeq_BASIS/PMDs/PMD_characteristics/hg19.Win1000__PMDs_merged.bed.distanceToClosest")
names(pos) <- c("wchr","wstart","wend","dist","pos","pchr","pstart","pend")

# do distance binning (distance to PMD) within a selected range
WIN <- 1000 
UPDOWN <- 50000 # only windows that are within 50 kb of a PMD
pos$dbin <- cut(pos$dist, breaks=seq(-UPDOWN, UPDOWN, by=WIN), labels=F)

# select the distance bins & counts that are within the selected range
l2 <- mclapply(l1, function(DF) {df<- DF ; df$dbin <- pos$dbin ; df[!is.na(df[,5]) & !is.na(df[,4]),4:5] }, 
    mc.cores=8)
l3 <- mclapply(l2, function(DF) {sapply(1:max(pos$dbin, na.rm=T), function(X) {mean(DF$V4[which(DF$dbin==X)])})}, 
    mc.cores=8)
distances <- sapply(1:max(pos$dbin, na.rm=T), function(X) {mean(pos$dist[!is.na(pos$dbin) & pos$dbin==X])})

# clean up/format names of the tracks
cleanNames <- function(DF) {
  str <- names(l3)
  for (i in 1:nrow(DF)) {
    str <- gsub(as.character(DF[i,1]), as.character(DF[i,2]), str, ignore.case=T)
  }
  str
}
names(l3) <- cleanNames(rbind(
  c("Ctcf"," CTCF "),
  c("Mcf7", "MCF7"),
  c("Imr90","IMR90"),
  c("Hmec","HMEC"),
  c("RepliSeq","RepliSeq "),
  c("MCF7S","MCF7 S"),
  c("MCF7G","MCF7 G"),
  c("imr90S","IMR90 S"),
  c("imr90G","IMR90 G"),
  c("Aln",""),
  c("StdAln",""),
  c("_"," "),
  c("Rep"," Rep"),
  c(" RepliSeq","repliSeq"),
  c("WaveSignal"," WaveSignal"),
  c(" DamID hg19","")
  )
)

l4 <- lapply(names(l3), function(X) {data.frame("dist"=distances, "value"=scale(l3[[X]]), 
    "assay"=rep(X, length(distances)))})
names(l4) <- as.character(unlist(lapply(l4, function(X) {X$assay[1]})))

# lineplot for each feature (to be formatted into a singe figure) #### Figure 2D ####
plotLines <- function(TERM1, TERM2, TERM3) {
  c1 <- grepl(TERM1, names(l4), ignore.case=T)
  c2 <- grepl(TERM2, names(l4), ignore.case=T)
  c3 <- grepl(TERM3, names(l4), ignore.case=T)
  sel <- c1 & c2 & c3
  df <- do.call("rbind", l4[sel])
  p <- ggplot(data=df, aes(dist/1000, value))
  p <- p + geom_line(aes(color=assay))
  colors <- c("blue","red","darkgreen", "purple","black","orange")
  p <- p + scale_color_manual(values=colors[1:length(levels(df$assay))])
  
  p <- p + theme_classic()
  p <- p + theme(panel.border=element_rect(fill=NA, size=0.8))
  p <- p + ylab("signal\n(z-score)")
  p <- p + xlab("distance to PMD border (kb)")
  ggsave(p, file=paste0("lineplot_PMDs_", TERM1, "_", TERM2, ".pdf"), scale=0.3, height=12, width=18)
}
plotLines("CTCF","rep1","std")
plotLines("RepliSeq","IMR90", "rep1")
plotLines("RepliSeq","IMR90 WaveSignal", "rep1")
plotLines("RepliSeq","MCF7", "rep1")
plotLines("RepliSeq","MCF7 WaveSignal", "rep1")
plotLines("domains","domains", "domains")
plotLines("loops","w5000", "w5000")
plotLines("Tig3","LaminB1", "LaminB1")


