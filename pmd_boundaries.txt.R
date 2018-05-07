library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)

pmd.gr <- get(load("~/BiSeq_BASIS/PMDs_all_samples/pmdmeth.basis.RData"))

# make genome tiles
tile.size <- 30000
chr.lengths <- seqlengths(Hsapiens)
chr.lengths <- chr.lengths[1:22]

tiles <- tileGenome(chr.lengths, tilewidth=tile.size, 
        cut.last.tile.in.chrom=T)

# randomize the PMDs, shuffle PMDs within each chromosome,
# but leave chromosomes untouched
randomizePmds <- function(GR) {
  d <- data.frame("start.min"=1+ceiling((width(GR)*0.5)),
      "end.max"=seqlengths(Hsapiens)[as.character(seqnames(GR))]-ceiling((width(GR)*0.5)))
  mids.ran <- apply(d, 1, function(x) {round(runif(1, x[1], x[2]),0)})
  starts.ran <- mids.ran-floor(width(GR)*0.5)
  ends.ran <- mids.ran+floor(width(GR)*0.5)
  gr.ran <- GRanges(seqnames=seqnames(GR), IRanges(start=starts.ran, end=ends.ran))
  gr.ran
}
pmd.gr.ran <- lapply(pmd.gr, randomizePmds)

# binary calls for no start/end (0) or at least one start/end (1)
binBounds <- function(GR) {ifelse(
        countOverlaps(tiles, flank(GR, width=1,start=TRUE)) >= 1 | 
        countOverlaps(tiles, flank(GR, width=1,start=FALSE)) >= 1, 1, 0)}

# binary calls for no overlap with at a PMD (0) or FULL overlap with a PMD (1)
binPmds <- function(GR) {ifelse(countOverlaps(tiles, GR, type="within") >= 1, 1, 0)} 

bound.count <- sapply(pmd.gr, binBounds, USE.NAMES=T)
bound.count.ran <- sapply(pmd.gr.ran, binBounds, USE.NAMES=T)

bound.sums <- apply(bound.count, 1, sum)
bound.sums.ran<- apply(bound.count.ran, 1, sum)

# save the bound.count and pmd.count tables for clustering on PMDs
bound.count.gr <- tiles
mcols(bound.count.gr) <- bound.count

# plot histograms
freqs <- 0:30
bound.freq <- sapply(freqs, function(X) {sum(bound.sums==X)})
bound.freq.ran <- sapply(freqs, function(X) {sum(bound.sums.ran==X)})
m.b <- rbind(bound.freq, bound.freq.ran)

# plot histograms of PMD boundary tiles
pdf(paste0("histograms_PMD_boundaries-tiles", tile.size, ".pdf")) #### Figure 2C ####
par(mfrow=c(2,1), mar=c(5,5,4,2))
barplot(m.b/1000, beside=T, col=c("red","black"), border=F, names.arg=freqs,
    xlab="frequency (#cases)",
    ylab=paste0("#genomic tiles (", tile.size/1000," kb)\nwith PMD boundaries (x1000)"), las=2)
legend("right", bty="n", legend=c("PMDs", "shuffled PMDs"),
    fill=c("red","black")) 
legend("topright", bty="n", legend="PMD boundaries")

dev.off()

