library(Gviz)
library(RColorBrewer)
library(gplots)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(data.table)

# get the WGBS data (GRanges object, with CpG positions as rows, and
# colnames PDxxxx = methylation value, only for CpGs with coverage >= 4)
load("~/BiSeq_BASIS/methcounts/meth_cov4.gr.RData")
hc <- get(load("tilemaps_PMDmeth_clustered_hclust.RData"))

makeScreenshot <- function(CHR, START, END) {
  gen <- 'hg19'
  
  # gene annotation
  genetrack <- UcscTrack(track="NCBI RefSeq", table="refGene",
      trackType="GeneRegionTrack",chromosome=CHR, genome=gen,
      rstart="exonStarts", rends="exonEnds", gene="name", symbol="name2",
      transcript="name", strand="strand", name="RefSeq Genes",
      feature="name2", showId=T, from=START, to=END,
      fill="darkblue", col.line="darkblue", col="darkblue", utr3="darkblue", utr5="darkblue", 
      protein_coding="darkblue", fontcolor.group="darkblue",
      background.title="white", fontcolor.title="black")

  # CpG islands
  cgi <- read.table("~/Datarepository/CpGislands_hg19_v20140606.full.table")[,1:3]
  colnames(cgi) <- c("chr","start","end")
  cgi <- makeGRangesFromDataFrame(cgi)
  cgitrack <- AnnotationTrack(cgi, name="CGI", start=START, end=END, chromosome=CHR, 
  	fill="darkgreen", col.line="darkgreen", col="white", stacking="dense",
  	background.title="white", fontcolor.title="black", rotation.title=0)

  # PMDs
  pmd <- makeGRangesFromDataFrame(read.table("~/BiSeq_BASIS/PMDs_all_samples/PMDs_10in30patients_noCent_merged.bed",
      col.names=c("chr","start","end")))
  pmd <- subsetByOverlaps(pmd, GRanges(seqnames=CHR, IRanges(start=START, end=END)))
  pmdtrack <- AnnotationTrack(pmd[,0], name="PMDs", start=START, end=END, chromosome=CHR,
        fill="pink", col.line="black", col="white", stacking="dense",
        background.title="white", fontcolor.title="black", rotation.title=0)

  # methylation
  meth.sel <- subsetByOverlaps(meth.gr,
        GRanges(seqnames=CHR, ranges=IRanges(start=START, end=END)))
  meth.sel <- meth.sel[,!colnames(mcols(meth.sel)) %in% c("HMEC","MCF7")]

  # use running median to smooth values
  runmed.k <- round((END-START)/5500, 0)
  runmed.k <- ifelse(runmed.k %% 2 == 0, runmed.k+1, runmed.k)
  smoothMethTrack <- function(NAME, BACK_COL) {
    meth.smooth <- meth.sel[,NAME]
    meth.smooth <- meth.smooth[!is.na(mcols(meth.smooth)[[NAME]])]
    mcols(meth.smooth)[[NAME]] <- runmed(mcols(meth.smooth)[[NAME]], k=runmed.k)
    DataTrack(range = meth.smooth, genome = gen, type = "l",
        name = NAME,
        chromosome = CHR, start=START, end=END,
        background.title="white", fontcolor.title="black", cex.title=0.5, col.line="darkblue", lwd=1.7, 
        background.panel=BACK_COL, showAxis=T, showTitle=T, rot.title=0)
   }
  # set background colors
  backgr.cols <- rep(c("lightgrey", "transparent"), length(colnames(mcols(meth.sel))))

  # use the hierarchical clustering order (done previously for the whole-genome/chromosome visualizations)
  methorder <- hc$labels[hc$order]
  methorder <- rev(methorder[!methorder %in% c("HMEC","MCF7","PD9590a")])

  methtracklist <- sapply(1:length(methorder), function(x) {
      smoothMethTrack(methorder[x], backgr.cols[x])}, USE.NAMES=T, simplify=F)

  # combine all the items to plot into one list
  plotlist <- c(list(genetrack, cgitrack,pmdtrack), methtracklist)
  
  pdfname <- paste0("lineplots_PMDmeth_clustered_", CHR, "_", START, "_", END, ".pdf" )
  pdf(pdfname, width=5, height=8)
  plotTracks(trackList=plotlist, 
      sizes=c(0.08, 0.03, 0.03, rep(0.6/length(methtracklist), length(methtracklist))), from = START, to = END)
  dev.off()
} 
makeScreenshot("chr11", 108046044, 110171456) #### Figure 1C ####
makeScreenshot("chr10", 128808505, 131541178) #### Figure 1C ####

# this code was also used for #### Supplemental Figure 5B ####
