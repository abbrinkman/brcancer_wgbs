library(Gviz)
library(RColorBrewer)
library(gplots)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggplot2)

# get the WGBS data (GRanges object, with CpG positions as rows, and
# colnames PDxxxx.T = total reads, PDxxxx.M = methylated reads)
meth <- get(load("~/BiSeq_BASIS/methcounts/meth.gr.RData"))

# for ordering the samples, use the clustering as in Figure 1ABC / Supplemental Figure 2AB
hc <- get(load("~/BiSeq_BASIS/PMDs_all_samples/PMDs_chrX/tilemaps_PMDmeth_clustered_hclust.RData"))
sample.order <- hc$labels[hc$order][hc$labels[hc$order] %in% gsub("\\.[MT]$","",colnames(mcols(meth)))]
sample.order <- sample.order[!sample.order %in% c("HMEC","MCF7", "PD9590a")]

# define function for creating the screenshot, 
# with accompanying methylation heatmaps (separate files) for all CGIs in used region
CGImeth_in_PMDs_screenshot <- function(CHR, START, END) {
  pdf.basename <- paste0(CHR, "_", START, "_", END)
  system(paste("mkdir -p", pdf.basename))
  gen <- 'hg19'
    
  # gene annotation
  #genetrack <- UcscTrack(track="RefSeq Genes", table="refGene",
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
  cgi <- subsetByOverlaps(cgi, GRanges(seqnames=CHR, IRanges(start=START, end=END)))
  cgi$id <- as.character(seq(1:length(cgi)))
  
  cgitrack <- AnnotationTrack(cgi, name="CGI", start=START, end=END, chromosome=CHR, genome=gen,
      fill="darkgreen", col.line="darkgreen", col="darkgreen", stacking="squish",
      background.title="white", fontcolor.title="black",
      showFeatureId=T, featureAnnotation="id", fontcolor.item="black")
  
  
  # PMDs
  pmds <- get(load("~/BiSeq_BASIS/PMDs_all_samples/pmdmeth.basis.RData"))
  pmds <- pmds[!names(pmds) %in% c("HMEC","MCF7", "PD9590a")]
  pmds.m <- reduce(unlist(GRangesList(lapply(pmds, function(x) {mcols(x) <- NULL ; x}))))
  pmds.m$PMDfreq <- apply(sapply(pmds, function(x) {ifelse(overlapsAny(pmds.m, x),1,0)}), 1, sum)
  pmds.m <- subsetByOverlaps(pmds.m, GRanges(seqnames=CHR, IRanges(start=START, end=END)))
  pmdtrack.heatmap <- DataTrack(pmds.m, type="heatmap", name="PMD freq",
        gradient=c("white", "red"),
        background.title="white", fontcolor.title="black")
  min.pmd.size <- 100000
  pmdtrack.separate <- lapply(rev(sample.order), function(x) {
      AnnotationTrack(pmds[[x]][width(pmds[[x]]) >= min.pmd.size], name=x, start=START, end=END,
      chromosome=CHR, fill="red", col.line="red", col="red", stacking="squish", background.title="white",
      fontcolor.title="black", cex.title=0.2, rotation.title=0)})
  
  # plot the annotation overview
  pdf(paste0(pdf.basename, "/", pdf.basename, "_overview_PMDs_separate.pdf"), height=7, width=5)
  plotTracks(trackList=c(pmdtrack.separate, list(genetrack, cgitrack)), from=START, to=END, 
      sizes=c(rep((1/length(pmds)),length(pmds)),1,3))
  dev.off()
  
  pdf(paste0(pdf.basename, "/", pdf.basename, "_overview_PMDs_heatmap.pdf"), height=5, width=5)
  plotTracks(trackList=list(pmdtrack.heatmap, genetrack, cgitrack), from=START, to=END, sizes=c(0.1, 0.1, 0.1))
  dev.off()
  
  # for each CGI in the specified region, make a heatmap of the methylation over the CGI
  cgi.list <- split(cgi, cgi$id)
  flank.length <- 1000
  
  for (i in names(cgi.list)) {
    message("CGI", i)
    pos <- cgi.list[[i]]
    seqlevels(pos) <- seqlevels(meth)
    pos.exp <- reduce(unlist(GRangesList(pos, flank(pos, flank.length, start=T), 
        flank(pos, flank.length, start=F))))

    # select CpGs belonging to the regions of interest 
    meth.sel <- subsetByOverlaps(meth, pos.exp)
    meth.sel <- meth.sel[,!c(
        grepl("HMEC", colnames(mcols(meth.sel))) | grepl(
        "MCF7", colnames(mcols(meth.sel))) | grepl("PD9590a",colnames(mcols(meth))))]

    # calculate methylation from methylated/covered reads    
    samples <- unique(gsub("\\.[TM]$","", colnames(mcols(meth.sel))))
    meth.sel.1 <- meth.sel
    for (j in samples) {
      mcols(meth.sel.1)[[j]] <- round(mcols(meth.sel)[[paste0(j, ".M")]]/mcols(meth.sel)[[paste0(j, ".T")]],2)
      mcols(meth.sel.1)[[paste0(j, ".M")]] <- NULL
      mcols(meth.sel.1)[[paste0(j, ".T")]] <- NULL
    }
    meth.sel.1 <- meth.sel.1[, sample.order] 
    pdf(paste0(pdf.basename, "/", pdf.basename, "_CGI", i, "_", start(pos), "_", end(pos), ".pdf"))
    heatmap(t(as.matrix(mcols(meth.sel.1))), scale="none", Rowv=NA, Colv=NA, 
        col=colorpanel(100, "yellow", "blue"), 
        ColSideColors=ifelse(overlapsAny(meth.sel.1, pos), "darkgreen", "white"), 
        labCol=rep("", length(meth.sel.1)), margins=c(30,5))
    dev.off()
  }
}
CGImeth_in_PMDs_screenshot("chr8",23233438, 25343262) #### Figure 3A ####


