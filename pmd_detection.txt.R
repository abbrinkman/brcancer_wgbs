library(FDb.UCSC.snp137common.hg19)
library(MethylSeekR)
library("BSgenome.Hsapiens.UCSC.hg19")
library(rtracklayer)
library(parallel)


# load SNP data
snp <- features(FDb.UCSC.snp137common.hg19)

# load cytoband data (for deselecting centromeres)
cytoband <- makeGRangesFromDataFrame(read.table("~/Datarepository/hg19_cytoBand.bed", col.names=c("chr","start","end","name","gieStain")), keep.extra.columns=T)
centro <- cytoband[mcols(cytoband)$gieStain %in% c("acen","gvar","stalk")]

# get the WGBS data (GRanges object, with CpG positions as rows, and
# colnames PDxxxx.T = total reads, PDxxxx.M = methylated reads)
load("~/BiSeq_BASIS/methcounts/meth.gr.RData")

# function for defining the PMDs

findPMDs <- function(X) { # X = sample name
  message(X)
  meth <- meth.gr[, grep(X, colnames(mcols(meth.gr)))]
  seqlevels(meth) <- seqlevelsInUse(meth)
  meth <- removeSNPs(meth, snp)
  pmds <- segmentPMDs(meth, chr.sel="chr20", num.cores=2, seqLengths=seqlengths(Hsapiens),
          pdfFilename=paste0(X, "_alphaDistribution_chr20.pdf"))
  rm(meth)
  pmds
}

# define the PMDs
basis.pmds <- sapply(names(meth.files), findPMDs, USE.NAMES=T, simplify=F)

# extract only the PMDs (remove non-PMDs) and remove centromeric PMDs
pmdsNoCent <- function(GR) {
  pmds.nocent <- GR[mcols(GR)$type=="PMD"]
  suppressWarnings(pmds.nocent <- setdiff(pmds.nocent, centro))
  pmds.nocent
}
basis.pmds.nocent <- lapply(basis.pmds, pmdsNoCent)

# export to BED format
system("~/scripts/fetchChromSizes hg19 > chromsizes_hg19")
exportPMDs <- function(NAME) {
  filename <- paste0("PMDs_", NAME, "_noCent.bed")
  export(basis.pmds.nocent[[NAME]], con=filename, format="BED")
  system(paste0("~/scripts/bedToBigBed ", filename, " chromsizes_hg19 ", gsub("bed$","bb", filename)))
}
sapply(names(basis.pmds.nocent), exportPMDs)
system("rm chromsizes_hg19")

# check WGBS tracks in UCSC
makeBigWig <- function(X, COV) {
  message(X)
  meth <- meth.gr[, grep(X, colnames(mcols(meth.gr)))]
  meth <- meth[!c(grepl("chrM", seqnames(meth)) | grepl("Un", seqnames(meth)) | grepl("random", seqnames(meth)))]
  seqlevels(meth) <- as.character(unique(seqnames(meth)))
  meth <- meth[mcols(meth)[[paste0(X, ".T")]] >= COV]
  mcols(meth)$score <- round(mcols(meth)[[paste0(X, ".M")]]/mcols(meth)[[paste0(X, ".T")]], 2)
  mcols(meth)$T <- NULL
  mcols(meth)$M <- NULL
  tmp.genome <- tempfile(tmpdir="/dev/shm")
  write.table(seqlengths(Hsapiens), file=tmp.genome, quote=F, col.names=F, sep="\t")
  tmp.wig <- tempfile(tmpdir="/dev/shm")
  seqlengths(meth) <- seqlengths(Hsapiens)[names(seqlengths(meth))]
  export(meth, con=tmp.wig, format="WIG")
  outname <- paste0(X, "_CpG_Meth_cov", COV, ".bw")
  system(paste("~/scripts/wigToBigWig", tmp.wig, tmp.genome, outname))
  unlink(tmp.genome)
  unlink(tmp.wig)
}
mclapply(unique(gsub("\\.[TM]$","", colnames(mcols(meth.gr)))), function(X) {makeBigWig(X, 4)}, mc.cores=2)


