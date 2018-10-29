library("BSgenome.Hsapiens.UCSC.hg19")
library(rtracklayer)
library(parallel)


# get the WGBS data (GRanges object, with CpG positions as rows, and
# colnames PDxxxx.T = total reads, PDxxxx.M = methylated reads)
load("~/BiSeq_BASIS/methcounts/meth.gr.RData")

makeBigWig <- function(X, COV) { # X = sample name, COV = minimal coverage (reads)
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


