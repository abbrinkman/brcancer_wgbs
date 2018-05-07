library(GenomicRanges)
library(ggplot2)
library(reshape2)
library(parallel)
library(data.table)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19)

meth <- get(load("~/BiSeq_BASIS/methcounts/meth.gr.RData"))
pmds <- get(load("~/BiSeq_BASIS/PMDs_all_samples/pmdmeth.basis.RData"))

# get the hg19 genome for constructing intergenic regions
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
seqlevels(hg19) <- seqlevelsInUse(hg19)

# get CGIs, promoters, genes
session <- browserSession()
genome(session) <- "hg19"
query <- ucscTableQuery(session, "CpG Islands")
cgi <- GRanges(track(query))
cgi$name <- NULL
cgi$elem <- "CGI"
cgi <- cgi[,"elem"]

shore.size <- 2000
cgi.shores <- cgi
start(cgi.shores) <- start(cgi.shores)-shore.size
end(cgi.shores) <- end(cgi.shores)+shore.size
cgi.shores$elem <- "CGI_shores"
cgi.shores <- cgi.shores[,"elem"]

query <- ucscTableQuery(session, "GENCODE Genes V19")
genes <- GRanges(track(query))

# set promoter size: size centered around transcription start site
prom.size <- 1000 
prom <- promoters(genes, upstream=prom.size/2, downstream=prom.size/2)
prom$elem <- "promoters"
prom <- prom[,"elem"]

bodies <- reduce(genes)
bodies$elem <- "gene_bodies"
bodies <- bodies[,"elem"]

# get chromHMM states from HMEC (to determine enhancer regions)
# http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmHmecHMM.bed.gz
chmm <- read.table("~/BiSeq_BASIS/PMDs_all_samples/LMR_DMR_CR_characterization/wgEncodeBroadHmmHmecHMM.bed.gz")
chmm <- chmm[, c(1:4,9)]
colnames(chmm) <- c("chr", "start", "end", "state", "color")
chmm <- makeGRangesFromDataFrame(chmm, keep.extra.columns=T)

# create an option to specifically use only 'strong enhancers' or all enhancers ('enhancers') (chromHMM state)
strong.enhancers <- TRUE
if (strong.enhancers) {
  enh <- chmm[grep("strong_enhancer", chmm$state, ignore.case=T)]
} else {
  enh <- chmm[grep("enhancer", chmm$state, ignore.case=T)]
}
enh$elem <- "HMEC_enhancers"
enh <-enh[,"elem"]

elements <- c(cgi, cgi.shores, prom, bodies, enh)
elements <- elements[!c(grepl("random", seqnames(elements)) | grepl("chrUn", 
    seqnames(elements)) | grepl("hap", seqnames(elements)) | grepl("chr[XYM]", seqnames(elements)))]
seqlevels(elements) <- seqlevelsInUse(elements)

intergenic <- setdiff(hg19, elements)
intergenic$elem <- "intergenic"
elements <- sort(c(elements, intergenic))

pat <- intersect(names(pmds), unique(gsub("\\.[TM]","", colnames(mcols(meth)))))
pat <- pat[!pat %in% c("HMEC","MCF7","PD9590a")]
min.cov <- 10 # set a CpG coverage of minimally 10
meth.cov <- meth
for (i in pat) { 
  message(i)
  mcols(meth.cov)[[i]] <- round(mcols(meth)[[paste0(i, ".M")]]/mcols(meth)[[paste0(i, ".T")]], 2)
  mcols(meth.cov)[[i]][mcols(meth)[[paste0(i, ".T")]] < min.cov] <- NaN
}
meth.cov <- meth.cov[, -grep(".[TM]", colnames(mcols(meth.cov)))]

cpg.pmd <- as.data.frame(sapply(pat, function(x) {overlapsAny(meth.cov, pmds[[x]])}))

cpg.elem <- sapply(unique(elements$elem), function(x) {
    overlapsAny(meth.cov, elements[elements$elem==x])}, USE.NAMES=T, simplify=F)

# due to memory restraints, remove some stuff
rm(meth, meth.gr)
gc()

getMethValues.1 <- function(pd, elem.name) {
  in.pmds = mcols(meth.cov)[[pd]][cpg.pmd[,pd] & cpg.elem[[elem.name]]]
  out.pmds = mcols(meth.cov)[[pd]][!cpg.pmd[,pd] & cpg.elem[[elem.name]]]
  in.pmds <- in.pmds[!is.na(in.pmds)]
  out.pmds <- out.pmds[!is.na(out.pmds)]
  # due to memory constrains, use only part of it (20% of CpGs is still OK on 48 Gb RAM):
  in.pmds <- sample(in.pmds, 0.2*length(in.pmds))
  out.pmds <- sample(out.pmds, 0.2*length(out.pmds))
  list("in"=in.pmds, "out"=out.pmds)
}

l2 <- sapply(pat, function(x) {
  message(x)
  sapply(unique(elements$elem), function(y) {
    message("   ", y)
    getMethValues.1(x, y)}, USE.NAMES=T, simplify=F)}, USE.NAMES=T, simplify=F)

d2 <- melt(l2)
colnames(d2) <- c("meth","PMD","element","patient")
d2$PMD <- factor(d2$PMD, levels=c("out","in"))
d2$element <- factor(d2$element,
    levels=c("promoters","CGI","CGI_shores","HMEC_enhancers", "gene_bodies","intergenic"))


p1 <- ggplot(d2, aes(element, meth)) +
      geom_boxplot(aes(fill=PMD), outlier.shape=NA) +
      facet_wrap(~PMD) +
      scale_fill_manual(values=c("in"="red", "out"="white")) +
      theme_classic() +
      theme(axis.text=element_text(color="black"),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks=element_line(color="black"))
if (strong.enhancers) {
  pdfname <- paste0("boxplots_meth_genomic_elements_perCpG_prom", prom.size, "_strongEnh_in_out_PMDs.pdf")
} else {
  pdfname <- paste0("boxplots_meth_genomic_elements_perCpG_prom", prom.size, "_allEnh_in_out_PMDs.pdf")
}
ggsave(p1, file=pdfname, height=4.5, width=7, scale=0.5) #### Figure 2H #### (using strong enhancers)

