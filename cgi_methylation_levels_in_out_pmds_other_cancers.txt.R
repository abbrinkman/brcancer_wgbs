library(GenomicRanges)
library(ggplot2)
library(gplots)
library(reshape2)

pmds <- get(load("~/BiSeq_BASIS/PMDs_other_celltypes/pmdmeth.select.RData"))
cgis <- get(load("~/BiSeq_BASIS/PMDs_other_celltypes/cgimeth.select.RData"))

pat <- intersect(names(pmds), names(cgis))


pmds <- pmds[pat]
cgis <- cgis[pat]

# select only the tumors for this
pmds <- pmds[grep("^tumor_", names(pmds))]
cgis <- cgis[grep("^tumor_", names(cgis))]

cgis <- lapply(cgis, function(x) {colnames(mcols(x)) <- "meth" ; x})
 
cgis <- sapply(names(cgis), function(x) {
     cgis[[x]]$PMD <- ifelse(overlapsAny(cgis[[x]], pmds[[x]]), "in", "out") ; cgis[[x]]}, 
     USE.NAMES=T, simplify=F)
 
cgis <- lapply(cgis, function(x) {x$meth.bin <- cut(x$meth, breaks=c(-0.01,0.2,0.4,0.6,0.8,1)) ; x})


# fix names
myMakeUnique <- function(x) {
  x <- make.unique(x)
  x <- gsub("([^0-9])$","\\1\\.0", x)
  x1 <- sapply(strsplit(x, "\\."), function(x) {x[1]})
  x2 <- as.numeric(sapply(strsplit(x, "\\."), function(x) {x[2]}))+1
  x <- paste(x1, x2, sep="_")
  x
}

cgis <- cgis[-grep("ThisStudy", names(cgis))]
names(cgis) <- myMakeUnique(sapply(strsplit(names(cgis), "\\."), function(x) {x[1]}))

# absolute counts
d1 <- melt(lapply(cgis, function(x) {table(x$meth.bin, x$PMD)}), varnames=c("meth", "PMD"), 
   value.name="count")
d1$tumor <- factor(d1$L1, levels=rev(unique(d1$L1)[order(unique(d1$L1))]))

# proportional counts
d2 <- melt(lapply(cgis, function(x) {
    meth.t <- table(x$PMD, x$meth.bin) ; 
    apply(meth.t, 1, function(y) {y/sum(y)})}), varnames=c("meth","PMD"), value.name="proportion")
d2$tumor <- factor(d1$L1, levels=rev(unique(d2$L1)[order(unique(d2$L1))]))


p1 <- ggplot(d2, aes(tumor, proportion)) + 
      geom_bar(stat="identity", aes(fill=meth)) + 
      coord_flip() +
      scale_fill_manual(values=colorpanel(length(unique(d2$meth)), "yellow", "blue")) +
      facet_grid(~PMD) +
      theme_classic() +
      theme(axis.text=element_text(colour="black"), axis.ticks=element_line(colour="black"))
ggsave(p1, 
    file="barplot_proportion_of_CGImeth_absolute_in_out_PMDs.pdf", width=10, height=20, scale=0.6) ##### Supplemental Figure 9B #####



