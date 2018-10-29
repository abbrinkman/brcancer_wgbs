library(GenomicRanges)
library(ggplot2)
library(gplots)
library(reshape2)

pmds <- c(get(load("~/BiSeq_BASIS/PMDs_all_samples/pmdmeth.tcga.RData")),
    get(load("~/BiSeq_BASIS/PMDs_all_samples/pmdmeth.lymph.RData")))
cgis <- c(get(load("~/BiSeq_BASIS/PMDs_all_samples/cgimeth.tcga.RData")), 
    get(load("~/BiSeq_BASIS/PMDs_all_samples/cgimeth.lymph.RData")))

names(pmds) <- make.unique(gsub("_.+$","", names(pmds)), sep="_")
names(cgis) <- make.unique(gsub("_.+$","", names(cgis)), sep="_")


pat <- intersect(names(pmds), names(cgis))

cgis <- cgis[pat]
pmds <- pmds[pat]

cgis <- lapply(cgis, function(x) {colnames(mcols(x)) <- "meth" ; x})

cgis <- sapply(names(cgis), function(x) {
    cgis[[x]]$PMD <- ifelse(overlapsAny(cgis[[x]], pmds[[x]]), "in", "out") ; cgis[[x]]}, 
    USE.NAMES=T, simplify=F)

cgis <- lapply(cgis, function(x) {x$meth.bin <- cut(x$meth, breaks=c(-0.01,0.2,0.4,0.6,0.8,1)) ; x})


# fix name ordering
pat.df <- as.data.frame(do.call(rbind,lapply(strsplit(pat, "_"), function(x) {
    if (length(x)==1) {x[2] <- "0" } ; x})))
rownames(pat.df) <- pat
pat.df$V2 <- as.numeric(as.character(pat.df$V2))+1
pat.df$newname <- paste(pat.df$V1, pat.df$V2, sep="_")
pat.df <- pat.df[with(pat.df, order(V1, V2)),]
pat.names <- pat.df$newname
names(pat.names) <- rownames(pat.df)

# proportional counts
d1 <- melt(lapply(cgis, function(x) {
    meth.t <- table(x$PMD, x$meth.bin) ; 
    apply(meth.t, 1, function(y) {y/sum(y)})}), varnames=c("meth","PMD"), value.name="proportion")
d1$tumor <- pat.df[as.character(d1$L1),]$newname
d1$tumor <- factor(d1$tumor, levels=rev(pat.df$newname))

p1 <- ggplot(d1, aes(tumor, proportion)) + 
      geom_bar(stat="identity", aes(fill=meth)) + 
      coord_flip() +
      scale_fill_manual(values=colorpanel(length(unique(d1$meth)), "yellow", "blue")) +
      facet_grid(~PMD) +
      theme_classic() +
      theme(axis.text=element_text(colour="black"), axis.ticks=element_line(colour="black"))
ggsave(p1, 
    file="barplot_proportion_of_CGImeth_absolute_in_out_PMDs.pdf", width=10, height=13, scale=0.6) #### Supplemental Figure 7B ####

