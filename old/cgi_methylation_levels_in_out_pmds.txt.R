library(GenomicRanges)
library(ggplot2)
library(gplots)
library(reshape2)

# load PMD and CGI location/methylation data
pmds <- get(load("~/BiSeq_BASIS/PMDs_all_samples/pmdmeth.basis.RData"))
cgis <- get(load("~/BiSeq_BASIS/PMDs_all_samples/cgimeth.basis.RData"))

pat <- intersect(names(pmds), names(cgis))
pat <- pat[! pat %in% c("HMEC","MCF7", "PD9590a")]

cgis <- cgis[pat]
pmds <- pmds[pat]

cgis <- lapply(cgis, function(x) {colnames(mcols(x)) <- "meth" ; x})

# determine which CGIs are inside or outside PMDs
cgis <- sapply(names(cgis), function(x) {
    cgis[[x]]$PMD <- ifelse(overlapsAny(cgis[[x]], pmds[[x]]), "in", "out") ; cgis[[x]]}, 
    USE.NAMES=T, simplify=F)

# bin the CGI methylation values into five bins
cgis <- lapply(cgis, function(x) {x$meth.bin <- cut(x$meth, breaks=c(-0.01,0.2,0.4,0.6,0.8,1)) ; x})

# plot the proportional counts of CGIs within each methylation bin, split over in/out PMD
d1 <- melt(lapply(cgis, function(x) {
    meth.t <- table(x$PMD, x$meth.bin) ; 
    apply(meth.t, 1, function(y) {y/sum(y)})}), varnames=c("meth","PMD"), value.name="proportion")

p1 <- ggplot(d1, aes(L1, proportion)) + 
      geom_bar(stat="identity", aes(fill=meth)) + 
      coord_flip() +
      scale_fill_manual(values=colorpanel(length(unique(d1$meth)), "yellow", "blue")) +
      facet_grid(~PMD) +
      theme_classic() +
      theme(axis.text=element_text(colour="black"), axis.ticks=element_line(colour="black"))
ggsave(p1, 
    file="barplot_proportion_of_CGImeth_absolute_in_out_PMDs.pdf", width=10, height=7, scale=0.6) #### Figure 3B ####

