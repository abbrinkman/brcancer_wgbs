library(rtracklayer)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(data.table)
library(ggplot2)
library(gplots)


# get the PMDs
pmds <- get(load("~/BiSeq_BASIS/PMDs_all_samples/pmdmeth.basis.RData"))
pmds <- pmds[!names(pmds) %in% c("HMEC","MCF7","PD9590a")]
pmds <- GRangesList(lapply(pmds, function(x) {mcols(x) <- NULL ; x}))
pmds.union <- reduce(unlist(pmds))

# get the CpG islands
session <- browserSession()
genome(session) <- "hg19"
cgi <- GRanges(track(ucscTableQuery(session, "CpG Islands")))
cgi <- cgi[seqnames(cgi) %in% seqlevels(pmds)]
seqlevels(cgi) <- seqlevelsInUse(cgi)

# determine how many of the CGIs are within the union of PMDs
cgi$in.pmd <- ifelse(overlapsAny(cgi, pmds.union, type="within"), "in", "out")

p <- ggplot(data.frame("PMD"=cgi$in.pmd), aes(PMD)) +
     geom_bar(aes(fill=PMD), color="black") +
     scale_fill_manual(values=c("out"="white","in"="red")) +
     theme_classic() + 
     theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"))
ggsave(p, file="CGIs_in_out_PMDs.pdf", height=3, width=2.8, scale=0.8) #### Supplemental Figure 4C ####


