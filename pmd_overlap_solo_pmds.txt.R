library(rtracklayer)
library(VennDiagram)
library(ggplot2)
library(reshape2)

pmds.solo <- get(load("~/BiSeq_BASIS/solo-WCGW/PMDs_solo-WCGW/pmdsolometh.basis.RData"))
pmds.solo <- pmds.solo[!names(pmds.solo) %in% c("HMEC", "MCF7")]

pmds <- get(load("~/BiSeq_BASIS/PMDs_all_samples/pmdmeth.basis.RData"))

pmds.solo <- lapply(pmds.solo, function(x) {x[!seqnames(x) %in% c("chrY","chrM")]})
pmds.solo <- lapply(pmds.solo, function(x) {colnames(mcols(x)) <- "meth" ; x})
pmds <- lapply(pmds, function(x) {x[!seqnames(x) %in% c("chrY","chrM")]})
pmds <- lapply(pmds, function(x) {colnames(mcols(x)) <- "meth" ; x})

# which fraction of the genome is covered in total (union)
pmds.solo.bases <- sum(width(reduce(unlist(GRangesList(pmds.solo)))))
pmds.bases <- sum(width(reduce(unlist(GRangesList(pmds)))))

pmds.solo.bases
# [1] 2267183395
pmds.bases
#  [1] 2264532147
pmds.solo.bases/3.2e09
# [1] 0.7084948
pmds.bases/3.2e09
# [1] 0.7076663


# determine overlaps
fractionOverlap <- function(x, y, x.name, y.name) { # x, y are GRanges objects
  x.y.frac <- sum(width(intersect(x,y)))/sum(width(union(x,y)))
  x.frac <- sum(width(setdiff(x, intersect(x,y))))/sum(width(union(x,y)))
  y.frac <- sum(width(setdiff(y, intersect(x,y))))/sum(width(union(x,y)))
  out <- round(c(x.frac, x.y.frac, y.frac), 2)
  names(out) <- c(x.name, paste0(x.name, " & ", y.name), y.name)
  out
}

fractionOverlap2Venn <- function(ov) {
  ov <- round(ov, 2)*100
  out <- list(x=seq(1, ov[1]+ov[2], by=1),
              y=seq(ov[1], 100, by=1))
  names(out) <- names(ov)[c(1,3)]
  out
}

# overlap of the two union tracks
ov.unions <- fractionOverlap(
  x=reduce(unlist(GRangesList(pmds.solo))), 
  y=reduce(unlist(GRangesList(pmds))),
  x.name="PMDs (solo-WCGW)",
  y.name="PMDs"
)
p1 <- venn.diagram(fractionOverlap2Venn(ov.unions), fill = c("red", "green"),
    alpha = c(0.3, 0.3), cex = 2, cat.fontface = 4, lty =1,
    filename=NULL, height=2000, width=2000)
pdf("venn_overlap_PMD_unions.pdf") ##### Supplemental Figure 4D #####
grid.newpage()
grid.draw(p1)
dev.off()
system("rm VennDiagram*.log")
