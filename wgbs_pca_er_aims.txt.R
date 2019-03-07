
load("~/BiSeq_BASIS/PCA/factominer/d4.pca.RData") # PCA output from factominer
load("~/BiSeq_BASIS/PCA/factominer/d4.dimdesc.RData") # dimdesc object from factominer

options(scipen=3)


# scatterplot of PC1 vs PC2
d2 <- as.data.frame(d4.pca$ind$coord[,1:2])
colnames(d2) <- gsub("Dim\\.","PC",colnames(d2))
d2$patient <- rownames(d2)
d2$ER <- d4.pca$call$X$ER.FPKM
d2$AIMS <- d4.pca$call$X$AIMS
d2$PMDmeth <- d4.pca$call$X$PMD_meth_wgbs
d2$mitotic_score <- d4.pca$call$X$mitotic_score
d2$ER <- factor(d2$ER, levels=c("positive", "negative"))
d2$AIMS <- factor(d2$AIMS, levels=c("Normal", "Basal", "Her2", "LumA", "LumB"))

pdf("boxplots_PC1_PC2_vs_PMDmeth_mitoticScoreER_AIMS.pdf") ##### Supplemental Figure 4A #####
par(mfcol=c(2,3))

plot(PC1 ~ PMDmeth, data=d2, 
	main=paste0("PMDmeth: p=",signif(d4.dimdesc$Dim.1$quanti["PMD_meth_wgbs","p.value"],2)))
plot(PC2 ~ PMDmeth, data=d2, 
	main=paste0("PMDmeth: p=",signif(d4.dimdesc$Dim.2$quanti["PMD_meth_wgbs","p.value"],2)))

boxplot(PC1 ~ ER, data=d2,
        main=paste0("ER: p=",signif(d4.dimdesc$Dim.1$quali["ER.FPKM","p.value"],2)))
boxplot(PC2 ~ ER, data=d2,
        main=paste0("ER: p=",signif(d4.dimdesc$Dim.2$quali["ER.FPKM","p.value"],2)))

boxplot(PC1 ~ mitotic_score, data=d2,
        main=paste0("mit. score: p=",signif(d4.dimdesc$Dim.1$quali["mitotic_score","p.value"],2)))
boxplot(PC2 ~ mitotic_score, data=d2,
        main=paste0("mit. score: p=",signif(d4.dimdesc$Dim.2$quali["mitotic_score","p.value"],2)))

boxplot(PC1 ~ AIMS, data=d2, 
	main=paste0("p=",signif(d4.dimdesc$Dim.1$quali["AIMS","p.value"],2)))
boxplot(PC2 ~ AIMS, data=d2, 
	main=paste0("p=",signif(d4.dimdesc$Dim.2$quali["AIMS","p.value"],2)))
dev.off()
