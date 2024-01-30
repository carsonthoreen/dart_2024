CountTable <- read.csv("analysis/SumBarcodes.counts.csv")
CountMatrix <- as.matrix(CountTable[,c(2:13, 26:31)])
row.names(CountMatrix) <- CountTable$X
sum <- colSums(CountMatrix)
Normalizing <- round(sapply(1:(length(sum)/2), function(x) sum(sum[((x * 2) - 2 + 1):(x * 2)]))/1000000)

CountMatrix <- rbind(CountMatrix, rep(Normalizing, each=2))
Normalizing_index <- NROW(CountMatrix)


library(DESeq2)
coldata <- data.frame(Sample = colnames(CountMatrix), Mod = factor(rep(c("Y", "U"), 9)), Type = c(rep("Monosome", 12), rep("IVT", 6)))
dds <- DESeqDataSetFromMatrix(countData = CountMatrix, colData = coldata, design = ~ Type + Mod + Type:Mod)
dds <- estimateSizeFactors(dds, controlGenes = Normalizing_index)
dds <- DESeq(dds)
res <- results(dds)



library(dplyr)
#res <- results(dds)
DEseqRes <- data.frame(res)
write.csv(DEseqRes, "analysis/DEseq2_mod_m1YvsUnmod_spikein.csv", quote=F)
DEseqRes <- DEseqRes[complete.cases(DEseqRes),]
mutateddf <- mutate(DEseqRes, significance=ifelse(DEseqRes$padj<0.01 & abs(DEseqRes$log2FoldChange) > 1, "padj < 0.01 and LFC > 1", "Not Sig")) #Will have different colors depending on significance
mutateddf$geneName <- sapply(row.names(mutateddf), function(x){strsplit(x, "\\.")[[1]][1]})

ggplot(mutateddf, aes(log2FoldChange, -log10(padj))) + #volcanoplot with log2Foldchange versus pvalue
  geom_point(aes(col=significance), alpha = 0.5) + #add points colored by significance
  scale_color_manual(values=c("grey", "red"))+
  theme_classic()+
  geom_vline(xintercept = c(-1,1), linetype = "dotted")+
  geom_hline(yintercept = 2, linetype = "dotted")+
  xlab("log2FoldChange m1Y vs Unmod")+ylab("Significance -log10(padj)")+
  theme(text = element_text(size = 10), legend.position = "none")
ggsave("analysis/Mod_volc.pdf", width = 6, height = 7)


