RRSdata <- read.csv("meanRRS.csv")
colnames(RRSdata)[1] <- "ID"
RRSdata$type <- sapply(RRSdata$ID, function(x) tail(strsplit(x, "\\.")[[1]], 1))
RRSdata$TranscriptID <- sapply(RRSdata$ID, function(x) strsplit(x, "\\.")[[1]][1])
RRSdata_EGFP <- RRSdata[RRSdata$type == "EGFPCDS",]
RRSdata_Endog <- RRSdata[RRSdata$type == "EndogCDS",]
RRSdata_Endog <- RRSdata_Endog[-grep("ENST00000260605.12.100.EndogCDS", RRSdata_Endog$ID),]
SeqName_intersect <- intersect(RRSdata_EGFP$TranscriptID, RRSdata_Endog$TranscriptID)

RRSdata_EGFP <- RRSdata_EGFP[RRSdata_EGFP$TranscriptID %in% SeqName_intersect,]
RRSdata_Endog <- RRSdata_Endog[RRSdata_Endog$TranscriptID %in% SeqName_intersect,]
EndogvsEGFP <- merge(RRSdata_EGFP[,c(2,3,5)], RRSdata_Endog[,c(2,3,5)], by = "TranscriptID")
library(ggpointdensity)
library(viridis)
ggplot(EndogvsEGFP, aes(RRS_Unmod.x, RRS_Unmod.y))+
  geom_pointdensity(size=0.1)+scale_color_viridis()+
  theme_classic()+scale_x_continuous(name="EGFPCDS log2RRS", limits=c(-6, 3))+
  scale_y_continuous(name="EndogCDS log2RRS", limits=c(-6, 3))+
  geom_smooth(se = F, color = "red", linewidth = 0.5, method = "lm")+
  theme(legend.position = "none", line = element_line(linewidth = 0.5))
ggsave("EGFPvsEndog_scatter.pdf", width = 3, height = 3)

CountTable <- read.csv("SumBarcodes.counts.csv")[,-c(14:25)]
colnames(CountTable)[1] <- "ID"
CountTable$type <- sapply(CountTable$ID, function(x) tail(strsplit(x, "\\.")[[1]], 1))
CountTable$TranscriptID <- sapply(CountTable$ID, function(x) strsplit(x, "\\.")[[1]][1])
CountTable_EGFP <- CountTable[CountTable$type == "EGFPCDS",]
CountTable_Endog <- CountTable[CountTable$type == "EndogCDS",]
CountTable_Endog <- CountTable_Endog[-grep("ENST00000260605.12.100.EndogCDS", CountTable_Endog$ID),]
SeqName_intersect <- intersect(CountTable_EGFP$TranscriptID, CountTable_Endog$TranscriptID)
CountTable_EGFP <- CountTable_EGFP[CountTable_EGFP$TranscriptID %in% SeqName_intersect,]
CountTable_Endog <- CountTable_Endog[CountTable_Endog$TranscriptID %in% SeqName_intersect,]


####Unmod
CountMatrix <- merge(CountTable_EGFP[,seq(3,21,2)], CountTable_Endog[,seq(3,21,2)], by = "TranscriptID")
row.names(CountMatrix) <- CountMatrix$TranscriptID
CountMatrix <- as.matrix(CountMatrix[,2:19])
sum <- round((colSums(CountMatrix[,1:9]) + colSums(CountMatrix[,10:18]))/100000)
CountMatrix <- rbind(CountMatrix, rep(sum,2))

library(DESeq2)
coldata <- data.frame(Sample = colnames(CountMatrix), RNA = rep(c("EGFP", "Endog"), each = 9), Type = rep(c(rep("Monosome", 6), rep("Input",3)), 2))
coldata$RNA <- as.factor(coldata$RNA)
coldata$Type <- as.factor(coldata$Type)


dds <- DESeqDataSetFromMatrix(countData = CountMatrix, colData = coldata, design = ~ Type + RNA + Type:RNA)
dds <- estimateSizeFactors(dds, controlGenes = NROW(CountMatrix))
dds <- DESeq(dds)
res <- results(dds)

DEseqRes <- data.frame(res)
write.csv(DEseqRes, "EndogvsEGFP_DEseq2_norm_unmod.csv")
DEseqRes <- DEseqRes[complete.cases(DEseqRes),]
mutateddf <- mutate(DEseqRes, significance=ifelse(DEseqRes$padj<0.01 & abs(DEseqRes$log2FoldChange) > 1, "padj < 0.01 and LFC > 1", "Not Sig")) #Will have different colors depending on significance
library(ggplot2)
ggplot(mutateddf, aes(log2FoldChange, -log10(padj))) + #volcanoplot with log2Foldchange versus pvalue
  geom_point(aes(col=significance), alpha = 0.5) + #add points colored by significance
  scale_color_manual(values=c("grey", "red"))+
  theme_classic()+
  geom_vline(xintercept = c(-1,1), linetype = "dotted")+
  geom_hline(yintercept = 2, linetype = "dotted")+
  xlab("log2FoldChange EndogvsEGFP")+ylab("Significance -log10(padj)")+
  theme(text = element_text(size = 10), legend.position = "none")
ggsave("EndogvsEGFP_DEseq_volc_norm.pdf", width = 6, height = 7)


####m1Y
CountMatrix <- merge(CountTable_EGFP[,c(seq(2,18,2), 21)], CountTable_Endog[,c(seq(2,18,2), 21)], by = "TranscriptID")
row.names(CountMatrix) <- CountMatrix$TranscriptID
CountMatrix <- as.matrix(CountMatrix[,2:19])
sum <- round((colSums(CountMatrix[,1:9]) + colSums(CountMatrix[,10:18]))/100000)
CountMatrix <- rbind(CountMatrix, rep(sum,2))

library(DESeq2)
coldata <- data.frame(Sample = colnames(CountMatrix), RNA = rep(c("EGFP", "Endog"), each = 9), Type = rep(c(rep("Monosome", 6), rep("Input",3)), 2))
coldata$RNA <- as.factor(coldata$RNA)
coldata$Type <- as.factor(coldata$Type)


dds <- DESeqDataSetFromMatrix(countData = CountMatrix, colData = coldata, design = ~ Type + RNA + Type:RNA)
dds <- estimateSizeFactors(dds, controlGenes = NROW(CountMatrix))
dds <- DESeq(dds)
res <- results(dds)

DEseqRes <- data.frame(res)
write.csv(DEseqRes, "EndogvsEGFP_DEseq2_norm_m1Y.csv")
DEseqRes <- DEseqRes[complete.cases(DEseqRes),]
mutateddf <- mutate(DEseqRes, significance=ifelse(DEseqRes$padj<0.01 & abs(DEseqRes$log2FoldChange) > 1, "padj < 0.01 and LFC > 1", "Not Sig")) #Will have different colors depending on significance
library(ggplot2)
ggplot(mutateddf, aes(log2FoldChange, -log10(padj))) + #volcanoplot with log2Foldchange versus pvalue
  geom_point(aes(col=significance), alpha = 0.5) + #add points colored by significance
  scale_color_manual(values=c("grey", "red"))+
  theme_classic()+
  geom_vline(xintercept = c(-1,1), linetype = "dotted")+
  geom_hline(yintercept = 2, linetype = "dotted")+
  xlab("log2FoldChange EndogvsEGFP")+ylab("Significance -log10(padj)")+
  theme(text = element_text(size = 10), legend.position = "none")
ggsave("EndogvsEGFP_DEseq_volc_norm_m1Y.pdf", width = 6, height = 7)
