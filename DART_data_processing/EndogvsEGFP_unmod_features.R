library(stringr)
library(ggplot2)

DEseqRes <- read.csv("EndogvsEGFP_DEseq2_norm_unmod.csv")
colnames(DEseqRes)[1] <- "TranscriptID"

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

DEseqRes <- merge(DEseqRes, CountTable_Endog[,c(1,21)], by = "TranscriptID")
DEseqRes$relativeExp <- "No_sig"
DEseqRes$relativeExp <- replace(DEseqRes$relativeExp, DEseqRes$log2FoldChange > 1 & DEseqRes$padj < 0.01, "Higher")
DEseqRes$relativeExp <- replace(DEseqRes$relativeExp, DEseqRes$log2FoldChange < -1 & DEseqRes$padj < 0.01, "Lower")


SeqFeatures <- read.csv("PoolSeqFeatures.csv")

EGFP_Sequence <- SeqFeatures[SeqFeatures$ID %in% CountTable_EGFP$ID,]
Endog_Sequence <- SeqFeatures[SeqFeatures$ID %in% CountTable_Endog$ID,]
EGFP_Sequence$TranscriptID <- sapply(EGFP_Sequence$ID, function(x) strsplit(x, "\\.")[[1]][1])
Endog_Sequence$TranscriptID <- sapply(Endog_Sequence$ID, function(x) strsplit(x, "\\.")[[1]][1])
EGFP_Sequence$CDS <- str_sub(EGFP_Sequence$Sequence, start = EGFP_Sequence$Length+1)
Endog_Sequence$CDS <- str_sub(Endog_Sequence$Sequence, start = Endog_Sequence$Length+1)

EndogvsEGFP <- merge(Endog_Sequence[,c(6,8,22,23)], EGFP_Sequence[,c(6,8,22,23)], by="TranscriptID")

Nucleotide = c("A", "T", "G", "C", "ATG", "GTG", "CTG")
deltaCount_list <- data.frame(TranscriptID = EndogvsEGFP$TranscriptID)
for(i in 1:7){
  nucleotide = Nucleotide[i]
  deltaCont = (str_count(EndogvsEGFP$CDS.x, pattern = nucleotide) - str_count(EndogvsEGFP$CDS.y, pattern = nucleotide))
  deltaPercent = deltaCont / nchar(EndogvsEGFP$CDS.x) * 100
  deltaCount_list <- data.frame(deltaCount_list, deltaCont, deltaPercent)
}
colnames(deltaCount_list)[-1] <- paste("delta", rep(Nucleotide, each=2), "CDS", c("Count", "Percent"), sep="_")
deltaCount_list$delta_GC_CDS_count <- deltaCount_list$delta_C_CDS_Count + deltaCount_list$delta_G_CDS_Count
deltaCount_list$delta_GC_CDS_Percent <- deltaCount_list$delta_C_CDS_Percent + deltaCount_list$delta_G_CDS_Percent
deltaCount_list$delta_MFE <- EndogvsEGFP$MFE.x - EndogvsEGFP$MFE.y
deltaCount_list$delta_Kozak <- EndogvsEGFP$Kozak_score.x - EndogvsEGFP$Kozak_score.y

DEseqRes <- merge(DEseqRes, deltaCount_list, by = "TranscriptID")

write.csv(DEseqRes, "EndogvsEGFP_sig_nucleotideContent.csv", quote = F, row.names = F)


pval_list <- vector()
for(i in 10:27){
  pval <- wilcox.test(DEseqRes[DEseqRes$relativeExp == "Higher", i], 
                      DEseqRes[DEseqRes$relativeExp == "Lower", i])$p.value
  pval_list <- c(pval_list, pval)
}
pval_list <- data.frame(padj = p.adjust(pval_list, n=length(pval_list), method = "BH"))
pval_list$Feature <- colnames(DEseqRes)[10:27]
pval_list$Feature <- factor(pval_list$Feature, levels = pval_list$Feature[order(pval_list$padj, decreasing = F)])
ggplot(pval_list, aes(-log10(padj), Feature))+
  geom_bar(stat = "identity")+theme_classic()+geom_vline(xintercept = 2, linetype = "dotted", color="red")
ggsave("EndogvsEGFP_feature_sig.pdf")
write.csv(pval_list, "EndogvsEGFP_feature_sig.csv", quote = F, row.names = F)

colnames(DEseqRes) <- str_replace_all(colnames(DEseqRes), pattern = "T_", replacement = "U_")
DEseqRes$relativeExp <- factor(DEseqRes$relativeExp, levels = c("Higher", "No_sig", "Lower"))
#DEseqRes$delta_CC <- (str_count(Endog_Sequence$CDS, pattern = "CC") - str_count(EGFP_Sequence$CDS, pattern = "CC"))

write.csv(DEseqRes, "EndogvsEGFP_DEseq_table_unmod.csv", quote = F, row.names = F)
