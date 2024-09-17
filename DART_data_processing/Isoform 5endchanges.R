library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(ggpointdensity)


IsoformTable <- read.table("Rebuttal/RRS_isoform.tsv", sep="\t", header = T)
IsoformTable$Seq1 <- str_sub(IsoformTable$Seq1, start = 4)
IsoformTable$Seq2 <- str_sub(IsoformTable$Seq2, start = 4)

SeqAlign <- function(Seq1, Seq2){
  match_pattern <- c(attr(adist(Seq1, Seq2, counts = TRUE), "trafos"))
  pattern_locate <- str_locate_all(match_pattern, "M+")[[1]]
  return(NROW(pattern_locate) == 1 & max(pattern_locate) == max(nchar(Seq1), nchar(Seq2)))
}
IsoformTable <- IsoformTable[,-c(10,15,16)]
IsoformTable$Endextension <- mapply(SeqAlign, IsoformTable$Seq1, IsoformTable$Seq2)
IsoformTable$ExtentedSeq <- ifelse(IsoformTable$Length1 > IsoformTable$Length2,
                                        str_sub(IsoformTable$Seq1, end = IsoformTable$LengthDiff),
                                        str_sub(IsoformTable$Seq2, end = IsoformTable$LengthDiff))

IsoformTable$LongerIsBetter <- "No Sig"
IsoformTable$LongerIsBetter[IsoformTable$padj < 0.05 & IsoformTable$deltaRRS * (IsoformTable$Length1 - IsoformTable$Length2) > 0] <- "Extension is better"
IsoformTable$LongerIsBetter[IsoformTable$padj < 0.05 & IsoformTable$deltaRRS * (IsoformTable$Length1 - IsoformTable$Length2) < 0] <- "Extension is worse"

IsoformTable$deltaCC_count <- str_count(IsoformTable$Seq1, "CC") - str_count(IsoformTable$Seq2, "CC")
IsoformTable$deltaCC_count[IsoformTable$deltaCC_count >= 10] <- 10
IsoformTable$deltaCC_count[IsoformTable$deltaCC_count <= -10] <- -10
IsoformTable$deltaCC_count <- factor(IsoformTable$deltaCC_count, levels = -10:10, labels = c("≤-10", -9:9, "≥10"))
IsoformTable$deltaUU_count <- str_count(IsoformTable$Seq1, "TT") - str_count(IsoformTable$Seq2, "TT")
IsoformTable$deltaUU_count[IsoformTable$deltaUU_count >= 10] <- 10
IsoformTable$deltaUU_count[IsoformTable$deltaUU_count <= -10] <- -10
IsoformTable$deltaUU_count <- factor(IsoformTable$deltaUU_count, levels = -10:10, labels = c("≤-10", -9:9, "≥10"))


IsoformTable_sub <- IsoformTable %>% filter(Length1 < 100 & Length2 < 100)
ggplot(IsoformTable_sub, aes(deltaCC_count, deltaRRS))+
  geom_boxplot()+theme_classic()+xlab("differences in CC motif counts")+
  geom_hline(yintercept = 0, linetype = "dotted")
ggsave("Rebuttal/Isoform_deltaCC.pdf", width = 5, height = 3)
ggplot(IsoformTable_sub, aes(deltaUU_count, deltaRRS))+
  geom_boxplot()+theme_classic()+xlab("differences in UU motif counts")+
  geom_hline(yintercept = 0, linetype = "dotted")
ggsave("Rebuttal/Isoform_deltaUU.pdf", width = 5, height = 3)


blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )
IsoformTable_5ext_summary <- IsoformTable %>% group_by(Endextension) %>% summarise(Count = n())
ggplot(IsoformTable_5ext_summary, aes(x = "", y = Count, fill = Endextension)) +
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+blank_theme +
  theme(axis.text.x=element_blank())
ggsave("Rebuttal/Isoform_extention.pdf")
IsoformTable_5ext <- IsoformTable[IsoformTable$Endextension,]
IsoformTable_5ext_sig_summary <- IsoformTable_5ext %>% group_by(padj < 0.05) %>% summarise(Count = n())
ggplot(IsoformTable_5ext_sig_summary, aes(x = "", y = Count, fill = `padj < 0.05`)) +
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+blank_theme +
  theme(axis.text.x=element_blank())+scale_fill_manual(values=c("#E69F00", "#56B4E9"))
ggsave("Rebuttal/Isoform_ext_sig.pdf")

IsoformTable_5ext <- IsoformTable[IsoformTable$Endextension,]
write.table(IsoformTable_5ext, "Rebuttal/Isoform_5ext.tsv", quote = F, row.names = F, sep = "\t")

UTR_mfe <- read.csv("Rebuttal/meanRRS_UTRMFE.csv")
UTR_mfe_EGFP <- UTR_mfe[UTR_mfe$Type == "EGFPCDS",]
UTR_mfe_EGFP$UTR <- str_sub(UTR_mfe_EGFP$UTR, start = 4)
MFE_table <- data.frame(IsoformTable_5ext[,c(1,8,9)], deltaMFE = 0)
for(i in 1:NROW(MFE_table)){
  MFE_table$deltaMFE[i] = UTR_mfe_EGFP$MFE_UTR[UTR_mfe_EGFP$UTR == MFE_table$Seq1[i]] - UTR_mfe_EGFP$MFE_UTR[UTR_mfe_EGFP$UTR == MFE_table$Seq2[i]]
}
IsoformTable_5ext <- merge(IsoformTable_5ext, MFE_table, by = c("geneID", "Seq1", "Seq2"))

IsoformTable_5ext$LongerIsBetter <- factor(IsoformTable_5ext$LongerIsBetter, levels = c("Extension is better", "No Sig", "Extension is worse"))
library(seqinr)
for(i in grep(pattern = "better", IsoformTable_5ext$LongerIsBetter)){
  if(IsoformTable_5ext$LengthDiff[i] >= 8) write.fasta(sequences = IsoformTable_5ext$ExtentedSeq[i], names = paste(IsoformTable_5ext$isoform1[i], IsoformTable_5ext$isoform2[i], sep = "_"), file.out = "DREME/ExtImprove.fasta", open = "a")
}
for(i in grep(pattern = "worse", IsoformTable_5ext$LongerIsBetter)){
  if(IsoformTable_5ext$LengthDiff[i] >= 8) write.fasta(sequences = IsoformTable_5ext$ExtentedSeq[i], names = paste(IsoformTable_5ext$isoform1[i], IsoformTable_5ext$isoform2[i], sep = "_"), file.out = "DREME/ExtDecrease.fasta", open = "a")
}
for(i in 1:NROW(IsoformTable_5ext)){
  if(IsoformTable_5ext$LengthDiff[i] >= 8) write.fasta(sequences = IsoformTable_5ext$ExtentedSeq[i], names = paste(IsoformTable_5ext$isoform1[i], IsoformTable_5ext$isoform2[i], sep = "_"), file.out = "DREME/ExtAll.fasta", open = "a")
}
IsoformTable_5ext$LengthDiff <- IsoformTable_5ext$Length1 - IsoformTable_5ext$Length2
ggplot(IsoformTable_5ext, aes(Length1 - Length2, deltaRRS))+
  geom_pointdensity(size=0.1)+theme_classic()+theme(legend.position = "none")+xlab("deltaLength")+
  scale_color_viridis()+geom_smooth(method = "lm", se = F, color = "red")
ggsave("Rebuttal/Extetion_Length.pdf", width = 4, height = 3)
IsoformTable_5ext$LengthDiff <- IsoformTable_5ext$Length1 - IsoformTable_5ext$Length2
ggplot(IsoformTable_5ext[abs(IsoformTable_5ext$LengthDiff) <= 50,], aes(Length1 - Length2, deltaRRS))+
  geom_pointdensity(size=0.1)+theme_classic()+theme(legend.position = "none")+xlab("deltaLength")+
  scale_color_viridis()+geom_smooth(method = "lm", se = F, color = "red")
ggsave("Rebuttal/Extetion_Length<50.pdf", width = 4, height = 3)

IsoformTable_5ext$deltaGC <- str_count(pattern = "[G|C]", IsoformTable_5ext$Seq1)/IsoformTable_5ext$Length1*100 - 
  str_count(pattern = "[G|C]", IsoformTable_5ext$Seq2)/IsoformTable_5ext$Length2*100
ggplot(IsoformTable_5ext, aes(deltaGC, deltaRRS))+
  geom_pointdensity()+theme_classic()+theme(legend.position = "none")+xlab("deltaLength")+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)
ggsave("Rebuttal/Extetion_GC.pdf", width = 4, height = 3)
ggplot(IsoformTable_5ext, aes(deltaMFE, deltaRRS))+
  geom_pointdensity()+theme_classic()+theme(legend.position = "none")+xlab("deltaMFE (kcal/mol)")+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+scale_x_reverse()
ggsave("Rebuttal/Extetion_MFE.pdf", width = 4, height = 3)

IsoformTable_5ext_within50 <- IsoformTable_5ext[abs(IsoformTable_5ext$LengthDiff) < 50,]
ggplot(IsoformTable_5ext_within50, aes(Length1 - Length2, deltaRRS))+
  geom_pointdensity()+theme_classic()+theme(legend.position = "none")+xlab("deltaLength")+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)
ggsave("Rebuttal/Extetion_Length_<50.pdf", width = 4, height = 3)

IsoformTable_5ext <- IsoformTable_5ext[rowMaxs(as.matrix(IsoformTable_5ext[,12:13])) < 100,]
IsoformTable_5ext$Extent_GC <- (str_count(IsoformTable_5ext$ExtentedSeq, "G") + str_count(IsoformTable_5ext$ExtentedSeq, "C"))/IsoformTable_5ext$LengthDiff * 100
IsoformTable_5ext$Extent_CC <- str_count(IsoformTable_5ext$ExtentedSeq, "CCC*")
IsoformTable_5ext$Extent_CC <- ifelse(IsoformTable_5ext$Length1 < IsoformTable_5ext$Length2, 0 - IsoformTable_5ext$Extent_CC, IsoformTable_5ext$Extent_CC)


IsoformTable_5ext$Extent_CWCC <- str_count(IsoformTable_5ext$ExtentedSeq, "C[A|T]CC")
IsoformTable_5ext$Extent_CWCC <- ifelse(IsoformTable_5ext$Length1 < IsoformTable_5ext$Length2, 0 - IsoformTable_5ext$Extent_CWCC, IsoformTable_5ext$Extent_CWCC)


ggplot(IsoformTable_5ext, aes(LongerIsBetter, Extent_GC))+
  geom_boxplot()+theme_classic()+theme(axis.title.x = element_blank())+
  ylab("GC % in extension")
ggsave("Extetion_GC_cont.pdf", width = 3, height = 3)
ggplot(IsoformTable_5ext, aes(LongerIsBetter, LengthDiff))+
  geom_boxplot()+theme_classic()+theme(axis.title.x = element_blank())+
  ylab("Length of the extension")
ggsave("Extetion_Length.pdf", width = 3, height = 3)




ggplot(IsoformTable_5ext, aes(factor(abs(Extent_CWCC)), deltaRRS))+
  geom_boxplot()+theme_classic()
ggsave("Extetion_CWCC_cont.pdf", width = 4, height = 3)



#Endchanges <- IsoformTable[IsoformTable$Endextension,c(1,3,8:10,12,13)]
Endchanges <- IsoformTable_5ext[, c(1,3,8:10,12,13)]
#Endchanges <- Endchanges[Endchanges$LengthDiff <= 20 & Endchanges$LengthDiff / rowMaxs(as.matrix(Endchanges[,6:7])) < 0.2,]
#Endchanges <- Endchanges[Endchanges$LengthDiff <= 20,]
deltaContent <- data.frame(ID = Endchanges$geneID)
for(i in 1:NROW(deltaContent)){
  Seq1 = str_sub(Endchanges$Seq1[i], end = 6)
  Seq2 = str_sub(Endchanges$Seq2[i], end = 6)
  deltaA = str_count(Seq1, "A") - str_count(Seq2, "A")
  deltaU = str_count(Seq1, "T") - str_count(Seq2, "T")
  deltaG = str_count(Seq1, "G") - str_count(Seq2, "G")
  deltaC = str_count(Seq1, "C") - str_count(Seq2, "C")
  deltaContent[i, 2:5] = c(deltaA, deltaU, deltaG, deltaC)
  #if (Endchanges$Length1[i] - Endchanges$Length2[i] > 0) {deltaContent[i, 2:5] = -deltaContent[i, 2:5]}
}
colnames(deltaContent)[-1] <- c("A", "U", "G", "C")
Endchanges <- cbind(Endchanges, deltaContent[,-1])

Endchanges_gather <- Endchanges[,c(2,8:11)] %>% gather(Nucleotide, Count, -deltaRRS)
write.csv(Endchanges, "Isoform_endchanges_6.csv", quote = F, row.names = F)

Endchanges_gather$Count[Endchanges_gather$Count <= -2] <- -2
Endchanges_gather$Count[Endchanges_gather$Count >= 2] <- 2
Endchanges_gather$Count <- factor(Endchanges_gather$Count, levels = -2:2, labels = c("≤-2", -1:1, "≥2"))
Endchanges_gather <- Endchanges_gather %>% group_by(Nucleotide, Count) %>% summarise(meandeltaRRS = mean(deltaRRS), se = sd(deltaRRS)/sqrt(n()))
library(ggplot2)
ggplot(Endchanges_gather, aes(Count,meandeltaRRS, fill=Nucleotide))+
  geom_bar(stat="identity", position=position_dodge(0.8), width = 0.8)+theme_classic()+xlab("Changes in nucleotide counts")+
  geom_errorbar(aes(ymin = meandeltaRRS - se, ymax = meandeltaRRS + se), position = position_dodge(0.8), width = 0.4)
ggsave("Rebuttal/Isoform_5endprox_ntchanges.pdf", width = 6, height = 3)

for(i in c("A", "U", "G", "C")){
  Subdata <- Endchanges_gather[Endchanges_gather$Nucleotide == i,]
  ggplot(Subdata, aes(Count,meandeltaRRS, fill=Nucleotide))+
    geom_bar(stat="identity", position=position_dodge(0.8), width = 0.8)+theme_classic()+
    geom_errorbar(aes(ymin = meandeltaRRS - se, ymax = meandeltaRRS + se), position = position_dodge(0.8), width = 0.4)+
    ylim(-0.5, 0.5)
  ggsave(paste0("Rebuttal/Isoform_5endprox_ntchanges_", i, "_.pdf"), width = 3, height = 3)
}

