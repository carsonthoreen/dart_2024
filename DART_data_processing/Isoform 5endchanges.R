library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(ggpointdensity)


IsoformTable <- read.table("analysis/RRS_isoform.tsv", sep="\t", header = T)
IsoformTable$Seq1 <- str_sub(IsoformTable$Seq1, start = 4)
IsoformTable$Seq2 <- str_sub(IsoformTable$Seq2, start = 4)

SeqAlign <- function(Seq1, Seq2){
  match_pattern <- c(attr(adist(Seq1, Seq2, counts = TRUE), "trafos"))
  pattern_locate <- str_locate_all(match_pattern, "M+")[[1]] #Select the isoform pairs with only 5' extension.
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
ggsave("analysis/Isoform_extention.pdf")
IsoformTable_5ext <- IsoformTable[IsoformTable$Endextension,]
IsoformTable_5ext_sig_summary <- IsoformTable_5ext %>% group_by(padj < 0.05) %>% summarise(Count = n())
ggplot(IsoformTable_5ext_sig_summary, aes(x = "", y = Count, fill = `padj < 0.05`)) +
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+blank_theme +
  theme(axis.text.x=element_blank())+scale_fill_manual(values=c("#E69F00", "#56B4E9"))
ggsave("analysis/Isoform_ext_sig.pdf")

IsoformTable_5ext <- IsoformTable[IsoformTable$Endextension,]
write.table(IsoformTable_5ext, "analysis/Isoform_5ext.tsv", quote = F, row.names = F, sep = "\t")


#Endchanges <- IsoformTable[IsoformTable$Endextension,c(1,3,8:10,12,13)]
Endchanges <- IsoformTable_5ext[, c(1,3,8:10,12,13)]
deltaContent <- data.frame(ID = Endchanges$geneID)
for(i in 1:NROW(deltaContent)){
  Seq1 = str_sub(Endchanges$Seq1[i], end = 6)
  Seq2 = str_sub(Endchanges$Seq2[i], end = 6)
  deltaA = str_count(Seq1, "A") - str_count(Seq2, "A")
  deltaU = str_count(Seq1, "T") - str_count(Seq2, "T")
  deltaG = str_count(Seq1, "G") - str_count(Seq2, "G")
  deltaC = str_count(Seq1, "C") - str_count(Seq2, "C")
  deltaContent[i, 2:5] = c(deltaA, deltaU, deltaG, deltaC)
}
colnames(deltaContent)[-1] <- c("A", "U", "G", "C")
Endchanges <- cbind(Endchanges, deltaContent[,-1])

Endchanges_gather <- Endchanges[,c(2,8:11)] %>% gather(Nucleotide, Count, -deltaRRS)
write.csv(Endchanges, "analysis/Isoform_endchanges_6.csv", quote = F, row.names = F)

Endchanges_gather$Count[Endchanges_gather$Count <= -2] <- -2
Endchanges_gather$Count[Endchanges_gather$Count >= 2] <- 2
Endchanges_gather$Count <- factor(Endchanges_gather$Count, levels = -2:2, labels = c("≤-2", -1:1, "≥2"))
Endchanges_gather <- Endchanges_gather %>% group_by(Nucleotide, Count) %>% summarise(meandeltaRRS = mean(deltaRRS), se = sd(deltaRRS)/sqrt(n()))
library(ggplot2)
for(i in c("A", "U", "G", "C")){
  Subdata <- Endchanges_gather[Endchanges_gather$Nucleotide == i,]
  ggplot(Subdata, aes(Count,meandeltaRRS, fill=Nucleotide))+
    geom_bar(stat="identity", position=position_dodge(0.8), width = 0.8)+theme_classic()+
    geom_errorbar(aes(ymin = meandeltaRRS - se, ymax = meandeltaRRS + se), position = position_dodge(0.8), width = 0.4)+
    ylim(-0.5, 0.5)
  ggsave(paste0("analysis/Isoform_5endprox_ntchanges_", i, "_.pdf"), width = 3, height = 3)
}

