DataRRS <- read.csv("analysis/meanRRSwithSequence.csv")
DataRRS$TransID <- gsub("\\..*$", "", DataRRS$ID)

Data_parent <- DataRRS[grep("parent", DataRRS$ID),]
library(dplyr)
Data_parent <- Data_parent %>% group_by(TransID) %>% summarise(RRS_m1Y = mean(RRS_m1Y), RRS_Unmod = mean(RRS_Unmod), UTR = unique(UTR))

Data_win <- DataRRS[grep("win", DataRRS$ID),]
Data_win$type <- sapply(Data_win$ID, function(x) as.numeric(str_remove(strsplit(x, "\\.")[[1]][2], pattern = "win")))
Data_win <- Data_win[Data_win$TransID %in% Data_parent$TransID,]

for(i in 1:NROW(Data_win)){
  transID = Data_win$TransID[i]
  location = Data_win$type[i]
  Data_win$delSeq[i] = str_sub(Data_parent$UTR[Data_parent$TransID == transID],location+4, location+9)
  Data_win$gainSeq[i] = str_sub(Data_win$Sequence[i], start = location+1, end = location+6)
  Data_win$delRRS_m1Y[i] = Data_win$RRS_m1Y[i]-Data_parent$RRS_m1Y[Data_parent$TransID == transID]
  Data_win$delRRS_Unmod[i] = Data_win$RRS_Unmod[i]-Data_parent$RRS_Unmod[Data_parent$TransID == transID]
}
write.csv(Data_win, "analysis/Deletion_delSeq.csv", quote = F, row.names = F)

TetromerList <- expand.grid(rep(list(c("A", "T", "G", "C")), 4))
TetromerList <- Reduce(paste0, TetromerList)

Tetromer <- data.frame(Seq = TetromerList, pval_m1Y = 1, pval_Unmod = 1)
for(i in 1:NROW(Tetromer)){
  Motif = Tetromer$Seq[i]
  DelTable <- na.omit(Data_win[str_count(Data_win$delSeq, pattern = Motif)>0,])
  Tetromer$Count[i] <- NROW(DelTable)
  Tetromer$meanRRS_m1Y[i] <- mean(DelTable$delRRS_m1Y)
  Tetromer$meanRRS_Unmod[i] <- mean(DelTable$delRRS_Unmod)
  if(NROW(DelTable) > 10) {
    Tetromer$pval_m1Y[i] <- wilcox.test(DelTable$delRRS_m1Y)$p.value
    Tetromer$pval_Unmod[i] <- wilcox.test(DelTable$delRRS_Unmod)$p.value
  }
}
Tetromer <- Tetromer[Tetromer$Count >= 30,]
Tetromer$padj_m1Y <- p.adjust(Tetromer$pval_m1Y, method = "fdr", n = NROW(Tetromer))
Tetromer$padj_Unmod <- p.adjust(Tetromer$pval_Unmod, method = "fdr", n = NROW(Tetromer))
Tetromer$Seq <- str_replace_all(Tetromer$Seq, pattern = "T", replacement = "U")
write.csv(Tetromer, "analysis/Del_tetromer.csv")