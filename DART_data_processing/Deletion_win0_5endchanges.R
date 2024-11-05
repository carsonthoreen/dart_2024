RRSdata <- read.csv("meanRRSwithSequence.csv")
RRSdata$TransID <- sapply(RRSdata$ID, function(x) strsplit(x, "\\.")[[1]][1])
RRSdata_parent <- RRSdata[grep("parent", RRSdata$ID),]
RRSdata_parent <- RRSdata_parent %>% group_by(TransID) %>% summarise(RRS_m1Y = mean(RRS_m1Y), RRS_Unmod = mean(RRS_Unmod), UTR = unique(UTR))

RRSdata_win <- RRSdata[grep("win", RRSdata$ID),]
RRSdata_win$type <- sapply(RRSdata_win$ID, function(x) as.numeric(str_remove(strsplit(x, "\\.")[[1]][2], pattern = "win")))
RRSdata_win <- RRSdata_win[RRSdata_win$TransID %in% RRSdata_parent$TransID,]
RRSdata_win <- RRSdata_win[RRSdata_win$type== 0,]

for(i in 1:NROW(RRSdata_win)){
  transID = RRSdata_win$TransID[i]
  location = RRSdata_win$type[i]
  RRSdata_win$delSeq[i] = str_sub(RRSdata_parent$UTR[RRSdata_parent$TransID == transID],location+4, location+9)
  RRSdata_win$delRRS_m1Y[i] = RRSdata_win$RRS_m1Y[i]-RRSdata_parent$RRS_m1Y[RRSdata_parent$TransID == transID]
  RRSdata_win$delRRS_Unmod[i] = RRSdata_win$RRS_Unmod[i]-RRSdata_parent$RRS_Unmod[RRSdata_parent$TransID == transID]
  RRSdata_win$deltaA[i] = str_count(pattern = "A", str_sub(RRSdata_win$UTR[i],location+4, location+9)) - str_count(pattern = "A", RRSdata_win$delSeq[i])
  RRSdata_win$deltaT[i] = str_count(pattern = "T", str_sub(RRSdata_win$UTR[i],location+4, location+9)) - str_count(pattern = "T", RRSdata_win$delSeq[i])
  RRSdata_win$deltaG[i] = str_count(pattern = "G", str_sub(RRSdata_win$UTR[i],location+4, location+9)) - str_count(pattern = "G", RRSdata_win$delSeq[i])
  RRSdata_win$deltaC[i] = str_count(pattern = "C", str_sub(RRSdata_win$UTR[i],location+4, location+9)) - str_count(pattern = "C", RRSdata_win$delSeq[i])
}
write.csv(RRSdata_win, "Win0_delta.csv", quote = F, row.names = F)

DeltaRRS <- RRSdata_win[,11:16] %>% gather(Nucleotide, Count, -c(delRRS_m1Y, delRRS_Unmod))
DeltaRRS <- DeltaRRS %>% mutate(Count=ifelse(Count < -2, -2, Count))
DeltaRRS <- DeltaRRS %>% mutate(Count=ifelse(Count > 2, 2, Count))

sres <- DeltaRRS %>% group_by(Nucleotide, Count) %>% summarise(mRRS_Unmod = mean(delRRS_Unmod), seRRS_Unmod=sd(delRRS_Unmod)/sqrt(n()), 
                                                               mRRS_m1Y = mean(delRRS_m1Y), seRRS_m1Y=sd(delRRS_m1Y)/sqrt(n()),
                                                               n = n())

sres$Nucleotide <- str_sub(string = sres$Nucleotide, start = 6)
for( nt in c('A','T','C','G') ) { 
  
  plot.data <- sres %>% filter(Nucleotide==nt)
  print(plot.data)
  p1 <- ggplot( plot.data, aes(x=Count, y=mRRS_Unmod) ) + geom_bar(stat='identity') +
    theme_classic() +
    geom_errorbar( aes(ymin=mRRS_Unmod-seRRS_Unmod, ymax=mRRS_Unmod+seRRS_Unmod, width=0.2) ) +
    ylim(-0.6,0.6)+theme(text = element_text(size = 8))
  ggsave(p1, file=paste0(nt,'_Unmod_barplot.pdf'), width = 1, height = 2.2)
  p2 <- ggplot( plot.data, aes(x=Count, y=mRRS_m1Y) ) + geom_bar(stat='identity') +
    theme_classic() +
    geom_errorbar( aes(ymin=mRRS_m1Y-seRRS_m1Y, ymax=mRRS_m1Y+seRRS_m1Y, width=0.2) ) +
    ylim(-0.6,0.8)+theme(text = element_text(size = 8))
  ggsave(p2, file=paste0(nt,'_m1Y_barplot.pdf'), width = 1, height = 2.2)
  
}




