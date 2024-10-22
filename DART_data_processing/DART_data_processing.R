Files <- list.files(pattern = ".txt")
Firstfile <- read.table(Files[1], header = F)

CountTable <- data.frame(row.names = Firstfile$V1)

for(i in Files){
  CurrFile <- read.table(i, header = F)[,8]
  CountTable <- data.frame(CountTable, CurrFile)
}

Filenames <- c(paste(rep(rep(c("m1y","Unmod"), each = 4),15), c(1:4), 
                     c(rep(rep(paste0("tube",1:6),each=8),2), rep(paste0("tube",1:3),each=8)),
                     c(rep(c("80S", "Input"), each=48), rep("IVT", 24)), sep = "-"))
colnames(CountTable) <- Filenames
FLUCnames <- CountTable[grep(row.names(CountTable), pattern = "10.FLUC"),]
CountTable <- as.matrix(CountTable[!row.names(CountTable) %in% row.names(FLUCnames),])


library(matrixStats)
#Combining counts from all sublibraries for TPM calculation
SumBarcodes <- sapply(1:(ncol(CountTable)/4), function(x) rowSums(CountTable[, ((x * 4) - 4 + 1):(x * 4)]))
colnames(SumBarcodes) <- paste(rep(c("m1Y", "Unmod"),15), c(rep(1:6, each = 2), rep(1:6, each = 2), rep(1:3, each = 2)), 
                               c(rep(c("80S", "Input"),each = 12), rep("IVT", 6)), sep = "-")
sum <- colSums(SumBarcodes)
sum <- sapply(1:(length(sum)/2), function(x) sum(sum[((x * 2) - 2 + 1):(x * 2)]))
SumBarcodes <- SumBarcodes[rowMins(SumBarcodes) > 10,]

write.csv(SumBarcodes, "analysis/SumBarcodes.counts.csv", quote = F)

library(dplyr)
RPMTable <- t(t(SumBarcodes)/rep(sum, each=2)*1000000)

IVT <- (RPMTable[,25:26]+RPMTable[,27:28]+RPMTable[,29:30])/3
RRS <- log2(RPMTable[,1:12]/rep(IVT, 6))
Stability <- log2(RPMTable[,13:24]/rep(IVT, 6))

write.csv(RRS, "analysis/RRS.csv", quote = F)
write.csv(Stability, "analysis/Stability.csv", quote = F)

meanRRS <- data.frame(RRS_m1Y = rowMeans(RRS[,seq(1,12,2)]), RRS_Unmod = rowMeans(RRS[,seq(2,12,2)]))
meanStab <-data.frame(Stability_m1Y = rowMeans(Stability[,seq(1,12,2)]), Stability_Unmod = rowMeans(Stability[,seq(2,12,2)]))

write.csv(meanRRS, "analysis/meanRRS.csv", quote = F)
write.csv(meanStab, "analysis/meanStability.csv", quote = F)



