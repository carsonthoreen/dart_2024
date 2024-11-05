library(Biostrings)
library(stringr)
args <- commandArgs(T)
K <- args[1]

#K=4

SequenceFile <- readDNAStringSet("Reference/New_Var_reduce_3.fa")
Sequences <- data.frame(ID = names(SequenceFile), Sequence = paste(subseq(SequenceFile, start =20, end = 282)))
Sequences$UTR <- sapply(base::strsplit(Sequences$Sequence, "ATG"), "[", 1)
Sequences$UTR <- str_sub(Sequences$UTR, start = 4)
Sequences$Length <- nchar(Sequences$UTR)-3
Sequences <- Sequences[Sequences$Length >= 10,]

KmerFreq_UTR <- t(sapply(Sequences$UTR, function(x){oligonucleotideFrequency(DNAString(x),K)}))
Sequences_UTR <- cbind(Sequences, KmerFreq_UTR)
write.csv(Sequences_UTR, "PoolKmerFreq.csv", quote = F, row.names = F)

KmerFreq_FL <- t(sapply(Sequences$Sequence, function(x){oligonucleotideFrequency(DNAString(x),K)}))
Sequences_FL <- cbind(Sequences, KmerFreq_FL)
write.csv(Sequences_FL, "PoolKmerFreq_FL.csv", quote = F, row.names = F)
