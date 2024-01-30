KAdata_mean <- read.csv("KAdata_mean.csv")
colnames(KAdata_mean)[1] <- "ID"
KAdata_mean$Type <- sapply(KAdata_mean$ID, function(x) tail(strsplit(x, "\\.")[[1]], 1))
Random <- KAdata_mean[KAdata_mean$Type %in% c("GFP", "FLUC", "Spike"),]
ggplot(Random, aes(x=RRS_Unmod, color=Type))+
  stat_ecdf()+ylab("Propotion")
ggplot(Random, aes(x=RRS_m1Y, color=Type))+
  stat_ecdf()+ylab("Propotion")
ggplot(Random, aes(x=RRS_Y, color=Type))+
  stat_ecdf()+ylab("Propotion")
ggplot(Random, aes(x=RRS_m6A, color=Type))+
  stat_ecdf()+ylab("Propotion")

Random$TransID <- gsub("\\..*$", "", Random$ID)
Random_intersect <- Random %>% group_by(TransID) %>% filter(n() == 3)
GFPTop10 <- Random_intersect[Random_intersect$Type == "GFP",] %>% group_by(Type) %>% slice_max(prop = 0.1, order = RRS_Unmod)
GFPBottom10 <- Random_intersect[Random_intersect$Type == "GFP",] %>% group_by(Type) %>% slice_min(prop = 0.1, order = RRS_Unmod)
SpikeTop10 <- Random_intersect[Random_intersect$Type == "Spike",] %>% group_by(Type) %>% slice_max(prop = 0.1, order = RRS_Unmod)
SpikeBottom10 <- Random_intersect[Random_intersect$Type == "Spike",] %>% group_by(Type) %>% slice_min(prop = 0.1, order = RRS_Unmod)
FLUCTop10 <- Random_intersect[Random_intersect$Type == "FLUC",] %>% group_by(Type) %>% slice_max(prop = 0.1, order = RRS_Unmod)
FLUCBottom10 <- Random_intersect[Random_intersect$Type == "FLUC",] %>% group_by(Type) %>% slice_min(prop = 0.1, order = RRS_Unmod)

TopVenn <- list(GFPTop10$TransID, FLUCTop10$TransID, SpikeTop10$TransID)
BottomVenn <- list(GFPBottom10$TransID, FLUCBottom10$TransID, SpikeBottom10$TransID)
library(ggvenn)
names(TopVenn) <- c("GFP", "FLUC", "Spike")
names(BottomVenn) <- c("GFP", "FLUC", "Spike")

ggvenn(TopVenn)
ggsave("KA_CDS_Venn_top.pdf")
ggvenn(BottomVenn)
ggsave("KA_CDS_Venn_bottom.pdf")

library(tidyr)
Random_intersect_spread <- Random_intersect[,c(2,6,7)] %>% spread(key = Type, value = RRS_Unmod)
ggplot(Random_intersect_spread, aes(GFP, Spike))+
  geom_pointdensity()+scale_color_viridis()+theme_classic()+
  geom_smooth(method = "lm", color="red", se=F)+
  theme(legend.position = "none",text = element_text(size = 8))
ggsave("Random_SpikevsGFP_scatter.pdf", width = 2.8, height = 2.8)
ggplot(Random_intersect_spread, aes(GFP, FLUC))+
  geom_pointdensity()+scale_color_viridis()+theme_classic()+
  geom_smooth(method = "lm", color="red", se=F)+
  theme(legend.position = "none",text = element_text(size = 8))
ggsave("Random_FLUCvsGFP_scatter.pdf", width = 2.8, height = 2.8)
ggplot(Random_intersect_spread, aes(FLUC, Spike))+
  geom_pointdensity()+scale_color_viridis()+theme_classic()+
  geom_smooth(method = "lm", color="red", se=F)+
  theme(legend.position = "none",text = element_text(size = 8))
ggsave("Random_SpikevsFLUC_scatter.pdf", width = 2.8, height = 2.8)

Kozak <- read.table("Kozak_New_longshort.txt", header = F, sep = "\t")
colnames(Kozak) <- c("ID", "Kozak_score")
colnames(Random)[1] <- "ID"
Random <- merge(Random, Kozak, by = "ID")
library(tidyr)
Random <- Random %>% gather(mod, RRS, -c(ID, Type, Kozak_score))
Random_Unmod <- Random[Random$mod == "RRS_Unmod",]
ggplot(Random_Unmod, aes(x=Kozak_score, y=RRS, color = Type))+
         geom_smooth(method = "lm")
