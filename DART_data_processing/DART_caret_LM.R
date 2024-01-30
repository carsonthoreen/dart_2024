argv <- commandArgs(T)
CountRRS <- read.csv(argv[1])
dir <- getwd()
dir.create(paste0(dir, "/LinerModel/"))


#CountRRS <- read.csv("RFfiles/ML_unmod_table.csv/ML_unmod_table_v2.csv")
colnames(CountRRS)[2] <- "RRS"
CountRRS <- CountRRS[CountRRS$RRS != -Inf,]
#CountRRS <- CountRRS[grep("GFP", CountRRS$ID),]

ML.data = CountRRS[,-1]
library(dplyr)
ML.data <- select_if(ML.data, is.numeric)


library(caret)
filtered_ML <- ML.data
set.seed(100)
samp <- createDataPartition(filtered_ML$RRS, p=0.8, list = F)
train_set <- filtered_ML[samp,]
test_set <- filtered_ML[-samp,]

Model_lm <- caret::train(RRS~., data = train_set, method = "lm", 
                             trControl = tr)
Model_prediction <- predict(Model_lm, test_set)
Predict_score <- postResample(Model_prediction, test_set$RRS)
write.table(Model_lm$results, paste0(dir, "/LinerModel/ModelResults.txt"))
write.table(Predict_score, paste0(dir, "/LinerModel/Predict_score.txt"))
saveRDS(Model_lm, paste0(dir, "/LinerModel/LinerModel.rds"))

library(ggplot2)
library(ggpointdensity)
library(viridis)
Prediction <- data.frame(Observation = test_set$RRS, Prediction = Model_prediction)
ggplot(Prediction, aes(Observation, Prediction))+
  geom_pointdensity(size=0.1)+scale_color_viridis()+theme_classic()+
  geom_smooth(method = "lm", color="red", se=F)+
  xlim(-6,5)+ylim(-6,5)+theme(legend.position = "none", text = element_text(size = 8))
ggsave("analysis/LinearModel_prediction1.pdf", width = 2.8, height = 2.8)

