library(devtools)
library(plyr)
library(ggpubr)
library(dplyr)
library(ggsignif)
library(ggplot2)
library(reshape2)
library(scales)


methods <-c('XGB', 'RF', 'LR','LR.GB', 'meta_MLP')

tissues <- c('cortex', 'cerebellum', 'brain_other')


final_df <-  read.csv( "../data/figureS9/tissues_pred_proba/Y_pred_prob_cerebellum_causal.csv", row.names = 1)
final_df <-  final_df[0:0]

for (tissue in tissues) {
  
  
  file_name <- sub("%s", tissue, "../data/figureS9/tissues_pred_proba/Y_pred_prob_%s_causal.csv")
  
  mydata <- read.csv(file_name, row.names = 1)
  
  mydata <- mydata[1:6]
  cols <- colnames(mydata)
  cols <- paste(cols, tissue, sep = "_")
  colnames(mydata) <- cols

  final_df <- cbind(final_df, mydata)


}

write.csv(final_df, file = "../data/figureS9/final_trace_scores.csv")
  
  
  
  





