library(devtools)
library(plyr)
library(ggpubr)
library(dplyr)
library(ggsignif)
library(ggplot2)
library(reshape2)
library(scales)


results <- read.csv('../data/figureS3/summary_tissues_xgb_auc.csv')


results$Tissue <- factor(results$Tissue, levels = c("whole_blood","whole_brain", "heart", "liver",
                                                        "muscle_skeletal", "nerve_tibial","skin", "testis"))


p <- ggplot(results, aes(x=Tissue, y=XGB)) +
  geom_bar(
    aes(color = "darkgreen", fill = "floralwhite"),
    stat = "identity", position = position_dodge(0.8),
    width = 0.7
  ) +
  labs(y = "Average AUC\n") +
  theme_classic() + 
  theme(legend.position = "none") +
  theme(
    axis.title.x=element_blank(),
    axis.title.y =element_text(color="black", size = 30, face="bold"),
    axis.text.x = element_text(color="black", size = 30, face="bold", hjust = 1, angle = 45),
    axis.text.y = element_text(color="black", size = 30, face="bold")
  ) +
  theme(plot.margin=unit(c(1,1,1,1),"cm")) +
  scale_y_continuous(breaks = c(0.5, 0.6, 0.7, 0.8, 0.9), limits = c(0.5,0.9), oob = rescale_none) +
  
  scale_fill_manual(values = c("floralwhite", "floralwhite", "floralwhite", "floralwhite", "floralwhite", "darkred")) +
  scale_color_manual(values = c("darkgreen", "deepskyblue", "deeppink", "gold", "blueviolet", "darkred"))  +
  scale_x_discrete(labels=c(
    "whole_brain" = 'Brain',
    "heart" = 'Heart',
    "skin" = "Skin",
    "muscle_skeletal" = "Muscle",
    "testis" = "Testis",
    "whole_blood" = "Blood",
    "nerve_tibial" = "Nerve",
    "liver" = "Liver")) +
  geom_errorbar(aes(ymin=XGB-SD, ymax=XGB+SD), width=.2,
                position=position_dodge(.9))


p



png(filename="../figureS3/panelA/all_tissues_xgb_auc_barplot.png", width = 1000, height = 700)
p
dev.off()

