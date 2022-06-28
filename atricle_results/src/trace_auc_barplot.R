library(devtools)
library(plyr)
library(ggpubr)
library(dplyr)
library(ggsignif)
library(ggplot2)
library(reshape2)
library(scales)


results <- read.csv('../data/figure4/methods_comparison_tissues_auc.csv')

names(results) <- gsub(x = names(results), pattern = "_causal", replacement = "")  
results <- melt(results)

results$variable <- factor(results$variable, levels = c("whole_blood","whole_brain", "heart", "liver",
                                                        "muscle_skeletal", "nerve_tibial","skin", "testis"))
results <- results[complete.cases(results), ]


p <- ggplot(data=results, aes(x=variable, y=value, fill=X)) + 
  geom_bar(stat = "identity", position=position_dodge(0.8),
           width = 0.7) + 
  theme_classic() + 
  theme(legend.title = element_text(color="black", size = 30, face="bold"),
        legend.text = element_text(color="black", size = 20, face="bold")) +
  labs(y = "Average AUC\n", fill = "") +
  scale_x_discrete(labels=c(
    "whole_brain" = 'Brain',
    "heart" = 'Heart',
    "skin" = "Skin",
    "muscle_skeletal" = "Muscle",
    "testis" = "Testis",
    "whole_blood" = "Blood",
    "nerve_tibial" = "Nerve",
    "liver" = "Liver")) +
  theme(
    axis.title.x=element_blank(),
    axis.title.y =element_text(color="black", size = 30, face="bold"),
    axis.text.x = element_text(color="black", size = 30, face="bold", hjust = 1, angle = 45),
    axis.text.y = element_text(color="black", size = 30, face="bold")
  ) +
  theme(plot.margin=unit(c(1,1,1,1),"cm")) +
  scale_y_continuous(breaks = c(0.5, 0.6, 0.7, 0.8, 0.9), limits = c(0.5,0.9), oob = rescale_none) +
  scale_fill_manual(values = c("grey", "darkred")) +
  scale_color_manual(values = c("grey", "darkred")) 

p


png(filename="../figure4/panelC/auc_trace_barplot_pbrit.png", width = 1000, height = 700)
p
dev.off()

###################################################################################################################


results <- read.csv('../data/figureS6/panelB/methods_comparison_tissues_pr_score.csv')

names(results) <- gsub(x = names(results), pattern = "_causal", replacement = "")  
results <- melt(results)

results$variable <- factor(results$variable, levels = c("whole_blood","whole_brain", "heart", "liver",
                                                        "muscle_skeletal", "nerve_tibial","skin", "testis"))
results <- results[complete.cases(results), ]



num_genes = 18927
num_causal_liver = 62
num_causal_nerve = 73
num_causal_blood = 80
num_causal_tetis = 115
num_causal_muscle = 128
num_causal_skin = 147
num_causal_heart = 154
num_causal_brain = 498

p <- ggplot(data=results, aes(x=variable, y=value, fill=X)) + 
  geom_bar(stat = "identity", position=position_dodge(0.8),
           width = 0.7) + 
  theme_classic() + 
  theme(legend.title = element_text(color="black", size = 30, face="bold"),
        legend.text = element_text(color="black", size = 20, face="bold")) +
  labs(y = "Summarized Precision-Recall\n", fill = "") +
  scale_x_discrete(labels=c(
    "whole_brain" = 'Brain',
    "heart" = 'Heart',
    "skin" = "Skin",
    "muscle_skeletal" = "Muscle",
    "testis" = "Testis",
    "whole_blood" = "Blood",
    "nerve_tibial" = "Nerve",
    "liver" = "Liver")) +
  theme(
    axis.title.x=element_blank(),
    axis.title.y =element_text(color="black", size = 30, face="bold"),
    axis.text.x = element_text(color="black", size = 30, face="bold", hjust = 1, angle = 45),
    axis.text.y = element_text(color="black", size = 30, face="bold")
  ) +
  theme(plot.margin=unit(c(1,1,1,1),"cm")) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4), limits = c(0, 0.4), oob = rescale_none) +
  scale_fill_manual(values = c("grey", "darkred")) +
  scale_color_manual(values = c("grey", "darkred"))  + 
  geom_segment(aes(x = 0.5, y = num_causal_blood/(num_genes-num_causal_blood), 
                   xend = 1.5, yend = num_causal_blood/(num_genes-num_causal_blood)), size=1.5) +
  geom_segment(aes(x = 1.5, y = num_causal_brain/(num_genes-num_causal_brain), 
                   xend = 2.5, yend = num_causal_brain/(num_genes-num_causal_brain)), size=1.5) +
  geom_segment(aes(x = 2.5, y = num_causal_heart/(num_genes-num_causal_heart), 
                   xend = 3.5, yend = num_causal_heart/(num_genes-num_causal_heart)), size=1.5) +
  geom_segment(aes(x = 3.5, y = num_causal_liver/(num_genes-num_causal_liver), 
                   xend = 4.5, yend = num_causal_liver/(num_genes-num_causal_liver)), size=1.5) +
  geom_segment(aes(x = 4.5, y = num_causal_muscle/(num_genes-num_causal_muscle), 
                   xend = 5.5, yend = num_causal_muscle/(num_genes-num_causal_muscle)), size=1.5) +
  geom_segment(aes(x = 5.5, y = num_causal_nerve/(num_genes-num_causal_nerve), 
                   xend = 6.5, yend = num_causal_nerve/(num_genes-num_causal_nerve)), size=1.5) +
  geom_segment(aes(x = 6.5, y = num_causal_skin/(num_genes-num_causal_skin), 
                   xend = 7.5, yend = num_causal_skin/(num_genes-num_causal_skin)), size=1.5) +
  geom_segment(aes(x = 7.5, y = num_causal_tetis/(num_genes-num_causal_tetis), 
                   xend = 8.5, yend = num_causal_tetis/(num_genes-num_causal_tetis)), size=1.5)

p


png(filename="../figureS6/panelB/pr_score_trace_barplot_pbrit.png", width = 1000, height = 700)
p
dev.off()
