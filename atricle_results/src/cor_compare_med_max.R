library(devtools)
library(plyr)
library(ggpubr)
library(dplyr)
library(ggsignif)
library(reshape2)
library(scales)
library(ggplot2)
library(ggridges)
library("FactoMineR")
library("factoextra")
library('hexbin')



my_dataset <- read.csv('../dataset/df_complete_dataset.csv', 
                       header = TRUE, row.names = 1, sep = ",")


###################################################################################################################


p <- ggplot(my_dataset, aes(x=Heart...Atrial.Appendage_median_tipa_pathways, y=Heart...Atrial.Appendage_max_tipa_pathways)) +
  geom_point()+
  geom_smooth(method=lm) +
  theme_classic() + 
  stat_cor(method = "pearson", label.x = 0.5, label.y = 6) +
  theme(
    axis.title.x = element_text(color="black", size=30, face="bold"),
    axis.title.y = element_text(color="black", size=30, face="bold"),
    axis.text.x = element_text(color="black", size=20),
    axis.text.y = element_text(color="black", size=20)
  ) + 
  theme(legend.position="none") +
  labs(x="\nHeart - atrial appendage median TIPA score", y = "Heart - atrial appendage max TIPA score\n") + 
  theme(plot.margin=unit(c(1,1,1,1),"cm"))

p



png(filename="../results/cor_compare_med_max_tipa.png", width = 700, height = 700)
p
dev.off()



p <- ggplot(my_dataset, aes(x=Heart_._Atrial_Appendage_diff_net_med, y=Heart_._Atrial_Appendage_diff_net_max)) +
  geom_point()+
  geom_smooth(method=lm) +
  theme_classic() + 
  stat_cor(method = "pearson", label.x = 0.1, label.y = 0) +
  theme(
    axis.title.x = element_text(color="black", size=30, face="bold"),
    axis.title.y = element_text(color="black", size=30, face="bold"),
    axis.text.x = element_text(color="black", size=20),
    axis.text.y = element_text(color="black", size=20)
  ) + 
  theme(legend.position="none") +
  labs(x="\nHeart - atrial appendage\nmedian differential PPI score ", y = "Heart - atrial appendage\n max differential PPI score \n") + 
  theme(plot.margin=unit(c(1,1,1,1),"cm"))

p



png(filename="../results/cor_compare_med_max_dif_ppi.png", width = 700, height = 700)
p
dev.off()





p <- ggplot(my_dataset, aes(x=Heart...Atrial.Appendage_num_elevated_interactors, y=Heart...Atrial.Appendage_num_elevated_interactors_dif_median)) +
  geom_point()+
  geom_smooth(method=lm) +
  theme_classic() + 
  stat_cor(method = "pearson", label.x = 10, label.y = 80) +
  theme(
    axis.title.x = element_text(color="black", size=30, face="bold"),
    axis.title.y = element_text(color="black", size=30, face="bold"),
    axis.text.x = element_text(color="black", size=20),
    axis.text.y = element_text(color="black", size=20)
  ) +  
  theme(legend.position="none") +
  labs(x="\nHeart - atrial appendage\nnum elevated interactions", 
       y = "Heart - atrial appendage\nnum elevated interactions dif median\n") +
theme(plot.margin=unit(c(1,1,1,1),"cm"))

p



png(filename="../results/cor_compare_dif_elevated.png", width = 700, height = 700)
p
dev.off()



p <- ggplot(my_dataset, aes(x=Heart_._Atrial_Appendage_diff_net_med, y=Heart...Atrial.Appendage_preferential_expression)) +
  geom_point()+
  geom_smooth(method=lm) +
  theme_classic() + 
  stat_cor(method = "pearson", label.x = 0, label.y = 20000) +
  theme(
    axis.title.x = element_text(color="black", size=30, face="bold"),
    axis.title.y = element_text(color="black", size=30, face="bold"),
    axis.text.x = element_text(color="black", size=20),
    axis.text.y = element_text(color="black", size=20)
  ) +  
  theme(legend.position="none") +
  labs(x="\nHeart - atrial appendage\nmedian differential PPI score", 
       y = "Heart - atrial appendage\nprefferential expression\n") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"))

p



png(filename="../results/cor_compare_difnet_pref_expr.png", width = 700, height = 700)
p
dev.off()






p <- ggplot(my_dataset, aes(x=Heart_._Atrial_Appendage_diff_net_med, y=Heart...Atrial.Appendage_num_specific_interactions)) +
  geom_point()+
  geom_smooth(method=lm) +
  theme_classic() + 
  stat_cor(method = "pearson", label.x = 0, label.y = 130) +
  theme(
    axis.title.x = element_text(color="black", size=30, face="bold"),
    axis.title.y = element_text(color="black", size=30, face="bold"),
    axis.text.x = element_text(color="black", size=20),
    axis.text.y = element_text(color="black", size=20)
  ) +  
  theme(legend.position="none") +
  labs(x="\nHeart - atrial appendage\nmedian differential PPI score", 
       y = "Heart - atrial appendage\nnum specific interactions\n") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"))

p



png(filename="../results/cor_compare_difnet_num_specific_inter.png", width = 700, height = 700)
p
dev.off()






my_dataset$Heart...Atrial.Appendage_preferential_expression
my_dataset$Heart...Atrial.Appendage_num_specific_interactions
