library(devtools)
library(plyr)
library(ggpubr)
library(dplyr)
library(ggsignif)
library(ggplot2)
library(reshape2)
library(scales)
library(easyGgplot2)


mydata <- read.csv('../data/figure2/panelA/pref_expression_membrane depolarization during atrial cardiac muscle cell action potential.csv', 
                   header = TRUE, sep = ",", row.names = 1)


colnames(mydata)


features <- c('Whole.Blood',
              'Brain1',
              "Brain0",
              'Heart...Atrial.Appendage',
              "Heart...Left.Ventricle",
              'Liver', 
              'Muscle...Skeletal',
              'Nerve...Tibial',
              'Skin...Not.Sun.Exposed..Suprapubic.',
              'Testis' 
)



depolarization_dataset <- mydata[, features]


mydata_melted <- melt(depolarization_dataset)


p <- ggplot(data=mydata_melted, aes(x=variable, y=value, fill=variable)) + 
  geom_violin(color="blueviolet") +
  geom_boxplot(alpha = 0.0, size = 0.6, color="blueviolet") +
  theme_classic() + 
  theme(legend.position="none") +
  theme(
    axis.title.x=element_blank(),
    axis.title.y = element_text(color="black", size=30, face="bold"),
    axis.text.x = element_text(color="black", size = 30, face="bold", hjust = 1, angle = 45),
    axis.text.y = element_text(color="black", size = 30, face="bold")
  ) +
  theme(plot.margin=unit(c(1,1,1,1),"cm")) +
  ylab("Membrane depolarization\npreferential expression\n") +
  
  scale_x_discrete(labels=c(
    "Brain0" = 'Brain - Cerebellum',
    "Brain1" = 'Brain - Cortex',
    "Heart...Atrial.Appendage" = 'Heart - Atrial appendage',
    "Heart...Left.Ventricle" = 'Heart - Left vntricle',
    "Skin...Not.Sun.Exposed..Suprapubic." = "Skin",
    "Muscle...Skeletal" = "Muscle",
    "Testis" = "Testis",
    "WholeBlood" = "Blood",
    "Nerve...Tibial" = "Nerve",
    "Liver" = "Liver")) + 
  scale_y_continuous(breaks = c(-8, -4, 0, 4, 8, 12, 16), limits = c(-8, 16)) + 
  scale_fill_manual(values = c("white", "white", "white", "white", "white",
                               "white", "blueviolet", "white", "white", "white"))

p


png(filename='../figure2/panelA/membrane_depolarization_example.png', width = 700, height = 700)
p
dev.off()


############################################################################################################




mydata <- read.csv('../data/figure2/panelB/pref_expression_muscle_filament_sliding.csv', 
                   header = TRUE, sep = ",", row.names = 1)


colnames(mydata)


features <- c('WholeBlood',
              'Brain1',
              "Brain0",
              'Heart.AtrialAppendage',
              "Heart.LeftVentricle",
              'Liver', 
              'Muscle.Skeletal',
              'Nerve.Tibial',
              'Skin.NotSunExposed.Suprapubic.',
              'Testis' 
)



filament_dataset <- mydata[, features]


mydata_melted <- melt(filament_dataset)


p <- ggplot(data=mydata_melted, aes(x=variable, y=value, fill=variable)) + 
  geom_violin(color="blueviolet") +
  geom_boxplot(alpha = 0.0, size = 0.6, color="blueviolet") +
  theme_classic() + 
  theme(legend.position="none") +
  theme(
    axis.title.x=element_blank(),
    axis.title.y = element_text(color="black", size=30, face="bold"),
    axis.text.x = element_text(color="black", size = 30, face="bold", hjust = 1, angle = 45),
    axis.text.y = element_text(color="black", size = 30, face="bold")
  ) +
  theme(plot.margin=unit(c(1,1,1,1),"cm")) +
  ylab("Muscle filament sliding\npreferential expression\n") +
  
  scale_x_discrete(labels=c(
    "Brain0" = 'Brain - Cerebellum',
    "Brain1" = 'Brain - Cortex',
    "Heart.AtrialAppendage" = 'Heart - Atrial appendage',
    "Heart.LeftVentricle" = 'Heart - Left vntricle',
    "Skin.NotSunExposed.Suprapubic." = "Skin",
    "Muscle.Skeletal" = "Muscle",
    "Testis" = "Testis",
    "WholeBlood" = "Blood",
    "Nerve.Tibial" = "Nerve",
    "Liver" = "Liver")) + 
  scale_y_continuous(breaks = c(-8, -4, 0, 4, 8, 12, 16), limits = c(-8, 16)) + 
  scale_fill_manual(values = c("white", "white", "white", "white", "white",
                               "white", "blueviolet", "white", "white", "white"))

p


png(filename='../figure2/panelB/DMD_muscle_filament_sliding_example.png', width = 700, height = 700)
p
dev.off()


