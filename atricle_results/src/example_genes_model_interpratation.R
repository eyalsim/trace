library(devtools)
library(plyr)
library(ggpubr)
library(dplyr)
library(ggsignif)
library(ggplot2)
library(reshape2)
library(scales)



mydata <- read.csv('../dataset/df_complete_dataset.csv', 
                   header = TRUE, row.names = 1, sep = ",")


###################################################################################################################

# CACNA1C - Heart

gene_ensg_name <- 'ENSG00000151067'


feature_name <- 'diff_proc_act_med'


features <- c('Whole.Blood_feature',
              "Brain1_feature",
              'Brain0_feature',
              'Heart...Atrial.Appendage_feature',
              "Heart...Left.Ventricle_feature",
              'Liver_feature', 
              'Muscle...Skeletal_feature',
              'Nerve...Tibial_feature',
              'Skin...Not.Sun.Exposed..Suprapubic._feature',
              'Testis_feature' 
)

features <- gsub(x = features, pattern = "feature", replacement = feature_name)  

mydata_filtered <- mydata %>% select(features)
mydata_filtered <- mydata_filtered[gene_ensg_name, ]

mydata_filtered <- melt(mydata_filtered)


p <- ggplot(data=mydata_filtered, aes(x=variable, y=value)) + 
  geom_bar(stat = "identity", fill=c("white", "white", "white", "blueviolet", "blueviolet",
                                     "white", "white", "white", "white", "white"),
           color=c("blueviolet", "blueviolet", "blueviolet", "blueviolet", "blueviolet",
                   "blueviolet", "blueviolet", "blueviolet", "blueviolet", "blueviolet")) +
  theme_classic() + 
  theme(legend.position="none") +
  theme(
    axis.title.x=element_blank(),
    axis.title.y = element_text(color="black", size=30, face="bold"),
    axis.text.x = element_text(color="black", size = 30, face="bold", hjust = 1, angle = 45),
    axis.text.y = element_text(color="black", size = 30, face="bold")
  ) +
  theme(plot.margin=unit(c(1,1,1,1),"cm")) +
  ylab("CACNA1C differential\nprocess activity\n") +
  
  scale_x_discrete(labels=c(
    "Brain1_diff_proc_act_med" = 'Brain - Cerebellum',
    "Brain0_diff_proc_act_med" = 'Brain - Cortex',
    "Heart...Atrial.Appendage_diff_proc_act_med" = 'Heart - Atrial appendage',
    "Heart...Left.Ventricle_diff_proc_act_med" = 'Heart - Left vntricle',
    "Skin...Not.Sun.Exposed..Suprapubic._diff_proc_act_med" = "Skin",
    "Muscle...Skeletal_diff_proc_act_med" = "Muscle",
    "Testis_diff_proc_act_med" = "Testis",
    "Whole.Blood_diff_proc_act_med" = "Blood",
    "Nerve...Tibial_diff_proc_act_med" = "Nerve",
    "Liver_diff_proc_act_med" = "Liver"))

p


file_name <- sub("%s1", feature_name, "../figure2/panelA/example_barplot_%s1_%s2.png")
file_name <- sub("%s2", gene_ensg_name, file_name)

png(filename=file_name, width = 700, height = 700)
p
dev.off()



###################################################################################################################


# DMD - muscle

gene_ensg_name <- 'ENSG00000198947'

feature_name <- 'diff_proc_act_max'



features <- c('Whole.Blood_feature',
              "Brain1_feature",
              'Brain0_feature',
              'Heart...Atrial.Appendage_feature',
              "Heart...Left.Ventricle_feature",
              'Liver_feature', 
              'Muscle...Skeletal_feature',
              'Nerve...Tibial_feature',
              'Skin...Not.Sun.Exposed..Suprapubic._feature',
              'Testis_feature' 
)

features <- gsub(x = features, pattern = "feature", replacement = feature_name)  

mydata_filtered <- mydata %>% select(features)
mydata_filtered <- mydata_filtered[gene_ensg_name, ]



mydata_filtered <- melt(mydata_filtered)


p <- ggplot(data=mydata_filtered, aes(x=variable, y=value)) + 
  geom_bar(stat = "identity", fill=c("white", "white", "white", "white", "white",
                                     "white", "blueviolet", "white", "white", "white"),
           color=c("blueviolet", "blueviolet", "blueviolet", "blueviolet", "blueviolet",
                   "blueviolet", "blueviolet", "blueviolet", "blueviolet", "blueviolet")) +
  theme_classic() + 
  theme(legend.position="none") +
  theme(
    axis.title.x=element_blank(),
    axis.title.y = element_text(color="black", size=30, face="bold"),
    axis.text.x = element_text(color="black", size = 30, face="bold", hjust = 1, angle = 45),
    axis.text.y = element_text(color="black", size = 30, face="bold")
  ) +
  theme(plot.margin=unit(c(1,1,1,1),"cm")) +
  ylab("DMD differential\nprocess activity\n") +
  
  scale_x_discrete(labels=c(
    "Brain1_diff_proc_act_max" = 'Brain - Cerebellum',
    "Brain0_diff_proc_act_max" = 'Brain - Cortex',
    "Heart...Atrial.Appendage_diff_proc_act_max" = 'Heart - Atrial appendage',
    "Heart...Left.Ventricle_diff_proc_act_max" = 'Heart - Left vntricle',
    "Skin...Not.Sun.Exposed..Suprapubic._diff_proc_act_max" = "Skin",
    "Muscle...Skeletal_diff_proc_act_max" = "Muscle",
    "Testis_diff_proc_act_max" = "Testis",
    "Whole.Blood_diff_proc_act_max" = "Blood",
    "Nerve...Tibial_diff_proc_act_max" = "Nerve",
    "Liver_diff_proc_act_max" = "Liver"))

p


file_name <- sub("%s1", feature_name, "../figure2/panelB/example_barplot_%s1_%s2.png")
file_name <- sub("%s2", gene_ensg_name, file_name)

png(filename=file_name, width = 700, height = 700)
p
dev.off()


