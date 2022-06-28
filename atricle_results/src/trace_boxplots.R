library(devtools)
library(plyr)
library(ggpubr)
library(dplyr)
library(ggsignif)
library(ggplot2)
library(reshape2)
library(scales)


method <- 'meta_MLP'



tissue <- 'whole_brain'
tissue <- 'heart'
tissue <- 'muscle_skeletal'
tissue <- 'skin'
tissue <- 'liver'
tissue <- 'nerve_tibial'
tissue <- 'testis'
tissue <- 'whole_blood'

tissue <- 'cortex'
tissue <- 'cerebellum'




file_name <- "../../data/tissues_pred_proba/causality_probabilities_%s_causal.csv"
file_name <- sub("%s", tissue, file_name)

mydata <- read.csv(file_name, row.names = 1)
mydata$Class <- factor(mydata$Class, levels = c("Non-causal","Causal in other tissues", "Causal in tissue"))


p <- ggboxplot(mydata, x = "Class", y = method,
               color = "Class",
               add = "jitter", xlab="\n", ylab = "TRACE\n") +
  # scale_color_manual(values=c("cornflowerblue", "firebrick3")) + 
  scale_color_manual(values=c( "#999999", "#56B4E9", "red3")) +
  
  theme(legend.position="NA") +
  theme(
    axis.title.x = element_text(color="black", size=40, face="bold"),
    axis.title.y = element_text(color="black", size=40, face="bold"),
    axis.text.x = element_text(size = 40, face="bold", angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 40, face="bold")
  )  + 
  scale_x_discrete(labels=c("Non-causal" = "Non-causal", "Causal in tissue" = "Causal in\ntissue",
                            "Causal in other tissues" = "Causal")) +
  theme(plot.margin=unit(c(1,1,1,1),"cm")) +
  
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10), limits = c(0,14))


p



p + stat_compare_means(comparisons=list(c("Non-causal","Causal in tissue"), c("Causal in other tissues","Causal in tissue")), 
                       label.x = 1.5, size=2, textsize = 10)



p + geom_signif(comparisons = list(c("Non-causal","Causal in tissue"), c("Causal in other tissues","Causal in tissue")), 
                map_signif_level=FALSE, size = 2, textsize = 10, y_position=c(13, 11))


p<- p + geom_signif(comparisons = list(c("Non-causal","Causal in tissue"), c("Causal in other tissues","Causal in tissue")), 
                    map_signif_level=TRUE, size = 2, textsize = 10, y_position=c(13, 11))


file_name <- sub("%s", tissue, "../figureS7/panelA/probab_compare_boxplot_%s.png")
png(filename=file_name, width = 800, height = 700)
p
dev.off()


p<- p + geom_signif(comparisons = list(c("Non-causal","Causal in tissue"), c("Causal in other tissues","Causal in tissue")), 
                    map_signif_level=TRUE, size = 2, textsize = 10, y_position=c(13, 11))


