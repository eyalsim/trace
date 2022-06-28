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


pref_expression <- my_dataset[ , grepl("_preferential_expression" , names(my_dataset))]

pref_expression['max_pref'] <- apply(pref_expression, 1, max)

not_pref_expressed <- pref_expression[pref_expression$max_pref < 2, ]

not_pref_expressed <- not_pref_expressed[complete.cases(not_pref_expressed), ]



tissue <- 'heart'
tissue <- 'whole_brain'
tissue <- 'muscle_skeletal'
tissue <- 'skin'
tissue <- 'liver'
tissue <- 'whole_blood'
tissue <- 'testis'
tissue <- 'nerve_tibial'


file_name <- "../data/causality_probabilities_%s_causal.csv"
file_name <- sub("%s", tissue, file_name)
mydata <- read.csv(file_name, row.names = 1)
mydata$Class <- factor(mydata$Class, levels = c("Non-causal","Causal in other tissues", "Causal in tissue"))

mydata_not_pref_expressed <- subset(mydata, rownames(mydata) %in% rownames(not_pref_expressed))



# Density ridgeline plots


##############################################################
##############################################################
########   if log scale   ####################################

mydata_not_pref_expressed$meta_MLP <- mydata_not_pref_expressed$meta_MLP + 0.01

##############################################################
##############################################################


p <- ggplot(mydata_not_pref_expressed, aes(x = meta_MLP, y = Class)) +
  geom_density_ridges(aes(fill = Class))+
  scale_fill_manual(values=c( "#999999", "#56B4E9", "red3"))+
  theme_classic()  + 
  theme(legend.text = element_text(color="black", size=18),
        legend.title = element_text(color="white"),
        legend.key.size = unit(3, 'lines')) +  
  theme(
    axis.title.x = element_text(color="black", size=40, face="bold"),
    axis.title.y = element_text(color="black", size=40, face="bold"),
    axis.text.x = element_text(color="black", size = 40, face="bold", hjust = 0.5),
    axis.text.y = element_text(color="black", size = 40, face="bold")
  )  +
  scale_y_discrete(labels=c("Non-causal" = "Non-causal ", "Causal in tissue" = "Causal in \ntissue ",
                            "Causal in other tissues" = "Causal ")) +
  theme(legend.position="none")+
  labs(y = "\n", x = "\nTRACE")  + 
  theme(plot.margin=unit(c(1,1,1,1),"cm")) +
  
  # scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10), limits = c(-1,11))
  
  ##############################################################
##############################################################
########   if log scale   ####################################

scale_x_continuous(trans = log2_trans(),
                   breaks = trans_breaks("log2", function(x) 2^x),
                   labels = trans_format("log2", math_format(2^.x)),
                   limits = c(-1,11))

##############################################################
##############################################################

p


file_name1 <- sub("%s", tissue, "../results/not_pref_expressed_trace_ridge_density_%s.png")

png(filename=file_name1, width = 1000, height = 700)
p
dev.off()





p <- ggboxplot(mydata_not_pref_expressed, x = "Class", y = method,
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
  
  # scale_y_continuous(trans = log2_trans(),
  #                    breaks = trans_breaks("log2", function(x) 2^x),
  #                    labels = trans_format("log2", math_format(2^.x)))
  
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10), limits = c(0,14))


p + geom_signif(comparisons = list(c("Non-causal","Causal in tissue"), c("Causal in other tissues","Causal in tissue")), 
                map_signif_level=FALSE, size = 2, textsize = 10, y_position=c(13, 11))


file_name <- sub("%s", tissue, "../results/not_pref_expressed_trace_boxplot_%s.png")

png(filename=file_name, width = 800, height = 700)
p + geom_signif(comparisons = list(c("Non-causal","Causal in tissue"), c("Causal in other tissues","Causal in tissue")), 
                map_signif_level=FALSE, size = 2, textsize = 10, y_position=c(13, 11))
dev.off()


