library(ggpubr)
library(ggplot2)
library(ggridges)
library("FactoMineR")
library("factoextra")
library('hexbin')



tissue <- 'heart'
tissue <- 'whole_brain'
tissue <- 'muscle_skeletal'
tissue <- 'skin'
tissue <- 'liver'
tissue <- 'whole_blood'
tissue <- 'testis'
tissue <- 'nerve_tibial'


file_name <- sub("%s", tissue, "../../data/tissues_pred_proba/causality_probabilities_%s_causal.csv")
mydata <- read.csv(file_name, row.names = 1)
mydata$Class <- factor(mydata$Class, levels = c("Non-causal","Causal in other tissues", "Causal in tissue"))



# Density ridgeline plots


##############################################################
##############################################################
########   if log scale   ####################################

mydata$meta_MLP <- mydata$meta_MLP + 0.01

##############################################################
##############################################################


p <- ggplot(mydata, aes(x = meta_MLP, y = Class)) +
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


file_name1 <- sub("%s", tissue, "../figureS7/panelB/trace_ridge_density_%s.png")

png(filename=file_name1, width = 1000, height = 700)
p
dev.off()



