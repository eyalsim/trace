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



expression_and_pref <- my_dataset[ , grepl("_expression" , names(my_dataset))]
pref_expression <- my_dataset[ , grepl("_pref_expression" , names(my_dataset))]

only_expression <- 
  expression_and_pref[ , ! names(expression_and_pref) %in% names(pref_expression)]


only_expression['expressed_ratio'] <-  rowSums(only_expression > 7) / rowSums(only_expression > -1)



widely_expressed <- only_expression[only_expression$expressed_ratio > 0.8, ]

widely_expressed <- widely_expressed[complete.cases(widely_expressed), ]




tissue <- 'heart'
tissue <- 'whole_brain'
tissue <- 'muscle_skeletal'
tissue <- 'skin'
tissue <- 'liver'
tissue <- 'whole_blood'
tissue <- 'testis'
tissue <- 'nerve_tibial'


method <- 'meta_MLP'


file_name <- "../../data/tissues_pred_proba/causality_probabilities_%s_causal.csv"
file_name <- sub("%s", tissue, file_name)
mydata <- read.csv(file_name, row.names = 1)
mydata$Class <- factor(mydata$Class, levels = c("Non-causal","Causal in other tissues", "Causal in tissue"))

mydata_widely_expressed <- subset(mydata, rownames(mydata) %in% rownames(widely_expressed))



# Density ridgeline plots


##############################################################
##############################################################
########   if log scale   ####################################

mydata_widely_expressed$meta_MLP <- mydata_widely_expressed$meta_MLP + 0.01

##############################################################
##############################################################


p <- ggplot(mydata_widely_expressed, aes(x = meta_MLP, y = Class)) +
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


file_name1 <- sub("%s", tissue, "../figureS8/widely_expressed_trace_ridge_density_%s.png")

png(filename=file_name1, width = 1000, height = 700)
p
dev.off()





p <- ggboxplot(mydata_widely_expressed, x = "Class", y = method,
               color = "Class",
               add = "jitter", xlab="\n", ylab = "TRACE\n") +
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


p + geom_signif(comparisons = list(c("Non-causal","Causal in tissue"), c("Causal in other tissues","Causal in tissue")), 
                map_signif_level=FALSE, size = 2, textsize = 10, y_position=c(13, 11))


file_name <- sub("%s", tissue, "../results/widely_expressed_trace_boxplot_%s.png")

png(filename=file_name, width = 800, height = 700)
p + geom_signif(comparisons = list(c("Non-causal","Causal in tissue"), c("Causal in other tissues","Causal in tissue")), 
                map_signif_level=FALSE, size = 2, textsize = 10, y_position=c(13, 11))
dev.off()




causal <- mydata_widely_expressed[which(mydata_widely_expressed$Class=='Causal in tissue'),]$meta_MLP
non_causal <- mydata_widely_expressed[which(mydata_widely_expressed$Class=='Non-causal'),]$meta_MLP
causal_other_tissues <- mydata_widely_expressed[which(mydata_widely_expressed$Class=='Causal in other tissues'), ]$meta_MLP

t1<-wilcox.test(x=causal, 
                y=non_causal, 
                na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)

print(t1$statistic)
print(t1$p.value)

t2<-wilcox.test(x=causal, 
                y=causal_other_tissues, 
                na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)

print(t2$statistic)
print(t2$p.value)



df <- data.frame(compared = c("causal-non_causal", "causal-causal_other_tissues"),
                 statistic = c(t1$statistic, t2$statistic),
                 p_value = c(t1$p.value, t2$p.value)
)

file_name <- sub("%s", tissue, "../figureS8/widely_expressed_trace_mw_%s.csv")
write.csv(df, file_name, row.names = FALSE)

