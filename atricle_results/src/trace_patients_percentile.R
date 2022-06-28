
library(ggpubr)
library(ggplot2)


mydata <- read.csv('../data/figure5/table_s4_june22.csv', 
                   header = TRUE, sep = ",")


mydata$trace_percentile <- mydata$TRACE_model_rank/mydata$Num.genes.TRACE * 100
mydata$gado_percentile <- mydata$GADO_rank/mydata$Num.genes.TRACE * 100
mydata$expression_percentile <- mydata$Expression_rank_candidate/mydata$Num.genes.TRACE * 100
mydata$plof_percentile <- mydata$plof_rank/mydata$Num.genes.TRACE * 100
mydata$missense_percentile <- mydata$missense_rank/mydata$Num.genes.TRACE * 100



mydata <- mydata[order(mydata$trace_percentile),]





p <- ggplot(mydata, aes(x=reorder(Patient_ID, -trace_percentile), y=trace_percentile,
                        color='darkred')) +
  geom_linerange(aes(ymin = 100, ymax = trace_percentile), size = 2) +
  scale_y_reverse(limits = c(100, -1), expand = c(0, 0.0),
                  breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  #  geom_point(size = 0.5) +
  # scale_color_gradient(low="darkred", high="pink") + 
  theme_classic() +
  theme(
    axis.title.x = element_text(color="black", size=30, face="bold"),
    axis.title.y = element_text(color="black", size=30, face="bold"),
    # axis.text.x = element_blank(),
    axis.text.x = element_text(color="black", size=10, angle=90, hjust = 1),
    axis.text.y = element_text(size = 30, face="bold", color="black") 
  ) + 
  theme(legend.position="none") +
  labs(x="\nPatient-personalized models", y = "Tissue-associated disease gene\nrank top %\n") +
  # scale_x_discrete(expand = c(0.0, 0.0)) +
  geom_hline(yintercept = 50, linetype = 1, color = "black") +
  
  geom_hline(yintercept = 10, linetype = 3, color = "black", size = 1) +
  geom_hline(yintercept = 25, linetype = 2, color = "black", size = 1) +
  geom_hline(yintercept = 90, linetype = 3, color = "black", size = 1) +
  geom_hline(yintercept = 75, linetype = 2, color = "black", size = 1) +
  
  theme(plot.margin=unit(c(1,1,1,1),"cm"))

p



png(filename="../figure5/panelC/trace_patients_percentiles_bottom.png", width = 700, height = 700)
p
dev.off()




# percentage of genes ranked above median
length(which(mydata$trace_percentile <= 50)) / length(mydata$trace_percentile)
length(which(mydata$gado_percentile <= 50)) / length(mydata$gado_percentile)
length(which(mydata$expression_percentile <= 50)) / length(mydata$expression_percentile)
length(which(mydata$plof_percentile <= 50)) / length(mydata$plof_percentile)
length(which(mydata$missense_percentile <= 50)) / length(mydata$missense_percentile)


# percentage of genes ranked top quartile
length(which(mydata$trace_percentile <= 25)) / length(mydata$trace_percentile)
length(which(mydata$gado_percentile <= 25)) / length(mydata$gado_percentile)
length(which(mydata$expression_percentile <= 25)) / length(mydata$expression_percentile)
length(which(mydata$plof_percentile <= 25)) / length(mydata$plof_percentile)
length(which(mydata$missense_percentile <= 25)) / length(mydata$missense_percentile)


# percentage of genes ranked top 10%
length(which(mydata$trace_percentile <= 10)) / length(mydata$trace_percentile)
length(which(mydata$gado_percentile <= 10)) / length(mydata$gado_percentile)
length(which(mydata$expression_percentile <= 10)) / length(mydata$expression_percentile)
length(which(mydata$plof_percentile <= 10)) / length(mydata$plof_percentile)
length(which(mydata$missense_percentile <= 10)) / length(mydata$missense_percentile)



# brain ranking subset
total_brain <- mydata[which(mydata$Tissue == 'Whole_Brain'),]
length(which(total_brain$trace_percentile < 25))
length(total_brain$trace_percentile)
length(which(total_brain$trace_percentile < 25)) / length(total_brain$trace_percentile)


# skin ranking subset
total_skin <- mydata[which(mydata$Tissue == 'Skin'),]
length(which(total_skin$trace_percentile < 25))
length(total_skin$trace_percentile)
length(which(total_skin$trace_percentile < 25)) / length(total_skin$trace_percentile)
