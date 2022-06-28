
library(ggpubr)
library(ggplot2)
library(reshape2)
library(devtools)
library(plyr)
library(dplyr)
library(ggsignif)
library(scales)



mydata <- read.csv('../data/figure5/table_s4_june22.csv', 
                   header = TRUE, sep = ",")

total <- mydata[order(mydata$TRACE_model_rank),]



# trace_median
median(total$TRACE_model_rank)
median(total$GADO_rank)
median(total$Expression_rank_candidate)
median(total$plof_rank)
median(total$missense_rank)



mydata_new <- total[, c("Patient_ID", "plof_rank", "missense_rank", "Expression_rank_candidate", "GADO_rank", "TRACE_model_rank")]

ks <- ks.test(total$TRACE_model_rank, total$GADO_rank, alternative = "greater")
ks

mydata_new <- melt(mydata_new, id=c("Patient_ID"))



p <- ggboxplot(mydata_new, x = "variable", y = "value",
               outlier.shape = NA,
               color = "variable",
               add = "jitter", 
               xlab="\n", ylab = "Tissue-associated disease gene rank") +

  theme(legend.position="NA") +
  theme(
    axis.title.x = element_text(color="black", size=30, face="bold"),
    axis.title.y = element_text(color="black", size=30, face="bold"),
    axis.text.x = element_text(size = 30, face="bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 30, face="bold")
  )  + 
  scale_x_discrete(labels=c("GADO_rank" = "GADO", 
                            "TRACE_model_rank" = "TRACE",
                            "Expression_rank_candidate" = "Expression",
                            "plof_rank" = "pLoF",
                            "missense_rank" = "Missense"
                            )
                   ) +
  theme(plot.margin=unit(c(1,1,1,1),"cm")) +

  scale_y_reverse(limits = c(1200, -500), expand = c(0, 0.0),
                  breaks = c(1, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200)
                   ) + 
  geom_violin(color = "darkgray", fill='NA', trim = TRUE) + 
  stat_compare_means(comparisons = list( c("TRACE_model_rank", "plof_rank"),
                                         c("TRACE_model_rank", "missense_rank"),
                                         c("TRACE_model_rank", "Expression_rank_candidate"), 
                                         c("TRACE_model_rank", "GADO_rank")
                                          ),
                     paired = TRUE, label = "p.format", method = "wilcox.test", 
                     label.y = c(50, 150, 250, 350),
                     method.args = list(alternative = "greater", correct = FALSE, exact=FALSE)
                     
                     )

p


file_name <- "../figure5/panelD/methods_compare_boxplot_cor.png"

png(filename=file_name, width = 700, height = 1000)
p
dev.off()


# Adjust p-values (BH)
p.adjust(c(0.00031, 0.01, 0.0097, 0.0082), "BH")


