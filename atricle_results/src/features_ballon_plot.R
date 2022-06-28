library(ggplot2)
library(ggpubr)
library("RColorBrewer")


mydata <- read.csv('../data/figure3/df_for_feature_plot.csv',
                   header = TRUE, row.names = 1, sep = ",")


my_cols <- c("#F0F921FF", "#FCA636FF", 'firebrick', "#B12A90FF")
my_cols <-  colorRampPalette(brewer.pal(9,'Greens'))(100)
my_cols <-  colorRampPalette(brewer.pal(9,'Reds'))(100)
my_cols <-  colorRampPalette(brewer.pal(9,'YlOrRd'))(100)

p <- ggballoonplot(mydata, fill = "value", size.range = c(1, 30))+
  scale_fill_gradientn(colors = my_cols) + 
  # scale_fill_viridis_c(option="viridis") +
  theme_classic() +
  labs(y = "Feature type\n", x = "\nTissue")  + 
  theme(
    axis.title.x = element_text(color="black", size=30, face="bold"),
    axis.title.y = element_text(color="black", size=30, face="bold"),
    axis.text.x = element_text(color="black", size = 30, face="bold", angle = 45, hjust = 1),
    axis.text.y = element_text(color="black", size = 30, face="bold", angle = 0, hjust = 1)
    ) +
   theme(
    legend.title = element_blank(),
    legend.text = element_text(color="black", size=30, face="bold")
  ) +
      theme(plot.margin=unit(c(1,1,1,1),"cm")
            ) 

p


png(filename="../figure3/panelC/features_ballon_plot.png", 
    width = 1800, height = 1000)
p

dev.off()


##############################################################################

mydata <- read.csv('../data/figure3/df_for_tissue_origin_plot.csv',
                   header = TRUE, row.names = 1, sep = ",")

mydata <- rbind(mydata, 1- c(colSums(mydata)))
rownames(mydata)[rownames(mydata) == "9"] <- 'Other'


my_cols <-  colorRampPalette(brewer.pal(9,'YlOrRd'))(100)
# my_cols <-  colorRampPalette(rev(brewer.pal(9,'RdYlGn')))(100)
# my_cols <-  colorRampPalette(brewer.pal(9,'YlGn'))(100)

p <- ggballoonplot(mydata, fill = "value", size.range = c(1, 30))+
  scale_fill_gradientn(colors = my_cols) + 
  # scale_fill_viridis_c(option="viridis") +
  theme_classic() +
  labs(y = "Tissue origin of features\n", x = "\nTissue")  + 
  theme(
    axis.title.x = element_text(color="black", size=30, face="bold"),
    axis.title.y = element_text(color="black", size=30, face="bold"),
    axis.text.x = element_text(color="black", size = 30, face="bold", angle = 45, hjust = 1),
    axis.text.y = element_text(color="black", size = 30, face="bold", angle = 0, hjust = 1)
  ) +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(color="black", size=30, face="bold")
  ) +
  theme(plot.margin=unit(c(1,1,1,1),"cm")
  ) 

p


png(filename="../figure3/panelB/tissue_origin_ballon_plot.png", 
    width = 1400, height = 1000)
p

dev.off()

