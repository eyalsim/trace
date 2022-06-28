library(ggplot2)
library(dplyr)
library(reshape2)


complete_dataset <- read.csv("../dataset/df_complete_dataset_ready_adapted_no_missing_values.csv", row.names = 1)


# Select column whose name ends with "causal"
complete_dataset_labels <- complete_dataset[, grep("causal", colnames(complete_dataset))]

brain_dataset_labels <- complete_dataset_labels[c('cortex_causal',
                                                  'cerebellum_causal',
                                                  'basal_ganglia_causal', 
                                                  'spinal_cord_causal',
                                                  'hypothalamus_causal',
                                                  'amygdala_causal',
                                                  'brain_other_causal',
                                                  'whole_brain_causal',
                                                  'whole_brain_lower_causal')]

complete_dataset_labels <- complete_dataset_labels[, !(names(complete_dataset_labels) %in% c('Heart_causal',
                                                                                             'brain_causal',
                                                                                             'brain.0_causal',
                                                                                             'brain.1_causal',
                                                                                             'brain.2_causal',
                                                                                             'brain_not_specific_causal',
                                                                                             'brain.not_specific_causal',
                                                                                             'cortex_causal',
                                                                                             'cerebellum_causal',
                                                                                             'heart_atrial_appendage_causal',
                                                                                             'heart_left_ventricle_causal',
                                                                                             'basal_ganglia_causal', 
                                                                                             'spinal_cord_causal',
                                                                                             'hypothalamus_causal',
                                                                                             'amygdala_causal',
                                                                                             'brain_other_causal',
                                                                                             'whole_brain_lower_causal'))]


chosen_tissues <- c('whole_brain_causal', 'heart_causal', 'muscle_skeletal_causal', 'skin_causal', 
                    'liver_causal', 'testis_causal', 'whole_blood_causal', 'nerve_tibial_causal'
)


complete_dataset_chosen_labels <- complete_dataset_labels[chosen_tissues]
complete_dataset_labels <- complete_dataset_labels[rowSums(complete_dataset_labels) > 0, ]


brain_dataset_labels <- brain_dataset_labels[rowSums(brain_dataset_labels) > 0, ]



complete_dataset_chosen_labels_no_other <- complete_dataset_chosen_labels
complete_dataset_chosen_labels_no_other <- complete_dataset_chosen_labels_no_other[rowSums(complete_dataset_chosen_labels_no_other) > 0, ]



complete_dataset_other_labels <- complete_dataset_labels[, !(names(complete_dataset_labels) %in% chosen_tissues)]

complete_dataset_other_labels$'other' <- rowSums(complete_dataset_other_labels)
complete_dataset_other_labels$'other'[complete_dataset_other_labels$'other' != 0] <- 1
complete_dataset_other_labels <- complete_dataset_other_labels[rowSums(complete_dataset_other_labels) > 0, ]


complete_dataset_other_labels <- complete_dataset_other_labels['other']
complete_dataset_other_labels <- melt(complete_dataset_other_labels)


complete_dataset_chosen_labels <- melt(complete_dataset_chosen_labels)
complete_dataset_chosen_labels <- subset(complete_dataset_chosen_labels, value == 1) 

complete_dataset_chosen_labels <- rbind(complete_dataset_chosen_labels, complete_dataset_other_labels)

complete_dataset_chosen_labels$variable <- factor(complete_dataset_chosen_labels$variable, 
                                                  levels=c(
                                                    "whole_brain_causal",
                                                    "other",
                                                    "heart_causal",
                                                    "skin_causal",
                                                    "muscle_skeletal_causal",
                                                    "testis_causal",
                                                    "whole_blood_causal",
                                                    "nerve_tibial_causal",
                                                    "liver_causal"))


p <- ggplot(complete_dataset_chosen_labels, aes(x=variable)) +
  geom_bar(stat="count", fill = 'coral3')  + 
  theme_classic()  +
  coord_flip() +
  theme(axis.text = element_text(angle = 0, size = 40, face="bold", colour = 'black')) +
  theme(axis.title = element_text(size=40, face="bold")) +
  labs(x = "Tissues\n", y = "\n# Tissue-associated disease genes") +   
  # geom_hline(yintercept=100) +
  geom_text(stat='count', aes(label=..count..), size=14, hjust=0) +
  scale_x_discrete(labels=c(
                            "whole_brain_causal" = 'Brain',
                            "heart_causal" = 'Heart',
                            "skin_causal" = "Skin",
                            "muscle_skeletal_causal" = "Muscle",
                            "testis_causal" = "Testis",
                            "whole_blood_causal" = "Blood",
                            "nerve_tibial_causal" = "Nerve",
                            "liver_causal" = "Liver",
                            "other" = "Other")) +
  
  scale_y_continuous(limits = c(0, 580), expand = c(0, 0.0)) +
  
  theme(plot.margin=unit(c(1,1,1,1),"cm"))
p  


png(filename="../figure1/panelB/dataset_barplot.png", width = 1000, height = 700)
p
dev.off()



only_cortex <- subset(brain_dataset_labels, cortex_causal == 1 &
                        cerebellum_causal == 0) 

only_cerebellum <- subset(brain_dataset_labels, cortex_causal == 0 &
                            cerebellum_causal == 1) 

cortex_and_cerebellum <- subset(brain_dataset_labels, cortex_causal == 1 &
                                  cerebellum_causal == 1) 

other_brain_regions <- subset(brain_dataset_labels, cortex_causal == 0 &
                                cerebellum_causal == 0 & 
                                whole_brain_lower_causal == 1) 


lower_confidence_brain <- subset(brain_dataset_labels, whole_brain_lower_causal == 1 &
                                   whole_brain_causal == 0) 

high_confidence_brain <- subset(brain_dataset_labels, whole_brain_causal == 1) 





brain_df <- data.frame(
  group = c("Cortex", "Cerebellum", "Cortex and Cerebellum", 'Other brain regions'),
  value = c(length(rownames(only_cortex)),
            length(rownames(only_cerebellum)),
            length(rownames(cortex_and_cerebellum)),
            length(rownames(other_brain_regions))
  )
)

# Add label position
brain_df <- brain_df %>%
  arrange(desc(group)) %>%
  mutate(lab.ypos = cumsum(value) - 0.5*value)

mycols <- c("#0073C2FF", "#EFC000FF", "aquamarine3", "#868686FF")


bp <- ggplot(brain_df, aes(x=2, y=value, fill=group))+
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0) +
  geom_text(aes(y = lab.ypos, label = value), color = "white", size=14)+
  scale_fill_manual(values = mycols) +
  theme_void() +
  xlim(0.5, 2.5) +
  theme(legend.text=element_text(size=30),
        legend.title = element_blank(),
        legend.position = "right")
bp

png(filename="../figureS9/panelA/brain_donut_chart.png", width = 700, height = 350)
bp
dev.off()



brain_df <- data.frame(
  group = c("High confidence", "Medium confidence"),
  value = c(length(rownames(high_confidence_brain)),
            length(rownames(lower_confidence_brain))
  )
)

# Add label position
brain_df <- brain_df %>%
  arrange(desc(group)) %>%
  mutate(lab.ypos = cumsum(value) - 0.5*value)

mycols <- c("firebrick4", "firebrick3")


bp <- ggplot(brain_df, aes(x=2, y=value, fill=group))+
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0) +
  geom_text(aes(y = lab.ypos, label = value), color = "white", size=14)+
  scale_fill_manual(values = mycols) +
  theme_void() +
  xlim(0.5, 2.5) +
  theme(legend.text=element_text(size=30),
        legend.title = element_blank(),
        legend.position = "left")
bp

png(filename="../figureS9/panelA/brain_donut_confidence_chart.png", width = 700, height = 350)
bp
dev.off()







