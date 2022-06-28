library(devtools)
library(plyr)
library(ggpubr)
library(dplyr)
library(ggsignif)
library(ggplot2)
library(reshape2)
library(scales)


all_omim_as_causal <- TRUE

method <- 'meta_MLP'


tissues <- c('whole_brain_causal', 'heart_causal', 'muscle_skeletal_causal', 'skin_causal',
             'liver_causal', 'nerve_tibial_causal', 'testis_causal', 'whole_blood_causal',
             'cortex_causal', 'cerebellum_causal', 'brain_other_causal')


genemap2 <- read.csv("../../data/genemap2.txt", sep='\t', skip = 3)
genemap2 <- genemap2[genemap2$Ensembl.Gene.ID != '',]
genemap2 <- genemap2[grepl("(3)", genemap2[["Phenotypes"]]), ]

omim_causal_genes <- unique(genemap2['Ensembl.Gene.ID'])

write.csv(omim_causal_genes, file = "../../data/omim_causal_genes.csv")



for (tissue in tissues) {
  
  
  file_name <- sub("%s", tissue, "../../data/tissues_pred_proba/Y_pred_prob_%s.csv")

  mydata <- read.csv(file_name, row.names = 1)
  
  
  mydata['pred'] <- mydata[method]
  
  
  mydata$Class <- as.character(mydata$Class)
  mydata$Class[mydata$Class=='1']<-'Causal in tissue'
  mydata$Class[mydata$Class==0]<-'Non-causal'
  
  if (all_omim_as_causal) {

    causal_other_tissues <- mydata[mydata$Class=='Non-causal',]
    
    causal_other_tissues <- subset(causal_other_tissues, (rownames(causal_other_tissues) %in% genemap2$Ensembl.Gene.ID))
    
  }
  
  
  mydata['Class'][row.names(causal_other_tissues),] <- 'Causal in other tissues'
  
  mydata$Class <- factor(mydata$Class, levels = c("Non-causal","Causal in other tissues", "Causal in tissue"))
  
  
  file_name <- "../../data/tissues_pred_proba/causality_probabilities_%s.csv"
  file_name <- sub("%s", tissue, file_name)
  
  write.csv(mydata, file = file_name)
  
  
  
}




