
library("IHW")
library("DESeq2")
library("ggplot2")
library(reshape2)
library(apeglm)
library(tidyverse)
library(Hmisc)

TFs =read.csv2('TF_enrichment_1020.txt', sep = "\t")
gene_names = read.csv2("TF_protein_metadata.txt", sep = "\t")

colnames(gene_names)[3] = "TF"

total <- inner_join(TFs,gene_names , by="TF") # join is a great way to merge data


total$p.value = as.numeric(total$p.value)

total = total %>% filter(p.value < 0.00005)

#now will plot 1) enrichment of each group, 2) abundance of these genes, 3) size??

write.table(total,"sigTFdata.txt", sep = "\t")
