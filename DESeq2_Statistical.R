#DESeq2 differential gene expression analysis
#Assuming you ran Kallisto_to_DESeq2_prepare.R
# Currently not a function
library("IHW")
library("DESeq2")
library("ggplot2")
library(reshape2)
library(apeglm)
library(tidyverse)
library(Hmisc)
setwd("/Users/gregbrittingham/Dropbox (Personal)/Scripts/Nacho_RNAseq_2020_final/Ignacio_RNAseq_drop_snf5_del")

dds = readRDS('load_files/dds.rda') #make sure to load correct data file here!! I originally had 31 files, but we cut down to 27 to get the highest quality data. 
sampledata = readRDS('load_files/sampledata.rda')
t2g = readRDS('load_files/t2g.rda')

alpha = 0.05
lfcThreshold = 1



# Create Wald Test contrasts ----
# log2fc below not neccesary as these here provide log2fc!!

WT_Starve_v_Dex_pH_7 = results(dds, contrast=c("Condition", "WT_Stv_pH7","WT_Dex"), alpha = alpha, lfcThreshold = lfcThreshold) #+ = up in stv vs dex for WT

    names = match(rownames(WT_Starve_v_Dex_pH_7),t2g[,'ens_gene'])
    ordered_anno = as.data.frame(t2g[names,])
    
    WT_Starve_v_Dex_pH_7$ext_gene = ordered_anno$ext_gene
    WT_Starve_v_Dex_pH_7$description = ordered_anno$description # add gene names and descriptions
    
    WT_Starve_v_Dex_pH_7_sig <- subset(WT_Starve_v_Dex_pH_7, padj < alpha, log2FoldChange > lfcThreshold)
    
#
WT_Starve_v_Dex = results(dds, contrast=c("Condition","WT_Stv_pH5","WT_Dex") , alpha = alpha, lfcThreshold = lfcThreshold)

    names = match(rownames(WT_Starve_v_Dex),t2g[,'ens_gene'])
    ordered_anno = as.data.frame(t2g[names,])
    
    WT_Starve_v_Dex$ext_gene = ordered_anno$ext_gene
    WT_Starve_v_Dex$description = ordered_anno$description # add gene names and descriptions
    
    WT_Starve_v_Dex_sig <- subset(WT_Starve_v_Dex, padj < alpha, log2FoldChange > lfcThreshold)

#
dQ_Starve_v_Dex = results(dds, contrast=c("Condition","dQ_snf5_Stv_pH5", "dQ_snf5_Dex"), alpha = alpha, lfcThreshold = lfcThreshold)
    names = match(rownames(dQ_Starve_v_Dex),t2g[,'ens_gene'])
    ordered_anno = as.data.frame(t2g[names,])
    
    dQ_Starve_v_Dex$ext_gene = ordered_anno$ext_gene
    dQ_Starve_v_Dex$description = ordered_anno$description # add gene names and descriptions
    
    dQ_Starve_v_Dex_sig <- subset(dQ_Starve_v_Dex, padj < alpha, log2FoldChange > lfcThreshold)

# Null_Starve_v_Dex = results(dds, contrast=c("Condition","Null_Stv_pH5","Null_Dex"), alpha = alpha, lfcThreshold = lfcThreshold)
# 
#   names = match(rownames(Null_Starve_v_Dex),t2g[,'ens_gene'])
#   ordered_anno = as.data.frame(t2g[names,])
#   
#   Null_Starve_v_Dex$ext_gene = ordered_anno$ext_gene
#   Null_Starve_v_Dex$description = ordered_anno$description # add gene names and descriptions
#   
#   Null_Starve_v_Dex_sig <- subset(Null_Starve_v_Dex, padj < alpha, log2FoldChange > lfcThreshold)
# 

HtoA_Starve_v_Dex = results(dds, contrast=c("Condition","HtoA_snf5_Stv_pH5", "HtoA_snf5_Dex"), alpha = alpha, lfcThreshold = lfcThreshold)

  names = match(rownames(HtoA_Starve_v_Dex),t2g[,'ens_gene'])
  ordered_anno = as.data.frame(t2g[names,])
  
  HtoA_Starve_v_Dex$ext_gene = ordered_anno$ext_gene
  HtoA_Starve_v_Dex$description = ordered_anno$description # add gene names and descriptions
  
  HtoA_Starve_v_Dex_sig <- subset(HtoA_Starve_v_Dex, padj < alpha, log2FoldChange > lfcThreshold)
  
dQ_v_WT_Dex = results(dds, contrast = c("Condition","dQ_snf5_Dex","WT_Dex"), alpha = alpha, lfcThreshold = lfcThreshold)
  names = match(rownames(dQ_v_WT_Dex),t2g[,'ens_gene'])
  ordered_anno = as.data.frame(t2g[names,])
  
  dQ_v_WT_Dex$ext_gene = ordered_anno$ext_gene
  dQ_v_WT_Dex$description = ordered_anno$description # add gene names and descriptions
  
  dQ_v_WT_Dex_sig = subset(dQ_v_WT_Dex, padj < alpha, log2FoldChange > lfcThreshold)

dQ_v_WT_Starve = results(dds, contrast = c("Condition","dQ_snf5_Stv_pH5","WT_Stv_pH5"), alpha = alpha, lfcThreshold = lfcThreshold)
  names = match(rownames(dQ_v_WT_Starve),t2g[,'ens_gene'])
  ordered_anno = as.data.frame(t2g[names,])
  
  dQ_v_WT_Starve$ext_gene = ordered_anno$ext_gene
  dQ_v_WT_Starve$description = ordered_anno$description # add gene names and descriptions
  
  dQ_v_WT_Starve_sig = subset(dQ_v_WT_Starve, padj < alpha, log2FoldChange > lfcThreshold)
  

#Results path
results_path = '/Users/gbrittingham/Documents/R/Scripts/Ignacio_RNAseq_drop_snf5_del/'
# Writing the list to CSV ----
write.table(as.data.frame(WT_Starve_v_Dex_pH_7_sig), file = paste(results_path,"WT_Starve_v_Dex_pH_7.txt", sep = ""), sep ="\t")
write.table(as.data.frame(WT_Starve_v_Dex_sig), file= paste(results_path, "WT_Starve_v_Dex.txt", sep = ""), sep ="\t")
write.table(as.data.frame(dQ_Starve_v_Dex_sig), file= paste(results_path, "dQ_Starve_v_Dex.txt", sep = ""), sep ="\t")
# write.table(as.data.frame(Null_Starve_v_Dex_sig), file= paste(results_path, "Null_Starve_v_Dex.txt", sep = ""), sep ="\t")
write.table(as.data.frame(HtoA_Starve_v_Dex_sig), file= paste(results_path, "HtoA_Starve_v_Dex.txt", sep = ""), sep ="\t")
write.table(as.data.frame(dQ_v_WT_Dex_sig), file= paste(results_path, "dQ_v_WT_Dex.txt", sep = ""), sep ="\t")
write.table(as.data.frame(dQ_v_WT_Starve_sig), file= paste(results_path, "dQ_v_WT_Starve.txt", sep = ""), sep ="\t")

# With filtering and with lfcShrinking (needed for volcano plot) ----

# log2fc depreciated ------------------------------------------------------


WT_Starve_v_Dex_pH_7 = lfcShrink(dds, contrast=c("Condition", "WT_Stv_pH7","WT_Dex"), coef = 2, type = "apeglm") #+ = up in stv vs dex for WT

names = match(rownames(WT_Starve_v_Dex_pH_7),t2g[,'ens_gene'])
ordered_anno = as.data.frame(t2g[names,])
  
  WT_Starve_v_Dex_pH_7$ext_gene = ordered_anno$ext_gene
  WT_Starve_v_Dex_pH_7$description = ordered_anno$description # add gene names and descriptions
  
  WT_Starve_v_Dex_pH_7_Volcano <- subset(WT_Starve_v_Dex_pH_7, padj < alpha, log2FoldChange > lfcThreshold)

#
WT_Starve_v_Dex = lfcShrink(dds, contrast=c("Condition","WT_Stv_pH5","WT_Dex") , coef = 2, type = "apeglm")

  names = match(rownames(WT_Starve_v_Dex),t2g[,'ens_gene'])
  ordered_anno = as.data.frame(t2g[names,])
  
  WT_Starve_v_Dex$ext_gene = ordered_anno$ext_gene
  WT_Starve_v_Dex$description = ordered_anno$description # add gene names and descriptions
  
  WT_Starve_v_Dex_Volcano <- subset(WT_Starve_v_Dex, padj < alpha, log2FoldChange > lfcThreshold)

#
dQ_Starve_v_Dex = lfcShrink(dds, contrast=c("Condition","dQ_snf5_Stv_pH5", "dQ_snf5_Dex"), coef = 2, type = "apeglm")
  names = match(rownames(dQ_Starve_v_Dex),t2g[,'ens_gene'])
  ordered_anno = as.data.frame(t2g[names,])
  
  dQ_Starve_v_Dex$ext_gene = ordered_anno$ext_gene
  dQ_Starve_v_Dex$description = ordered_anno$description # add gene names and descriptions
  
  dQ_Starve_v_Dex_Volcano <- subset(dQ_Starve_v_Dex, padj < alpha, log2FoldChange > lfcThreshold)

Null_Starve_v_Dex = lfcShrink(dds, contrast=c("Condition","Null_Stv_pH5","Null_Dex"), coef = 2, type = "apeglm")

  names = match(rownames(Null_Starve_v_Dex),t2g[,'ens_gene'])
  ordered_anno = as.data.frame(t2g[names,])
  
  Null_Starve_v_Dex$ext_gene = ordered_anno$ext_gene
  Null_Starve_v_Dex$description = ordered_anno$description # add gene names and descriptions
  
  Null_Starve_v_Dex_Volcano <- subset(Null_Starve_v_Dex, padj < alpha, log2FoldChange > lfcThreshold)


HtoA_Starve_v_Dex = lfcShrink(dds, contrast=c("Condition","HtoA_snf5_Stv_pH5", "HtoA_snf5_Dex"), coef = 2, type = "apeglm")

  names = match(rownames(HtoA_Starve_v_Dex),t2g[,'ens_gene'])
  ordered_anno = as.data.frame(t2g[names,])
  
  HtoA_Starve_v_Dex$ext_gene = ordered_anno$ext_gene
  HtoA_Starve_v_Dex$description = ordered_anno$description # add gene names and descriptions
  
  HtoA_Starve_v_Dex_Volcano <- subset(HtoA_Starve_v_Dex, padj < alpha, log2FoldChange > lfcThreshold)

dQ_v_WT_Dex = lfcShrink(dds, contrast=c("Condition","dQ_snf5_Dex","WT_Dex"), coef = 2, type = "apeglm")
  
  names = match(rownames(dQ_v_WT_Starve),t2g[,'ens_gene'])
  ordered_anno = as.data.frame(t2g[names,])
  
  dQ_v_WT_Dex$ext_gene = ordered_anno$ext_gene
  dQ_v_WT_Dex$description = ordered_anno$description # add gene names and descriptions
  
  dQ_v_WT_Dex_volcano <- subset(dQ_v_WT_Dex, padj < alpha, log2FoldChange > lfcThreshold)  
 
  
  #comparing within starvation - log2fc depciated ####----- 
dQ_v_WT_Starve = lfcShrink(dds, contrast=c("Condition","WT_Stv_pH5", "dQ_snf5_Stv_pH5"), coef = 2, type = "apeglm")
  
  names = match(rownames(dQ_v_WT_Starve),t2g[,'ens_gene'])
  ordered_anno = as.data.frame(t2g[names,])
  
  dQ_v_WT_Starve$ext_gene = ordered_anno$ext_gene
  dQ_v_WT_Starve$description = ordered_anno$description # add gene names and descriptions
  
  dQ_v_WT_Starve_volcano <- subset(dQ_v_WT_Starve, padj < alpha, log2FoldChange > lfcThreshold)
  
dQ_v_WT_Starve_pH7 = lfcShrink(dds, contrast=c("Condition","WT_Stv_pH7", "dQ_snf5_Stv_pH5"), coef = 2, type = "apeglm")
  
  names = match(rownames(dQ_v_WT_Starve_pH7),t2g[,'ens_gene'])
  ordered_anno = as.data.frame(t2g[names,])
  
  dQ_v_WT_Starve_pH7$ext_gene = ordered_anno$ext_gene
  dQ_v_WT_Starve_pH7$description = ordered_anno$description # add gene names and descriptions
  
  dQ_v_WT_Starve_pH7_volcano <- subset(dQ_v_WT_Starve_pH7, padj < alpha, log2FoldChange > lfcThreshold)
  
dQ_v_HtoA_Starve = lfcShrink(dds, contrast=c("Condition","HtoA_snf5_Stv_pH5", "dQ_snf5_Stv_pH5"), coef = 2, type = "apeglm"E)
  
  names = match(rownames(dQ_v_HtoA_Starve),t2g[,'ens_gene'])
  ordered_anno = as.data.frame(t2g[names,])
  
  dQ_v_HtoA_Starve$ext_gene = ordered_anno$ext_gene
  dQ_v_HtoA_Starve$description = ordered_anno$description # add gene names and descriptions
  
  dQ_v_WT_Starve_pH7_sig <- subset(dQ_v_HtoA_Starve, padj < alpha, log2FoldChange > lfcThreshold)
  
# dQ_v_Null_Starve = lfcShrink(dds, contrast=c("Condition","Null_Stv_pH5", "dQ_snf5_Stv_pH5"), coef = 2, type = "apeglm")
#   
#   names = match(rownames(dQ_v_Null_Starve),t2g[,'ens_gene'])
#   ordered_anno = as.data.frame(t2g[names,])
#   
#   dQ_v_Null_Starve$ext_gene = ordered_anno$ext_gene
#   dQ_v_Null_Starve$description = ordered_anno$description # add gene names and descriptions
#   
#   dQ_v_WT_Starve_pH7_sig <- subset(dQ_v_Null_Starve, padj < alpha, log2FoldChange > lfcThreshold)
#   
# WT_v_Null_Starve = lfcShrink(dds, contrast=c("Condition","WT_Stv_pH5", "Null_Stv_pH5"), coef = 2, type = "apeglm")
#   
#   names = match(rownames(WT_v_Null_Starve),t2g[,'ens_gene'])
#   ordered_anno = as.data.frame(t2g[names,])
#   
#   WT_v_Null_Starve$ext_gene = ordered_anno$ext_gene
#   WT_v_Null_Starve$description = ordered_anno$description # add gene names and descriptions
#   
#   dQ_v_WT_Starve_pH7_sig <- subset(WT_v_Null_Starve, padj < alpha, log2FoldChange > lfcThreshold)
  
WT_pH4_v_7_Starve = lfcShrink(dds, contrast=c("Condition","WT_Stv_pH5","WT_Stv_pH7") , coef = 2, type = "apeglm")
  
  names = match(rownames(WT_pH4_v_7_Starve),t2g[,'ens_gene'])
  ordered_anno = as.data.frame(t2g[names,])
  
  WT_pH4_v_7_Starve$ext_gene = ordered_anno$ext_gene
  WT_pH4_v_7_Starve$description = ordered_anno$description # add gene names and descriptions
  
  WT_pH4_v_7_Starve_Volcano <- subset(WT_pH4_v_7_Starve, padj < alpha, log2FoldChange > lfcThreshold) 

# ---- Create DE genes log2FC list for use with Heatmap ----
  
  
  #need to add the following from individual!!
  #---- turning DE genes into Log2FC values for plotting - use with other scripts---- 
  # WT_Starve_v_Dex_df = as.data.frame(WT_Starve_v_Dex)
  # dQ_Starve_v_Dex_df = as.data.frame(dQ_Starve_v_Dex)
  # HtoA_Starve_v_Dex_df = as.data.frame(HtoA_Starve_v_Dex)
  # Null_Starve_v_Dex_df = as.data.frame(Null_Starve_v_Dex)
  # WT_Starve_v_Dex_pH_7_df = as.data.frame(WT_Starve_v_Dex_pH_7)
  # 
  # # Should set a criteria for abs expression level here maybe? - Will revisit
  # a = as.data.frame(WT_Starve_v_Dex_df$log2FoldChange)
  # rownames(a) = paste(rownames(WT_Starve_v_Dex_df), '_', WT_Starve_v_Dex_df$ext_gene, sep = "")
  # colnames(a) = 'WT'
  # 
  # b =as.data.frame(dQ_Starve_v_Dex_df$log2FoldChange)
  # rownames(b) = paste(rownames(dQ_Starve_v_Dex_df), '_', dQ_Starve_v_Dex_df$ext_gene, sep = "")
  # colnames(b) = 'dQ'
  # 
  # c =as.data.frame(HtoA_Starve_v_Dex_df$log2FoldChange)
  # rownames(c) = paste(rownames(HtoA_Starve_v_Dex_df), '_', HtoA_Starve_v_Dex_df$ext_gene, sep = "")
  # colnames(c) = 'HtoA'
  # 
  # d =as.data.frame( Null_Starve_v_Dex_df$log2FoldChange)
  # rownames(d) = paste(rownames( Null_Starve_v_Dex_df), '_',  Null_Starve_v_Dex_df$ext_gene, sep = "")
  # colnames(d) = 'Null'
  # 
  # e =as.data.frame(  WT_Starve_v_Dex_pH_7_df$log2FoldChange)
  # rownames(e) = paste(rownames(  WT_Starve_v_Dex_pH_7_df), '_',   WT_Starve_v_Dex_pH_7_df$ext_gene, sep = "")
  # colnames(e) = 'WT_pH7'
  # 
  # DE_Genes_in_Starvation = cbind(a,b,c,d,e)
  # saveRDS(DE_Genes_in_Starvation, "DE_Genes_in_Starvation_combined_master.rda")
  # 
  # 
  # Without the '_' and common gene name ----
  
  
  # Calculating number of genes induced or repressed
  a
  WT_Starve_v_Dex_df = as.data.frame(WT_Starve_v_Dex)
  dQ_Starve_v_Dex_df = as.data.frame(dQ_Starve_v_Dex)
  HtoA_Starve_v_Dex_df = as.data.frame(HtoA_Starve_v_Dex)
 # Null_Starve_v_Dex_df = as.data.frame(Null_Starve_v_Dex)
  WT_Starve_v_Dex_pH_7_df = as.data.frame(WT_Starve_v_Dex_pH_7)
  
  a = as.data.frame(WT_Starve_v_Dex_df$log2FoldChange)
  rownames(a) = rownames(WT_Starve_v_Dex_df)
  colnames(a) = 'WT'
  
  b =as.data.frame(dQ_Starve_v_Dex_df$log2FoldChange)
  rownames(b) = rownames(dQ_Starve_v_Dex_df)
  colnames(b) = 'dQ'
  
  c =as.data.frame(HtoA_Starve_v_Dex_df$log2FoldChange)
  rownames(c) = rownames(HtoA_Starve_v_Dex_df)
  colnames(c) = 'HtoA'
  
  # d =as.data.frame( Null_Starve_v_Dex_df$log2FoldChange)
  # rownames(d) = rownames( Null_Starve_v_Dex_df)
  # colnames(d) = 'Null'
  
  e =as.data.frame(  WT_Starve_v_Dex_pH_7_df$log2FoldChange)
  rownames(e) = rownames(  WT_Starve_v_Dex_pH_7_df)
  colnames(e) = 'WT_pH7'
  
  #DE_Genes_in_Starvation = cbind(a,b,c,d,e)
  DE_Genes_in_Starvation = cbind(a,b,c,e)
  saveRDS(DE_Genes_in_Starvation, "DE_Genes_in_Starvation_combined_master.rda")
  
  # gn.most.sign <- rownames(res.DESeq2)[1]
  # gn.most.diff.val <- counts(dds.norm, normalized=T)[gn.most.sign,]
  # barplot(gn.most.diff.val, col=col.pheno.selected, main=gn.most.sign, las=2, cex.names=0.5)
  
  #just adding in gene naems ect.
  res_table = counts(dds, normalized=T) #gene counts results table
  colnames(res_table) =  dds$sample
  res_ordered = res_table[,c(1,10,23,2,11,24,3,12,25,4,13,26,5,14,19,27,6,15,20,28,7,16,21,29,8,17,22,30,9,18,31)]
  #for when I dropped a bunch of samples
  
  res_ordered = res_table[,c(1,7,19,2,8,20,3,9,21,4,10,22,5,11,23,6,15,27,12,16,24,13,17,25,14,18,26)] #for when I dropped a bunch of samples

  
  DEstarveGenes <- WT_Starve_v_Dex_sig[order(WT_Starve_v_Dex_sig$pvalue),] # can use these with heatmap

  dQ_genes_expr = dQ_Starve_v_Dex_sig[order(dQ_Starve_v_Dex_sig$pvalue),]
  
  
  # 
  # sum(a[,1])
  # sum(a[,2])
  # sum(a[,3])
  # sum(a[,4])
  # sum(a[,5])
  # 
  ### -------
  #making fold repression

  fold_repression_prep =  a[a < -2,]
  fold_repression_prep = as.data.frame(fold_repression_prep)
  rownames(fold_repression_prep) = subset(rownames(a), a < -2)
  
  fold_repression = DE_Genes_in_Starvation[rownames(fold_repression_prep),]
  fold_repression = data.matrix(fold_repression, rownames.force = NA)
  
  write.table(fold_repression, "Fold_repression_in_starvation.txt")
  
  fold_induction_prep =  a[a > 2,]
  fold_induction_prep = as.data.frame(fold_induction_prep)
  rownames(fold_induction_prep) = subset(rownames(a), a > 2)
  
  fold_induction = DE_Genes_in_Starvation[rownames(fold_induction_prep),]
  fold_induction = data.matrix(fold_induction, rownames.force = NA)
  
  write.table(fold_induction, "fold_induction_n_in_starvation_new.txt")
  
  #-----
  #read in induction and repression
  
  fold_repression = read.table("Fold_repression_in_starvation.txt")
  fold_induction = read.table("fold_induction_n_in_starvation_new.txt")

  
  melt_fold_rep = melt(as.matrix(fold_repression))
  colnames(melt_fold_rep) =  c("Gene","Condition", "FC")
  
  melt_fold_ind = melt(as.matrix(fold_induction))
  colnames(melt_fold_ind) =  c("Gene","Condition", "FC")
  
  melt_fold_rep = melt_fold_rep %>% add_column(State = "Repressed")
  melt_fold_ind = melt_fold_ind %>% add_column(State = "Induced")
  
 tuna_melt =  bind_rows(melt_fold_ind,melt_fold_rep)
  
 dodge <- position_dodge(width = 0.4)
#- Jitter plot
 palette = c("red","blue")
 
 ggplot(tuna_melt, aes( x = factor(Condition), y = FC, color = State)) +  geom_boxplot(notch = F, dodge.with = 0 ) + geom_jitter(cex = 2, position = position_jitterdodge(jitter.width = 0.3)) + xlab("Condition") +
   ylab("Average Log2 Fold Change") + ylim(-10,10) + theme_classic(base_size = 18)
 
 
 
 
 
 
 #----- OLD
 
 
 
 
 
 
 
 
 
 
 

  ggplot(melt_fold_rep, aes(factor(Condition), Average_Fold_Repression)) + geom_boxplot(aes(fill = Condition)) + geom_boxplot(width=.1, outlier.colour=NA, position = dodge) + xlab("Condition") + 
    ylab("Average Log2 Fold Repression") + ylim(-10,5) + scale_fill_manual(values=rep_palette) + theme_bw()
  ## ---- fold induction NEW

 
  # glu WT, dQ, HtoA c("black","#7fc97f", "#2065d4"),
  #stv WT, dQ, HtoA, pH7 c("#b2730d","red","#f4ce0c","pink", "purple")
  
  dodge <- position_dodge(width = 0.4)
  #2065d4
  
  ggplot(melt_fold_ind, aes(factor(Condition), Average_Fold_Induction)) + geom_boxplot(aes(fill = Condition)) + geom_boxplot(width=.1, outlier.colour=NA, position = dodge) + xlab("Condition") +
    ylab("Average Log2 Fold Induction")  + ylim(-5,10) +  scale_fill_manual(values=ind_palette)
  
  ### -----------------what happens if we now take dQ as standard

# dQ ind rep --------------------------------------------------------------

  
  fold_repression_prep =  b[b < -2,]
  fold_repression_prep = as.data.frame(fold_repression_prep)
  rownames(fold_repression_prep) = subset(rownames(b), b < -2)
  
  fold_repression = DE_Genes_in_Starvation[rownames(fold_repression_prep),]
  fold_repression = data.matrix(fold_repression, rownames.force = NA)
  
  write.table(fold_repression, "Fold_repression_in_starvation.txt")
  
  melt_fold_rep = melt(as.matrix(fold_repression))
  colnames(melt_fold_rep) =  c("Gene","Condition", "Average_Fold_Repression")
  
  
  
  dodge <- position_dodge(width = 0.4)
  
  
  ggplot(melt_fold_rep, aes(factor(Condition), Average_Fold_Repression)) + geom_boxplot(aes(fill = Condition)) + geom_boxplot(width=.1, outlier.colour=NA, position = dodge) + xlab("Condition") +
    ylab("Average Log2 Fold Repression") + ylim(-10,5)
  ## ---- fold induction NEW
  fold_induction_prep =  b[b > 2,]
  fold_induction_prep = as.data.frame(fold_induction_prep)
  rownames(fold_induction_prep) = subset(rownames(b), b > 2)
  
  fold_induction = DE_Genes_in_Starvation[rownames(fold_induction_prep),]
  fold_induction = data.matrix(fold_induction, rownames.force = NA)
  
  write.table(fold_induction, "fold_induction_n_in_starvation_new.txt")
  
  melt_fold_ind = melt(as.matrix(fold_induction))
  colnames(melt_fold_ind) =  c("Gene","Condition", "Average_Fold_Induction")
  
  
  
  dodge <- position_dodge(width = 0.4)
  
  
  ggplot(melt_fold_ind, aes(factor(Condition), Average_Fold_Induction)) + geom_boxplot(aes(fill = Condition)) + geom_boxplot(width=.1, outlier.colour=NA, position = dodge) + xlab("Condition") +
    ylab("Average Log2 Fold Induction")  + ylim(-5,15)
  
  
  
  ### ---- a
  ## fold induction OLD!!!
  
  fold_induction = read.table("Fold_induction_starvation_response_genes.txt")
  #going to add section to make violin plots of fold induction
  #fold_induction_fac = cbind(fold_induction[,1],fold_induction[,2], fold_induction[,3], fold_induction[,4], fold_induction[,5])
  
   melt_fold = melt(fold_induction)
   colnames(melt_fold) = c("Condition", "Average_Fold_Induction")
  #fold_induction_fac = as.factor(fold_induction_fac)
   dodge <- positiofn_dodge(width = 0.4)
   ggplot(melt_fold, aes(factor(Condition), Average_Fold_Induction)) + geom_violin(aes(fill = Condition)) + geom_boxplot(width=.1, outlier.colour=NA, position = dodge) + scale_y_log10() + xlab("Condition") +
     ylab("Average Fold Induction") 
   
 
  
  a = (fold_induction > 5)
  WT_5_names = rownames(fold_induction[a[,1],])
  WT_7_names = rownames(fold_induction[a[,2],])
  dQ_names = rownames(fold_induction[a[,3],])
  HtoA_names = rownames(fold_induction[a[,4],])
  Null_names = rownames(fold_induction[a[,5],])

write.table(WT_5_names, "WT_5_names.txt", sep = "\t")
write.table(WT_7_names, "WT_7_names.txt", sep = "\t")
write.table(dQ_names, "dQ_Names.txt", sep = "\t")
write.table(HtoA_names,"HtoA_names.txt", sep = "\t")
write.table(Null_names, "Null_names.txt", sep = "\t")


fold_induction_ratio = cbind(fold_induction[,2]/fold_induction[,1], fold_induction[,3]/fold_induction[,1], fold_induction[,4]/fold_induction[,1], fold_induction[,5]/fold_induction[,1])

fold_induction_ratio[is.na(fold_induction_ratio)] = 0



fold_induction_ratio[(fold_induction_ratio) == Inf] = 0
colnames(fold_induction_ratio) =  colnames(fold_induction[,2:5])
rownames(fold_induction_ratio) = rownames(fold_induction)

subset_induction = fold_induction_ratio > 0.85

WT_pH5_cutoff = rownames(subset_induction)
pH7_cutoff = rownames(subset_induction)[subset_induction[,1]]
dQ_cutoff = rownames(subset_induction)[subset_induction[,2]]
HtoA_cutoff = rownames(subset_induction)[subset_induction[,3]]
Null_cutoff = rownames(subset_induction)[subset_induction[,4]]

89
19
5
11
23

v_left = read.table('/Users/HoltLab/Dropbox/ADH2_Switch_snf_RNAseq_12.28.17/volcano_leftside.txt', sep = "\t")
v_right = read.table('/Users/HoltLab/Dropbox/ADH2_Switch_snf_RNAseq_12.28.17/volcano_rightside.txt', sep = "\t")

v_left = as.character(v_left[,1])
v_right = as.character(v_right[,1])

t2g_sub_left = t2g$ext_gene %in% v_left
t2g_sub_right = t2g$ext_gene %in% v_right

left_genes = t2g[t2g_sub_left,]
right_genes = t2g[t2g_sub_right,]

write.table(left_genes, '/Users/HoltLab/Dropbox/ADH2_Switch_snf_RNAseq_12.28.17/volcano_annot_left.txt', sep = "\t")
write.table(right_genes, '/Users/HoltLab/Dropbox/ADH2_Switch_snf_RNAseq_12.28.17/volcano_annot_right.txt', sep = "\t")