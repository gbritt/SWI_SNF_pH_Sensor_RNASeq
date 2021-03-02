#DESeq2_Heatmaps_ect
#Will need to load in .csv files from DESeq2_Statistical and datastructure from DESeq2_prepare



DESeq2Clusterplot = function(filepath='~/Dropbox (Personal)/Scripts/Nacho_RNAseq_2020_final/Ignacio_RNAseq_drop_snf5_del/text_files/WT_Starve_v_Dex.txt', filename='/WT_Starve_v_Dex', numbergenes = 40){
   # Heatmap cutting = https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/
    library("pheatmap")
   # should think about using "copmlex heatmap in future perhaps?"
    library("DESeq2")
    library("RColorBrewer")
   library("tidyverse")
   
    dds = readRDS('load_files/dds.rda') #from wherever you saved the data structure (should automatically save in prepare)----
    sampledata = readRDS('load_files/sampledata.rda')# sampledata - these both assume that these are in the current folder - can add to function w/ path if not the case
    t2g = readRDS('load_files/t2g.rda')
    # create filenames
    filename = read.table(filepath, sep = "\t", header = TRUE)
    filenameup = paste(filename,"upinstarve","_",numbergenes, sep = "") 
    filenamedown = paste(filename,"downinstarve","_",numbergenes, sep = "") 
    
    filename_dQ = read.table("~/Dropbox (Personal)/Scripts/Nacho_RNAseq_2020_final/Ignacio_RNAseq_drop_snf5_del/text_files/dQ_Starve_v_Dex.txt", sep = "\t", header = TRUE)
    filenameup_dQ = paste(filename_dQ,"upinstarve","_",numbergenes, sep = "") 
    filenamedown_dQ = paste(filename_dQ,"downinstarve","_",numbergenes, sep = "")
    
    
    filename_HtoA = read.table("~/Dropbox (Personal)/Scripts/Nacho_RNAseq_2020_final/Ignacio_RNAseq_drop_snf5_del/text_files/HtoA_Starve_v_Dex.txt", sep = "\t", header = TRUE)
    filenameup_HtoA = paste(filename_HtoA,"upinstarve","_",numbergenes, sep = "") 
    filenamedown_HtoA = paste(filename_HtoA,"downinstarve","_",numbergenes, sep = "") 
   
    
    filename_pH7 = read.table("~/Dropbox (Personal)/Scripts/Nacho_RNAseq_2020_final/Ignacio_RNAseq_drop_snf5_del/text_files/WT_Starve_v_Dex_pH_7.txt", sep = "\t", header = TRUE)
    filenameup_pH7 = paste(filename_pH7,"upinstarve","_",numbergenes, sep = "") 
    filenamedown_pH7 = paste(filename_pH7,"downinstarve","_",numbergenes, sep = "") 
    
    filename_dQ_WT_S = read.table("~/Dropbox (Personal)/Scripts/Nacho_RNAseq_2020_final/Ignacio_RNAseq_drop_snf5_del/text_files/dQ_v_WT_Starve.txt", sep = "\t", header = TRUE)
    filenameup_dQ_WT_S = paste(filename_dQ_WT_S,"upinstarve","_",numbergenes, sep = "") 
    filenamedown_dQ_WT_S = paste(filename_dQ_WT_S,"downinstarve","_",numbergenes, sep = "") 
    
    filename_dQ_WT_D = read.table("~/Dropbox (Personal)/Scripts/Nacho_RNAseq_2020_final/Ignacio_RNAseq_drop_snf5_del/text_files/dQ_v_WT_Dex.txt", sep = "\t", header = TRUE)
    filenameup_dQ_WT_D = paste(filename_dQ_WT_D,"upinstarve","_",numbergenes, sep = "") 
    ilenamedown_dQ_WT_D = paste(filename_dQ_WT_D,"downinstarve","_",numbergenes, sep = "") 
    
    
    #filter SigGene lists into up and downregulated genes and cut by number desired ----
    filenamedown = filename[order(filename$log2FoldChange),][1:numbergenes,]
    filenameup = filename[order(filename$log2FoldChange, decreasing = TRUE),][1:numbergenes,]
    
    filenamedown_dQ = filename_dQ[order(filename_dQ$log2FoldChange),][1:numbergenes,]
    filenameup_dQ = filename_dQ[order(filename_dQ$log2FoldChange, decreasing = TRUE),][1:numbergenes,]
    
    filenamedown_HtoA = filename_HtoA[order(filename_HtoA$log2FoldChange),][1:numbergenes,]
    filenameup_HtoA = filename_HtoA[order(filename_HtoA$log2FoldChange, decreasing = TRUE),][1:numbergenes,]
    
    filenamedown_pH7 = filename_pH7 [order(filename_pH7 $log2FoldChange),][1:numbergenes,]
    filenameup_pH7 = filename_pH7 [order(filename_pH7 $log2FoldChange, decreasing = TRUE),][1:numbergenes,]
    
    filenamedown_dQ_WT_S = filename_dQ_WT_S[order(filename_dQ_WT_S$log2FoldChange),][1:numbergenes,]
    filenameup_dQ_WT_S = filename_dQ_WT_S[order(filename_dQ_WT_S$log2FoldChange, decreasing = TRUE),][1:numbergenes,]
    
    filenamedown_dQ_WT_D = filename_dQ_WT_D[order(filename_dQ_WT_D$log2FoldChange),][1:numbergenes,]
    filenameup_dQ_WT_D = filename_dQ_WT_D[order(filename_dQ_WT_D$log2FoldChange, decreasing = TRUE),][1:numbergenes,]
    
    
    
    #create master list of all genes you selected ----
    select = c(rownames(filenameup), rownames(filenamedown),rownames(filenamedown_dQ),rownames(filenameup_dQ),rownames(filenamedown_HtoA),rownames(filenameup_HtoA),rownames(filenamedown_pH7),rownames(filenameup_pH7))#,rownames(filenamedown_dQ_WT_S),rownames(filenameup_dQ_WT_S),rownames(filenamedown_dQ_WT_D),rownames(filenameup_dQ_WT_D)
   
    select = unique(select)
    #create transformation to have a more reasonable distribution of data----
    ntd <- normTransform(dds) # can probably play with different normalizations here 
   #ntd = vst(dds)
    #ntd = rlog(dds)
   
   cluster = assay(ntd)[select,]
   cluster = cluster - rowMeans(cluster)# making matrix with gene names and gene counts
   colnames(cluster)= as.character(sampledata$sample) #giving colnames to clustered object using the 'sampledata' variable created in the prepare script
   
   # Adding common names to ens_ids for genes ----
   shared_results <- t2g[t2g[,2] %in% rownames(cluster),] # Select rownames in the 'gene mapping' matrix that are also in our clustering matrix
   a = as.matrix(shared_results)
   ordered_genes_mapping <- a[ order(a[,1]), ]# ordering genes mapping so that the rownames in cluster to be alphabetical
   cluster <- cluster[ order(row.names(cluster)), ] # order genes in cluster to be alphabetical (will be important below)
   
   common_and_ens_ids = matrix()  # create matrix to hold a pasted name containing common and ens_id names for a gene
   
   
   for(i in 1:length(ordered_genes_mapping[,1])){
       common_and_ens_ids[i] = paste(ordered_genes_mapping[i,3], "_",ordered_genes_mapping[i,1], sep = "")
   }# paste together each common name with the ens_id (will make heatmap easier to read)
  
   rownames(cluster) = common_and_ens_ids # now our cluster rownames contain common names if available
   
   #create annotations for top of heatmap ----
   df <- as.data.frame(colData(dds)[,c("Condition","GrowthCondition")])
   rownames(df) <- colnames(cluster)
   
   
   quantile_breaks <- function(cluster, n = 10) {
      breaks <- quantile(cluster, probs = seq(0, 1, length.out = n))
      breaks[!duplicated(breaks)]
   }
   
   mat_breaks <- quantile_breaks(cluster, n = 100)
   breaksList = seq(-3,3, by = .1)
   #my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 99)
   
   cluster_ordered = cluster[,c(1,6,15,2,7,16,3,8,17,4,9,18,10,13,19,11,14,20,5,12,21)]
   
  Heatmap = pheatmap(cluster_ordered,  border_color = NA,cluster_rows=TRUE, show_rownames=T,
            cluster_cols=F, show_colnames = T, annotation_col = df,  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),  breaks = breaksList, drop_levels = F, cutree_rows = 16)
  Heatmap
  heatmaprows = rownames(cluster_ordered[Heatmap$tree_row[["order"]],])
  #scaling by row!
  #https://rockefelleruniversity.github.io/RU_VisualizingGenomicsData/viz_course/Answers/Viz_part3_answers.html
  
  saveRDS(Heatmap, "Heatmap.RDA")
   write.table(heatmaprows, "Genes_heatmap_Garb.txt", sep = "\t",col.names = F, row.names = F)
}


 



