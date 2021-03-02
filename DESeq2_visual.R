#visually exploring sample relationships, in  transform counts for computing distances or making plots
# For visual analysis we are transforming so:
#assuming you prepped the dds using Kallisto_to_DESeq2_prepare.R
library('pheatmap')
library('RColorBrewer')
library('DESeq2')

dds = readRDS('load_files/dds.rda')
sampledata = readRDS('load_files/sampledata.rda')
vsd <- vst(dds)
vsd$Condition = factor(vsd$Condition, levels =c("WT_Dex","dQ_snf5_Dex","HtoA_snf5_Dex","WT_Stv_pH5","dQ_snf5_Stv_pH5","HtoA_snf5_Stv_pH5","WT_Stv_pH7"))

General_palette = c("black","#00a746", "#0065b8","#bd5112","red","#ffb800", "#652a96")
res <- results(dds,contrast=c("GrowthCondition","Stv","Dex"), cooksCutoff=FALSE, independentFiltering=FALSE)

#resLFC <- lfcShrink(dds, coef=2) for shrinkage in plotting heatmaps and PCA
source('scripts/plotPCAWithSampleNames.R')
z = plotPCAWithSampleNames(vsd, intgroup=c("Condition"))
z
#can probably find a way to look at other principal components
#can also do regular plotPCA for prettier graph

#sample - sample variances in heatmap form ----

dds = readRDS('dds.rda')
sampledata = readRDS('sampledata.rda')
vsd <- vst(dds)
# sample to sample distances ----
sampleDists <- dist(t(assay(vsd)))

library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Condition,  sep="-")
colnames(sampleDistMatrix) <- paste(vsd$Condition,  sep="-")
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors, rownames = T)

# Sample to sample scatter plots ----
scat_matrix = assay(vsd)
colnames(scat_matrix) <- paste(vsd$Condition,  sep="-")
lev = levels(vsd$Condition)

avg_matrix = matrix(,nrow = nrow(scat_matrix), ncol = length(lev))

for(column in 1:length(lev)){
  
  sub = colnames(scat_matrix) == lev[column]
  tmp = scat_matrix[,sub]
  tmp1 = as.matrix(rowSums(tmp)/ncol(tmp))
  avg_matrix[,column] = tmp1
  colnames(avg_matrix) = lev
}

plotFun <- function(x,y){ 
  dns <- densCols(x,y); 
  points(x,y, col=dns, pch=".", panel.first=grid());  
  #  abline(a=0, b=1, col="brown")
}

epsion = 1
nb.pairs <- 9
pairs(log2(avg_matrix[,sample(ncol(avg_matrix), nb.pairs)] + epsilon), 
      panel=plotFun, lower.panel = NULL)



## venn diagram

library('gplots')
lists = c(WT_5_names,WT_7_names,dQ_names,HtoA_names,Null_names)
venndiagram = venn(lists)

lists = list(WT_pH5_cutoff,pH7_cutoff,dQ_cutoff,HtoA_cutoff,Null_cutoff)
