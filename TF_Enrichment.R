#TF enrichment tool!! https://bioconductor.org/packages/release/bioc/vignettes/RcisTarget/inst/doc/RcisTarget.html
DEseq2_TF_enrichment = function(filepath='~/Dropbox (Personal)/Scripts/Nacho_RNAseq_2020_final/Ignacio_RNAseq_drop_snf5_del/text_files/WT_Starve_v_Dex.txt', filename='/WT_Starve_v_Dex', numbergenes = 100){
  library("RcisTarget")
  library("DESeq2")
  library("RColorBrewer")
  library("tidyverse")
  
  
  dds = readRDS('load_files/dds.rda') #from wherever you saved the data structure (should automatically save in prepare)----
  sampledata = readRDS('load_files/sampledata.rda')# snampledata - these both assume that these are in the current folder - can add to function w/ path if not the case
  t2g = readRDS('load_files/t2g.rda')
  # create filenames
  DE_Genes = read.table(filepath, sep = "\t", header = TRUE)
  
  
  DE_UP = rownames(DE_Genes[DE_Genes$log2FoldChange > 1,])
  DE_Down = rownames(DE_Genes[DE_Genes$log2FoldChange < 1,])

  
  data(motifAnnotations_hgnc)
  motifRankings <- importRankings("~/databases/hg19-tss-centered-10kb-7species.mc9nr.feather")
  
  # should create gene set then make a gene set collection so that we can do TF enrichment on all genes
  #
  # # Motif enrichment analysis:
  # motifEnrichmentTable_wGenes <- cisTarget(geneLists, motifRankings,
  #                                          motifAnnot=motifAnnotations_hgnc)
  # Motif enrichment analysis:
  motifEnrichmentTable_wGenes <- cisTarget(DE_UP, motifRankings,
                                           motifAnnot=motifAnnotations_hgnc)
  motifEnrichmentTable_wGenes <- addSignificantGenes(motifEnrichmentTable, 
                                                     geneSets=geneLists,
                                                     rankings=motifRankings, 
                                                     nCores=8,
                                                     method="iCisTarget")
  
  
  print('done')
}
  