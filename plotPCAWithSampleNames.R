plotPCAWithSampleNames = function(x, intgroup="condition", ntop=100)
{
 
  library(RColorBrewer)
  library(genefilter)
  library(lattice)
  
  # pca
  rv = rowVars(assay(x))
  select = order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(x)[select,]))
  
  # proportion of variance
  variance = pca$sdev^2 / sum(pca$sdev^2)
  variance = round(variance, 3) * 100
  
  # sample names
  #names = colnames(x)
  names = as.character(x$sample)
  
  # factor of groups
  fac = factor(apply(as.data.frame(colData(x)[, intgroup, drop=FALSE]), 1, paste, collapse=" : "))
  #fac = factor(fac, levels = c("WT_Dex","dQ_snf5_Dex","HtoA_snf5_Dex","WT_Stv_pH5","dQ_snf5_Stv_pH5","HtoA_snf5_Stv_pH5","WT_Stv_pH7"))
  #General_palette = c("black","#00a746", "#0065b8","#bd5112","red","#ffb800", "#652a96")
  # colors
  if( nlevels(fac) >= 10 )
    colors = rainbow(nlevels(fac))
  else if( nlevels(fac) >= 3 )
    colors =  brewer.pal(nlevels(fac), "Set1") #General_palette
  else
    colors = c( "dodgerblue3", "firebrick3" )
  
  # plot
  xyplot(
    PC2 ~ PC1, groups=fac, data=as.data.frame(pca$x), pch=16, cex=1.5,
    aspect = "fill",
    col = colors,
    xlab = list(paste("PC1 (", variance[1], "%)", sep=""), cex=0.8),
    ylab = list(paste("PC2 (", variance[2], "%)", sep=""), cex=0.8),
    panel = function(x, y, ...) {
      panel.xyplot(x, y, ...);
      ltext(x=x, y=y,labels=names,  pos=1, offset=0.8, cex=0.7) #
    },
    main = draw.key(
      key = list(
        rect = list(col = colors),
        text = list(levels(fac)),
        rep = FALSE
      )
    )
  )
}