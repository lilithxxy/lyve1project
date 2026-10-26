# Customised PCA plot for the bulk samples (prcomp on top-variable genes).
library(dplyr)
library(ggplot2)
library(ggalt)
library(ggforce)

# Helper to build the PCA data frame ----
getPCs = function(object, intgroup="condition",
                                  ntop=500, returnData=FALSE, pcsToUse=1:2)
{
  message(paste0("using ntop=",ntop," top features by variance"))
  
  # calculate the variance for each gene
  rv <- rowVars(assay(object))
  
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))
  
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
  
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=":"))
  } else {
    colData(object)[[intgroup]]
  }
  
  # assembly the data for the plot
  pcs <- paste0("PC", pcsToUse)
  d <- data.frame(V1=pca$x[,pcsToUse[1]],
                  V2=pca$x[,pcsToUse[2]],
                  group=group, intgroup.df, name=colnames(object))
  colnames(d)[1:2] <- pcs
  
  if (returnData) {
    attr(d, "percentVar") <- percentVar[pcsToUse]
    return(d)
  }
}
# Run PCA ----
pc.df <- getPCs(vsdata, intgroup = 'condition', ntop = 500, returnData = T)


rv <- rowVars(assay(vsdata))

# select the ntop genes by variance
select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]

# perform a PCA on the data in assay(x) for the selected genes
pca <- prcomp(t(assay(vsdata)[select,]))

# the contribution to the total variance for each component
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )





mycolours <- c("red", "blue", "#016124", "#e36002", "purple")
names(mycolours) <- levels(pc.df$group)
colScale <- scale_colour_manual(name = "group",values = mycolours)
# Plot ----
ggplot(data = pc.df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3) +
  geom_mark_ellipse(expand = 0.05,aes(fill=group)) + 
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  colScale +
  theme_classic() + 
  ggtitle('PCA Plot for Bulk Samples') 
   coord_fixed()


