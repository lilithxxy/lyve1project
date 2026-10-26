#####slingshot
library(Seurat)
library(cowplot)
library(slingshot)
library(RColorBrewer)
library(scales)
data <- readRDS('scored_final.rds')
sds <- slingshot(Embeddings(data, "umap"), clusterLabels = data$seurat_clusters, 
                 start.clus = 0, stretch = 0)
sce <- as.SingleCellExperiment(mdix)
sds <- slingshot(sce, clusterLabels = 'seurat_clusters',reducedDim = "UMAP",
                      allow.breaks = FALSE, start.clus = '3')
# Lineages and pseudotime ----
lnes <- getLineages(reducedDim(sds,"UMAP"),
                    sds$seurat_clusters, start.clus = '3')

lnes <- getLineages(reducedDim(sds,"PCA"),
                    sds$seurat_clusters, start.clus = '0')
my_color <- c('grey', 'purple', 
              'blue', 'yellow', 
              'red', 'green')
my_color <- c('red', '#bd8215', '#a3c208', 'green', '#05ebe3', '#059eeb', 'purple', '#ce07e8')

names(my_color) <- unique(as.character(sds$seurat_clusters))
plot(reducedDims(sds)$UMAP, col = my_color[as.character(sds$seurat_clusters)], 
     pch=16, 
     asp = 1)
legend("bottomleft",legend = names(my_color[levels(sds$seurat_clusters)]),  
       fill = my_color[levels(sds$seurat_clusters)])
lines(SlingshotDataSet(lnes), lwd=2, type = 'lineages', col = c("black"))
##
sce$slingPseudotime_1 <- sds$slingPseudotime_1
sce$slingPseudotime_2 <- sds$slingPseudotime_2
sce$slingPseudotime_3 <- sds$slingPseudotime_3
slingshot_df <- as.data.frame(colData(sce))
library(ggplot2)
library(ggbeeswarm)

ggplot(slingshot_df, aes(x = slingPseudotime_1, y = seurat_clusters, 
                         colour = seurat_clusters)) +
  geom_quasirandom(groupOnX = FALSE) + theme_classic() +
  xlab("First Slingshot pseudotime") + ylab("cell type") +
  ggtitle("Cells ordered by Slingshot pseudotime")+scale_colour_manual(values = my_color)
#
