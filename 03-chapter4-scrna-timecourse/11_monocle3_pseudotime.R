# Monocle3 pseudotime on the merged Seurat object via as.cell_data_set.
library(Seurat)
library(dplyr)
library(SeuratDisk)
library(tidyverse)
merged <- readRDS('final.rds')
merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged, selection.method = 'vst')
merged <- ScaleData(merged)
merged <- RunPCA(merged)
merged <- FindNeighbors(merged, dims = 1:10)
merged <- FindClusters(merged, resolution = 0.6)
merged<- RunUMAP(merged, dims = 1:10)
DimPlot(merged, reduction = 'pca', label = T, label.size = 5)
library(Nebulosa)
library(scCustomize)
genes <- c('C1QA', 'C1QC', 'APOE','CD68','F13A1', 'MRC1', 'LYVE1', 'FCGR3A', 'CD14', 'S100A8')
p2 <- FeaturePlot(rna.data, features = genes)
p1 <- DimPlot(rna.data, label = T)
p1 + p2
DotPlot_scCustom(seurat_object = rna.data, features = genes, colors_use = viridis_plasma_dark_high)
library(monocle3)
library(SeuratWrappers)
cds <- as.cell_data_set(merged)
fData(cds)$gene_short_name <- rownames(fData(cds))
counts(cds)
recreate.partition <- c(rep(1,length(cds@colData@rownames)))
names(recreate.partition) <- cds@colData@rownames
recreate.partition <- as.factor(recreate.partition)
cds@clusters$UMAP$partitions <- recreate.partition
Idents(merged)
list_cluster <- merged@active.ident
cds@clusters$UMAP$clusters <- list_cluster
cds@int_colData@listData$reducedDims$UMAP <- merged@reductions$umap@cell.embeddings

cds <- cluster_cells(cds)
cds <- learn_graph(cds, use_partition = F)

plot_cells(cds,reduction_method = 'UMAP',
           color_cells_by = 'cluster',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 5)

cds <- order_cells(cds, reduction_method = 'UMAP')
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=T,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
p1 + p2
plot_cells(cds, genes=c('CCR2','CCL2', 'CCL3', 'CCL4', 'IL4R', 'IL6', 'CXCL8', 'TNF',
                                   'IL1B', 'IL3RA', 'CSF1R', 'HS3ST2', 'SELENOP', 'MRC1', 'IGF1', 'MMP9',
                                   'F13A1', 'TREM2', 'APOE'),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)

###########addmodulescore#######
l1list <- list(l1genes_filtered)
mhclist <- list(mhcgenes_filtered)

merged <- AddModuleScore(
  object = merged,
  features = l1list,
  name = "l1score")

merged <- AddModuleScore(
  object = merged,
  features = mhclist,
  name = "mhcscore")
FeatureScatter(merged, feature1 = 'l1score1', feature2 = 'mhcscore1')
FeaturePlot(merged, features = 'l1score')
# Create Plots
FeatureScatter_scCustom(seurat_object = merged, feature1 = "l1score1", feature2 = "mhcscore1",
             num_columns = 2, pt.size = 1)
FeaturePlot(merged, features = 'mhcscore1')
VlnPlot(merged, features = c('l1score1', 'mhcscore1'), group.by = 'seurat_clusters')
###########

# 0, immature; 1,2, intermediate; 4,7LYVE1hi; 6, postLYVE1hi; 5, MDIX; 3, MDI
cluster.ids.new <- c('immature', 'intermediate', 'intermediate', 'MDI', 'LYVE1hi',
                     'MDIX', 'postLYVE1hi', 'LYVE1hi')
names(cluster.ids.new) <- levels(merged)
merged <- RenameIdents(merged, cluster.ids.new)
merged$cell.type <- Idents(merged)
Idents(merged)
mm <- FindAllMarkers(merged)
top_markers <- Extract_Top_Markers(marker_dataframe = mm, num_genes = 7, named_vector = FALSE,
                                   make_unique = TRUE)

Clustered_DotPlot(seurat_object = merged, features = top_markers)



FeaturePlot(merged, features = 
              c('LYVE1', 'C1QC', 'CD14', 'FCGR3A', 'F13A1', 'CD68', 'APOE', 'C1QA', 'S100A8'), reduction = 'pca')
