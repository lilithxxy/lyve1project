# HTODemux hashtag demultiplexing (CLR, positive.quantile 0.99), QC and merge of the timecourse samples.
library(Seurat)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)   

# Demultiplex the first sample ----
both <- Read10X(data.dir = "~/Desktop/bioinfo/singlecell_human/tmp/R1_DMo_D0M/outs/filtered_feature_bc_matrix")
seu.obj = CreateSeuratObject(counts = both$'Gene Expression')

#Normalize RNA and find and scale variable feature
seu.obj <- NormalizeData(seu.obj)
seu.obj <- FindVariableFeatures(seu.obj, selection.method = "mean.var.plot")
seu.obj <- ScaleData(seu.obj, features = VariableFeatures(seu.obj))

#Normalize HTO data
seu.obj[['HTO']] = CreateAssayObject(counts = both$'Antibody Capture')
seu.obj <- NormalizeData(seu.obj, assay = "HTO", normalization.method = "CLR")
seu.obj <- HTODemux(seu.obj, assay = "HTO", positive.quantile = 0.99)

View(seu.obj@meta.data)

filtered <- subset(seu.obj, subset = HTO_classification != 'Negative' & HTO_classification.global == 'Singlet')

View(filtered@meta.data)
obj.list <- SplitObject(filtered, split.by = "HTO_classification")

#decontamination
library(decontX)
# Demultiplex the second sample ----
# Note: decontX is loaded but never called; no ambient-RNA correction applied.
#Normalize RNA and find and scale variable feature
counts1 <- Read10X(data.dir = "~/Desktop/bioinfo/singlecell_human/tmp/R1_D5MDI_D5A/outs/filtered_feature_bc_matrix")
unfiltered = CreateSeuratObject(counts = counts1$'Gene Expression')
unfiltered <- NormalizeData(unfiltered)
unfiltered <- FindVariableFeatures(unfiltered, selection.method = "mean.var.plot")
unfiltered <- ScaleData(unfiltered, features = VariableFeatures(unfiltered))

#Normalize HTO data
unfiltered[['HTO']] = CreateAssayObject(counts = counts1$'Antibody Capture')
unfiltered <- NormalizeData(unfiltered, assay = "HTO", normalization.method = "CLR")
unfiltered <- HTODemux(unfiltered, assay = "HTO", positive.quantile = 0.99)

View(unfiltered@meta.data)
raw.seu <- subset(unfiltered, subset = HTO_classification == 'Negative')
filtered <- subset(unfiltered, subset = HTO_classification != 'Negative' & HTO_classification.global == 'Singlet')

obj.list <- SplitObject(filtered, split.by = "HTO_classification")

DefaultAssay(day5mdix) <- 'RNA'
day5mdix[["percent.mt"]] <- PercentageFeatureSet(day5mdix, pattern = "^MT-")
VlnPlot(day5mdi, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)



# Cell-level QC filters ----
filtered <- subset(filtered, subset = nCount_RNA<50000 & nFeature_RNA > 500 & nFeature_RNA < 8000 & percent.mt < 15)

plot1 <- FeatureScatter(filtered, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

saveRDS




###
for (i in obj.list) {
  hash_id <- i@meta.data$hash.ID[1] 
  file_name <- paste0(hash_id, ".rds")
  saveRDS(i, file = file_name)
}

seu.obj <- readRDS('240316-Donor2-MoMac-Day2-MDI.rds')
day4mdi[["percent.mt"]] <- PercentageFeatureSet(day4mdi, pattern = "^MT-")
VlnPlot(day4mdi, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(seu.obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seu.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

seu.obj.filtered <- subset(seu.obj, subset = nCount_RNA<60000 & nFeature_RNA > 300 & nFeature_RNA < 8000 & percent.mt < 25)
saveRDS(seu.obj.filtered, 'day1mdix')
#####
seu.obj <- readRDS('Momac_Day0_M.rds')
seu.obj@meta.data$timepoint <- 'Day0'
seu.obj@meta.data$condition <- 'M'

saveRDS(seu.obj, file = 'Momac_Day0_M.rds')

# Merge per-sample objects ----
#####
rds.path <- "~/Desktop/bioinfo/sc_R/unfiltered/all/"  
rds_files <- list.files(path = rds.path, pattern = "\\.rds$", full.names = TRUE)
obj.list <- list()
for (file in rds_files) {
  seurat_obj <- readRDS(file)  
  timepoint <- gsub(".*day(\\d+).*", "\\1", basename(file))
  seurat_obj$timepoint <- paste0("day", timepoint)
  obj.list[[file]] <- seurat_obj 
}
merged <- merge(obj.list[[1]], obj.list[2:6])

DefaultAssay(merged) <- 'RNA'
merged[["HTO"]] = NULL
# merged[["RNA"]] <- split(merged[["RNA"]], f = merged$timepoint)
merged$condition <- 'MDIX'
rna.list = SplitObject(merged, split.by = "hash.ID")
for (i in 1:length(rna.list)) {
  rna.list[[i]] = FindVariableFeatures(rna.list[[i]])
}
rna.list = SelectIntegrationFeatures(rna.list, nfeatures = 5000)
#merged = subset(merged, features = grep("MT-", rownames(merged), invert = T, value = T))
merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged, selection.method = 'vst')
merged <- ScaleData(merged, features = rna.list)
merged <- RunPCA(merged, features = rna.list)
merged <- FindNeighbors(merged, dims = 1:20)
merged <- FindClusters(merged, resolution = 0.6)
merged <- RunUMAP(merged, dims = 1:20)

DimPlot(merged, reduction = 'umap', group.by = 'hash.ID', label = F)
DimPlot(merged, reduction = 'umap', group.by = 'seurat_clusters', label = T, split.by = 'condition')
FeaturePlot(merged, reduction = 'umap', features = c('CD14', 'S100A8', 'FCGR3A'))
FeaturePlot(allwithscores, reduction = 'pca', 
            features = c('NUPR1', 'LYVE1', 'CROCC', 
                         'FABP4', 'ADGRL2', 'AC108134.2', 'TMIGD3', 'PMP22', 'CCL22', 'RYR1', 
                         'C1QA', 'ARMH1', 'IFIT2', 'Z93241.1', 'AL513283.1', 'F13A1', 'OSTN',
                         'ABCA6' ,
                         'SH2D4A', 'MPPED2', 'CCL13', 'HIST1H4E', 'CSGALNACT1', 
                         'SCN9A', 'CCDC141', 'CCL18', 'FBLN5'),
             ncol = 5)
p3 <- DimPlot(merged, group.by = 'hash.ID')
p1 + p3
p1 + p2
merged <- IntegrateLayers(
  object = merged, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)
merged <- FindNeighbors(merged, reduction = "integrated.cca", dims = 1:30)
merged <- FindClusters(merged, resolution = 1, cluster.name = "cca_clusters")
merged <- RunUMAP(merged, reduction = "integrated.cca", dims = 1:30)
DimPlot(merged, reduction = 'umap', group.by = 'seurat_clusters', label = T)
DimPlot(merged, reduction = 'umap', group.by = 'hash.ID', label = T)
FeaturePlot(merged, reduction = 'umap', features = c('C1QA', 'C1QC', 'S100A8', 'LYVE1', 'F13A1', 'FCGR3A', 'CD14',
                                                               'SELENOP', 'MRC1'))
FeaturePlot(allwithscores, reduction = 'umap', features = c('C1QA', 'C1QC', 'S100A8', 'LYVE1', 'F13A1', 'FCGR3A', 'CD14',
                                                     'SELENOP', 'MRC1', 'APOE'))
VlnPlot(merged, features = c('LYVE1', 'F13A1', 'FOLR2', 'MRC1', 'C1QA'), group.by = 'seurat_clusters')
p1 + p2
merged <- JoinLayers(merged)
rmerged
markers <- FindMarkers(merged, ident.1 = '240316-Donor2-MoMac-Day5-MDI', ident.2 = '240316-Donor2-MoMac-Day5-MDIX', group.by = 'hash.ID')
DimPlot(merged, reduction = 'umap', group.by = 'hash.ID')
p1 <- DimPlot(merged, reduction = 'umap', group.by = 'seurat_clusters', label = T)
p2 <- DimPlot(merged, group.by = 'timepoint')
DimPlot(merged, group.by = 'condition')
p1 + p2 
FeaturePlot(merged, features = c('LYVE1'), min.cutoff = "q10", max.cutoff = "q90")
VlnPlot(allwithscores, features = 'LYVE1', group.by = 'seurat_clusters')

FeaturePlot(merged, features = c('HLA-DPA1'), 
            split.by = "hash.ID", ncol = 3)


mm = Seurat::FindMarkers(merged, ident.1 = c(7,8), ident.2 = 0)
l1pos <- top500_invitrol1$...1
l1neg <- last500_invitrol1$...1
l1pos_present <- l1pos[l1pos %in% rownames(merged)]
l1neg_present <- l1neg[l1neg %in% rownames(merged)]
l1pos_present <- list(l1pos_present)
l1neg_present <- list(l1neg_present)

merged <- AddModuleScore(
  object = merged,
  features = l1pos_present,
  name = "l1pos"
)

merged <- AddModuleScore(
  object = merged,
  features = l1neg_present,
  name = "l1neg"
)
FeaturePlot(merged, features = c('l1pos1', 'l1neg1'), split.by = 'condition')



library(ggplot2)
library(dplyr)
library(magrittr)
library(patchwork)
library(viridis)
library(Seurat)
library(scCustomize)

DotPlot_scCustom(seurat_object = merged, features = c('LYVE1', 'F13A1', 'FOLR2', 'MRC1', 'CD68'), colors_use = viridis_plasma_dark_high,
                 group.by = 'hash.ID')

Clustered_DotPlot(seurat_object = merged, features =  c('LYVE1', 'F13A1', 'FOLR2', 'MRC1'), group.by = 'hash.ID')
Clustered_DotPlot(seurat_object = allwithscores, features =  c('LYVE1', 'F13A1', 'FOLR2', 'MRC1'), group.by = 'seurat_clusters')



mdiwithdcinfo <- CreateDimReducObject(mdiwithdcinfo, key = 'DC_', )
# Find markers and limit to those expressed in greater than 75% of target population
all_markers <- FindAllMarkers(object = merged) %>%
  Add_Pct_Diff() %>%
  filter(pct_diff > 0.6)

top_markers <- Extract_Top_Markers(marker_dataframe = all_markers, num_genes = 7, named_vector = FALSE,
                                   make_unique = TRUE)

Clustered_DotPlot(seurat_object = merged, features = top_markers)

# pseudobulk cells only by cell typ
stats <- Cluster_Stats_All_Samples(merged, group_by_var = 'cell.type')
mm <- FindMarkers(merged, ident.1 = 1, ident.2 = 3)
#######
dc1 <- mdiwithdcinfo$mdi_dc1
dc2 <- mdiwithdcinfo$mdi_dc2

# Create a matrix for the new embeddings with row names
new_embedding <- cbind(DC1 = dc1, DC2 = dc2)
rownames(new_embedding) <- colnames(mdiwithdcinfo)
mdiwithdcinfo[['diffusionmap']] <- CreateDimReducObject(embeddings = new_embedding, key = "DC_", assay = DefaultAssay(mdiwithdcinfo))
###### heatmap######

markers <- FindAllMarkers(allwithscores, only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(allwithscores, features = top10$gene) + NoLegend()
##### 
FeaturePlot(allwithscores, features = c('CD34', 'KIT', 'MYB', 'AZU1', 'MPO', 'MS4A3',
                                        'LYZ', 'CCR2', 'HLA-DRA', 'CSF1R', 'MEF2C', 'CD14'))
f <- FeaturePlot(mdixwithdc, features = c('CCR2','CCL2', 'CCL3', 'CCL4', 'IL4R', 'IL6', 'CXCL8', 'TNF',
                                        'IL1B', 'IL3RA', 'CSF1R', 'HS3ST2', 'SELENOP', 'MRC1', 'IGF1', 'MMP9',
                                        'F13A1', 'TREM2', 'APOE'), reduction = 'diffusionmap') +xlim(-0.05, 0.05) +  ylim(-0.05, 0.05)
 
print(f)

FeaturePlot(allwithscores, features = c('CCR2','CCL2', 'CCL3', 'CCL4', 'IL4R', 'IL6', 'CXCL8', 'TNF',
                                        'IL1B', 'IL3RA', 'CSF1R', 'HS3ST2', 'SELENOP', 'MRC1', 'IGF1', 'MMP9'))

FeaturePlot(allwithscores, features = c('FOLR2'))
DoHeatmap(allwithscores, features = c('CCR2','CCL2', 'CCL3', 'CCL4', 'IL4R', 'IL6', 'CXCL8', 'TNF',
                                      'IL1B', 'IL3RA', 'CSF1R', 'HS3ST2', 'SELENOP', 'MRC1', 'IGF1', 'MMP9',
                                      'F13A1', 'TREM2', 'APOE')) + NoLegend()





#####

allwithscores <- CellCycleScoring(allwithscores, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, 
                                  set.ident = TRUE)
allwithscores <- RunPCA(allwithscores, features = c(cc.genes$s.genes, cc.genes$g2m.genes))
DimPlot(allwithscores, split.by = 'Phase', group.by = 'condition')
# view cell cycle scores and phase assignments
head(marrow[[]])



