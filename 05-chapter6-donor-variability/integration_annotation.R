# Merge, QC and manually annotate the donor PBMC dataset.
setwd("/home/yc/Desktop/monos")
library(Seurat)
library(dplyr)
library(ggplot2)
path <- "/home/yc/Desktop/monos"
pbmc_files <- list.files(path, pattern = "PBMC.*\\.rds$", full.names = TRUE)
pbmc_list <- lapply(pbmc_files, readRDS)
names(pbmc_list) <- basename(pbmc_files)
for (i in seq_along(pbmc_list)) {
  pbmc_list[[i]]@assays$HTO <- NULL
  pbmc_list[[i]]@assays$RNA@layers$data <- NULL
  pbmc_list[[i]]@assays$RNA@layers$scale.data <- NULL
}
integrated <- merge(pbmc_list[[1]], y = pbmc_list[2:4])
integrated <- JoinLayers(integrated)
remove(pbmc_list)
integrated[["percent.mt"]] <- PercentageFeatureSet(integrated, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
##
plot1 <- FeatureScatter(integrated, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(integrated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#33538 features across 15123 samples within 1 assay 
integrated <- subset(integrated, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & nCount_RNA < 80000 & percent.mt < 10)
#33538 features across 9315 samples within 1 assay 
table(integrated$hash.ID)
DefaultAssay(integrated) <- "RNA"
integrated <- NormalizeData(integrated, verbose = FALSE)
data("cc.genes.updated.2019", package = "Seurat")
s.genes   <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
s.genes   <- s.genes[s.genes %in% rownames(integrated)]
g2m.genes <- g2m.genes[g2m.genes %in% rownames(integrated)]

integrated <- CellCycleScoring(
  object        = integrated,
  s.features    = s.genes,
  g2m.features  = g2m.genes,
  set.ident     = TRUE
)
FeatureScatter(integrated, feature1 = "S.Score", feature2 = "G2M.Score")
plot1 + plot2
# Note: vars.to.regress left empty; default applies and nothing is regressed.
integrated <- SCTransform(
  integrated,
  verbose = FALSE,
  vars.to.regress = 
)

DefaultAssay(integrated) <- "SCT"

integrated <- RunPCA(integrated, npcs = 50, verbose = FALSE)
integrated <- RunUMAP(integrated, dims = 1:30, verbose = FALSE)
integrated <- FindNeighbors(integrated, dims = 1:30, verbose = FALSE)
integrated <- FindClusters(integrated, resolution = 0.4, verbose = FALSE)
DimPlot(integrated, reduction = "umap", label = TRUE, group.by = "Phase",
        pt.size = 1) 
FeaturePlot(integrated, features = c(
  'CD14', 'FCGR3A', 'S100A8', 'LYZ'
))

donor_map <- c(
  "240106-Donor1-PBMC-L1-15" = "donor_15",
  "240120-Donor1-PBMC-L1-64" = "donor_64",
  "240229-Donor1-PBMC-L1"    = "donor_6",
  "240309-Donor2-PBMC-L1"    = "donor_90"
)
integrated$donor_name <- unname(donor_map[integrated$hash.ID])
# Map donor_name to efficacy
efficacy_map <- c(
  "donor_15" = "low",
  "donor_6"  = "low",
  "donor_64" = "high",
  "donor_90" = "high"
)

integrated$efficacy <- unname(efficacy_map[integrated$donor_name])

FeaturePlot(integrated, features =c(
  'CD3G', 'CCR7', 'IL7R', 'S100A4', 
  'CD4', 'CD8A', 
  'MS4A1', 'CD19',
  'CD14','S100A8', 'LYZ',
  'FCGR3A', 
  'FCER1A', 'CST3',
  'GNLY', 'NKG7', 'NCAM1',
  'PPBP'
))
DefaultAssay(integrated) <- 'RNA'
DefaultAssay(integrated) <- 'SCT'
VlnPlot(integrated, features = c(
  'CD3G', 'CCR7', 'IL7R', 'S100A4', 
  'CD4', 'CD8A', 
  'MS4A1', 'CD19',
  'CD14','S100A8', 'LYZ',
  'FCGR3A',  'C1QA', 'C1QC',
  'FCER1A', 'CST3',
  'GNLY', 'NKG7', 'NCAM1', 'KLRD1', 'KLRF1',
  'PPBP'
))
### 0: CD14+ Mono; 1: CD56+ NK; 2: CD4 + T; 3 : CD56 + NK; 4 : CD4 + T; 5 : CD8+ T; 
### 6: B; 7: CD16+ Mono; 8: CD8+ T; 9: B; 10: CD56- NK;
### 11: CD8+ T; 12: DC; 13: DC; 14: Debris; 15: Platelet
new.cluster.ids <- c(
  "CD14+ Mono",     # 0
  "CD56+ NK",       # 1
  "CD4+ T",         # 2
  "CD56+ NK",       # 3
  "CD4+ T",         # 4
  "CD8+ T",         # 5
  "B",              # 6
  "CD16+ Mono",     # 7
  "CD8+ T",         # 8
  "B",              # 9
  "CD56- NK",       #10
  "CD8+ T",         #11
  "DC",             #12
  "DC",             #13
  "Debris",         #14
  "Platelet"        #15
)
Idents(integrated) <- 'seurat_clusters'
names(new.cluster.ids) <- levels(integrated)
integrated <- RenameIdents(integrated, new.cluster.ids)
integrated$cell.type <- Idents(integrated)

pbmc_markers <- list(
  "Naive CD4+ T"   = c("CD3G", "CCR7", "IL7R"),
  "Memory CD4+ T"  = c("S100A4"),
  "CD8+ T"         = c("CD8A"),
  "CD14+ Mono"     = c("CD14", "LYZ", "S100A8"),
  "FCGR3A+ Mono"   = c("FCGR3A", "MS4A7"),
  "B cell"         = c("MS4A1", "CD19"),
  "NK cell"        = c("GNLY", "NKG7", "NCAM1", "KLRD1", "KLRF1"),
  "Dendritic cell" = c("FCER1A", "CST3"),
  "Platelet"       = c("PPBP")
)
library(scCustomize)
Clustered_DotPlot(seurat_object = integrated, features = c(
  'CD3G', 'CCR7', 'IL7R', 'S100A4', 
  'CD4', 'CD8A', 
  'MS4A1', 'CD19',
  'CD14','S100A8', 'LYZ',
  'FCGR3A', 
  'FCER1A', 'CST3',
  'GNLY', 'NKG7', 'NCAM1', 'KLRD1', 'KLRF1',
  'PPBP'
), k = 7)

FeaturePlot_scCustom(seurat_object = integrated, features = c(
  'CD14','S100A8', 'LYZ', 'C1QA', 'C1QC',
  'FCGR3A'
))
DimPlot(integrated)



DefaultAssay(integrated) <- 'SCT'
DoHeatmap(integrated, features = c(
  'CD3G', 'CCR7', 'IL7R', 'S100A4', 
  'CD4', 'CD8A', 
  'MS4A1', 'CD19',
  'CD14','S100A8', 'LYZ',
  'FCGR3A', 
  'FCER1A', 'CST3',
  'GNLY', 'NKG7', 'NCAM1', 'KLRD1', 'KLRF1',
  'PPBP'
), group.by = 'cell.type') + NoLegend()

DimPlot(integrated, reduction = "umap", label = TRUE, group.by = "seurat_clusters",
        pt.size = 3) 
DimPlot(integrated, reduction = "umap", label = TRUE, group.by = "seurat_clusters", split.by = 'efficacy',
        pt.size = 3) 

saveRDS(integrated, file = "integrated.rds")

FeaturePlot(integrated, features = c("LAIR2"), split.by = "efficacy")

FeaturePlot(integrated, features = c("LAIR2"), split.by = "efficacy")
DefaultAssay(integrated) <- 'RNA'
VlnPlot_scCustom(seurat_object = integrated, features = "LAIR2", split.by = 'donor_name')
VlnPlot_scCustom(seurat_object = integrated, features = "LAIR2", split.by = 'efficacy')
