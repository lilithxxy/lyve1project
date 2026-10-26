# SCTransform merge of donor PBMCs with cell-cycle scoring and marker analysis.
setwd("~/Desktop/scanpy/monocytes")
library(Seurat)
library(dplyr)
library(ggplot2)
path <- "~/Desktop/bioinfo/sc_R/unfiltered/toy"
pbmc_files <- list.files(path, pattern = ".rds$", full.names = TRUE)
pbmc_list <- lapply(pbmc_files, readRDS)
names(pbmc_list) <- basename(pbmc_files)
for (i in seq_along(pbmc_list)) {
  pbmc_list[[i]]@assays$HTO <- NULL
  pbmc_list[[i]]@assays$RNA@layers$data <- NULL
  pbmc_list[[i]]@assays$RNA@layers$scale.data <- NULL
}
integrated <- merge(pbmc_list[[1]], y = pbmc_list[2:5])
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
integrated <- SCTransform(
  integrated,
  verbose = FALSE
)

DefaultAssay(integrated) <- "SCT"

integrated <- RunPCA(integrated, npcs = 50, verbose = FALSE)
integrated <- RunUMAP(integrated, dims = 1:30, verbose = FALSE)
integrated <- FindNeighbors(integrated, dims = 1:30, verbose = FALSE)
integrated <- FindClusters(integrated, resolution = 0.4, verbose = FALSE)

donor_map <- c(
  "240106-Donor1-PBMC-L1-15" = "donor_15",
  "240120-Donor1-PBMC-L1-64" = "donor_64",
  "240229-Donor1-PBMC-L1"    = "donor_6",
  "240309-Donor2-PBMC-L1"    = "donor_90",
  "240316-Donor2-MoMac-Day-0-M" = 'day0M'
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

pbmc_markers <- list(
  "Naive CD4+ T"   = c("CCR7"),
  "CD14+ Mono"     = c("CD14", "LYZ"),
  "Memory CD4+ T"  = c("IL7R", "S100A4"),
  "B cell"         = c("MS4A1"),
  "CD8+ T"         = c("CD8A"),
  "FCGR3A+ Mono"   = c("FCGR3A", "MS4A7"),
  "NK cell"        = c("GNLY", "NKG7"),
  "Dendritic cell" = c("FCER1A", "CST3"),
  "Platelet"       = c("PPBP")
)

DotPlot(integrated, features = c('SLC1A4', 'SLC1A5'), cols = c("blue", "red")
        , col.max = 1.5, col.min = -1, group.by = 'donor_name') 
DimPlot(integrated, reduction = "umap", label = TRUE, group.by = "seurat_clusters") 
p <- DimPlot(integrated, reduction = "umap", label = TRUE, group.by = "seurat_clusters",
        split.by = "efficacy") 
DimPlot(integrated, reduction = "umap", label = TRUE, group.by = "Phase",
        pt.size = 1) 
DimPlot(integrated, reduction = "umap", label = TRUE, group.by = "percent.mt",
        pt.size = 1) 
p <- DimPlot(
  integrated,
  reduction = "umap",
  label = TRUE,
  group.by = "Phase",
  pt.size = 1
)

ggsave(
  filename = "umap.png",
  plot = p,
  width = 16,      
  height = 8,     
  dpi = 300       
)
pbmc.markers <- FindAllMarkers(integrated, only.pos = TRUE)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)


pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(integrated, features = top10$gene, group.by = 'seurat_clusters') + NoLegend()


saveRDS(integrated, file = "~/Desktop/scanpy/monos_final/monos_withm.rds")

FeaturePlot_scCustom(integrated, features = c('SLC1A4', 'SLC1A5', 'LDLR'), split.by = "donor_name")
DefaultAssay(integrated) <- 'RNA'
VlnPlot(integrated, features =  c('SLC1A4', 'SLC1A5', 'LDLR'
                                  ), group.by = "donor_name")
ggsave("ggsave1.png", plot = p, width = 40, height = 40, units = "cm", dpi = 300)
library(scCustomize)
FeaturePlot_scCustom(integrated, features = c('CD14', 'C1QA', 'C1QC', 'LYVE1',
                                              'FCGR3A', 'SLC1A4', 'SLC1A5', 'LDLR'))
                     