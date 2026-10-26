# Project bulk and single-cell signature scores onto the MoMacVerse reference with Seurat AddModuleScore.
setwd('/Users/yanjun/Desktop/scanpy/8_figures/score_genes_momacverse') 
library(Seurat) 
library(scCustomize)
library(ggplot2)
path <- '/Users/yanjun/Desktop/scanpy/8_figures/score_genes_momacverse'  
momac <- readRDS(file.path(path, 'MoMac_VERSE.RDS')) 

DefaultAssay(momac) <- 'RNA'
length(rownames(momac))

files_path <- list.files(path = path, pattern = "top.*\\.csv$", full.names = F)
files <- list()
for (i in seq_along(files_path)) {
  files[[i]] <- read.csv(files_path[i], col.names = F)
  colnames(files[[i]] ) <- 'genes'
}
names(files) <- basename(files_path)

hla <- c('HLA-A', 'HLA-DRA', 'CD68', 'PSAP')

momac <- AddModuleScore(object = momac, features = list(files[[1]]$genes), name = 'bulktop50features', 
               ctrl = 5, replace = TRUE)

momac <- AddModuleScore(object = momac, features = list(files[[2]]$genes), name = 'hlafeatures', 
                        ctrl = 5, replace = TRUE)

momac <- AddModuleScore(object = momac, features = list(files[[3]]$genes), name = 'eraslanfeatures', 
                        ctrl = 5, replace = TRUE)

momac <- AddModuleScore(object = momac, features = list(files[[4]]$genes), name = 'mac4top50features', 
                        ctrl = 5, replace = TRUE)

momac <- AddModuleScore(object = momac, features = list(files[[5]]$genes[1:50]), name = 'mac7top50features', 
                        ctrl = 5, replace = TRUE)

png(filename = 'momacverse_dimplot.png', width = 3000, height = 3600, res = 600)
p1 <- DimPlot_scCustom(momac, group.by = 'Clusters', pt.size = 0.5, label = TRUE, label.size = 4) +
  theme(legend.position = "none")
print(p1)
dev.off()

png(filename = 'mac4top50features_dimplot.png', width = 3000, height = 3600, res = 600)
f1 <- FeaturePlot_scCustom(momac, features = 'mac4top50features1',
                           pt.size = 0.5, label = FALSE,
                           colors_use = c('Blue', 'White', 'Red'))+
  theme(legend.position = "none")
print(f1)
dev.off()

png(filename = 'eraslantop50features1.png', width = 1600, height = 1800, res = 300)
f1 <- FeaturePlot_scCustom(momac, features = 'eraslantop50features1', 
                           colors_use = c('Blue', 'White', 'Red'),pt.size = 2)
print(f1)
dev.off()

png(filename = 'mac4top50features1.png', width = 1600, height = 1800, res = 300)
f1 <- FeaturePlot_scCustom(momac, features = 'mac4top50features1', 
                           colors_use = c('Blue', 'White', 'Red'),pt.size = 2)
print(f1)
dev.off()
#####heatmap####
feature_mac4top50 <- "mac4top50features1"
feature_eraslantop50 <- "eraslantop50features1"
feature_bulktop50 <- "bulktop50features1"
cluster_averages <- momac@meta.data %>%
  group_by(Clusters) %>%
  summarize(
    avg_mac4top50 = mean(get(feature_mac4top50), na.rm = TRUE),
    avg_eraslantop50 = mean(get(feature_eraslantop50), na.rm = TRUE),
    avg_bulktop50 = mean(get(feature_bulktop50), na.rm = TRUE)
  )


data_to_plot <- cluster_averages[, 2:4] 
rownames(data_to_plot) <- make.names(cluster_averages$Clusters, unique = TRUE)  
heatmap_filename <- "heatmap_momac_mac4top100_eraslantop50_bulktop100.png"
library(pheatmap)
library(ComplexHeatmap)
library(colorRamp2)
data_to_plot <- as.matrix(data_to_plot)
scaled_data <- scale(data_to_plot, center = TRUE, scale = TRUE)  
png(filename = 'heatmap_momacverse.png', width = 1600, height = 1800, res = 300)
p1 <- Heatmap(scaled_data, 
              row_title = NULL, 
              column_title = NULL,
              width = ncol(scaled_data) * unit(0.5, "cm"), 
              height = nrow(scaled_data) * unit(0.5, "cm"),
              col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
              heatmap_legend_param = list(
                title = "Scaled Expression Level",        
                title_gp = gpar(fontsize = 12)
              )
)
print(p1)  
dev.off() 

unique(momacverse$Tissue)
list <- SplitObject(momacverse, split.by = 'Tissue')

for (i in names(list)) {
  list[[i]] <- AddModuleScore(object = list[[i]], features = list(files[[3]]$genes), 
                              name = paste0('mac4top100in_', i), ctrl = 5, replace = TRUE)
  
  list[[i]] <- AddModuleScore(object = list[[i]], features = list(files[[2]]$genes), 
                              name = paste0('eraslantop50in_', i), ctrl = 5, replace = TRUE)
  
  list[[i]] <- AddModuleScore(object = list[[i]], features = list(files[[1]]$genes), 
                              name = paste0('bulktop100in_', i), ctrl = 5, replace = TRUE)
}
library(ggplot2)
for (i in seq_along(list)) {
  feature_mac4top100 <- paste0("mac4top100in_", names(list)[i], "1")
  feature_eraslantop50 <- paste0("eraslantop50in_", names(list)[i], "1")
  feature_bulktop100 <- paste0("bulktop100in_", names(list)[i], "1")
  plot_mac4top100 <- FeaturePlot_scCustom(list[[i]], features = feature_mac4top100)
  plot_eraslantop50 <- FeaturePlot_scCustom(list[[i]], features = feature_eraslantop50)
  plot_bulktop100 <- FeaturePlot_scCustom(list[[i]], features = feature_bulktop100)
  plot1 <- DimPlot_scCustom(list[[i]], group.by = 'Clusters')
  ggsave(paste0("mac4top100features_", names(list)[i], ".png"), plot = plot_mac4top100, dpi = 300)
  ggsave(paste0("eraslantop50features_", names(list)[i], ".png"), plot = plot_eraslantop50, dpi = 300)
  ggsave(paste0("bulktop100features_", names(list)[i], ".png"), plot = plot_bulktop100, dpi = 300)
  ggsave(paste0("DimPlot_", names(list)[i], ".png"), plot = plot1, width = 6, height = 3, dpi = 300)
}

#### heatmap: x is tissues; y is signatures of cluster 1 to 7; each block is a signature####
library(dplyr)
cluster_averages_list <- list()

for (i in seq_along(list)) {
  feature_mac4top100 <- paste0("mac4top100in_", names(list)[i], "1")
  feature_eraslantop50 <- paste0("eraslantop50in_", names(list)[i], "1")
  feature_bulktop100 <- paste0("bulktop100in_", names(list)[i], "1")
  cluster_averages <- list[[i]]@meta.data %>%
    group_by(Clusters) %>%
    summarize(
      avg_mac4top100 = mean(get(feature_mac4top100), na.rm = TRUE),
      avg_eraslantop50 = mean(get(feature_eraslantop50), na.rm = TRUE),
      avg_bulktop100 = mean(get(feature_bulktop100), na.rm = TRUE)
    )
  cluster_averages_list[[names(list)[i]]] <- cluster_averages
}
library(pheatmap)
for (i in names(cluster_averages_list)) {
  data_to_plot <- cluster_averages_list[[i]][, 2:4]
  rownames(data_to_plot) <- make.names(cluster_averages_list[[i]]$Clusters, unique = TRUE)
  heatmap_filename <- paste0("heatmap_", i, "_mac4top100_eraslantop50_bulktop100.png")
  pheatmap(data_to_plot, 
           scale = "column",          
           main = paste0("Scaled Signature Scores in ", i),  
           clustering_distance_rows = "euclidean", 
           clustering_method = "complete",
           show_rownames = TRUE,   
           show_colnames = TRUE,   
           color = colorRampPalette(c("blue", "white", "red"))(100), 
           filename = heatmap_filename
  )
}






















######heatmap#####
####kevin scaled by row####
DefaultAssay(momac) <- 'RNA'
Idents(momac) <- 'RNA'
pbulk <- AverageExpression(object = momac, slot = 'data',
                  group.by = 'Tissue',
                  features = rownames(momac))
pbulk <- as.matrix(pbulk$RNA)
deg <- read.csv(file.path(path, 'deg.csv'), header = TRUE)
dfs <- list()
total_cols <- ncol(deg)
group_size <- 5
for (i in seq(1, total_cols, by = group_size)) {
  subset_df <- deg[, i:(i + group_size - 1)]
  dfs[[length(dfs) + 1]] <- subset_df
}

top100 <- unlist(lapply(dfs, function(df) head(df[[1]], 10)))
pbulk <- pbulk[rownames(pbulk) %in% top100, ]
pbulk <- pbulk[match(top100, rownames(pbulk)), , drop = FALSE]
pheatmap::pheatmap(pbulk, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = TRUE,
                   show_colnames = TRUE,
                   scale = 'row')
pheatmap::pheatmap(pbulk, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE,
                   show_colnames = TRUE,
                   scale = 'column')
#####










hla <- c("NKG7", "CNTN4", "DHCR24", "APOE", "NQO1", "TRAF3IP3", "KCNN4", "JOSD2",
  "KCNJ15", "FTSJ2", "EGLN3", "GRIP1", "AC108938.5", "ORMDL2", "ACTN1", "SGTA",
  "MFNG", "SMARCD2", "FAT1", "NFKBIE", "DPAGT1", "CHTF18", "NDUFA1", "PLPP7",
  "TMEM119", "DNAH17", "KIF24", "TRABD2B", "CDA", "MORN1", "TMEM156", "MTRF1",
  "ZNF668", "TMEM14C", "APH1A", "HAMP", "BTF3", "PPM1K", "PRELID1", "C1orf131",
  "CD300LF", "AC010287.1", "CYC1", "MTRR", "SUPT4H1", "ATP5G2", "IDI1", "XPO5",
  "ZNHIT1", "TMEM203")











