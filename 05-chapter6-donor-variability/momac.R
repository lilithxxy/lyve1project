# Merge donor PBMC objects, QC, cluster and test high- vs low-yield donors.
library(Seurat)
path <- "~/Desktop/bioinfo/sc_R/unfiltered/"
pbmc_files <- list.files(path, pattern = "PBMC.*\\.rds$", full.names = TRUE)
pbmc_list <- lapply(pbmc_files, readRDS)
names(pbmc_list) <- basename(pbmc_files)
for (i in seq_along(pbmc_list)) {
  pbmc_list[[i]]@assays$HTO <- NULL
  pbmc_list[[i]]@assays$RNA@layers$data <- NULL
  pbmc_list[[i]]@assays$RNA@layers$scale.data <- NULL
}
integrated <- merge(pbmc_list[[1]], y = pbmc_list[2:4])
integrated
remove(pbmc_list)
integrated[["percent.mt"]] <- PercentageFeatureSet(integrated, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
##
plot1 <- FeatureScatter(integrated, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(integrated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#33538 features across 15123 samples within 1 assay 
integrated <- subset(integrated, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & nCount_RNA < 80000 & percent.mt < 5)
#33538 features across 9315 samples within 1 assay 
table(integrated$hash.ID)
###
integrated <- NormalizeData(integrated, normalization.method = "LogNormalize", scale.factor = 10000)
integrated <- FindVariableFeatures(integrated, selection.method = "vst", nfeatures = 3000)
integrated <- ScaleData(integrated, features = rownames(integrated))
integrated <- RunPCA(integrated, features = VariableFeatures(object = integrated))
integrated <- RunUMAP(integrated, dims = 1:30)
integrated <- FindNeighbors(integrated, dims = 1:30)
integrated <- FindClusters(integrated, resolution = 0.5)

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
integrated <- JoinLayers(integrated)
#saveRDS(integrated, file = "~/Desktop/bioinfo/sc_R/unfiltered/integrated.rds")
integrated <- NormalizeData(integrated)
# Find DE features between CD16 Mono and CD1 Mono
Idents(integrated) <- "efficacy"
de.markers <- FindMarkers(integrated, ident.1 = "high", ident.2 = "low")
# view results
head(de.markers, n = 10)
##

pseudo_integrated <- AggregateExpression(integrated, assays = "RNA", 
                                         return.seurat = F, 
                                         group.by = c("efficacy", "donor_name"))
head(pseudo_integrated$RNA)

library(edgeR)

# 1. Prepare count matrix
counts <- as.matrix(pseudo_integrated$RNA)

# 2. Sample group labels
group <- factor(c("high", "high", "low", "low"), levels = c("low", "high"))
design <- model.matrix(~ group)
samples <- colnames(counts)
sample_info <- data.frame(row.names = samples, group = group)

# 3. Create DGEList
dge <- DGEList(counts = counts, group = sample_info$group)

# 4. Filter and normalize
keep <- filterByExpr(dge)
dge <- dge[keep,, keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

# 5. Estimate dispersion and fit model
dge <- estimateDisp(dge, design, robust = TRUE)
fit <- glmFit(dge, design)
colnames(fit$design)
# 6. Perform differential expression test
lrt <- glmLRT(fit, coef = 2)

de_results <- topTags(lrt, n = Inf)

head(de_results$table, n = 20)
df <- de_results$table
##
de.markers$dif <- abs(de.markers$pct.1 - de.markers$pct.2)
de.markers <- de.markers[order(-de.markers$dif), ]
##edgeR
features <- c(
  "CXCL9", "CADM2", "ARHGEF10", "LAIR2", "RNF17", "TDRG1",
  "HRASLS2", "IGLC3", "PLEKHS1", "SYCP2L", "CNTNAP4",
  "IGHG2", "TRDV2", "IGHG4", "P2RX6", "LRFN2",
  "HLA-DRB5", "FAM222A-AS1", "COL13A1", "FLNC-AS1",
  "ANGPTL6", "TBX15", "MYOM2", "IFI27", "ZMAT4",
  "FBXO39", "G0S2", "MTTP", "TRPC3", "HFM1"
)
VlnPlot(integrated, features = features, group.by = 'donor_name')
##wilcox
features <- c(
  "MKX", "ARHGEF10", "RNF17", "LAIR2", "CADM2", "HRASLS2", "TDRG1",
  "SYCP2L", "IGHG2", "TRDV2", "IGHG4", "LRFN2", "TRHDE", "COL13A1",
  "CNTNAP4", "MYOM2", "RGS6", "ANGPTL6", "ZMAT4",
  "FAM222A-AS1", "TRPC3", "HLA-DRB5", "IGLC3", "PHETA2", "MAFIP",
  "PARD6G", "ANXA10", "CCDC192", "G0S2", "OLFM1", "LIPC", "NINL", "CCL3L1"
)
VlnPlot(integrated, features = features, group.by = 'donor_name')

VlnPlot(integrated, features = c('LAIR2', 'LAIR1' ), group.by = 'donor_name')
library(EnhancedVolcano)
