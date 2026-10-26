# MuSiC deconvolution of the bulk profile using the Day-5 single-cell reference.
library(MuSiC)
library(SingleCellExperiment)

sce <- zellkonverter::readH5AD("/Users/yanjun/Desktop/scanpy/decon/day5mdix.h5ad", verbose = TRUE)

# Load bulk matrix ----
bulk_data <- read.csv("/Users/yanjun/Desktop/scanpy/decon/mdix.csv", row.names=1)

# Ensure genes match between the reference and bulk data
genes <- intersect(rownames(sce), rownames(bulk_data))
sce <- sce[genes, ]
bulk_data <- bulk_data[genes, ]

# Load necessary package
library(SingleCellExperiment)

num_samples <- ncol(sce)
set.seed(123)  
sample_ids <- sample(c("Sample1", "Sample2"), num_samples, replace = TRUE)
colData(sce)$sampleID <- sample_ids
head(colData(sce))

bulk_data <- as.matrix(bulk_data)

# Deconvolve with MuSiC ----
# Produce the first step information
musicRes = music_prop(bulk.mtx = bulk_data, sce, clusters = 'leiden_res_0_8', samples = 'sampleID',
                      verbose = F)

est_prop_weighted <- musicRes$Est.prop.weighted
var_prop <- musicRes$Var.prop
library(pheatmap)
# Plot estimated proportions ----
# Heatmap for Est.prop.weighted
pheatmap(
  est_prop_weighted,
  main = "Heatmap of Estimated Proportion Weighted",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  color = colorRampPalette(c("white", "blue"))(50)
)

# Heatmap for Var.prop
pheatmap(
  var_prop,
  main = "Heatmap of Variance Proportion",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  color = colorRampPalette(c("white", "red"))(50)
)
