# Convert per-donor Seurat objects to h5ad for the Python side.
library(Seurat)

path <- "~/Desktop/bioinfo/sc_R/unfiltered/"
pbmc_files <- list.files(path, pattern = "PBMC.*\\.rds$", full.names = TRUE)
pbmc_list <- lapply(pbmc_files, readRDS)
names(pbmc_list) <- basename(pbmc_files)
for (i in seq_along(pbmc_list)) {
  pbmc_list[[i]]@assays$HTO <- NULL
}
library(sceasy)
library(reticulate)
library(SeuratDisk)
library(SingleCellExperiment)

library(sceasy)
library(SingleCellExperiment)

output_dir <- getwd()  # current working directory

for (i in seq_along(pbmc_list)) {
  seurat_obj <- pbmc_list[[i]]
  sample_name <- sub("\\.rds$", "", names(pbmc_list)[i])
  sce_obj <- as.SingleCellExperiment(seurat_obj)
  
  out_file <- file.path(output_dir, paste0(sample_name, ".h5ad"))
  message("Saving to: ", out_file)
  
  sceasy::convertFormat(sce_obj, from = "sce", to = "anndata", outFile = out_file)
}