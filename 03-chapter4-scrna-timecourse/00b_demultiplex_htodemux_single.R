# Compact single-run HTODemux demultiplexing and merge to final.h5ad.
library(Seurat)
# Demultiplex ----
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

filtered <- subset(seu.obj, subset = HTO_classification != 'Negative' & HTO_classification.global == 'Singlet')
obj.list <- SplitObject(filtered, split.by = "HTO_classification")
# Merge per-day objects ----
# saveRDS here
for (file in rds_files) {
  seurat_obj <- readRDS(file)  
  timepoint <- gsub(".*day(\\d+).*", "\\1", basename(file))
  seurat_obj$timepoint <- paste0("day", timepoint)
  obj.list[[file]] <- seurat_obj 
}
merged <- merge(obj.list[[1]], obj.list[2:11])
DefaultAssay(merged) <- 'RNA'
merged[["HTO"]] = NULL
library(sceasy)
sce <- as.Seurat(merged)
# Convert to h5ad ----
sceasy::convertFormat(final, from="sce", to="anndata",
                      outFile='final.h5ad')