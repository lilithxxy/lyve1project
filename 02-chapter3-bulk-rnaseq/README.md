# 02-chapter3-bulk-rnaseq — thesis Chapter 3

Bulk RNA-seq of FACS-sorted/polarised macrophage populations
(classical monocytes; in-vitro LYVE-1⁺, LYVE-1⁻ and M1-like).

## Scripts (run order)

### 1. `01_deseq2_analysis.R` (source: `scripts/chapter2/FinalVersionForBulk_thesis.R`)

Full exploratory DESeq2 analysis. Steps, in script order: read count matrix →
drop three `ubc_LYVE-1+` columns → filter rows (rowSums > 20) → build
condition factor (`Classical_monocytes`, `in_vitro_LYVE1_pos`, `M1_like`,
`in_vitro_LYVE1_neg`) → `DESeqDataSetFromMatrix(design = ~condition)` → DESeq
→ vst → results for
`condition_in_vitro_LYVE1_pos_vs_in_vitro_LYVE1_neg` with apeglm LFC shrinkage
→ PCA, MA plot, hierarchical clustering (dendextend/factoextra), ComplexHeatmap
distance/gene heatmaps (rlog, z-score), scree plot, EnhancedVolcano →
gseGO/enrichGO/KEGG (enrichKEGG + gseKEGG), chord diagrams, gene-list checks.

- **Working directory:** the script begins with
  `setwd("~/Desktop/scanpy/bulk_for_paper/everything")` — change to your data
  location.
- **Inputs (expected in that directory):**
  - `tablefordeseqfinal.csv` — raw count matrix, genes × samples. Sample
    columns: `Classical_monocytes_1..7`, `ubc_LYVE-1+_1..3`,
    `in_vitro_LYVE-1+_1..4`, `in_vitro_M1-like_1..4`, `in_vitro_LYVE-1-_1..4`.
  - `batch12withgenelength.csv` — per-gene annotation with gene length.
  - `tpmtable.csv` — TPM matrix (same samples).
- **Outputs:** plots to screen/plot window (no `ggsave`); one fgsea chunk
  contains a placeholder gene set (`...`) and is not executable as-is.
- **Libraries:** tibble, DESeq2, dplyr, ggplot2, tidyr, readr, tidyverse,
  dendextend, factoextra, apeglm, ComplexHeatmap, circlize, pheatmap,
  EnhancedVolcano, clusterProfiler, fgsea, org.Hs.eg.db, enrichplot, DOSE,
  reshape2, R.utils.

### 2. `02_decoupler_tf_activity.R` (source: `scripts/chapter4/decoupleR_for_bulk.R`)

Transcription-factor activity (decoupleR ulm on the CollecTRI network from
OmnipathR) for the LYVE1pos-vs-LYVE1neg contrast.

- **Note:** does **not** re-read files; it continues from the objects
  `dds`, `resLFC`, `vsdata` created by `01_deseq2_analysis.R` — run both in the
  same R session.
- **Inputs:** CollecTRI network via `get_collectri(organism="human")`.
- **Outputs:** `tf_up_LYVE1pos.pdf`, `tf_up_LYVE1pos.csv`,
  `tf_up_LYVE1pos_heatmap.pdf`.
- **Libraries:** decoupleR, OmnipathR, dplyr, tibble, tidyr, ggplot2, pheatmap.

### 3. `03_momacverse_projection.R` (source: `scanpy/8_figures/score_genes_momacverse/momacverse.R`)

Reverse projection of the single-cell/bulk signatures onto MoMacVerse (Mulder
et al. 2021) — thesis Fig 4.4.13 / Methods 9.1. Loads `MoMac_VERSE.RDS`, reads
five `top*.csv` gene lists (bulk top-50, HLA set, Eraslan set, Mac4 top-50, Mac7
top-50), scores each with Seurat `AddModuleScore(ctrl = 5)`, plots DimPlots /
FeaturePlots per signature, per-cluster signature-average heatmaps
(ComplexHeatmap / pheatmap), repeats scoring per `Tissue` after `SplitObject`
with per-tissue heatmaps, and ends with a pseudobulk DE heatmap.

- **Inputs:** `MoMac_VERSE.RDS`, `top*.csv` gene lists, `deg.csv` (all in the
  script's `setwd` folder).
- **Outputs:** `momacverse_dimplot.png`, `mac4top50features_dimplot.png`,
  `*features1.png`, `heatmap_momacverse.png`, per-tissue `*_<Tissue>.png` and
  `heatmap_<Tissue>_*.png`.
- **Libraries:** Seurat, scCustomize, ggplot2, dplyr, pheatmap, ComplexHeatmap,
  colorRamp2.

### 4. `04_music_deconvolution.R` (source: `scanpy/decon/music.R`)

MuSiC deconvolution (Methods 9.3): converts the Day-5 MDIX single-cell
reference `day5mdix.h5ad` to a SingleCellExperiment (`zellkonverter`), reads
the bulk matrix `mdix.csv`, intersects genes, assigns pseudo-sample IDs, runs
`music_prop(bulk.mtx, sce, clusters = 'leiden_res_0_8', samples = 'sampleID')`
and plots `Est.prop.weighted` and `Var.prop` as pheatmaps.

- **Inputs:** `decon/day5mdix.h5ad`, `decon/mdix.csv` (bulk counts, genes × samples).
- **Outputs:** two heatmaps (plot window).

### 5. `05_celldistinguisher_deconvolution.R` (source: `scanpy/decon/CellDistinguisher.R`)

Marker-free deconvolution of the bulk table with CellDistinguisher
(`gecd_CellDistinguisher`, `numCellClasses = 4`, then
`gecd_DeconvolutionByDistinguishers`) → sample compositions and cell-class
signatures. Not named in the thesis Methods.

- **Inputs:** `tablefordeseqfinal` object from `01_deseq2_analysis.R` (columns 6, 7, 8, 23).

### `helpers/plotpca.R` (source: `scanpy/bulk_for_paper/paper/codes/plotpca.R`)

Customised PCA for the bulk samples: `getPCs()` (top-500 genes by variance,
`prcomp(t(assay(vsdata)))`, percent-variance labels) and a ggplot with
`geom_mark_ellipse` per group and a fixed five-colour palette. Expects the
`vsdata` object from `01_deseq2_analysis.R` in the session.

## Note on script versions

Several near-identical earlier versions of the bulk analysis exist at source —
an interactive R console history and an older, larger copy of the DESeq2
script (`scanpy/bulk_for_paper/everything/FinalVersionForBulk_thesis.R`,
2026-05-25) — none copied. `01_deseq2_analysis.R` is the single final script.
