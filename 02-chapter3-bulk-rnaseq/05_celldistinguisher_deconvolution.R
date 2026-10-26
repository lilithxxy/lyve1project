# CellDistinguisher marker-free deconvolution of the bulk count table.
library(CellDistinguisher)
df <- tablefordeseqfinal[, c(6,7,8,23)]
rownames(df) <- make.names(tablefordeseqfinal$...1, unique = T)
df <- as.matrix(df)
# Find distinguishers ----
distinguishers <- gecd_CellDistinguisher(
  df,
  genesymb = NULL,  # Use gene symbols if available, or set to NULL
  numCellClasses = 4,  # Adjust based on your estimated number of cell types
  minDistinguisherAlternatives = 1,
  maxDistinguisherAlternatives = 100,
  minAlternativesLengthsNormalized = 0.5,
  expressionQuantileForFilter = 0.999,
  expressionConcentrationRatio = 0.333,
  verbose = 0
)
# Deconvolve by distinguishers ----
deconvolution <- gecd_DeconvolutionByDistinguishers(
  df,
  distinguishers$bestDistinguishers,
  nonNegativeOnly = TRUE,
  convexSolution = TRUE,
  verbose = 0
)
sampleCompositions <- t(deconvolution$sampleCompositions)

# Cell-type-specific expression profiles
celltypeSignatures <- deconvolution$cellSubclassSignatures
