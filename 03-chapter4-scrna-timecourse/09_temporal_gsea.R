# GO, GSEA and KEGG loop over per-cluster single-cell DE results.
setwd("~/Desktop/PQE_thesis/scripts/chapter4/timeporal_gsea")
library(readr)
sc_deg <- read_csv("Desktop/PQE_thesis/scripts/chapter5/sc_deg.csv")
dfs <- list()
total_cols <- ncol(sc_deg)
group_size <- 5
for (i in seq(1, total_cols, by = group_size)) {
  subset_df <- sc_deg[, i:(i + group_size - 1)]
  dfs[[length(dfs) + 1]] <- subset_df
}
setwd("~/Desktop/PQE_thesis/scripts/chapter5/sc_deg")
for (i in seq_along(dfs)) {
  dfs[[i]] <- as.data.frame(dfs[[i]])
  write.csv(dfs[[i]], file = paste0('heart_sc', i-1, '.csv'))
}

library(clusterProfiler)
library(org.Hs.eg.db)

######bp and cc#######
cc_results <- list()
bp_results <- list()
gsea_results <- list()
for (i in seq_along(dfs)) {
  dfs[[i]] <- as.data.frame(dfs[[i]])
  LFC <- dfs[[i]][, 2]
  name <- dfs[[i]][, 1]
  names(LFC) <- name
  LFC <- na.omit(LFC)
  LFC <- sort(LFC, decreasing = TRUE)
  gse <- gseGO(
    geneList = LFC, 
    ont = 'BP', 
    keyType = "SYMBOL", 
    minGSSize = 3, 
    maxGSSize = 800, 
    pvalueCutoff = 0.05, 
    verbose = TRUE, 
    OrgDb = org.Hs.eg.db, 
    pAdjustMethod = "none"
  )
  gse.df <- as.data.frame(gse)
  gsea_results[[i]] <- gse
  bp_results[[i]] <- gse.df
}
####overrepresentation####
go_results <- list()
for (i in seq_along(dfs)) {
  dfs[[i]] <- as.data.frame(dfs[[i]])
  df <- dfs[[i]][1:1000, ]
  df <- df[df[, 4] < 0.05, ]
  ego <- enrichGO(gene         = df[,1 ],
                  OrgDb         = org.Hs.eg.db,
                  universe      = dfs[[1]][, 1],
                  keyType       = 'SYMBOL',
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05)
  go.df <- as.data.frame(ego)
  go_results[[i]] <- go.df
}
setwd("~/Desktop/PQE_thesis/scripts/chapter5/sc_go")
for (i in seq_along(go_results)) {
  go_results[[i]] <- as.data.frame(go_results[[i]])
  write.csv(go_results[[i]], file = paste0('heart_sc', i-1, '_go.csv'))
}
####
setwd("~/Desktop/PQE_thesis/scripts/chapter5/sc_gsea")
for (i in seq_along(bp_results)) {
  bp_results[[i]] <- as.data.frame(bp_results[[i]])
  write.csv(bp_results[[i]], file = paste0('heart_sc', i-1, '_gsea.csv'))
}
for (i in seq_along(gsea_results)) {
  saveRDS(gsea_results[[i]], file = paste0('heart_sc', i-1, '_gsea.rds'))
}
library(enrichplot)

gsea_mac2 <- readRDS("heart_sc2_gsea.rds")

# Plot wound healing + inflammatory response together
p <- gseaplot2(
  gsea_mac2,
  geneSetID = c("GO:0042060", "GO:0006954"),
  title = "Mac2 fibrosis-related GSEA",
  pvalue_table = TRUE
)

ggsave("mac2_fibrosis_gsea.png", p, width = 7, height = 6, dpi = 600)
####
## sc mac0: GO:0007015	actin filament organization; 
## sc mac1: GO:0019882 antigen processing and presentation; 
##
##
##
png(filename = 'mac5.png', width = 2000, height = 1000, res = 200)
dev.off()
enrichplot::gseaplot2(gsea_results[[1]], 
                      geneSetID = c('GO:0019882', 'GO:0007015', 'GO:0006119'), 
                      subplots = 1:3,
                      base_size = 12,
                      title = 'heart_sc_Mac0_GSEA')
enrichplot::gseaplot2(gsea_results[[2]], 
                      geneSetID = c('GO:0006119', 'GO:0002483', 'GO:0002483'), 
                      subplots = 1:3,
                      base_size = 12,
                      title = 'heart_sc_Mac1_GSEA')
enrichplot::gseaplot2(gsea_results[[3]], 
                      geneSetID = c('GO:0002526'), 
                      subplots = 1:3,
                      base_size = 12,
                      title = 'heart_sc_Mac2_GSEA')
enrichplot::gseaplot2(gsea_results[[4]], 
                      geneSetID = c('GO:0140014', 'GO:0000280'), 
                      subplots = 1:3,
                      base_size = 12,
                      title = 'heart_sc_Mac3_GSEA')
enrichplot::gseaplot2(gsea_results[[5]], 
                      geneSetID = c('GO:0098754', 'GO:0098849', 'GO:0098869'), 
                      subplots = 1:3,
                      base_size = 12,
                      title = 'heart_sc_Mac4_GSEA')

enrichplot::gseaplot2(gsea_results[[5]], 
                      geneSetID = c('GO:0097193', 'GO:0030574', 'GO:0080171'), 
                      subplots = 1:3,
                      base_size = 12,
                      title = 'heart_sc_Mac5_GSEA')



gsea_mac0 <- gsea_results[[1]]

####
go_results <- list()
for (i in seq_along(dfs)) {
  dfs[[i]] <- as.data.frame(dfs[[i]])
  LFC <- dfs[[i]][, 2]
  name <- dfs[[i]][, 1]
  names(LFC) <- name
  LFC <- na.omit(LFC)
  LFC <- sort(LFC, decreasing = TRUE)
  ego <- enrichGO(gene         = cluster_4$Mac4_names[1:1000],
                  OrgDb         = org.Hs.eg.db,
                  keyType       = 'SYMBOL',
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05)
  go.df <- as.data.frame(ego)
  go_results[[i]] <- gse.df
}
for (i in seq_along(cc_results)) {
  cc_results[[i]] <- as.data.frame(cc_results[[i]])
  write.csv(cc_results[[i]], file = paste0('cluster_', i-1, '_cc.csv'))
}


#####filtered bp and cc######

for (i in seq_along(dfs)) {
  dfs[[i]] <- dfs[[i]][!(dfs[[i]][, 1] %in% c(ribo, mito)), ]
}

cc_results <- list()
bp_results <- list()
for (i in seq_along(dfs)) {
  dfs[[i]] <- as.data.frame(dfs[[i]])
  LFC <- dfs[[i]][, 2]
  name <- dfs[[i]][, 1]
  names(LFC) <- name
  LFC <- na.omit(LFC)
  LFC <- sort(LFC, decreasing = TRUE)
  gse <- gseGO(
    geneList = LFC, 
    ont = 'BP', 
    keyType = "SYMBOL", 
    minGSSize = 3, 
    maxGSSize = 800, 
    pvalueCutoff = 0.05, 
    verbose = TRUE, 
    OrgDb = org.Hs.eg.db, 
    pAdjustMethod = "none"
  )
  gse.df <- as.data.frame(gse)
  bp_results[[i]] <- gse.df
}
for (i in seq_along(bp_results)) {
  bp_results[[i]] <- as.data.frame(bp_results[[i]])
  write.csv(bp_results[[i]], file = paste0('cluster_', i-1, '_bp_filtered.csv'))
}


for (i in seq_along(dfs)) {
  dfs[[i]] <- as.data.frame(dfs[[i]])
  LFC <- dfs[[i]][, 2]
  name <- dfs[[i]][, 1]
  names(LFC) <- name
  LFC <- na.omit(LFC)
  LFC <- sort(LFC, decreasing = TRUE)
  gse <- gseGO(
    geneList = LFC, 
    ont = 'CC', 
    keyType = "SYMBOL", 
    minGSSize = 3, 
    maxGSSize = 800, 
    pvalueCutoff = 0.05, 
    verbose = TRUE, 
    OrgDb = org.Hs.eg.db, 
    pAdjustMethod = "none"
  )
  gse.df <- as.data.frame(gse)
  cc_results[[i]] <- gse.df
}
for (i in seq_along(cc_results)) {
  cc_results[[i]] <- as.data.frame(cc_results[[i]])
  write.csv(cc_results[[i]], file = paste0('cluster_', i-1, '_cc_filtered.csv'))
}

####kegg####
kegg_results <- list()  

for (i in seq_along(dfs)) {
  dfs[[i]] <- as.data.frame(dfs[[i]])
  LFC <- dfs[[i]][, 2]
  name <- dfs[[i]][, 1]
  names(LFC) <- name
  LFC <- na.omit(LFC)
  LFC <- sort(LFC, decreasing = TRUE)
  kegg.df <- bitr(
    names(LFC), 
    fromType = 'SYMBOL', 
    toType = 'UNIPROT', 
    OrgDb = org.Hs.eg.db
  )
  
  
  mapped_LFC <- LFC[names(LFC) %in% kegg.df$SYMBOL]
  names(mapped_LFC) <- kegg.df$UNIPROT[match(names(mapped_LFC), kegg.df$SYMBOL)]
  
  gse <- gseKEGG(
    geneList = mapped_LFC,
    organism = 'hsa',
    nPerm = 10000,
    minGSSize = 3,
    maxGSSize = 800,
    pvalueCutoff = 0.05,
    pAdjustMethod = "none",
    keyType = "uniprot"
  )
  
  gse.df <- as.data.frame(gse)
  kegg_results[[i]] <- gse.df
}

for (i in seq_along(kegg_results)) {
  kegg_results[[i]] <- as.data.frame(kegg_results[[i]])
  write.csv(kegg_results[[i]], file = paste0('cluster_', i, '_kegg.csv'))
}

##### plots #####

library(dplyr)
library(tidyr)
gsea_plot <- function(data, title) {
  ggplot(data, aes(x = zScore, y = Description, fill = p.adjust)) +
    geom_bar(stat = "identity") +
    geom_label(aes(label = geneID, x = 0), 
               size = 5, hjust = 0, vjust = 0.5, fontface = "bold") +
    labs(x = "zScores", y = "Description") +
    ggtitle(title) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 15, colour = 'black', face = 'bold'),
      plot.title = element_text(face = "bold"),
      panel.background = element_rect(fill = "white", colour = NA), # White panel background
      plot.background = element_rect(fill = "white", colour = NA),  # White plot background
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank()   # Remove minor grid lines
    ) +
    scale_fill_gradient(low = "red", high = "pink")
}
library(ggplot2)
prepare_core_enrichment <- function(data) {
  data %>%
    group_by(Description) %>%
    slice_max(order_by = FoldEnrichment, n = 10) %>%  
    mutate(geneID = sapply(strsplit(geneID, "/"), 
                           function(x) paste(head(x, 5), collapse = ", "))) %>%
    ungroup()
}
# cluster0 actin filaments; endocytosis; cell-matrix
cluster_0_go <- read_csv("cluster_0_go.csv")
filtered_0 <- cluster_0_go %>%
  dplyr::filter(stringr::str_detect(pattern = "actin|endocytosis|vesicle", Description)) %>%
  group_by(Description) %>%
  add_count() %>%
  dplyr::arrange(dplyr::desc(n), dplyr::desc(Description)) %>%
  mutate(Description= forcats::fct_inorder(Description))
filtered_0 <- prepare_core_enrichment(filtered_0)
gsea_plot(filtered_0, "cluster0")

### cluster1 cardiac; muscle; MHC; endocytosis down
cluster_1_go <- read_csv("cluster_1_go.csv")
filtered_1 <- cluster_1_go %>%
  dplyr::filter(stringr::str_detect(pattern = "cardi|muscle|MHC", Description)) %>%
  group_by(Description) %>%
  add_count() %>%
  dplyr::arrange(dplyr::desc(n), dplyr::desc(Description)) %>%
  mutate(Description= forcats::fct_inorder(Description))
filtered_1 <- prepare_core_enrichment(filtered_1)
gsea_plot(filtered_1, "cluster1")

### cluster2 cardiac; muscle; endocytosis up
cluster_2_go <- read_csv("cluster_2_go.csv")
filtered_2 <- cluster_2_go %>%
  dplyr::filter(stringr::str_detect(pattern = "cardi|muscle|endocytosis", Description)) %>%
  group_by(Description) %>%
  add_count() %>%
  dplyr::arrange(dplyr::desc(n), dplyr::desc(Description)) %>%
  mutate(Description= forcats::fct_inorder(Description))
filtered_2 <- prepare_core_enrichment(filtered_2)
gsea_plot(filtered_2, "cluster2")
##cluster3 recruit t cell ; endocytosis
cluster_3_go <- read_csv("cluster_3_go.csv")
filtered_3 <- cluster_3_go %>%
  dplyr::filter(stringr::str_detect(pattern = "T cell|endocytosis", Description)) %>%
  group_by(Description) %>%
  add_count() %>%
  dplyr::arrange(dplyr::desc(n), dplyr::desc(Description)) %>%
  mutate(Description= forcats::fct_inorder(Description))
filtered_3 <- prepare_core_enrichment(filtered_3)
gsea_plot(filtered_3, "cluster3")
### cluster4 exocytosis; secretory; phagocytosis; endocytosis
cluster_4_go <- read_csv("cluster_4_go.csv")
filtered_4 <- cluster_4_go %>%
  dplyr::filter(stringr::str_detect(pattern = "exocytosis|secretory|phagocytosis|endocytosis", Description)) %>%
  group_by(Description) %>%
  add_count() %>%
  dplyr::arrange(dplyr::desc(n), dplyr::desc(Description)) %>%
  mutate(Description= forcats::fct_inorder(Description))
filtered_4 <- prepare_core_enrichment(filtered_4)
gsea_plot(filtered_4, "cluster4")
