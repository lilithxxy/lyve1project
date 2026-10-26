# =============================================================================
# TF activity inference with decoupleR
# Continues from: dds, resLFC, vsdata (LYVE1_pos vs LYVE1_neg contrast)
# =============================================================================
library(decoupleR)
library(OmnipathR)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(pheatmap)
# -----------------------------------------------------------------------------
# 1. Get CollecTRI prior (human; use organism = "mouse" if needed)
# -----------------------------------------------------------------------------
net <- get_collectri(organism = "human", split_complexes = FALSE)
# -----------------------------------------------------------------------------
# 2. Contrast-level TF activity (LYVE1+ vs LYVE1-)
deg <- as.data.frame(resLFC) %>%
  rownames_to_column("gene") %>%
  filter(!is.na(stat)) %>%
  column_to_rownames("gene")

mat_contrast <- deg[, "log2FoldChange", drop = FALSE]
colnames(mat_contrast) <- "LYVE1pos_vs_LYVE1neg"

tf_contrast <- run_ulm(
  mat = as.matrix(mat_contrast),
  net = net,
  .source = "source",
  .target = "target",
  .mor = "mor",
  minsize = 5
)

# Keep only UP-activated TFs in LYVE1+ (score > 0), sorted by score
tf_up <- tf_contrast %>%
  filter(statistic == "ulm", score > 0) %>%
  arrange(desc(score))

top_up <- tf_up %>% slice_head(n = 25)

ggplot(top_up, aes(x = reorder(source, score), y = score, fill = -log10(p_value))) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient(low = "#fcae91", high = "#a50f15") +
  labs(x = NULL, y = "TF activity score (ulm)",
       title = "Top 25 TFs UP-activated in LYVE1+ macrophages",
       fill = "-log10(p)") +
  theme_minimal()

ggsave("tf_up_LYVE1pos.pdf", width = 6, height = 7)
write.csv(tf_up, "tf_up_LYVE1pos.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 3. Sample-level heatmap restricted to those up-activated TFs
# -----------------------------------------------------------------------------
mat_samples <- assay(vsdata)

tf_samples <- run_ulm(
  mat = mat_samples,
  net = net,
  .source = "source",
  .target = "target",
  .mor = "mor",
  minsize = 5
)

up_tf_names <- top_up$source

tf_mat <- tf_samples %>%
  filter(statistic == "ulm", source %in% up_tf_names) %>%
  pivot_wider(id_cols = source, names_from = condition, values_from = score) %>%
  column_to_rownames("source") %>%
  as.matrix()

annot <- data.frame(condition = coldata$condition,
                    row.names = rownames(coldata))

pheatmap(tf_mat,
         scale = "row",
         annotation_col = annot,
         color = colorRampPalette(c("steelblue", "white", "firebrick"))(100),
         fontsize_row = 8,
         main = "Top up-activated TFs in LYVE1+ — sample-level activity",
         filename = "tf_up_LYVE1pos_heatmap.pdf",
         width = 9, height = 8)
