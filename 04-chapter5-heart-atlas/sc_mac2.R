# Mac2 fibrosis barplot and wound-healing GSEA plot.
library(dplyr)
library(ggplot2)
library(stringr)

go <- read.csv("/Users/yanjun/Desktop/PQE_thesis/scripts/chapter5/sc_go/heart_sc2_go.csv",
               row.names = 1)

# Filter for fibrosis-relevant terms
fibrosis_pattern <- paste(c("wound", "collagen", "extracellular matrix",
                            "inflammat", "tumor necrosis", "TNF",
                            "cytokine production", "angiogen", "fibroblast",
                            "remodel"), collapse = "|")

fib <- go %>%
  filter(str_detect(Description, regex(fibrosis_pattern, ignore_case = TRUE))) %>%
  # Filter out neuroinflammatory (off-topic for cardiac)
  filter(!str_detect(Description, regex("neuro", ignore_case = TRUE))) %>%
  arrange(p.adjust) %>%
  slice_head(n = 12) %>%
  mutate(
    neg_log_padj = -log10(p.adjust),
    top_genes = sapply(strsplit(geneID, "/"), function(g) {
      paste(head(g, 6), collapse = ", ")
    }),
    Description = str_wrap(Description, width = 40)
  ) %>%
  arrange(neg_log_padj) %>%
  mutate(Description = factor(Description, levels = Description))

p <- ggplot(fib, aes(x = neg_log_padj, y = Description, fill = RichFactor)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(
    aes(label = top_genes, x = 0.1),
    hjust = 0, vjust = 0.5,
    size = 2.8, color = "black", fontface = "italic"
  ) +
  scale_fill_gradient(low = "#fee5d9", high = "#a50f15", name = "RichFactor") +
  labs(
    x = expression(-log[10]~"(p.adjust)"),
    y = NULL,
    title = "Fibrosis-related GO terms in heart sc Mac2"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 9),
    plot.title = element_text(size = 11, face = "bold"),
    legend.position = "right"
  )

ggsave(
  "/Users/yanjun/Desktop/PQE_thesis/scripts/chapter5/sc_go/mac2_fibrosis_barplot.png",
  p, width = 8, height = 5, dpi = 600
)
print(p)
go <- read.csv("/Users/yanjun/Desktop/PQE_thesis/scripts/chapter5/sc_gsea/heart_sc2_gsea.csv",
               row.names = 1)

library(enrichplot)
p_enrich <- gseaplot2(
  gsea_mac2,
  geneSetID = "GO:0042060",  # wound healing
  title = "Wound healing — heart sc Mac2",
  color = "#a50f15",
  base_size = 11,
  pvalue_table = TRUE
)

ggsave(
  "/Users/yanjun/Desktop/PQE_thesis/scripts/chapter5/sc_go/mac2_wound_healing_gsea.png",
  p_enrich, width = 7, height = 5, dpi = 600
)