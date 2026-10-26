# Detox and fibrosis GO and GSEA barplots for Mac4 and Mac2.
library(dplyr)
library(ggplot2)
library(stringr)

# Read GO results
go <- read.csv("/Users/yanjun/Desktop/PQE_thesis/scripts/chapter5/sc_go/heart_sc4_go.csv",
               row.names = 1)

# Filter for endocytosis terms
endo <- go %>%
  filter(str_detect(Description, regex("detoxi", ignore_case = TRUE))) %>%
  mutate(
    neg_log_padj = -log10(p.adjust),
    top_genes = sapply(strsplit(geneID, "/"), function(g) {
      paste(head(g, 6), collapse = ", ")
    }),
    # Wrap description for cleaner labels
    Description = str_wrap(Description, width = 40)
  ) %>%
  arrange(neg_log_padj) %>%
  mutate(Description = factor(Description, levels = Description))

# Plot
p <- ggplot(endo, aes(x = neg_log_padj, y = Description, fill = RichFactor)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(
    aes(label = top_genes, x = 0.1),
    hjust = 0, vjust = 0.5,
    size = 2.8, color = "black", fontface = "italic"
  ) +
  scale_fill_gradient(low = "white", high = "red", name = "RichFactor") +
  labs(
    x = expression(-log[10]~"(p.adjust)"),
    y = NULL,
    title = "Up-regulated Detox GO Terms in heart_sc_Mac4"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )

ggsave(
  "/Users/yanjun/Desktop/PQE_thesis/Chapter5_figures/detox_barplot.png",
  p, width = 10, height = 4, dpi = 300
)
print(p)