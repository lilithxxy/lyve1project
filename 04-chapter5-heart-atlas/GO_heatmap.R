# Curated GO dotplot and heatmap across heart-atlas macrophage clusters.
library(ggplot2)
library(stringr)
library(dplyr)

# Load curated GO terms
sc_go_curated <- read_csv("~/Desktop/PQE_thesis/scripts/chapter5/sc_go_curated.csv")
library(ggplot2)
print(GO_curated)
unique(sc_go_curated$Description)
#####dotplot######
library(RColorBrewer)
filtered_go <- sc_go_curated %>%
  filter(RichFactor > 0, p.adjust < 1) %>%
  mutate(RichFactor = as.numeric(RichFactor) ,
   p.adjust <- as.numeric(p.adjust))
#
custom_order <- c(
     "antigen processing and presentation",
     "oxidative phosphorylation",
     "regulation of endocytosis" , 
     "receptor-mediated endocytosis",          
     "clathrin-dependent endocytosis"   ,  
     "regulation of inflammatory response",  
     "nuclear division" ,     
     "DNA replication",                     
     "detoxification"  , 
     "cellular detoxification",               
     "tumor necrosis factor production" ,
     "regulation of apoptotic signaling pathway"
   )
# Apply manual order (reverse to have top-to-bottom)
filtered_go$Description <- factor(filtered_go$Description, levels = rev(custom_order))

# Plot
ggplot(filtered_go, aes(x = Clusters, y = Description)) +
  geom_point(aes(size = RichFactor, color = p.adjust)) +
  scale_size(range = c(2, 6))+
  scale_color_gradientn(
    colors = c("red", "blue"),
    limits = c(0, 0.008),
    oob = scales::squish
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    legend.title = element_text(size = 11)
  ) +
  labs(
    size = "Rich Factor",
    color = "Adjusted p-value",
    x = "Clusters",
    y = "Biological Pathways",
    title = "heart sc-atlas\nGO Dotplot Across Clusters"
  )
ggsave("heart_go_dotplot.png", width = 10, height = 8, units = "in", dpi = 600)
#####heatmap#####
sc_go_curated$Description <- factor(sc_go_curated$Description, levels = rev(custom_order))
library(scales)  
sc_go_curated$log_p <- pmin(-log10(sc_go_curated$p.adjust), 4)  

ggplot(sc_go_curated, aes(x = Clusters, y = Description, fill = log_p)) +
  geom_tile(color = "black", size = 0.1) +  
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 2,  
                       name = "-log10(p.adjust)", 
                       limits = c(0, 4),  # Set maximum limit to 4
                       oob = squish) +  # Ensures out-of-bounds values are squished into range
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(face = "bold", size = 14)
  ) +
  labs(x = "Clusters", y = "GO Description", title = "GO Enrichment Heatmap")
ggplot(GO_curated, aes(x = Clusters, y = Description, fill = -log10(p.adjust))) +
  geom_tile(color = "black", size = 0.1) +  # Adds black border for visibility
  scale_fill_gradient2(low = "skyblue", mid = "white", high = "red", 
                       midpoint = 2,  # Manually set midpoint (~ p.adjust = 0.05)
                       name = "-log10(p.adjust)",
                       limits = c(0, 8),  
                       oob = squish) +  # Squishes out-of-bound values
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(face = "bold", size = 14)
  ) +
  labs(x = "Clusters", y = "GO Description", title = "GO Enrichment Heatmap")
