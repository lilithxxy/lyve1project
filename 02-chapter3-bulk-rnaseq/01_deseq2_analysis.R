# Bulk RNA-seq differential expression: DESeq2 on sorted macrophage populations with PCA, MA, volcano, heatmaps and GSEA.
library(tibble) 
library(DESeq2)
library(dplyr)
library(ggplot2)
library(tidyr)
setwd("~/Desktop/scanpy/bulk_for_paper/everything")
library(readr)
tablefordeseqfinal <- read_csv("tablefordeseqfinal.csv")
dds <- tablefordeseqfinal[,-1]


colnames(dds) <- c('Classical_monocytes_1', 
                   'ubc_LYVE-1+_1', 'ubc_LYVE-1+_2', 'ubc_LYVE-1+_3', 
                   'in_vitro_LYVE-1+_1', 'in_vitro_LYVE-1+_2', 'in_vitro_LYVE-1+_3',
                   'in_vitro_M1-like_1', 'in_vitro_M1-like_2', 'in_vitro_M1-like_3', 
                   'Classical_monocytes_2', 'Classical_monocytes_3', 
                   'in_vitro_LYVE-1-_1', 'in_vitro_LYVE-1-_2', 'in_vitro_LYVE-1-_3', 'in_vitro_LYVE-1-_4',
                   'in_vitro_M1-like_4', 'Classical_monocytes_4', 'Classical_monocytes_5', 'Classical_monocytes_6', 
                   'Classical_monocytes_7', 'in_vitro_LYVE-1+_4')
dds <- dds[, -c(2,3,4)]
row.names(dds) <- make.names(tablefordeseqfinal$...1, unique = T) 
backup <- dds            
dds <- as.matrix(dds)
dds <- dds[which(rowSums(dds) > 20), ]
#condition <- factor(c('Classical_monocytes', rep('ubc_LYVE1_pos', 3), rep('in_vitro_LYVE1_pos', 3),  rep('in_vitro_m1_like', 3), 'Classical_monocytes', 'Classical_monocytes', rep('in_vitro_LYVE1_neg', 4), 'in_vitro_m1_like', rep('Classical_monocytes', 4), 'in_vitro_LYVE1_pos'))
condition <- factor(c('Classical_monocytes', 
                      rep('in_vitro_LYVE1_pos', 3),  rep('M1_like', 3),
                      'Classical_monocytes', 'Classical_monocytes', rep('in_vitro_LYVE1_neg', 4), 'M1_like',
                       rep('Classical_monocytes', 4), 'in_vitro_LYVE1_pos'))
coldata <- data.frame(row.names = colnames(dds), condition)
dds <- DESeqDataSetFromMatrix(countData = dds, colData = coldata, design = ~condition)
dds <- DESeq(dds)
resultsNames(dds)
vsdata <- vst(dds, blind = F)
condition1 <- relevel(condition, 'in_vitro_LYVE1_neg')
dds$condition <- condition1
df <- plotPCA(vsdata, intgroup = 'condition', return = T)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds, name = "condition_in_vitro_LYVE1_pos_vs_in_vitro_LYVE1_neg")
library(apeglm)
resLFC <- lfcShrink(dds, coef=3 , type="apeglm")
plotMA(resLFC)
resLFC.df <- as.data.frame(resLFC)
ordered <- resLFC.df[order(resLFC.df$log2FoldChange, decreasing = T), ]
top1000 <- head(ordered, n = 1000)
##### ridgeplot, correlation heatmap, PCA #####
library(tidyverse)
vsdata <- vst(dds)
plotPCA(vsdata)
df <- t(assay(vsdata))
dist <- dist(df)
plot(hclust(dist))
########distance by using TPM #####
tpm <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}
a <- read.csv('batch12withgenelength.csv')
colnames(a)
b <- tpm(a[, c(3:12, 18, 22:27, 29:34)], a$GeneLength)
row.names(b) <- make.names(a$GeneName, unique = F)
colnames(b)
dists <- dist(t(b))
plot(hclust(dists), main = "Hierarchy")
tpm.res <- read.csv('tpmtable.csv')

###############################3

df %>% 
  dist() %>% 
  hclust() %>% 
  as.dendrogram() -> dend
par(mar=c(7,3,1,1))  # Increase bottom margin to have the complete label
plot(dend)
library(dendextend)
dend %>% 
  # Custom branches
  set("branches_col", "grey") %>% set("branches_lwd", 3) %>%
  # Custom labels
  set("labels_col", "orange") %>% set("labels_cex", 0.8) %>%
  plot()
par(mar = c(5, 5, 5, 10)) 
dend %>%
  set("labels_col", value = c("red", "blue", "green", "yellow"), k=4) %>%
  set("branches_k_color", value = c("red", "blue", "green", "yellow"), k = 4) %>%
  plot(horiz=TRUE, axes=FALSE) 
abline(v = 350, lty = 2) 
title('Hierachical Clustering for Bulk Samples')
#####dendro#####
library(factoextra)
fviz_dend2(x = dend, cex = 0.8, lwd = 0.8, k = 4,
          # Manually selected colors
          k_colors = c("red", "blue", "green", "black"),
          color_labels_by_k = T,
          rect = TRUE, 
          rect_border = c("red", "blue", "green", "black"), 
          rect_fill = TRUE, horiz = T,
          ggtheme = theme_grey())



plotMA(resLFC)
####dist heatmaps####
disMatrix <- as.matrix(dist)
rownames(disMatrix) <- paste(vsdata$condition, vsdata$type, sep="-")
colnames(disMatrix) <- condition
colors <- colorRamp2(c(0, 60 ,120), c("red", "white", "blue"))
disMatrix <- t(apply(disMatrix, 1, scale))
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
pheatmap(disMatrix)
Heatmap(disMatrix, 
        column_title_gp = gpar(fontsize = 10, fontface = "bold", fill = "pink"), 
        column_title = 'Correlation Among Samples \n by Euclidean Distance',
        row_split = 5, column_split = 5,
        show_row_names =F,
        show_column_names = F,
        column_names_rot = 60,
        row_dend_side = "right",
        width = unit(8, "cm"), height = unit(8, "cm"),
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("blue", "red","#016124",  "#e36002" , "purple")),
                                                            labels = c("classical \n monocytes", "in vitro \n M1-like", 'ubc \n LYVE-1+',  "in vitro \n LYVE-1-",  "in vitro \n LYVE-1+" 
                                                            ), 
                                                            labels_gp = gpar(col = "white", fontsize = 7, fontface = 'bold'))),
        left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = c("blue", "red","#016124",  "#e36002" , "purple")),
                                                         labels = c("classical \n monocytes", "in vitro \n M1-like", 'ubc \n LYVE-1+',  "in vitro \n LYVE-1-",  "in vitro \n LYVE-1+" 
                                                         ),
                                                            labels_gp = gpar(col = "white", fontsize = 7, fontface = 'bold'))),
        heatmap_legend_param = list(
          at = c(0, 100, 200, 300,400),
          title = "distance",
          legend_height = unit(4, "cm"),
          title_position = "lefttop-rot"))
lgd = Legend(labels = c("in vitro LYVE-1-", "classical monocytes", "in vitro LYVE-1+", "in vitro M1-like"),
             legend_gp = gpar(fill = c("green", "#c6eb34","#34ebab","#d9eb34" )))
draw(lgd, x = unit(1, "cm"), y = unit(1, "cm"), just = c("left", "bottom"))

#######genes heatmap####
top1000 <- resLFC.df[order(abs(resLFC.df$log2FoldChange), decreasing = TRUE), ]

top1000<- head(top1000, n = 1000)



genes <- c('LYVE1', 'FOLR2', 'F13A1', 'SELENOP', 'BNC2', 'TRPS1', 'SAMD4A', 'SNTB2', 'SLC18B1',
           'SLC2A13', 'PLPP3', 'ABCA1', 'STARD13', 'MYO1D', 'MRC1', 'DAAM2', 'ARHGEF12', 'SRGAP1',
           'FMN1', 'HRH1', 'TBC1D16')

genes <- c('F13A1', 'MRC1', 'SELENOP', 'LYVE1', 'PSAP', 'HLA.A', 'CD68',
            'HLA.DRA', 'VEGFC', 'VEGFA', 'VEGFB', 'NRP1')

genes <- c("F13A1", "FRMD4B", "RBPJ", "SLC9A9", "SEPP1", "LYVE1", "STAB1", "STARD13", "NAV2", "MPPED2", "MAN1A1", "PDE4D", "RNASE1", "RNF150", "MTSS1", "RGL1", "PLTP", "PMP22", "NRP1", "LGMN", "SCN9A", "MRC1", "MERTK", "CCDC141", "SIGLEC1", "COLEC12", "P2RY14", "AUTS2", "FGF13", "ME1", "CD163", "THRB", "EMB", "EDA", "CPM", "MS4A6A", "WWP1", "ITSN1", "DNM1", "DAAM2", "NPAS3", "MYO5A", "AP2A2", "MEF2C", "PDGFC", "EMP1", "PEAK1", "FMN1", "MAF", "TTN", "LRRC16A", "FCHSD2", "HDAC9", "WLS", "GNG2", "PTPRG", "MS4A4A", "MAFB", "LDLRAD4", "ARHGAP24", "ABCA6", "PDE7B", "CD163L1", "MUC6", "HRH1", "SESN1", "ABCA1", "ADAP2", "AP000304.12", "MS4A4E", "TRPS1", "BNC2", "EPS15", "CSGALNACT1", "NEURL2", "NAIP", "SLC40A1", "SLCO2B1", "IQGAP2", "RBMS1", "MAMDC2", "NPL", "NFIA", "C20orf194", "SYT1", "ZFP36L1", "ATP8B4", "RCAN1", "TBC1D14", "MED12L", "DIP2C", "ZFHX3", "TMEM2", "IL12RB2", "TCF4", "ITGA9", "DIRC3", "BCAT1", "SGMS1", "SLC18B1", "VSIG4", "WDFY3", "KHDRBS2", "EIF2B5", "HECTD2", "ST6GAL1", "ANKH", "CTTNBP2", "GABRB2", "PSD3", "GAS7", "FAM20A", "TGFBI", "TSHZ3", "FOLR2", "TPCN1", "ABCA9", "ARHGAP6", "CMKLR1", "EPS8", "FGFR1", "RCSD1", "HSF4", "CR1", "ADGRG6", "CEBPD", "CTDSPL", "RIN2", "BLVRB", "GPR34", "NISCH", "SIAH2", "APBA1", "STON2", "SCN1B", "ZDHHC14", "WDR36", "FAM13A", "PAPSS2", "EEF2K", "SNX9", "TDRD10", "LIMCH1", "LARGE", "DEPTOR", "RHOBTB1", "GPR155", "PLXND1", "DAAM1", "TNFRSF25", "EYA2", "SNCA", "PTK2", "BLNK", "RYR1", "SLC2A13", "ROR1", "TNS1", "NFATC2", "TSPAN4", "P2RY6", "PLEKHG5", "TNFRSF21", "CD209", "ETV5", "GAS6", "LIMD1", "FCGBP", "PLEKHA1", "MFHAS1", "RPS6KA2", "CFTR", "CD14", "SLC22A23", "SHMT1", "CREG1", "IGF1", "HPGDS", "ABCC5", "CLEC4G", "MYO7A", "AGAP1", "GLIS3", "SH3PXD2A", "SHE", "SLC36A1", "RCN3", "ANTXR1", "GNB1L", "GRAMD3", "PLEK", "CD28", "PRSS36", "PER3", "NOL3", "TGFBR3", "HSPA6", "CNRIP1", "NCKAP5", "IL2RA", "ARHGEF12", "FAM177B", "PLPP3", "IGSF21", "ELL", "P3H2", "RPH3AL", "ASAP2", "NRP2", "ICA1", "C4BPB", "MAP7", "ARL4C", "NINJ1", "CRYL1", "SEMA6A", "NAV3", "PRR5", "GARNL3", "FCGR2B", "CD200R1", "MAN1C1", "ZNF704", "PCDH12", "RIMKLB", "ACSM5", "VNN1", "ASTN2", "CCDC149", "MCC", "AKAP6", "A2M", "FUCA1", "ULK3", "GNG12", "EPN2", "TCN2", "DSCAML1", "CA8", "CCDC186", "MYLIP", "RAB3IL1", "HPN", "STOX2", "STARD8", "STXBP6", "LRRC2", "KANK2", "GFRA2", "CHID1", "TGFB2", "PRSS23", "RIC3", "SH3RF3", "GULP1", "ZNF812P", "LGI2", "HS3ST2", "SV2B", "SFXN2", "KLRC4-KLRK1", "LRP5L", "CXCL12", "KCNJ5", "TSPAN5", "KATNAL2", "TRIM47", "OSBPL6", "ABLIM3", "PKD1L3", "TNS2", "C5orf30", "TRAF5", "DAB2", "CHSY3", "DNAH2", "HTRA1", "TFAP2C", "SNX25", "GRAPL", "SCAMP5", "PDE2A", "SPTLC3", "KBTBD12", "TMEM37", "CCL8", "VIPR1", "ETV1", "GPRC5C", "TRAM2", "SOX5", "EPOR", "SDC3", "DCHS1", "CCL13", "FRK", "TRIM50", "WDR17", "MICU3", "PGBD1", "ALS2CL", "ABCD2", "RHPN1", "PYROXD2", "CALCR", "C9", "JPH1", "ARHGEF37", "TEKT3", "GATM", "VENTX", "EMP2", "PLAG1", "GSDMA", "LURAP1", "AASS", "SHROOM1", "PDGFB", "ECM1", "KITLG", "HHIPL1", "IGSF10", "C8orf46", "HOXB3", "SLC5A9", "FAT2", "LTBP2", "HTR1F", "WBSCR27", "COBLL1", "GIPC2", "SLC46A1", "CABP4", "CD247", "TMEM26", "LRP6", "SCML2", "PTPN21", "GSTA4", "DIRAS3", "SLC1A1", "QPRT", "RASL10A", "SPATA9", "TBKBP1", "KCNJ1", "CCDC136", "ADGRA3", "GPR12", "TDGF1", "EGFL7", "FOXRED2", "CH507-9B2.5", "ZKSCAN4", "ACSM4", "PWP2", "MEGF6", "DPYSL3", "SPATA25", "LRRC4", "SCTR", "FAM198A", "CRHBP", "CES3", "SULT1C4", "PTPRB", "RP11-65B7.2", "COLGALT2", "KRBA1", "CCDC152", "IRGM", "N6AMT1", "FCER2", "LBX2", "ZNF395", "HOXB7", "GATSL3", "FSTL3", "MAP1LC3B2", "MAPK12", "ARC", "ZNF585B", "TBX18", "EMID1", "PLEKHH3", "MRO", "CYTL1", "LTC4S", "SKIDA1", "ARHGEF16", "SLC22A3", "CRYBB1", "SSPO", "FAM212A", "CYP26B1", "DACT1", "SH2D6", "POM121L2", "FGF20", "HCRT", "CRB2", "SCARF2", "TGFB1I1", "TMPRSS11F", "TDRD1", "SYNDIG1L", "NEU4", "FSIP2", "ROBO4")



genes <- c('LYVE1', 'MRC1', 'F13A1', 'SELENOP', 'PTPRG', 'HDAC9', 'WWP1', 'CSGALNACT1',
           'RASSF4', 'MARCO', 'CD209', 'CLEC7A', 'IL1R2')
genes <- c('LYVE1', 'MRC1', 'F13A1', 'SELENOP', 'PTPRG', 'HDAC9', 'WWP1', 'CSGALNACT1',
           'RASSF4', 'MARCO', 'CD209', 'CLEC7A', 'IL1R2')
genes <- c('SLC1A4', 'SLC1A5', 'LDLR')
conditions_to_subset <- c('in_vitro_LYVE1_neg', 'Classical_monocytes')
genes_in_rlog <- genes
rlog.dds <- assay(rlog(dds))
mat2 <- rlog.dds[rownames(rlog.dds) %in% genes, 
                
                        coldata$condition %in% conditions_to_subset]
#genes <- c('LDLR', 'SLC1A4', 'SLC1A5')
#mat2 <- rlog.dds[rownames(rlog.dds) %in% genes, , drop = FALSE]
mat2
mat2.scaled <- t(apply(mat2, 1, scale)) #center and scale each column (Z-score) then transpose
colnames(mat2.scaled) <- colnames(mat2)
ann_col <- data.frame(
  Condition = condition
)
rownames(ann_col) <- colnames(mat2.scaled)
ann_colors$Condition <- c(
  "Classical_monocytes" = "red",
  "in_vitro_LYVE1_neg"  = "blue"
  #"in_vitro_LYVE1_pos"  = "#e36002",
  #"in_vitro_m1_like"    = "purple",
 # "ubc_LYVE1_pos"       = "cyan"
)
library(pheatmap)
pheatmap(
  mat2.scaled,
  cluster_rows = T,
  cluster_cols = T,
  annotation_col    = ann_col,
  annotation_colors = ann_colors,
  legend = TRUE,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  cellwidth  = 14.17,
  cellheight = 14.17,
  cutree_rows = 3,  # split rows into 3 groups
  cutree_cols = 2   # split cols into 2 groups
)
#num_positive <- colSums(mat2.scaled > 0)
#num_positive <- num_positive/409 *100
hlt <- c('LYVE1', 'F13A1', 'MRC1', 'MERTK', 'FOLR2')
matching_row <- which(rownames(mat2.scaled) %in% hlt)
col_fun = colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
ha <- rowAnnotation(foo = anno_mark(at = matching_row, labels = hlt, lines_gp = gpar(fill = 'blue')))
ha = rowAnnotation(foo = anno_density(mat2.scaled, type = "violin", 
                                      gp = gpar(fill = 1:30)))

ha = rowAnnotation(foo = anno_text(hlt, location = 0.5, just = "center"))

ha = rowAnnotation(foo = anno_histogram(mat2.scaled, gp = gpar(fill = 1:30)))
mat2.t <- t(mat2.scaled)
col_fun = colorRamp2(c(-2.5, 0, 2.5), c("#779906", "white", "#b34b0b"))
Heatmap(mat2.scaled,
        col = col_fun,
        column_split = factor(ann_col$Condition,
                              levels = c("in_vitro_m1_like", "in_vitro_LYVE1_neg",
                                         "Classical_monocytes", "ubc_LYVE1_pos",
                                         "in_vitro_LYVE1_pos")),
        column_title_gp = gpar(fontsize = 8, fontface = "bold", fill = "#caf2a7"),
        row_names_gp = gpar(fontface = 'bold'),
        show_column_dend = F,
        show_row_names = T,
        show_column_names = F,
        column_names_rot = 60,
        row_dend_side = "left",
        width = ncol(mat2.scaled)*unit(0.5, "cm"),
        height = nrow(mat2.scaled)*unit(0.5, "cm"),
        top_annotation = HeatmapAnnotation(
          Condition = ann_col$Condition,
          col = ann_colors,
          annotation_name_side = "left"
        ),
        heatmap_legend_param = list(
          at = c(-2.5, 0, 2.5),
          title = "z-score",
          legend_height = unit(4, "cm"),
          title_position = "lefttop-rot"))
#############eRASLAn's genes######
num_positive <- colSums(mat2.scaled > 0)
num_positive <- num_positive/409 * 100
data <- data.frame(Column = names(num_positive), Positive_Count = num_positive)
data <- data[order(data$Positive_Count, decreasing = T),]
data$Column <- factor(data$Column, levels = data$Column)
ggplot(data, aes(x = Column, y = Positive_Count, fill = Positive_Count)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(Positive_Count, 2)), vjust = -0.5, size = 3) + 
  labs(title = "Percentage of LYVE-1+ Markers Highly Expressed in Corresponding Samples",
       x = "Samples",
       y = "Pertage of LYVE-1+ Markers with Positive Z-score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_fill_gradient(low = "pink", high = "red")


ggplot(scree_plot, aes(x = pc_number, y = variance, fill = variance)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(variance, 2)), vjust = -0.5, size = 3) +  # Add labels
  labs(x = "Principal Component Number", y = "Variance") +
  ggtitle("Scree Plot for Principal Component Numbers") +
  scale_fill_gradient(low = "pink", high = "red")  # Set the range of fill colors




#####volcano
library(EnhancedVolcano)

genes_to_label <- c('C1QA', 'C1QB', 'C1QC', 'CD68', 'CD14', 'S100A8', 'S100A9',
                    'FCGR3A', 'LYZ', 'LYVE1', 'MRC1', 'FOLR2', 'SELENOP',
                    'F13A1', 'HLA.A', 'HLA.DRA')

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                selectLab = genes_to_label,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'grey50',
                boxedLabels = TRUE,
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2,
                labSize = 3.5,
                title = 'Volcano Plot',
                subtitle = NULL,
                caption = NULL)
##### fgsea
library(clusterProfiler)
library(fgsea)

# 1. ranked list from bulk DESeq2
ranked <- res$stat  # or log2FoldChange
names(ranked) <- rownames(res)
ranked <- sort(ranked, decreasing = TRUE)

# 2. gene set from sc top marker genes
# e.g. top 50 markers from your LYVE1+ sc cluster
gene_set <- list(LYVE1_signature = c("LYVE1", "MRC1", "F13A1", "FOLR2", ...))  # paste your marker genes

# 3. run fgsea
result <- fgsea(pathways = gene_set,
                stats    = ranked,
                minSize  = 10,
                maxSize  = 500)

# 4. plot
plotEnrichment(gene_set$LYVE1_signature, ranked) +
  ggtitle("LYVE1 Pos vs LYVE1 Neg")
#######scree####

rv <- rowVars(assay(vsdata))
## select the ntop genes by variance
select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
## perform a PCA on the data in assay(x) for the selected genes
pca <- prcomp(t(assay(vsdata)[select,]))
## the contribution to the total variance for each component
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
##plot the "percentVar"
scree_plot=data.frame(percentVar)
scree_plot[,2]<- c(1:22)
colnames(scree_plot)<-c("variance","pc_number")

library(ggplot2)

ggplot(scree_plot, aes(x = pc_number, y = variance, fill = variance)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(variance, 2)), vjust = -0.5, size = 3) +  # Add labels
  labs(x = "Principal Component Number", y = "Variance") +
  ggtitle("Scree Plot for Principal Component Numbers") +
  scale_fill_gradient(low = "pink", high = "red")  # Set the range of fill colors



ggplot(scree_plot, aes(x = pc_number, y = variance, group = 1)) +
  geom_line() +
  geom_point() +  # Add points for each PC
  geom_text(aes(label = round(variance, 2)), vjust = -0.5) +  # Add labels
  labs(x = "Principal Component Number", y = "Variance") +
  ggtitle("Scree Plot") +
  scale_x_continuous(breaks = seq(1, max(scree_plot$pc_number), by = 1))


#####Violinplot####
rlog.df <- as.data.frame(rlog.dds1)
vio.df <-rlog.df[c('LYVE1', 'FOLR2', 'SELENOP', 'F13A1'),] 
ggplot(timd4, aes(x = Condition, y = Number, fill = Condition)) +
  geom_violin() +
  geom_boxplot(width = 0.1, fill = "white", color = "black") +  # Add boxplot for additional information
  labs(x = "Condition", y = "Number") +
  ggtitle('Normalised TIMD4 Expression in different Conditions') +
  theme_minimal()  # Optional: Customize the theme if needed
                 

##### GSEA #####

sigs <- na.omit(resLFC.df)
sigs <- sigs[sigs$pvalue < 0.05, ]
sorted_sigs <- resLFC.df[order(resLFC.df$log2FoldChange, 
                          decreasing = TRUE), ]
#LFC need to be a named vector here
LFC <- sorted_sigs$log2FoldChange
name <- row.names(sorted_sigs)
names(LFC) <- name
LFC<-na.omit(LFC)
LFC = sort(LFC, decreasing = TRUE)
library(clusterProfiler)
library(org.Hs.eg.db)
gse <- gseGO(geneList=LFC, 
             ont = 'BP', 
             keyType = "SYMBOL", nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "none")
gsedf <- as.data.frame(gse)
gse1 <- pairwise_termsim(gse)
gse <- enrichplot::pairwise_termsim(gse)
a <- selected$Description
library(enrichplot)
emapplot(gse1, showCategory = c("postsynaptic endocytosis", "synaptic vesicle endocytosis", "presynaptic endocytosis",
                                                   "apoptotic cell clearance", "glycosaminoglycan metabolic process", "glycosaminoglycan catabolic process",
                                                   "clathrin-dependent endocytosis"), param_group_category = TRUE)

dotplot(gse, showCategory = c("postsynaptic endocytosis", "synaptic vesicle endocytosis", "presynaptic endocytosis",
                              "apoptotic cell clearance", "glycosaminoglycan metabolic process", "glycosaminoglycan catabolic process",
                              "clathrin-dependent endocytosis")) + ggtitle("dotplot for GSEA")

cnetplot(gse, showCategory =  c('synaptic vesicle endocytosis', 'cell adhesion mediated by integrin', 
                                'regulation of cell adhesion mediated by integrin', 'response to interleukin-4', 'cellular response to interleukin-4','apoptotic cell clearance',
                                'glycosaminoglycan catabolic process', 'hyaluronan catabolic process', 'folic acid transport', 'clathrin-dependent endocytosis'), 
         colorEdge = T, cirular = T)

treeplot(gse)
#enrichgo
ego <- enrichGO(gene = row.names(sorted_sigs),
universe      = LFC,
keyType = 'SYMBOL',
OrgDb         = org.Hs.eg.db,
ont           = "BP",
pAdjustMethod = "BH",
pvalueCutoff  = 0.01,
qvalueCutoff  = 0.05,
readable      = TRUE)
go.df <- as.data.frame(ego)
#####filter terms pathways ####
#MRC1|DNM1|AP2A2|CSF1R|TIMD4|STAB1|STAB2|CD163|CD209|PIEZO1|MERTK|AXL|TYRO3|S1PR1|LYVE1|FOLR2|CLTB|CLTC|TIMD4|
gseres <- gsedf %>%
  dplyr::filter(stringr::str_detect(pattern = "endocytosis|apoptotic|gly", Description)) %>%
  group_by(Description) %>%
  add_count() %>%
  dplyr::arrange(dplyr::desc(n), dplyr::desc(Description)) %>%
  mutate(Description= forcats::fct_inorder(Description))
selected_genes <- gseres[c(97,93,78,24,42,4),]
selected <- separate_rows(selected_genes, core_enrichment, sep = "/")

gseres <- gseres %>%
  dplyr::filter(stringr::str_detect(pattern = "MRC1|DNM1|AP2A2|CSF1R|TIMD4|STAB1|STAB2|CD163|CD209|PIEZO1|MERTK|AXL|TYRO3|S1PR1|LYVE1|FOLR2|CLTB|CLTC|TIMD4", core_enrichment)) %>%
  group_by(core_enrichment) %>%
  add_count() %>%
  dplyr::arrange(dplyr::desc(n), dplyr::desc(core_enrichment)) %>%
  mutate(core_enrichment= forcats::fct_inorder(core_enrichment))
gseres <- gseres[order(gseres$NES, decreasing = T), ]

gse.df <- selected_genes %>%
  mutate(core_enrichment = str_split(core_enrichment, "/")) %>%
  mutate(core_enrichment = sapply(core_enrichment, function(x) paste(head(x, 5), collapse = "/")))

ggplot(gse.df, aes(x = NES, y = Description, fill = p.adjust)) + 
  geom_bar(stat = "identity", aes(fill = p.adjust)) +  # Specify fill gradient here
  geom_label(aes(label = core_enrichment, x = 0), size = 5, hjust = 0, vjust = 0.5, fontface = "bold") +   
  labs(x = "Normalised Enrichment Scores", y = "Description") +
  ggtitle("Gene Set Enrichment Analysis \nPathways Up Regulated in in vitro LYVE-1+ Macrophages") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 15, colour = 'black', face = 'bold'), 
        plot.title = element_text(face = "bold")) +  # Making title bold
  scale_fill_gradient(low = "red", high = "pink")


ggplot(gseres, aes(x = NES, y = p.adjust, size = setSize, color = Description)) +
  geom_point() +
  scale_size_area(max_size = 5) +
  scale_colour_discrete()





cnetplot(gse, showCategory = c('positive regulation of endocytosis', 'regulation of mononuclear cell migration'), circular = TRUE, colorEdge = TRUE)
library(circlize)
library(ComplexHeatmap)

#######chorddiagram####
selected$core_enrichment
# order = order of genes; 
circos.initialize(sectors, c(0,1))

circos.rect(xleft, ybottom, xright, ytop, sector.index, track.index)
order1 = c('DNM1', 'PIEZO1', 'MRC1', 'MERTK', 'LYVE1', 'FOLR2', 'AP2A2', 'synaptic vesicle endocytosis', 'cell adhesion mediated by integrin', 
           'regulation of cell adhesion mediated by integrin', 'response to interleukin-4', 'cellular response to interleukin-4','apoptotic cell clearance',
           'glycosaminoglycan catabolic process', 'hyaluronan catabolic process', 'folic acid transport', 'clathrin-dependent endocytosis')
col1 = c('DNM1' = "#91b1e6", 'PIEZO1' = "#84dbae", 'MRC1' = "#82adb0", 'MERTK' = "#704c87", 'LYVE1' = "#7b74cf",
         'FOLR2' = '#a9ab6a', 'AP2A2'= 'Black',
         'synaptic vesicle endocytosis' = '#b4ff05', 'cell adhesion mediated by integrin' = '#ff8052', 
         'regulation of cell adhesion mediated by integrin' = '#ff00cc', 'response to interleukin-4' = '#05f1f5', 
         'cellular response to interleukin-4' = '#03ff03' ,'apoptotic cell clearance' = '#f5a002',
         'glycosaminoglycan catabolic process' = '#e37b8f', 'hyaluronan catabolic process' = '#fc4c4c', 
         'folic acid transport' = '#048ec4', 'clathrin-dependent endocytosis' = '#046110')
col_fun = colorRamp2(c(1.4, 2.1, 2.8), c("#ffffff",'#cf6d7a' ,"#FF0000"), transparency = 0.5)
col_fun1 = colorRamp2(c(0.005,0.02, 0.035), c("white",'#a9f5bc', "#2ae859"), transparency = 0.5)
col_fun2 <- colorRamp2(c(9, 50, 91), c("white",'#f2b679', "orange"), transparency = 0)
circos.par(start.degree = 90)
cdm_res <- chordDiagram(selected[, c(11,2,5)], annotationTrack ="grid", col = col_fun, 
             directional = 1,
             order = order1,
             grid.col = col1,
             direction.type = c("arrows"),
             link.arr.type = "big.arrow",
             link.sort = TRUE, preAllocateTracks = 1)
             
          


circos.clear()
#add the labels and axis

circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  #print labels 
  circos.text(x = mean(xlim), y = ylim[1] + 2.5, sector.name, facing = "clockwise",
                           niceFacing = T, adj = c(0.5,0), cex=0.8)  #bending.inside
  circos.axis(h = "bottom", labels.cex = 0.5, major.tick.length = 0.2, 
              sector.index = sector.name, track.index = 2)
   },
  bg.border = NA)


for (i in seq_along(cdm_res$rn)) {
  circos.text(x1, y1, sector.index = sector.name, labels = i, col = "white", 
              font = 2, cex = 0.5, adj = c(0.5, 0.5)) #bending.inside
}
circos.trackText(CELL_META$xcenter, CELL_META$ycenter,
                 labels = c('a', 'b', 'c'),
                 cex = 0.8, sectors = get.cell.meta.data("sector.index"), col = "#EEEEEE", font = 2, facing = "downward")

xlim = get.cell.meta.data("xlim", sector.index = selected$Description[1], track.index = 1)
x1 = xlim[1]
x2 = xlim[2]
sector.name = get.cell.meta.data("sector.index")
ylim = get.cell.meta.data("ylim", sector.index = selected$Description[1], track.index = 1)
y1 = ylim[1]
y2 = y1 + 0.6

xlim1 = get.cell.meta.data("xlim", sector.index = selected$core_enrichment[1], track.index = 1)
x3 = xlim1[1]
x4 = xlim1[2]
for (i in seq_along(cdm_res$rn)) {
  circos.rect(x1, y1, cdm_res$x1[i], y2, 
              col = col_fun1(selected$p.adjust[i]), 
              border = col_fun1(selected$p.adjust[i]),
              sector.index = selected$core_enrichment[i], track.index = 1)
  }

for (i in seq_along(cdm_res$cn)) {    
    circos.rect(x1, y1, cdm_res$x2[i], y2, 
                col = col_fun2(selected$setSize[i]), 
                border = col_fun2(selected$setSize[i]),
                sector.index = selected$Description[i], track.index = 1)
}

lgd1 = Legend(col_fun = col_fun, at = c(1,2,3,4), legend_gp = gpar(fill = 1:4), title = "NES")
lgd2 = Legend(col_fun = col_fun1, title = "p.adjust")
lgd3 = Legend(col_fun = col_fun2, at = c(16,32,48,64,92), legend_gp = gpar(fill = 5:9), title = 'setsize')
lgd4 = Legend(col_fun = col1, title = 'Genes and Pathways')


lgd = Legend(labels = order1, title = "Genes and Pathways", legend_gp = gpar(fill = col1),
             title_position = "lefttop")

pd = packLegend(lgd1, lgd2, lgd3, direction = "horizontal")



draw(lgd, x = unit(3, "cm"), y = unit(3, "cm"), just = c("left", "bottom"))
draw(pd, x = unit(3, "cm"), y = unit(3, "cm"), just = c("left", "bottom"))
circos.clear()




lgd_links = Legend(at = c(0,0.2,0.4,0.6), col_fun = col_fun, 
                   title_position = "topleft", title = "EnrichmentScores")

selected <- as.data.frame(selected)
pathways <- unique(selected$Description)
lgd_pathway = Legend(title = " Pathways", labels = pathways, legend_gp = gpar(fill = 1:15)
                     , ncol = 3, by_row = TRUE, title_position = "leftcenter", gap = unit(1, "cm"))

lgd_list_horizontal = packLegend(lgd_links, lgd_pathway,
                                 direction = "horizontal")

draw(lgd_list_horizontal, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))
circos.info()
circos.clear()
#
library(tidyr)
library(reshape2)
long_data <- separate_rows(selected_genes, core_enrichment, sep = "/")
long_data1 <- separate_rows(padj, core_enrichment, sep = "/")
long_data <- as.data.frame(long_data)
test <- long_data[c(6,7,14,31,33,37,44,52,59,117,133, 149,158,162,171,179),]
chordDiagram(test, transparency = 0.5)
rownames_to_extract <- c("GO:2000097", "GO:1990384", "GO:0002941", "GO:1901645", "GO:1901647", "GO:0001892", "GO:0097084", "GO:0060669", "GO:0060452", "GO:0061302", "GO:0060931")

require(DOSE)
#####KEGG #####
search_kegg_organism('hsa', by='kegg_code')
human <- search_kegg_organism('Homo sapiens', by='scientific_name')


kegg.df <- bitr(row.names(resLFC.df), fromType = 'SYMBOL', toType = 'UNIPROT'
                         , OrgDb = 'org.Hs.eg.db')
resLFC.df$SYMBOL <- row.names(resLFC.df)
sorted_sigs.kegg <- merge(resLFC.df, kegg.df, by = 'SYMBOL')

sorted_sigs.kegg <- sorted_sigs.kegg[order(sorted_sigs.kegg$log2FoldChange, 
                          decreasing = TRUE), ]
#LFC need to be a named vector here
LFC.kegg <- sorted_sigs.kegg$log2FoldChange
name <- sorted_sigs.kegg$UNIPROT
names(LFC.kegg) <- name
LFC.kegg<-na.omit(LFC.kegg)
LFC.kegg = sort(LFC.kegg, decreasing = TRUE)
options(clusterProfiler.download.method = "wget")
library(R.utils)
R.utils::setOption("clusterProfiler.download.method",'auto')
top1000 <- head(sorted_sigs.kegg, n = 1000)
kegg.enrich <- enrichKEGG(gene         = top1000$UNIPROT,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05,
                 keyType = "uniprot")
kegg.enrich.df <- as.data.frame(kegg.enrich)

gsekegg <- gseKEGG(geneList  = LFC.kegg,
                  organism     = 'hsa',
                  nPerm        = 10000,
                  minGSSize    = 3,
                  maxGSSize    = 800,
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "none",
                  keyType       = "uniprot")
gse.kegg.df <- as.data.frame(gsekegg)

gse.kegg.df <- gse.kegg.df %>%
  dplyr::filter(stringr::str_detect(pattern = "N-Glycan biosynthesis|Efferocytosis|Endocytosis", Description)) %>%
  group_by(Description) %>%
  add_count() %>%
  dplyr::arrange(dplyr::desc(n), dplyr::desc(Description)) %>%
  mutate(Description= forcats::fct_inorder(Description))


ggplot(gse.kegg.df, aes(x = NES, y = Description, fill = p.adjust)) + 
  geom_bar(stat = "identity", aes(fill = p.adjust)) +  # Specify fill gradient here
  #geom_label(aes(label = core_enrichment, x = 0), size = 5, hjust = 0, vjust = 0.5, fontface = "bold") +   
  labs(x = "Normalised Enrichment Scores", y = "Description") +
  ggtitle("KEGG \nPathways Up Regulated in in vitro LYVE-1+ Macrophages") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 15, colour = 'black', face = 'bold'), 
        plot.title = element_text(face = "bold")) +  # Making title bold
  scale_fill_gradient(low = "red", high = "pink")


















res.lyve1ubc <- resLFC.df
# endocytosis related#
all_genes_of_interest <- c("DNM1", "CANX", "AP2M1", "FCHO2", "DNM2", "ITSN1", "SNX9", "SH3GL1", "GAK", "AP2B1", 
                           "SYNJ1", "AP1G1", "VAMP4", "BIN1", "LRRK2", "PICALM", "ARF6", "SCRIB", "PIP5K1C", 
                           "SYT11", "ROCK1", "TOR1A", "STON2", "AP3M2", "PPP3CC", "DGKQ", "ITSN2", "CDK5",
                           "ARC", "AP2A2", "EPS15", "AP2M1", "ATAD1", "AP2B1", "VAC14", "NUMB", "HIP1", 
                           "SCRIB", "AP2A1", "CLTC", "SCYL2", "FNBP1L", "ITSN1", "BMP2K", "UBQLN2", 
                           "DAB2", "SIGLEC1", "GAK", "AP2B1", "PICALM", "DGKD", "PIK3CB", "PIP5K1C", 
                           "STON2", "AP2A1", "CCDC32", "ITSN2", "CLTB", "FCHO1", "FCHSD2", "NEU3", 
                           "UNC119", "CLTA", "LMBRD1")

# Remove duplicates if any
all_genes_of_interest <- unique(all_genes_of_interest)


genes_in_row_names <- all_genes_of_interest %in% row.names(top1000)
for (i in 1:length(all_genes_of_interest)) {
  print(paste(all_genes_of_interest[i], ": ", genes_in_row_names[i]))
}
#dnm1, canx,ARC,AP2A2,EPS15,AP2M1,FCHO2,CLTC
#migration
new_gene_list <- c("P2RY12", "CTTN", "NRP1", "NCK1", "SDCBP", "CD2AP",
                   "P2RY12", "CTTN", "NRP1", "CSPG4", "EPB41L5", "ANKS1A", "ATP5F1B", "NCK1", "SDCBP", "CD2AP", "ITGA11",
                   "P2RY12", "CD200R1", "SLAMF1", "CYP19A1", "PTK2", "CSF1R", "CMKLR1", "MMP14", "MTUS1", "C3AR1", "EMILIN1", "SLAMF8", "AKIRIN1", "MAPK1", "BCR",
                   "CD200R1", "IL27RA", "CRK", "KLRK1", "SLC8B1", "ADAM17", "FADD", "RIPK3", "AKT1", "ADTRP", "CRKL", "ADAM10", "MIA3", "MSN", "STK10", "ABL1", "NEDD9", "DOCK8", "PTK2B", "APP", "CORO1A", "CD99L2",
                   "P2RY12", "SLAMF1", "PTK2", "CSF1R", "CMKLR1", "MMP14", "C3AR1", "AKIRIN1", "MAPK1",
                   "P2RY12", "SAA1", "CD200R1", "SLAMF1", "CYP19A1", "PTK2", "CSF1R", "CMKLR1", "MMP14", "MTUS1", "C3AR1", "EMILIN1", "SLAMF8", "AKIRIN1", "MAPK1")

# Remove duplicates if any
new_gene_list <- unique(new_gene_list)
genes_in_row_names <- new_gene_list %in% row.names(top1000)
for (i in 1:length(all_genes_of_interest)) {
  print(paste(all_genes_of_interest[i], ": ", genes_in_row_names[i]))
}
#SNX9, AP1G1,VAMP4,BIN1,LRRK2,CLTB,AP2A2,ARC,CDK5,ITSN2,DGKQ,AP2M1,CANX,DNM1
