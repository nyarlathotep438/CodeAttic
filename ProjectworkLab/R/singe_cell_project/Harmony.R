#Import Package####
library(Seurat)
library(harmony)
library(dplyr)
library(readr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(reshape2)
source("./Script/My_Function.R")

# Import Data####
MBM <- readRDS("./data/MBM.RDS")
Mel_metastases <- readRDS("./data/Mel_metastases.RDS")

# Subset T cell####
T_MBM <- subset(MBM, subset = seurat_clusters %in% c(0,1,2))
T_Mel_metastases <- subset(Mel_metastases, subset = seurat_clusters %in% c(0,1,6))

#Save/Load Data####
save.image("./Harmony.RData")
load("./Harmony.RData")

# Harmony
# 1. Data merging and preprocessing----
# Add dataset identifier
T_MBM$orig.ident <- "MBM"
T_Mel_metastases$orig.ident <- "Mel_meta"
merged_T <- merge(T_MBM, T_Mel_metastases)

# Standardized Process
merged_T <- NormalizeData(merged_T) %>%
  FindVariableFeatures(nfeatures = 3000) %>%
  ScaleData(vars.to.regress = c("nCount_RNA")) %>% 
  RunPCA(npcs = 50, verbose = FALSE)

# 2. Harmony integrates batch effects----
merged_T <- RunHarmony(
  merged_T,
  group.by.vars = "orig.ident",
  theta = 2, # Adjust batch correction strength
  lambda = 1,
  plot_convergence = TRUE
)

# Check whether batch effects are eliminated
DimPlot(merged_T, reduction="harmony", group.by="orig.ident") + 
  ggtitle("Batch Integration Check")

# Verify that biological structure is preserved
DimPlot(merged_T, reduction="harmony", group.by="cell_type")

# 3. Dimensionality reduction and clustering----
merged_T <- merged_T %>%
  RunUMAP(reduction = "harmony", dims = 1:15) %>%
  FindNeighbors(reduction = "harmony", dims = 1:15) %>%
  FindClusters(resolution = 0.6)

# 4. Functional score calculation----
# Custom gene set
ifn_genes <- c("IFNG", "STAT1", "IRF1", "CXCL9") 
cyto_genes <- c("GZMB", "PRF1", "GNLY", "NKG7")
exhaust_genes <- c("PDCD1", "LAG3", "HAVCR2", "TIGIT")

merged_T <- AddModuleScore(
  merged_T,
  features = list(ifn_genes),
  name = "IFN_Gamma"
) %>% 
  AddModuleScore(
    features = list(cyto_genes),
    name = "Cytotoxicity"
  ) %>%
  AddModuleScore(
    features = list(exhaust_genes),
    name = "Exhaustion"
  )

# Extract rating data and merge with grouping information
scores <- FetchData(merged_T, 
                    vars = c("orig.ident", 
                             "IFN_Gamma1", 
                             "Cytotoxicity1", 
                             "Exhaustion1"))

# Long format conversion for easy batch processing
melted_scores <- melt(scores, id.vars = "orig.ident", 
                      variable.name = "Score", 
                      value.name = "Value")

# Statistical tests ----------------------------------------------------------------
# Batch Wilcoxon rank sum test
stats_res <- compare_means(
  Value ~ orig.ident, 
  data = melted_scores,
  group.by = "Score",
  method = "wilcox.test"
)
print(stats_res)

# Visual Comparison -----------------------------------------------------------------
# Generate box plot with significance markers
ggplot(melted_scores, aes(x = orig.ident, y = Value, fill = orig.ident)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  stat_compare_means(
    label = "p.signif", 
    method = "wilcox.test",
    comparisons = list(c("MBM", "Mel_meta"))
  ) +
  scale_fill_manual(values = c("MBM" = "#E64B35", "Mel_meta" = "#4DBBD5")) +
  facet_wrap(~Score, scales = "free_y", ncol = 3) +
  labs(x = "", y = "Module Score") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Violin plot
VlnPlot(merged_T, 
        features = c("IFN_Gamma1", "Cytotoxicity1", "Exhaustion1"),
        group.by = "orig.ident",
        pt.size = 0, 
        ncol = 3,
        cols = c("#E64B35", "#4DBBD5")) 

FeaturePlot(merged_T, 
            features = c("IFN_Gamma1", "Cytotoxicity1", "Exhaustion1"),
            split.by = "orig.ident",
            blend = FALSE,
            order = TRUE
)

#------------------------------------------------------------------#
# 5. CD8+ T cell subset extraction----
cd8_cells <- merged_T@meta.data %>%
  filter(
    grepl("CD8[+\\s].*\\bT\\b.*\\bcell[s]?\\b", 
          cell_type, 
          ignore.case = TRUE, 
          perl = TRUE)
  ) %>%
  rownames()

cd8_sub <- subset(merged_T, cells = cd8_cells)

remaining_cells <- setdiff(rownames(merged_T@meta.data), cd8_cells)
cd4_sub <- subset(merged_T, cells = remaining_cells)

saveRDS(merged_T,file = "./data/harmony/merged_T.RDS")
saveRDS(cd8_sub,file = "./data/harmony/cd8_sub.RDS")
saveRDS(cd4_sub,file = "./data/harmony/cd4_sub.RDS")

# 6. Integrated spatial activity contrast----
# Visualizing feature score differences
FeaturePlot(cd8_sub, 
            features = c("IFN_Gamma1", "Cytotoxicity1", "Exhaustion1"),
            split.by = "orig.ident",
            blend = FALSE,
            order = TRUE
)

FeaturePlot(cd4_sub, 
            features = c("IFN_Gamma1", "Cytotoxicity1", "Exhaustion1"),
            split.by = "orig.ident",
            blend = FALSE,
            order = TRUE
)

FeaturePlot(merged_T, 
            features = c("IFN_Gamma1", "Cytotoxicity1", "Exhaustion1"),
            split.by = "orig.ident",
            blend = FALSE,
            order = TRUE
)

DimPlot(object = merged_T, group.by = "cell_type", reduction = "umap")

# 7. Differential gene expression analysis----
# Activate integrated calibration data
DefaultAssay(cd8_sub) <- "RNA"
cd8_sub <- NormalizeData(cd8_sub)

# different analysis
de_markers <- FindMarkers(
  cd8_sub,
  ident.1 = "MBM",
  ident.2 = "Mel_meta",
  group.by = "orig.ident",
  logfc.threshold = 0.25,
  min.pct = 0.1,
  test.use = "wilcox" 
)

# Screening and annotation of key genes
top_de_genes <- de_markers %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5) %>%
  arrange(desc(avg_log2FC))

# 8. Key gene verification----
# Create a plotting data frame
volcano_data <- data.frame(
  gene = rownames(de_markers),
  log2FC = de_markers$avg_log2FC,
  p_val_adj = de_markers$p_val_adj,
  stringsAsFactors = FALSE
) %>% 
  mutate(
    significant = p_val_adj < 0.05 & abs(log2FC) > 0.5  # Mark significant genes
  )

# Draw Volcano plot
ggplot(volcano_data, aes(x = log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = significant), alpha = 0.6, size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +  # p-value threshold 
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey50") +   # log2FC threshold 
  scale_color_manual(values = c("grey", "red")) +                                  # Custom Colour
  geom_text_repel(                                                                 # label
    data = subset(volcano_data, significant), 
    aes(label = gene), 
    size = 3, 
    max.overlaps = 20  # Display up to 20 tags
  ) +
  labs(
    x = "Log2(Fold Change)", 
    y = "-Log10(Adj.P-value)", 
    title = "MBM vs Mel_meta"
  ) +
  theme_minimal(base_size = 12) +  
  theme(legend.position = "none") 
