#Import Package####
library(ggplot2)
library(dplyr)
library(patchwork)
library(Seurat)
library(harmony)
library(readr)
library(clusterProfiler)
library(org.Hs.eg.db)
source("./Script/My_Function.R")

#Import Data####
Mel <- readRDS("./data/All_Celltypes_seurat_Arnon_data_metadata_M_mel_subset_revised.RDS")

#Save/Load Data####
save.image("./Mel_Pathway_analyse.RData")
load("./Mel_Pathway_analyse.RData")

#Clustering with pre-set resolutions####
Mel <- RunPCA(Mel)
Mel <- FindNeighbors(object = Mel, dims = 1:15)
Mel <- FindClusters(object = Mel, resolution = 0.4, cluster.name = "TEST_res_0.4")
Mel <- RunUMAP(Mel,
               dims = 1:15,
               min.dist = 0.5,
               spread = 1.5)
Mel_cluster_markers <- FindAllMarkers(Mel,
                                      only.pos = TRUE,
                                      min.pct = 0.25,
                                      logfc.threshold = 0.25)

# Pathway analyse####
# Prepare work
top_markers <- Mel_cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 15, wt = avg_log2FC)

#Gene Expression compare
DotPlot(Mel, 
        features = marker_db$`T cell`, 
        group.by = "TEST_res_0.4") + 
  RotatedAxis()

FeaturePlot(Mel, features = marker_db$`B cell`)

# Mel_cluster0####
Mel_cluster0_markers <- subset(Mel_cluster_markers, cluster == 0 & p_val_adj < 0.05)

# Convert gene symbols to Entrez ID
entrez_ids <- bitr(Mel_cluster0_markers$gene, 
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)

# GO enrichment analysis
ego <- enrichGO(gene          = entrez_ids$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",  # Biological Process
                pAdjustMethod = "BH",    # Multiple testing correction method
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.2,
                readable      = TRUE)

# KEGG pathway analysis
kk <- enrichKEGG(gene         = entrez_ids$ENTREZID,
                 organism     = 'hsa',  # human
                 pvalueCutoff = 0.05)

# Visualization
#dotplot(ego, showCategory=20)  # GO entry dot plot
png("GO_plot.png", width=3500, height=3500, res=300)
dotplot(ego, showCategory=20) 
dev.off()

#dotplot(kk, showCategory=20)   # KEGG pathway dot plot
png("kk_plot.png", width=3500, height=3500, res=300)
dotplot(kk, showCategory=20) 
dev.off()

## Cluster0-15####
# Create a dedicated directory to store results
output_dir <- "Cluster_Enrichment_Results"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Cycle Analysis Framework
for (cluster_num in 0:15) {
  tryCatch({
    # ---- Step 1: Genetic Screening ----
    current_markers <- subset(Mel_cluster_markers, 
                              cluster == cluster_num & p_val_adj < 0.05)
    
    # Empty data check
    if (nrow(current_markers) == 0) {
      message(paste0("Cluster ", cluster_num, ": No significant marker genes, skip analysis"))
      next
    }
    
    # ---- Step 2: ID conversion (adding fault tolerance mechanism) ----
    entrez_ids <- bitr(unique(current_markers$gene), # Remove duplicates to avoid repeated calculations
                       fromType = "SYMBOL",
                       toType = "ENTREZID",
                       OrgDb = org.Hs.eg.db)
    
    # Valid conversion check
    if (nrow(entrez_ids) < 5) { # A minimum of 5 genes are required for enrichment analysis
      message(paste0("Cluster ", cluster_num, ": Insufficient valid Entrez IDs, skipping analysis"))
      next
    }
    
    # ---- Step 3: Enrichment Analysis ----
    # GO enrichment analyse
    ego <- enrichGO(gene = entrez_ids$ENTREZID,
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2,
                    minGSSize = 10, # Setting minimum gene set size
                    maxGSSize = 500)
    
    # KEGG pathway analyse
    kk <- enrichKEGG(gene = entrez_ids$ENTREZID,
                     organism = 'hsa',
                     pvalueCutoff = 0.05,
                     minGSSize = 10,
                     maxGSSize = 500)
    
    # ---- Step 4: Result storage ----
    # Save original data
    saveRDS(list(ego = ego, kk = kk),
            file = file.path(output_dir, 
                             paste0("cluster", cluster_num, "_results.rds")))
    
    # Text format output
    write.csv(ego@result,
              file.path(output_dir, paste0("GO_cluster", cluster_num, ".csv")),
              row.names = FALSE)
    
    write.csv(kk@result,
              file.path(output_dir, paste0("KEGG_cluster", cluster_num, ".csv")),
              row.names = FALSE)
    
    # ---- Step 5: Visual Optimization ----
    # GO dotplot
    if (nrow(ego) > 0) {
      png(file.path(output_dir, paste0("GO_cluster", cluster_num, ".png")),
          width = 3000, height = 2500, 
          res = 300)
      print(dotplot(ego, 
                    showCategory = 20,
                    font.size = 14,
                    label_format = 50) + 
              ggtitle(paste("Cluster", cluster_num, "GO Enrichment")))
      dev.off()
    }
    
    # KEGG dotplot
    if (nrow(kk) > 0) {
      png(file.path(output_dir, paste0("KEGG_cluster", cluster_num, ".png")),
          width = 2800, height = 2000,
          res = 300)
      print(dotplot(kk, 
                    showCategory = 15,
                    font.size = 12) + 
              ggtitle(paste("Cluster", cluster_num, "KEGG Pathways")))
      dev.off()
    }
    
    message(paste0("Cluster ", cluster_num, " Analysis Complete"))
    
  }, error = function(e) {
    message(paste0("Cluster ", cluster_num, " Analysis failed: ", e$message))
  })
}