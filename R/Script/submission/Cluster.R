#Import Package####
library(ggplot2)
library(dplyr)
library(patchwork)
library(Seurat)
library(harmony)
library(readr)
source("./Script/My_Function.R")

#Import Data####
Mel <- readRDS("./data/All_Celltypes_seurat_Arnon_data_metadata_M_mel_subset_revised.RDS")
MBM <- readRDS("./data/All_Celltypes_seurat_MBM_Christopher_m_epi_subset_revised.RDS")

#Save/Load Data####
save.image("./Cluster.RData")
load("./Cluster.RData")

#Data Analyse####
Mel_NA_Empty_Report <- Check_NA_Empty(Mel)
MBM_NA_Empty_Report <- Check_NA_Empty(MBM)

Mel_Frequency_Report <- Frequency_Count_All(Mel)
MBM_Frequency_Report <- Frequency_Count_All(MBM)

#Cluster####
#Preparation before clustering
Preparation_folders()

#Run PCA and choose the dimensions(Mel:1-15, MBM:1-12)
Mel <- RunPCA_PlotElbow(
  scdata = Mel,
  sample_name = "Mel"
)

MBM <- RunPCA_PlotElbow(
  scdata = MBM,
  sample_name = "MBM"
)

#Try different resolution and Visualization
Mel <- Try_different_resolutions(scdata = Mel,
                          resolutions = seq(0.1, 1.2, by = 0.1),
                          dims = 1:15)
MBM <- Try_different_resolutions(scdata = MBM,
                          resolutions = seq(0.1, 1.2, by = 0.1),
                          dims = 1:12)

#Run UMAP and create marker gene table to find the best resolution
Mel <- Get_marker_genes(scdata = Mel,
                 resolutions = c(0.1, 0.4, 1),
                 genes_amount = 15)

MBM <- Get_marker_genes(scdata = MBM,
                 resolutions = c(0.1, 0.2, 0.4, 0.6, 0.7, 1.2),
                 genes_amount = 15)

# Use the determined clustering results to perform formal clustering and draw heat map and feature plot
## Mel
Mel <- FindNeighbors(object = Mel, dims = 1:15)
Mel <- FindClusters(object = Mel, resolution = 0.4, cluster.name = "TEST_res_0.4")
Mel_cluster_markers <- FindAllMarkers(Mel,
                                      only.pos = TRUE,
                                      min.pct = 0.25,
                                      logfc.threshold = 0.25)

Get_marker_genes_all(scdata = Mel, resolutions = 0.4)

# Feature plot
ggsave("./plot/Mel/Mel_PTPRC_Feature.tiff", plot = FeaturePlot(Mel, features = c("PTPRC")), 
       width = 8, height = 6, dpi = 300)

ggsave("./plot/Mel/Mel_CD4_8A_Feature.tiff", plot = FeaturePlot(Mel, features = c("CD4","CD8A")), 
       width = 16, height = 6, dpi = 300)

# Universal exhaustion marker (applicable to both in situ and metastatic lesions)
Mel_exhaustion_core <- c("PDCD1",    # PD-1 (Immune Checkpoint Basics)
                     "CTLA4",    # CTLA-4 (early exhaustion marker)
                     "LAG3",     # LAG-3 (combination therapy target)
                     "TIGIT",    # TIGIT (Melanoma Emerging Target)
                     "HAVCR2")   # TIM-3 (metastasis-associated depletion)

# Metastasis-specific markers
Mel_metastasis_specific <- c("ENTPD1",  # CD39 (ATP metabolic depletion)
                         "CD38",    # Adenosine Signalling
                         "TOX")     # Exhaustion program regulator

# Visualization (in situ/metastasis comparison)
Mel_exhaustion_markers <- FeaturePlot(Mel, features = c(Mel_exhaustion_core, Mel_metastasis_specific), 
            split.by = "sample_primary_met",
            ncol = 4, 
            order = TRUE)

ggsave("./plot/Mel/Mel_exhaustion_markers_Feature.tiff",
       plot = Mel_exhaustion_markers, 
       width = 16, height = 36, dpi = 300)

# Heatmap
Mel_top5_markers <- Mel_cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

Mel_heatmap_plot <- DoHeatmap(
  object = Mel,
  features = Mel_top5_markers$gene,
  group.by = "ident",
  group.bar = TRUE,
  size = 3
) + 
  theme(axis.text.y = element_text(size = 8)) 

ggsave("./plot/Mel/Mel_marker_genes_heatmap.tiff", plot = Mel_heatmap_plot, 
       width = 12, height = 10, dpi = 300)

# MBM
MBM <- FindNeighbors(object = MBM, dims = 1:12)
MBM <- FindClusters(object = MBM, resolution = 0.7, cluster.name = "TEST_res_0.7")
MBM_cluster_markers <- FindAllMarkers(MBM,
                                      only.pos = TRUE,
                                      min.pct = 0.25,
                                      logfc.threshold = 0.25)

Get_marker_genes_all(scdata = MBM, resolutions = 0.7)

# Feature Plot
ggsave("./plot/MBM/MBM_PTPRC_Feature.tiff", plot = FeaturePlot(MBM, features = c("PTPRC")), 
       width = 8, height = 6, dpi = 300)

ggsave("./plot/MBM/MBM_CD4_8A_Feature.tiff", plot = FeaturePlot(MBM, features = c("CD4","CD8A")), 
       width = 16, height = 6, dpi = 300)

MBM_exhaustion_markers <- c("PDCD1",   # PD-1 (basic immune checkpoint marker)
                          "HAVCR2",  # TIM-3 (brain metastasis-specific exhaustion marker)
                          "LAG3",    # LAG-3 (synergistic inhibition with PD-1)
                          "TIGIT",   # TIGIT (Melanoma Immune Escape Key)
                          "ENTPD1",  # CD39 (Metabolic Depletion Core Molecule)
                          "TOX"      # Depletion of master transcription factors
)

ggsave("./plot/MBM/MBM_exhaustion_markers_Feature.tiff",
       plot = FeaturePlot(MBM, features = MBM_exhaustion_markers, ncol = 3), 
       width = 22, height = 12, dpi = 300)

# Heatmap
MBM_top5_markers <- MBM_cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

MBM_heatmap_plot <- DoHeatmap(
  object = MBM,
  features = MBM_top5_markers$gene,
  group.by = "ident",
  group.bar = TRUE,
  size = 3
) + 
  theme(axis.text.y = element_text(size = 8)) 

ggsave("./plot/MBM/MBM_marker_genes_heatmap.tiff", plot = MBM_heatmap_plot, 
       width = 12, height = 10, dpi = 300)
