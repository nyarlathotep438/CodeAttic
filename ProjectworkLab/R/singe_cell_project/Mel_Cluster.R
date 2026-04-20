#Import Package####
library(ggplot2)
library(dplyr)
library(patchwork)
library(Seurat)
library(harmony)
library(readr)
source("./Script/My_Function.R")

#Memory Set
options(future.globals.maxSize = 10 * 1024^3)

#Import Data####
Mel.data <- readRDS("./data/All_Celltypes_seurat_Arnon_data_metadata_M_mel_subset_revised.RDS")
Mel.treatment <- read_csv("./data/Mel_treatment.csv")

#Save/Load Data####
save.image("./Mel_Cluster.RData")
load("./Mel_Cluster.RData")

#Data Analyse####
#General analyse
Mel_NA_Empty_Report <- Check_NA_Empty(Mel.data)
Mel_Columns <- names(Mel.data@meta.data)
Mel_Frequency_Columns <- c("cell_types",
                          "treatment.group",
                          "Cohort",
                          "tumor",
                          "treatment",
                          "Immune_resistance",
                          "Immune_resistance.up",
                          "Immune_resistance.down",
                          "technology",
                          "patient",
                          "sex",
                          "sample_primary_met",
                          "site",
                          "genetic_hormonal_features",
                          "outcome",
                          "chemotherapy_exposed",
                          "chemotherapy_response",
                          "targeted_rx_exposed",
                          "targeted_rx_response",
                          "ICB_exposed",
                          "pre_post",
                          "enough_cells",
                          "Study_name",
                          "Primary_or_met")
Mel_Frequency_analyse <- Frequency_Count_List(scdata = Mel.data, columns = Mel_Frequency_Columns)
write_csv(Mel_Frequency_analyse$ICB_exposed, "ICB_exposed.csv")

##Cluster
#Determine suitable dimensions
print(Mel.data@reductions$pca, dims = 1:5, nfeatures = 5)
ElbowPlot(Mel.data, ndims = 50)
DimPlot(Mel.data, reduction = "pca", dims = c(1,2)) + NoLegend()

#Try different resolutions
Run_cluster_umap_and_plot(scdata = Mel.data, resolutions = seq(0.9, 1.2, by = 0.1))

Mel.data <- FindNeighbors(object = Mel.data, dims = 1:15)
Mel.data <- FindClusters(object = Mel.data,
                        resolution = 0.4,
                        cluster.name = "TEST_res_0.4")
Mel.data <- RunUMAP(Mel.data,
                   dims = 1:15,
                   min.dist = 0.5,
                   spread = 1.5)
Mel_Dimplot <- DimPlot(Mel.data,
        reduction = "umap",
        group.by = "ICB_exposed")
ggsave("Mel_ICB.png", plot = Mel_Dimplot, width = 12, height = 6, dpi = 300)


cluster_markers = FindAllMarkers(Mel.data,
                                 only.pos = TRUE,
                                 min.pct = 0.25,
                                 logfc.threshold = 0.25)
top_markers = cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

top_marks_table = top_markers %>%
  group_by(cluster) %>%
  summarise(genes = paste(gene, collapse = "/")) %>%
  ungroup()
write_csv(top_marks_table, "cluster_genes.csv")