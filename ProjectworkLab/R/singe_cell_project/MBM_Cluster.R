#Import Function####
library(ggplot2)
library(dplyr)
library(patchwork)
library(Seurat)
library(harmony)
library(readr)
library(scmap)
source("./Script/My_Function.R")

#Memory Set
options(future.globals.maxSize = 10 * 1024^3)

#Import Data####
MBM.data <- readRDS("./data/All_Celltypes_seurat_MBM_Christopher_m_epi_subset_revised.RDS")
MBM.treatment <- read_csv("./data/MBM_treatment.csv")

#Save/Load Data####
save.image("./MBM_Cluster.RData")
load("./MBM_Cluster.RData")

#Data Analyse####
##General Analyse
MBM_Columns <- names(MBM.data@meta.data)
MBM_NA_Empty_Report <- Check_NA_Empty(MBM.data)
MBM_Frequency_Report <- Frequency_Count_All(MBM.data)

##Cluster
#Determine suitable dimensions
print(MBM.data@reductions$pca, dims = 1:5, nfeatures = 5)
ElbowPlot(MBM.data, ndims = 50)
DimPlot(MBM.data, reduction = "pca", dims = c(1,2)) + NoLegend()

#Try different resolutions
Run_cluster_umap_and_plot(scdata = MBM.data,
                          resolutions = seq(0.1, 1.2, by = 0.1),
                          dims = 1:12)

#Try resolution for find marker gene
resolutions <- c(0.3, 0.5, 0.6, 1.1)

for (res in c(0.6, 1.1)){
  cluster_col_name <- paste0("TEST_res_", res)
  MBM.data = FindClusters(object = MBM.data,
                          resolution = res,
                          cluster.name = cluster_col_name)
  
  #Find marker gene
  cluster_markers = FindAllMarkers(MBM.data,
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
  
  file_name <- paste0("cluster_genes_", res, ".csv")
  write_csv(top_marks_table, file_name)
}


#Final cluster
MBM.data <- RunUMAP(MBM.data, dims = 1:12, min.dist = 0.5, spread = 1.5)
MBM_Dimplot <- DimPlot(MBM.data,
                      reduction = "umap",
                      group.by = "Pre/post ICI")
ggsave("MBM_Pre_post ICI.png", plot = MBM_Dimplot, width = 8, height = 6, dpi = 300)
