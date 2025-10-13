#Import Package####
library(dplyr)
library(readr)
library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)

#Import Data####
original_Mel <- readRDS("./data/All_Celltypes_seurat_Arnon_data_metadata_M_mel_subset_revised.RDS")
original_MBM <- readRDS("./data/All_Celltypes_seurat_MBM_Christopher_m_epi_subset_revised.RDS")

Mel_meta_cell_type <- read_csv("./data/Mel_meta_cell_type.csv")
MBM_cell_type <- read_csv("./data/MBM_cell_type.csv")

#Data Processing####
Mel_metastases <- subset(original_Mel, subset = sample_primary_met == "met")

Mel_metastases <- RunPCA(Mel_metastases)
Mel_metastases <- FindNeighbors(object = Mel_metastases, dims = 1:15)
Mel_metastases <- FindClusters(object = Mel_metastases, resolution = 0.3, cluster.name = "TEST_res_0.3")
Mel_metastases <- RunUMAP(Mel_metastases,
                          dims = 1:15,
                          min.dist = 0.5,
                          spread = 1.5)

MBM <- RunPCA(original_MBM)
MBM <- FindNeighbors(object = MBM, dims = 1:12)
MBM <- FindClusters(object = MBM, resolution = 0.7, cluster.name = "TEST_res_0.7")
MBM <- RunUMAP(MBM,
               dims = 1:12,
               min.dist = 0.5,
               spread = 1.5)

cell_type_dict <- setNames(Mel_meta_cell_type$`Cell type`, Mel_meta_cell_type$cluster)
Mel_metastases@meta.data$cell_type <- cell_type_dict[as.character(Mel_metastases@meta.data$seurat_clusters)]

cell_type_dict <- setNames(MBM_cell_type$`Cell type`, MBM_cell_type$cluster)
MBM@meta.data$cell_type <- cell_type_dict[as.character(MBM@meta.data$seurat_clusters)]

#Save Seurat project after processing####
saveRDS(MBM, file = "./data/MBM.RDS")
saveRDS(Mel_metastases, file = "./data/Mel_metastases.RDS")
