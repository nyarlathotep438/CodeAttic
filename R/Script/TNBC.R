#Import Package####
library(ggplot2)
library(dplyr)
library(patchwork)
library(Seurat)
library(harmony)

#Import Data####
TNBC.normalized = read.table("./data/All_Celltypes_Rshinydata_hscSEG_normalized_data/Bassez_TNBC_data.csv",
                             header = TRUE,
                             row.names = 1,
                             sep=",")
TNBC.pseudobulk = read.table("./data/All_Celltypes_Rshinydata_pseudobulk/Bassez_TNBC_data.csv",
                             header = TRUE,
                             row.names = 1,
                             sep=",")

TNBC.RDS = readRDS("./data/All_Celltypes_seurat_Bassez_counts_treatment_naive_neoadjuvant_epi_subset_revised.RDS")
TNBC.RDS = subset(TNBC.RDS, subset = Cancer_type == 'TNBC')

#Data clean####
meta.table = data.frame(TNBC.RDS@meta.data)
counts.table = data.frame(TNBC.RDS@assays$RNA$counts)

ggplot(meta.table, aes(x = nCount_RNA)) + geom_histogram()
ggplot(meta.table, aes(x = nFeature_RNA)) + geom_histogram()
ggplot(counts.table, aes(x = M21_3p2_Lym_A03_S3)) + geom_histogram()

VlnPlot(MBM.RDS, features = c("nCount_RNA"), pt.size = 0) + NoLegend()
VlnPlot(MBM.RDS, features = c("nCount_RNA","nFeature_RNA"), pt.size = 0) + NoLegend()
FeatureScatter(MBM_RDS, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
