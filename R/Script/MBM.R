#Import Package####
library(ggplot2)
library(dplyr)
library(patchwork)
library(Seurat)
library(harmony)

#Import Data####
MBM.normalized = read.table("./data/All_Celltypes_Rshinydata_hscSEG_normalized_data/Alvarez_Breckenridge_data.csv",
                            header = TRUE,
                            row.names = 1,
                            sep=",")
MBM.pseudobulk = read.table("./data/All_Celltypes_Rshinydata_pseudobulk/Alvarez_Breckenridge_data.csv",
                            header = TRUE,
                            row.names = 1,
                            sep=",")
original_MBM.RDS = readRDS("./data/All_Celltypes_seurat_MBM_Christopher_m_epi_subset_revised.RDS")

#Save/Load Data####
save.image("./MBM.RData")
load("./MBM.RData")

#Data clean####
#Data Statistic before Quality Control
#original_counts.table = data.frame(original_MBM.RDS@assays$RNA$counts)
original_meta.table = data.frame(original_MBM.RDS@meta.data)

ggplot(original_meta.table, aes(x = nCount_RNA)) +
  geom_histogram()

ggplot(original_meta.table, aes(x = nFeature_RNA)) + 
  geom_histogram()

#ggplot(counts.table, aes(x = M21_3p2_Lym_A03_S3)) + geom_histogram()
FeatureScatter(original_MBM.RDS, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

original_cell_stats = original_meta.table %>%
  count(cell_types, name = "Cell_Count") %>%
  arrange(desc(Cell_Count))

#Quality Control
original_MBM.RDS = PercentageFeatureSet(original_MBM.RDS, pattern = "^MT-", col.name = "percent_mito")
original_meta.table = data.frame(MBM.RDS@meta.data)

VlnPlot(MBM.RDS, features = c("nCount_RNA", "nFeature_RNA" ,"percent_mito"), pt.size = 0) + NoLegend()
VlnPlot(MBM.RDS, features = c("nCount_RNA"), pt.size = 0) + NoLegend()
VlnPlot(MBM.RDS, features = c("nFeature_RNA"), pt.size = 0) + NoLegend()
VlnPlot(MBM.RDS, features = c("percent_mito"), pt.size = 0) + NoLegend()

FeatureScatter(MBM.RDS, feature1 = "nCount_RNA", feature2 = "percent_mito")
FeatureScatter(MBM.RDS, feature1 = "nFeature_RNA", feature2 = "percent_mito")

MBM.RDS = subset(original_MBM.RDS,
                 subset = nFeature_RNA > 200 &
                   percent_mito < 10)

#Data Statistic after QC
QC_meta.table = data.frame(MBM.RDS@meta.data)

QC_cell_stats = QC_meta.table %>%
  count(cell_types, name = "Cell_Count") %>%
  arrange(desc(Cell_Count))

#Normalized####
##Classical Log Normalize Method
MBM.RDS = NormalizeData(MBM.RDS)
MBM.RDS = ScaleData(MBM.RDS)

##SCT Method
MBM.RDS = SCTransform(MBM.RDS, vars.to.regress = "percent_mito", verbose = TRUE)

MBM.RDS = FindVariableFeatures(MBM.RDS,
                               selection.method = "vst",
                               nfeatures = 2000)

#DefaultAssay(sc.data) = "RNA"
DefaultAssay(sc.data) = "SCT"

#PCA, UMAP & Clustering####
#PCA
MBM.RDS = RunPCA(MBM.RDS)

print(MBM.RDS@reductions$pca, dims = 1:5, nfeatures = 10)

VizDimLoadings(MBM.RDS, dims = 1:2, reduction = "pca", nfeatures = 10)
VizDimLoadings(MBM.RDS, dims = 1:4, reduction = "pca", nfeatures = 5)

DimHeatmap(MBM.RDS, dims = 1:5, cells = 500, balanced = TRUE)
DimHeatmap(MBM.RDS, dims = 1:2, cells = 100, balanced = TRUE)

DimPlot(MBM.RDS, reduction = "pca", dims = c(1,2)) + NoLegend()
DimPlot(MBM.RDS, reduction = "pca", dims = c(3,4)) + NoLegend()
DimPlot(MBM.RDS, reduction = "pca", dims = c(5,6)) + NoLegend()

#UMAP
ElbowPlot(MBM.RDS, ndims = 50)
MBM.RDS = RunUMAP(MBM.RDS,dims = 1:20, min.dist = 0.3, spread = 1)
DimPlot(MBM.RDS, reduction = "umap")

#Clustering
MBM.RDS = FindNeighbors(object = MBM.RDS, dims = 1:20)
MBM.RDS = FindClusters(object = MBM.RDS, cluster.name = "TEST1")

Clustered_meta.table = data.frame(MBM.RDS@meta.data)
DimPlot(MBM.RDS, reduction = "umap", label = TRUE)
