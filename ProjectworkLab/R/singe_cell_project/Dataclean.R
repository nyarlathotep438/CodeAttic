#Import Package####
library(ggplot2)
library(dplyr)
library(patchwork)
library(Seurat)
library(harmony)

#Import Data####
MBM.data = readRDS("./data/All_Celltypes_seurat_MBM_Christopher_m_epi_subset_revised.RDS")
Mel.data = readRDS("./data/All_Celltypes_seurat_Arnon_data_metadata_M_mel_subset_revised.RDS")

#Save/Load Data####
save.image("./Meta_Data.RData")
load("./Meta_Data.RData")

#Define function####
##Data Statistic & Quality Control
#Frequency count for any column
Frequency_Count = function(scdata, column){
  require(dplyr)
  meta.table = scdata@meta.data
  counts_table <- meta.table %>%
    dplyr::count({{ column }}, name = "Counts") %>%
    dplyr::arrange(dplyr::desc(Counts))
  return(counts_table)
}

#Cell Count table
Cells_counts = function(meta.table){
  require(dplyr)
  cells_counts_table <- meta.table %>%
    count(cell_types, name = "Cell_Count") %>%
    arrange(desc(Cell_Count))
  return(cells_counts_table)
}

#Donors Number Count
Donors_number = function(meta.table){
  require(dplyr)
  num_donors <- meta.table %>% 
    summarise(n_distinct(donor_id)) %>% 
    pull()
  return(num_donors)
}

#Histogram
Yield_histogram = function(meta.table){
  require(ggplot2)
  nCount_RNA_plot = ggplot(meta.table, aes(x = nCount_RNA)) +
    geom_histogram()
  nfeature_RNA_plot = ggplot(meta.table,aes(x = nFeature_RNA)) + 
    geom_histogram()
  histogram = list("nCount_RNA" = nCount_RNA_plot,
                   "nFeature_RNA" = nfeature_RNA_plot
                   )
  return(histogram)
}

#Feature Scatter
Yield_Feature_Scatter = function(scdata){
  require(Seurat)
  C_v_F_Scatter = FeatureScatter(scdata,
                                 feature1 = "nCount_RNA",
                                 feature2 = "nFeature_RNA")
  C_v_m_Scatter = FeatureScatter(scdata,
                                 feature1 = "nCount_RNA",
                                 feature2 = "percent_mito")
  F_v_m_Scatter = FeatureScatter(scdata,
                                 feature1 = "nCount_RNA",
                                 feature2 = "percent_mito")
  scatter = list("C_v_F_Scatter" = C_v_F_Scatter,
                 "C_v_m_Scatter" = C_v_m_Scatter,
                 "F_v_m_Scatter" = F_v_m_Scatter)
  return(scatter)
}

#Data Statistic Report
Data_summary = function(meta.table, scdata){
  Cells_counts = Cells_counts(meta.table)
  Donors_number = Donors_number(meta.table)
  histogram = Yield_histogram(meta.table)
  scatter = Yield_Feature_Scatter(scdata)
  Data_summary = list("Cells_counts" = Cells_counts,
                      "Donors_number" = Donors_number,
                      "histogram" = histogram,
                      "scatter" = scatter)
  return(Data_summary)
}

##Print unified format plot
Save_Plot= function(plot, plot.path, plot.height, plot.width){
  tiff(plot.path, width = plot.width, height = plot.height)
  print(plot)
  dev.off()
}

#Data clean####
#Add the column about the percentage of Mitochondria RNA
MBM.data = PercentageFeatureSet(MBM.data,
                                pattern = "^MT-",
                                col.name = "percent_mito")
Mel.data = PercentageFeatureSet(Mel.data,
                                pattern = "^MT-",
                                col.name = "percent_mito")

#Get Meta data table
MBM_meta.table = data.frame(MBM.data@meta.data)
Mel_meta.table = data.frame(Mel.data@meta.data)

#Compare the data
cols_MBM = colnames(MBM_meta.table)
cols_Mel = colnames(Mel_meta.table)
common_columns = intersect(cols_MBM, cols_Mel)

#Data describe report
MBM_Data_summary = Data_summary(MBM_meta.table, MBM.data)
Mel_Data_summary = Data_summary(Mel_meta.table, Mel.data)

#Normalized,Scaled & Find Variable gene####
##SCT method
MBM.data = SCTransform(MBM.data,
                       vars.to.regress = "percent_mito",
                       verbose = TRUE)
Mel.data = SCTransform(Mel.data,
                       vars.to.regress = "percent_mito",
                       verbose = TRUE)

#PCA, UMAP & Clustering####
#PCA
MBM.data = RunPCA(MBM.data)
Mel.data = RunPCA(Mel.data)

#UMAP
MBM.data = RunUMAP(MBM.data,
                   dims = 1:20,
                   min.dist = 0.3,
                   spread = 1)
Mel.data = RunUMAP(Mel.data,
                   dims = 1:20,
                   min.dist = 0.3,
                   spread = 1)

#Clustering
MBM.data = FindNeighbors(object = MBM.data, dims = 1:20)
MBM.data = FindClusters(object = MBM.data, cluster.name = "TEST1")

Mel.data = FindNeighbors(object = Mel.data, dims = 1:20)
Mel.data = FindClusters(object = Mel.data, cluster.name = "TEST1")
