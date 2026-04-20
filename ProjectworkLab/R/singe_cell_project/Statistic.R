##This script is used to count the meta-table data of the Seurat project
#Import Package####
library(dplyr)
library(patchwork)
library(Seurat)
library(harmony)

#Import Data####
MBM.data = readRDS("./data/All_Celltypes_seurat_MBM_Christopher_m_epi_subset_revised.RDS")
Mel.data = readRDS("./data/All_Celltypes_seurat_Arnon_data_metadata_M_mel_subset_revised.RDS")

Mel.treatment = read.table("./data/Mel_treatment.csv",
                       header = TRUE,
                       row.names = 1,
                       sep=",")
MBM.treatment = read.table("./data/MBM_treatment.csv",
                           header = TRUE,
                           row.names = 1,
                           sep=",")

#Save/Load Data####
save.image("./Statistic.RData")
load("./Statistic.RData")

#Define function####
#Frequency Count
Frequency_Count = function(scdata, column){
  require(dplyr)
  meta.table = scdata@meta.data
  frequency_counts_table <- meta.table %>%
    dplyr::count({{ column }}, name = "Counts") %>%
    dplyr::arrange(dplyr::desc(Counts))
  return(frequency_counts_table)
}

#Unique Value Count
Unique_Value_Count = function(scdata, column){
  require(dplyr)
  meta.table = scdata@meta.data
  num_unique <- meta.table %>%
    dplyr::summarise(n_distinct({{column}})) %>%
    pull()
  return(num_unique)
}

#Get Statistic result####
#Mel
Mel_Cell_Frequency_table = Frequency_Count(scdata = Mel.data,column =  cell_types)
Mel_donor_Frequency_table = Frequency_Count(scdata = Mel.data,column =  donor_id)
Mel_patient_Frequency_table = Frequency_Count(scdata = Mel.data,column =  patient)
Mel_orig_Frequency_table = Frequency_Count(scdata = Mel.data,column =  orig.ident)
Mel_extent_Frequency_table = Frequency_Count(scdata = Mel.data,column =  disease_extent)
Mel_Cohort_Frequency_table = Frequency_Count(scdata = Mel.data,column =  Cohort)
Mel_tumour_Frequency_table = Frequency_Count(scdata = Mel.data,column =  tumor)
Mel_LABELS_Frequency_table = Frequency_Count(scdata = Mel.data,column =  LABELS)
Mel_treatment_Frequency_table = Frequency_Count(scdata = Mel.data,column =  treatment)
Mel_sex_Frequency_table = Frequency_Count(scdata = Mel.data,column =  sex)
Mel_ICB_table = Frequency_Count(scdata = Mel.data,column =  ICB_exposed)

donor_in_dataset = Mel_donor_Frequency_table$donor_id
filtered_out_donors = setdiff(all_donors, donor_in_dataset)
Mel_columns = names(Mel.data@meta.data)

#MBM
MBM_Cell_Frequency_table = Frequency_Count(scdata = MBM.data,column =  cell_types)
#MBM_Cell_sex_table = Frequency_Count(scdata = MBM.data,column =  sex)
