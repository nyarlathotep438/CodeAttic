#Import Package####
library(ggplot2)
library(ggrepel)
library(dplyr)
library(patchwork)
library(Seurat)
library(harmony)

#Import Data####
MBM.data = readRDS("./data/All_Celltypes_seurat_MBM_Christopher_m_epi_subset_revised.RDS")
Mel.data = readRDS("./data/All_Celltypes_seurat_Arnon_data_metadata_M_mel_subset_revised.RDS")

#Save/Load Data####
save.image("./Plot.RData")
load("./Plot.RData")

#Define function####
Frequency_Count = function(scdata, column){
  require(dplyr)
  meta.table = scdata@meta.data
  counts_table <- meta.table %>%
    dplyr::count({{ column }}, name = "Counts") %>%
    dplyr::arrange(dplyr::desc(Counts))
  return(counts_table)
}

#Mel DimPlot####
#Replace the original name for plot
Mel.data@meta.data$cell_types <- gsub("\\?", "UNKOWN", Mel.data@meta.data$cell_types)
Mel.data@meta.data$cell_types <- gsub("\\T.CD4", "T.cell", Mel.data@meta.data$cell_types)
Mel.data@meta.data$cell_types <- gsub("\\T.CD8", "T.cell", Mel.data@meta.data$cell_types)
Mel.data@meta.data$cell_types <- factor(Mel.data@meta.data$cell_types)

#Change the original order for pot
Mel_current_levels <- levels(Mel.data@meta.data$cell_types)
Mel_new_levels <- c("T.cell",
                    "B.cell",
                    "Macrophage",
                    "CAF",
                    "Endo.",
                    "Mal",
                    "UNKOWN")
Mel.data@meta.data$cell_types <- factor(Mel.data@meta.data$cell_types, levels = Mel_new_levels)

#Statistic
#Frequency Count
Mel_Cell_Counts_table = Frequency_Count(scdata = Mel.data,column =  cell_types)


#Set Colours for plot
cell_colors <- c(
  "T.cell" = "#E64B35FF",
  "B.cell" = "#4DBBD5FF",
  "Macrophage" = "#00A087FF",
  "CAF" = "#925E9FFF",
  "Endo." = "#800000",
  "Mal" = "#7FFD4F",
  "UNKOWN" = "#000000"
  )

#Core Plot code
Mel_Dimplot = DimPlot(Mel.data,
        reduction = "umap",
        group.by = "cell_types") +
  scale_color_manual(values = cell_colors) +
  labs(title = "Cell Types") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        legend.text = element_text(size = 10))

#Save plot
ggsave("Mel_CellTypes.png", plot = Mel_Dimplot, width = 8, height = 6, dpi = 300)

#MBM DimPlot####
#Replace the original name for plot
MBM.data@meta.data$cell_types <- gsub("\\hypoxia-associated monocyte-derived cluster 1",
                                      "HAMC-1", MBM.data@meta.data$cell_types)
MBM.data@meta.data$cell_types <- gsub("\\hypoxia-associated monocyte-derived cluster 2",
                                      "HAMC-2", MBM.data@meta.data$cell_types)
MBM.data@meta.data$cell_types <- gsub("\\calprotectin-high neutrophils",
                                      "Calpro+ Neutro", MBM.data@meta.data$cell_types)
MBM.data@meta.data$cell_types <- factor(MBM.data@meta.data$cell_types)

#Change the original order for pot
MBM_current_levels <- levels(MBM.data@meta.data$cell_types)

#Core Plot code
MBM_Dimplot = DimPlot(MBM.data,
        reduction = "umap",
        group.by = "cell_types") +
  labs(title = "Cell Types") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        legend.text = element_text(size = 10))

#Save plot
ggsave("MBM_CellTypes.png", plot = MBM_Dimplot, width = 8, height = 6, dpi = 300)
