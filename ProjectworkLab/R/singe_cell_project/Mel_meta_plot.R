#Import Package####
library(ggplot2)
library(ggrepel)
library(dplyr)
library(patchwork)
library(Seurat)
library(readr)
library(RColorBrewer)
library(devEMF)
library(clusterProfiler)
library(org.Hs.eg.db)
source("./Script/My_Function.R")

#Import Data####
Mel <- readRDS("./data/All_Celltypes_seurat_Arnon_data_metadata_M_mel_subset_revised.RDS")
Mel_metastases <- subset(Mel, subset = sample_primary_met == "met")
Mel_meta_cell_type <- read_csv("./data/Mel_meta_cell_type.csv")
Mel_meta_marker_gene <- read_csv2("./marker_gene/Mel_metastases/cluster_genes_res0.3.csv")

Mel_meta_GO_activeT <- readRDS("./Mel_metastases_Cluster_Enrichment_Results/cluster6_results.rds")
Mel_meta_GO_CD4T <- readRDS("./Mel_metastases_Cluster_Enrichment_Results/cluster0_results.rds")
Mel_meta_GO_CD8T <- readRDS("./Mel_metastases_Cluster_Enrichment_Results/cluster1_results.rds")
Mel_meta_GO_M2 <- readRDS("./Mel_metastases_Cluster_Enrichment_Results/cluster5_results.rds")

#Save/Load Data####
save.image("./Mel_meta_plot.RData")
load("./Mel_meta_plot.RData")

write.csv(Mel_meta_cell_type, file = "./data/Mel_meta_cell_type.csv", row.names = FALSE)

#Statistic
OR_result <- Mel_metastases@meta.data %>% 
  filter(outcome == "OR") %>% 
  group_by(cell_type) %>% 
  summarise(count = n()) %>%
  arrange(desc(count))

R_result <- Mel_metastases@meta.data %>% 
  filter(outcome == "R") %>% 
  group_by(cell_type) %>% 
  summarise(count = n()) %>%
  arrange(desc(count))

UT_result <- Mel_metastases@meta.data %>% 
  filter(outcome == "UT") %>% 
  group_by(cell_type) %>% 
  summarise(count = n()) %>%
  arrange(desc(count))

#Clustering with pre-set resolutions####
Mel_metastases <- RunPCA(Mel_metastases)
Mel_metastases <- FindNeighbors(object = Mel_metastases, dims = 1:15)
Mel_metastases <- FindClusters(object = Mel_metastases, resolution = 0.3, cluster.name = "TEST_res_0.3")
Mel_metastases <- RunUMAP(Mel_metastases,
               dims = 1:15,
               min.dist = 0.5,
               spread = 1.5)

Mel_metastases_cluster_markers <- FindAllMarkers(Mel_metastases,
                                      only.pos = TRUE,
                                      min.pct = 0.25,
                                      logfc.threshold = 0.25)

#Cell type annotation
cell_type_dict <- setNames(Mel_meta_cell_type$`Cell type`, Mel_meta_cell_type$cluster)
Mel_metastases@meta.data$cell_type <- cell_type_dict[as.character(Mel_metastases@meta.data$seurat_clusters)]

# UMAP Plot####
# Define the target clustering order
target_order <- c(0,1,6,2,16,5,12,3,4,7,8,9,10,11,13,14,15)

# ---- Step 1: Create sorting mapping rules ----
# Extract the correspondence between the current cluster and the cell type
cluster_celltype <- unique(Mel_metastases@meta.data[, c("seurat_clusters", "cell_type")])
# Sort by target order
ordered_celltypes <- cluster_celltype[
  order(match(cluster_celltype$seurat_clusters, target_order)),
  "cell_type"
] %>% as.character()

# ---- Step 2: Adjust metadata order ----
Mel_metastases@meta.data$cell_type <- factor(
  Mel_metastases@meta.data$cell_type,
  levels = ordered_celltypes  
)

# ---- Step 3: Create an ordered colour vector ----
# Original colour vector
original_colors <- c("#4169E1", "#1E90FF", "#000080", "#FFD700", "#FF8C00",
                     "#32CD32", "#8B4513", "#6B8E23", "#BA55D3", "#A52A2A",
                     "#8B0000", "#4B0082", "#800080", "#2F4F4F", "#B22222",
                     "#696969","#654321")

# Create a new mapping data frame to ensure the correct colour correspondence
color_mapping <- data.frame(
  cell_type = ordered_celltypes,
  color = original_colors[1:length(target_order)] 
)
custom_ordered_colors <- as.character(color_mapping$color)

# ---- Step 4: Verify mapping correctness ----
print(data.frame(
  Cluster = target_order,
  CellType = ordered_celltypes,
  AssignedColor = custom_ordered_colors
))

# ---- Step 5: Drawing ----
group <- c("pre_post", "cell_type", "Primary_or_met")

pre_post_counts <- table(Mel_metastases@meta.data$pre_post)
print(pre_post_counts)

p <- DimPlot(
  object = Mel_metastases,
  
  group.by = group[2],
  cols = custom_ordered_colors,  # Choose to use the costumed colors
  label = FALSE,
  pt.size = 0.6,
  reduction = "umap"
) + 
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle("") +
  theme_bw(base_size = 12) +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_blank(),   
    plot.background = element_blank(),
    panel.grid = element_blank()      
  ) +
  guides(color = guide_legend(override.aes = list(size = 4)))

p_no_legend <- p + theme(legend.position = "none")

p_legend_bottom <- p + 
  theme(
    legend.position = "bottom",  
    legend.justification = "center",  
    legend.box.spacing = unit(0.3, "cm")  
  ) +
  guides(color = guide_legend(
    override.aes = list(size = 4),
    nrow = 9,  
    title.position = "top"  
  ))

print(p)
print(p_no_legend)
print(p_legend_bottom)

svg("UMAP_Mel_meta_post_pre.svg")
print(p)
dev.off()

svg("UMAP_Mel_meta_cell_type_noLegend.svg")
print(p_no_legend)
dev.off()

svg("UMAP_Mel_meta_cell_type_Legend_bottom.svg", width = 10, height = 10)
print(p_legend_bottom)
dev.off()

# Dot Plot####
Mel_metastases$TEST_res_0.3 <- factor(
  Mel_metastases$TEST_res_0.3,
  levels = c("0","1","6","2","16","5","12","3","4","7","8","9","10","11","13","14","15") 
)

dot_p <- DotPlot(Mel_metastases, 
        features = marker_db$`T cell`, 
        group.by = "TEST_res_0.3") + 
  RotatedAxis() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

svg("UMAP_Mel_meta_dot.svg")
print(dot_p)
dev.off()

emf("UMAP_Mel_meta_dot.emf")
print(dot_p)
dev.off()

# Heatmap####
# 提取目标亚群
t_cells <- subset(Mel_metastases, idents = c(0,1,6))

# 设置直观Cluster命名
t_cells$named_clusters <- factor(
  t_cells$seurat_clusters, 
  labels = c("CD8+ T cell", "CD4+ T cell", "Cycling CD8+ T cell")
)

# activation
activation_genes <- c("CD69", "IL2RA", "MKI67", "IFNG", "GZMB", "PRF1")
heatmap_active_p <- DoHeatmap(t_cells,
          features = activation_genes,
          group.by = "named_clusters",
          slot = "scale.data",    # 使用标准化数据
          disp.min = -2,          # 截断Z-score范围
          disp.max = 2,
          angle = 45,             # 调整标签角度
          size = 4) + 
  scale_fill_gradientn(colors = c("#2E22EA","#F9F8F6","#FD1222"))

# exhaustion
exhaustion_genes <- c("PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT")
DoHeatmap(t_cells,
          features = exhaustion_genes,
          group.by = "named_clusters",
          cells = order(t_cells$named_clusters), # 按cluster排序
          size = 4,
          disp.min = -1.5,
          disp.max = 2.5) +
  scale_fill_gradientn(colors = c("#2E22EA","#F9F8F6","#FD1222"))


# 按细胞类型拆分成三个子对象
split_groups <- split(t_cells$group, t_cells$named_clusters)

# 生成每个细胞类型的对比图
heatmap_list <- lapply(names(split_groups), function(cell_type) {
  # 提取当前细胞类型的子集
  sub_obj <- subset(t_cells, named_clusters == cell_type)
  
  # 排列顺序：确保Pre在前，Post在后
  ordered_cells <- Cells(sub_obj)[order(sub_obj$pre_post, decreasing = FALSE)]
  
  # 双面板图布局（激活标记左，耗竭标记右）
  p1 <- DoHeatmap(sub_obj,
                  features = activation_genes,
                  cells = ordered_cells,
                  group.by = "pre_post",
                  slot = "scale.data",
                  angle = 30,                 # 更小的标签角度
                  size = 3.5,                 # 缩小标签尺寸
                  disp.min = -2, 
                  disp.max = 2) +
    scale_fill_gradientn(colors = c("#2E22EA", "#F9F8F6", "#FD1222")) +
    ggtitle(paste(cell_type, "Activation")) 
  
  p2 <- DoHeatmap(sub_obj,
                  features = exhaustion_genes,
                  cells = ordered_cells,
                  group.by = "pre_post",
                  slot = "scale.data",
                  angle = 30,
                  size = 3.5,
                  disp.min = -1.5, 
                  disp.max = 2.5) +
    scale_fill_gradientn(colors = c("#2E22EA", "#F9F8F6", "#FD1222")) +
    ggtitle(paste(cell_type, "Exhaustion"))
  
  # 组合双面板
  p1 + p2 + plot_layout(widths = c(1, 0.8))  # 调整面板宽度比例
})

# 独立展示每个子图（以CD8+为例） 
heatmap_list[[1]]  # CD8+ T细胞
heatmap_list[[2]]  # CD4+ T细胞
heatmap_list[[3]]  # Cycling CD8+




sub_obj <- subset(t_cells, named_clusters == cell_type)
VlnPlot(sub_obj, 
        features = c("CD69", "LAG3"),
        split.by = "pre_post",
        pt.size = 0, 
        ncol = 2) +
  geom_violin(scale = "width", adjust = 1.2) +
  theme(legend.position = "none")

svg("heatmap_Mel_meta_active.svg")
print(heatmap_active_p)
dev.off()

emf("heatmap_Mel_meta_active.emf")
print(heatmap_active_p)
dev.off()

# GO
dotplot(Mel_meta_GO_M2$ego, 
        showCategory = 20,          # 显示前20个条目
        color = "p.adjust",         # 颜色映射调整p值
        label_format = 60) +        # 条目名称换行长度
  ggtitle("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # x轴标签倾斜45度




# 提取前6类中心坐标
centroid_df <- FetchData(Mel_metastases, vars = c("umap_1", "umap_2", group[2])) |>
  group_by(across(all_of(group[2]))) |>
  summarise(
    x = median(umap_1)+ 2.5,
    y = median(umap_2) + 3,  
    .groups = "drop"
  ) |>
  slice_head(n = 6)

p <- DimPlot(
  object = Mel_metastases,
  
  group.by = group[2],
  cols = custom_ordered_colors,  # Choose to use the costumed colors
  label = FALSE,
  pt.size = 0.6,
  reduction = "umap"
) + 
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle("") +
  theme_bw(base_size = 12) +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_blank(),   
    plot.background = element_blank(),
    panel.grid = element_blank()      
  ) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  geom_text(
    data = centroid_df,
    aes(x = x, y = y, label = !!sym(group[2])),  # 动态引用分组列名
    size = 3.5,
    color = "black",
    fontface = "bold"
  )

p_no_legend <- p + theme(legend.position = "none")

p_legend_bottom <- p + 
  theme(
    legend.position = "bottom",  
    legend.justification = "center",  
    legend.box.spacing = unit(0.3, "cm")  
  ) +
  guides(color = guide_legend(
    override.aes = list(size = 4),
    nrow = 9,  
    title.position = "top"  
  ))