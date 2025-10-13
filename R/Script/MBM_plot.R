#Import Package####
library(ggplot2)
library(ggrepel)
library(dplyr)
library(patchwork)
library(Seurat)
library(readr)
library(RColorBrewer)
library(devEMF)
library(tidyr)
library(clusterProfiler)
library(org.Hs.eg.db)
source("./Script/My_Function.R")

#Import Data####
MBM <- readRDS("./data/MBM.RDS")

#Save/Load Data####
save.image("./MBM_plot.RData")
load("./MBM_plot.RData")

write.csv(MBM_cell_type, file = "./data/MBM_cell_type2.csv", row.names = FALSE)
MBM_cell_type <- read_csv("./data/MBM_cell_type2.csv")

#Cluster
MBM <- RunPCA(MBM)
MBM <- FindNeighbors(object = MBM, dims = 1:12)
MBM <- FindClusters(object = MBM, resolution = 0.4, cluster.name = "TEST_res_0.4")
MBM <- RunUMAP(MBM,
               dims = 1:12,
               min.dist = 0.5,
               spread = 1.5)

MBM_cluster_markers <- FindAllMarkers(MBM,
                                      only.pos = TRUE,
                                      min.pct = 0.25,
                                      logfc.threshold = 0.25)


DimPlot(object = MBM, reduction = "umap")

top15_genes <- MBM_cluster_markers %>% 
  group_by(cluster) %>%
  dplyr::slice_head(n = 15) %>% 
  split(.$cluster) %>% 
  lapply(function(x) x$gene)

#Get_marker_genes_all(scdata = MBM, resolutions = 0.4)

# UMAP Plot####
# Define the target clustering order
target_order <- c(0,1,11,13,14,20,3,10,4,19,2,5,6,7,8,9,12,15,16,17,18)

#Cell type annotation
cell_type_dict <- setNames(MBM_cell_type$`Cell type`, MBM_cell_type$cluster)
MBM@meta.data$cell_type <- cell_type_dict[as.character(MBM@meta.data$seurat_clusters)]

# ---- Step 1: Create sorting mapping rules ----
# Extract the correspondence between the current cluster and the cell type
cluster_celltype <- unique(MBM@meta.data[, c("seurat_clusters", "cell_type")])
# Sort by target order
ordered_celltypes <- cluster_celltype[
  order(match(cluster_celltype$seurat_clusters, target_order)),
  "cell_type"
] %>% as.character()

# ---- Step 2: Adjust metadata order ----
MBM@meta.data$cell_type <- factor(
  MBM@meta.data$cell_type,
  levels = ordered_celltypes  
)

# ---- Step 3: Create an ordered colour vector ----
# Original colour vector
original_colors <- c("#4169E1", "#1E90FF", "#4682B4", "#000080", "#FFD700", 
                     "#FF8C00", "#FF7F55", "#CC5500", "#F4A460", "#D2691E", 
                     "#8B0000", "#800080", "#FF4500", "#B22222", "#8B0010",
                     "#4B0082", "#9400D3", "#8A2BE2", "#9370DB", "#6A5ACD",
                     "#2E8B57")

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
group <- c("Pre/post ICI", "cell_type")

p <- DimPlot(
  object = MBM,
  group.by = group[2],
  
  cols = custom_ordered_colors, # Choose to use the costumed colors
  
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
    legend.text = element_text(size = 9),
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
    nrow = 11,  
    title.position = "top"  
  ))


svg("UMAP_MBM_Legend_bottom.svg", width = 10, height = 11)
print(p_legend_bottom)
dev.off()

# Dot Plot####
MBM$TEST_res_0.4 <- factor(
  MBM$TEST_res_0.4,
  levels = c(0,1,11,13,14,20,3,10,4,19,2,5,6,7,8,9,12,15,16,17,18) 
)

dot_p <- DotPlot(MBM, 
                 features = marker_db$`T cell`, 
                 group.by = "TEST_res_0.4") + 
  RotatedAxis() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


# Heatmap
# 提取目标亚群
t_cells <- subset(MBM, idents = c(0,13,1,11))

cd8_cells <- subset(MBM, idents = c(0))
cyc_cd8_cells <- subset(MBM, idents = c(13))
cd4_cells <- subset(MBM, idents = c(1))
act_cd4_cells <- subset(MBM, idents = c(11))

# 设置更直观的cluster命名（可选）
t_cells$named_clusters <- factor(
  t_cells$seurat_clusters, 
  labels = c("CD8+ T cell", "CD4+ T cell", "Activated CD4+ T cell","Cycling CD8+ T cell")
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

# 分组热图（pre_post为分组变量）
cd8h <- DoHeatmap(cd8_cells,
          features = activation_genes,
          group.by = "pre_post",          # 按pre_post分组
          angle = 45,                     # 调整列标签倾斜角度
          size = 4,                       # 基因名字体大小
          label = TRUE) +
  ggtitle("CD8+ Tcell") + 
  scale_fill_gradientn(colors = c("#2E22EA","#F9F8F6","#FD1222"))

cd4h <- DoHeatmap(cd4_cells,
          features = activation_genes,
          group.by = "pre_post",          # 按pre_post分组
          angle = 45,                     # 调整列标签倾斜角度
          size = 4,                       # 基因名字体大小
          label = TRUE) +
  ggtitle("CD4+ Tcell") + 
  scale_fill_gradientn(colors = c("#2E22EA","#F9F8F6","#FD1222"))

cyccd8h <- DoHeatmap(cyc_cd8_cells,
          features = activation_genes,
          group.by = "pre_post",          # 按pre_post分组
          angle = 45,                     # 调整列标签倾斜角度
          size = 4,                       # 基因名字体大小
          label = TRUE) +
  ggtitle("Cycling CD8+ Tcell") + 
  scale_fill_gradientn(colors = c("#2E22EA","#F9F8F6","#FD1222"))

actcd4h <- DoHeatmap(act_cd4_cells,
          features = activation_genes,
          group.by = "pre_post",          # 按pre_post分组
          angle = 45,                     # 调整列标签倾斜角度
          size = 4,                       # 基因名字体大小
          label = TRUE) +
  ggtitle("Actived CD4+ Tcell") + 
  scale_fill_gradientn(colors = c("#2E22EA","#F9F8F6","#FD1222"))

combined_plot <- 
  (cd8h | cd4h) /  
  (cyccd8h | actcd4h)    

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

DoHeatmap(t_cells,
          features = exhaustion_genes,
          group.by = "pre_post",
          cells = order(t_cells$named_clusters), # 按cluster排序
          size = 4,
          disp.min = -1.5,
          disp.max = 2.5) +
  scale_fill_gradientn(colors = c("#2E22EA","#F9F8F6","#FD1222"))

cd8h <- DoHeatmap(cd8_cells,
          features = exhaustion_genes,
          group.by = "pre_post",          # 按pre_post分组
          angle = 45,                     # 调整列标签倾斜角度
          size = 4,                       # 基因名字体大小
          label = TRUE) +
  ggtitle("CD8+ Tcell") + 
  scale_fill_gradientn(colors = c("#2E22EA","#F9F8F6","#FD1222"))

cd4h <- DoHeatmap(cd4_cells,
          features = exhaustion_genes,
          group.by = "pre_post",          # 按pre_post分组
          angle = 45,                     # 调整列标签倾斜角度
          size = 4,                       # 基因名字体大小
          label = TRUE) +
  ggtitle("CD4+ Tcell") + 
  scale_fill_gradientn(colors = c("#2E22EA","#F9F8F6","#FD1222"))

cyccd8h <- DoHeatmap(cyc_cd8_cells,
          features = exhaustion_genes,
          group.by = "pre_post",          # 按pre_post分组
          angle = 45,                     # 调整列标签倾斜角度
          size = 4,                       # 基因名字体大小
          label = TRUE) +
  ggtitle("Cycling CD8+ Tcell") + 
  scale_fill_gradientn(colors = c("#2E22EA","#F9F8F6","#FD1222"))

actcd4h <- DoHeatmap(act_cd4_cells,
          features = exhaustion_genes,
          group.by = "pre_post",          # 按pre_post分组
          angle = 45,                     # 调整列标签倾斜角度
          size = 4,                       # 基因名字体大小
          label = TRUE) +
  ggtitle("Actived CD4+ Tcell") + 
  scale_fill_gradientn(colors = c("#2E22EA","#F9F8F6","#FD1222"))

combined_plot <- 
  (cd8h | cd4h) /  
  (cyccd8h | actcd4h)  



# UMAP note immunity
# 提取前7类中心坐标
centroid_df <- FetchData(MBM, vars = c("umap_1", "umap_2", group[2])) |>
  group_by(across(all_of(group[2]))) |>
  summarise(
    x = median(umap_1)+ 2.5,
    y = median(umap_2) + 3,  
    .groups = "drop"
  ) |>
  slice_head(n = 8)

p <- DimPlot(
  object = MBM,
  
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
    nrow = 11,  
    title.position = "top"  
  ))





# 对细胞占比进行分析：MBM
# 定义核心分析组别：过滤无效数据（保留所有转移瘤）
analysis_groups <- MBM@meta.data %>%
  filter(
    Primary_or_met == "Metastatic"  # 显式过滤确保数据纯净
  ) %>%
  mutate(
    group = case_when(
      pre_post == "Pre" ~ "Untreated",  # Pre组
      TRUE ~ "Treated"                  # 其他自动归为Post组
    )
  )

# 生成统计表（频数+百分比）
celltype_stats <- analysis_groups %>%
  group_by(group, cell_type) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(percentage = round(n / sum(n) * 100, 1)) %>%
  pivot_wider(
    names_from = group, 
    values_from = c(n, percentage),
    names_sort = TRUE  # 按分组名称排序
  ) 

# 输出表格（保持格式一致性）
knitr::kable(celltype_stats, caption = "Cell Type Distribution in Metastatic Groups")

# 数据准备：合并颜色映射（保持原色方案）
plot_data <- analysis_groups %>%
  count(group, cell_type) %>%
  group_by(group) %>%
  mutate(prop = n / sum(n)) %>%
  left_join(color_mapping, by = "cell_type")

# 设置分组的显式排序（重要：保证图形分组顺序）
plot_data <- plot_data %>%
  mutate(
    group = factor(
      group,
      levels = c("Untreated", "Treated")  # 新排序
    ),
    cell_type = factor(cell_type, levels = ordered_celltypes)  # 保持原有类型顺序
  )

# 调整后的可视化（简化分组后的版式优化）
ggplot(plot_data, 
       aes(x = group, 
           y = prop, 
           fill = cell_type)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  geom_text(
    aes(label = ifelse(prop * 100 >= 1, paste0(round(prop*100, 1), "%"), "")),
    position = position_stack(vjust = 0.5, reverse = TRUE),
    color = "white", 
    size = 3
  ) +
  scale_fill_manual(
    values = setNames(color_mapping$color, color_mapping$cell_type),
    breaks = ordered_celltypes,
    guide = guide_legend(reverse = TRUE)
  ) +
  labs(x = "Treatment Status", y = "Proportion") +  # 调整坐标轴标签
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  # 强化分组标签
    legend.position = "right"
  )



ggplot(plot_data, 
       aes(x = group, 
           y = prop, 
           fill = cell_type)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  geom_text(
    aes(label = ifelse(prop * 100 >= 1, paste0(round(prop*100, 1), "%"), "")),
    position = position_stack(vjust = 0.5, reverse = TRUE),
    color = "white", 
    size = 3
  ) +
  scale_fill_manual(
    values = setNames(color_mapping$color, color_mapping$cell_type),
    breaks = ordered_celltypes,
    guide = guide_legend(  # 关键修改部分
      reverse = TRUE,
      ncol = 1,            # 强制单列排列
      title = NULL         # 移除图例标题
    )
  ) +
  labs(x = "Treatment Status", y = "Proportion") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "right",
    legend.text = element_text(size = 8),  # 调整图例字体大小
    legend.key.size = unit(0.4, "cm"),        # 缩小图例色块大小
    legend.spacing.y = unit(0.1, "cm")        # 减小条目间距
  )




# 对比第1簇和第11簇(CD4+和活性CD4+)
# 提取簇1和11的细胞子集（需确认ID是否正确）
cd4_subset <- subset(MBM, 
                     subset = seurat_clusters %in% c(1, 11))

# 标注分组信息（将簇编号转化为细胞状态）
cd4_subset$cell_state <- ifelse(
  cd4_subset$seurat_clusters == 1, 
  "CD4+ T cell", 
  "Activated CD4+ T cell"
)

# 重置默认分组变量
Idents(cd4_subset) <- "cell_state"


# 绘制小提琴图
ISG15_vil_p <- VlnPlot(cd4_subset, 
             features = c("ISG15"), 
             pt.size = 0,    # 隐藏数据点
             ncol = 2,       # 两列布局
             split.by = NULL) + 
  scale_fill_manual(values = c("#1E90FF", "#4682B4")) +  # 自定义颜色
  theme_classic() +
  labs(title = "ISG15", 
       x = "", 
       y = "Expression Level") +
  theme(
    axis.text.x = element_blank()
  )


RSAD2_vil_p <- VlnPlot(cd4_subset, 
                 features = c("RSAD2"), 
                 pt.size = 0,    # 隐藏数据点
                 ncol = 2,       # 两列布局
                 split.by = NULL) + 
  scale_fill_manual(values = c("#1E90FF", "#4682B4")) +  # 自定义颜色
  theme_classic() +
  labs(title = "RSAD2", 
       x = "", 
       y = "Expression Level") +
  theme(
    axis.text.x = element_blank()
  )

IFIT1_vil_p <- VlnPlot(cd4_subset, 
                         features = c("IFIT1"), 
                         pt.size = 0,    # 隐藏数据点
                         ncol = 2,       # 两列布局
                         split.by = NULL) + 
  scale_fill_manual(values = c("#1E90FF", "#4682B4")) +  # 自定义颜色
  theme_classic() +
  labs(title = "IFIT1", 
       x = "", 
       y = "Expression Level")+
  theme(
    axis.text.x = element_blank()
  )

OAS1_vil_p <- VlnPlot(cd4_subset, 
                       features = c("OAS1"), 
                       pt.size = 0,    # 隐藏数据点
                       ncol = 2,       # 两列布局
                       split.by = NULL) + 
  scale_fill_manual(values = c("#1E90FF", "#4682B4")) +  # 自定义颜色
  theme_classic() +
  labs(title = "OAS1", 
       x = "", 
       y = "Expression Level")+
  theme(
    axis.text.x = element_blank()
  )

combined_plot <- 
  (ISG15_vil_p | RSAD2_vil_p) /  
  (IFIT1_vil_p | OAS1_vil_p)     




