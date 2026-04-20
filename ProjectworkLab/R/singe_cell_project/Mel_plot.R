#Import Package####
library(ggplot2)
library(ggrepel)
library(dplyr)
library(patchwork)
library(Seurat)
library(harmony)
library(readr)
library(devEMF)
library(tidyr)
library(clusterProfiler)
library(org.Hs.eg.db)
source("./Script/My_Function.R")

#Import Data####
Mel <- readRDS("./data/All_Celltypes_seurat_Arnon_data_metadata_M_mel_subset_revised.RDS")
Mel_cell_type <- read_csv("./data/Mel_cluster_cell.csv")
Mel_GO_activeT <- readRDS("./Mel_Cluster_Enrichment_Results/cluster7_results.rds")
Mel_GO_CD4T <- readRDS("./Mel_Cluster_Enrichment_Results/cluster0_results.rds")
Mel_GO_CD8T <- readRDS("./Mel_Cluster_Enrichment_Results/cluster1_results.rds")

#Save/Load Data####
save.image("./Mel_plot.RData")
load("./Mel_plot.RData")

write.csv(Mel_cell_type, file = "./data/Mel_cluster_cell.csv", row.names = FALSE)

#Statistic####
count_table <- Mel@meta.data %>% 
  group_by(Primary_or_met, pre_post) %>% 
  summarise(Cell_Count = n(), .groups = "drop") %>% 
  mutate(
    Percentage = round(Cell_Count / sum(Cell_Count) * 100, 1),
    .by = Primary_or_met
  )

# 对细胞占比进行分析
# 定义核心分析组别：过滤无效数据
analysis_groups <- Mel@meta.data %>%
  filter(
    (Primary_or_met == "Primary" & pre_post == "Pre") |  # 原发未治组
      (Primary_or_met == "Metastatic")                     # 转移瘤（含pre/post）
  ) %>%
  mutate(
    group = case_when(
      Primary_or_met == "Primary" ~ "Primary (Untreated)",
      Primary_or_met == "Metastatic" & pre_post == "Pre" ~ "Metastatic (Untreated)",
      TRUE ~ "Metastatic (Treated)"
    )
  )

# 生成统计表（频数+百分比）
celltype_stats <- analysis_groups %>%
  group_by(group, cell_type) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(percentage = round(n / sum(n) * 100, 1)) %>%
  pivot_wider(names_from = group, values_from = c(n, percentage)) 

# 输出可发表的格式化表格（R Markdown适用）
knitr::kable(celltype_stats, caption = "Cell Type Distribution Across Clinical Groups")

# 数据准备：合并颜色映射
plot_data <- analysis_groups %>%
  count(group, cell_type) %>%
  group_by(group) %>%
  mutate(prop = n / sum(n)) %>%
  left_join(color_mapping, by = "cell_type")  # 合并颜色信息

# 设置细胞类型的显示顺序
plot_data <- plot_data %>%
  mutate(
    cell_type = factor(cell_type, levels = ordered_celltypes) 
  )

# 调整分组顺序的因子水平
plot_data <- plot_data %>%
  mutate(
    group = factor(
      group,
      levels = c("Primary (Untreated)", 
                 "Metastatic (Untreated)",
                 "Metastatic (Treated)")
    )
  )

# 调整后的可视化代码
ggplot(plot_data, 
       aes(x = group, 
           y = prop, 
           fill = cell_type)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  geom_text(
    aes(label = ifelse(prop * 100 >= 1, paste0(round(prop*100, 1), "%"), "")),  # 只显示≥1%的标签
    position = position_stack(vjust = 0.5, reverse = TRUE),
    color = "white", 
    size = 3
  ) +
  scale_fill_manual(
    values = setNames(color_mapping$color, color_mapping$cell_type),
    breaks = ordered_celltypes,
    guide = guide_legend(reverse = TRUE)
  ) +
  labs(x = "Clinical Group", y = "Proportion", title = "") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  )

# 统计分析：多组卡方检验
cont_table <- table(analysis_groups$group, analysis_groups$cell_type)
chi_test <- chisq.test(cont_table)

# 生成组间差异矩阵（标准化残差）
adjusted_residuals <- chi_test$stdres %>% 
  as.data.frame() %>% 
  mutate(significance = ifelse(abs(Freq) > 2.58, "p<0.01", 
                               ifelse(abs(Freq) > 1.96, "p<0.05", "")))



#Clustering with pre-set resolutions####
Mel <- FindNeighbors(object = Mel, dims = 1:15)
Mel <- FindClusters(object = Mel, resolution = 0.4, cluster.name = "TEST_res_0.4")
Mel <- RunUMAP(Mel,
                    dims = 1:15,
                    min.dist = 0.5,
                    spread = 1.5)

Mel_cluster_markers <- FindAllMarkers(Mel,
                                      only.pos = TRUE,
                                      min.pct = 0.25,
                                      logfc.threshold = 0.25)

#Cell type annotation
cell_type_dict <- setNames(Mel_cell_type$`cell type`, Mel_cell_type$cluster)
Mel@meta.data$cell_type <- cell_type_dict[as.character(Mel@meta.data$seurat_clusters)]

# UMAP Plot
# Define the target clustering order
target_order <- c(0,1,7,5,2,15,4,6,10,3,8,9,11,12,13,14)

# ---- Step 1: Create sorting mapping rules ----
# Extract the correspondence between the current cluster and the cell type
cluster_celltype <- unique(Mel@meta.data[, c("seurat_clusters", "cell_type")])
# Sort by target order
ordered_celltypes <- cluster_celltype[
  order(match(cluster_celltype$seurat_clusters, target_order)),
  "cell_type"
] %>% as.character()

# ---- Step 2: Adjust metadata order ----
Mel@meta.data$cell_type <- factor(
  Mel@meta.data$cell_type,
  levels = ordered_celltypes  
)

# ---- Step 3: Create an ordered colour vector ----
# Original colour vector
original_colors <- c("#4169E1", "#1E90FF", "#000080", "#32CD32", "#FFD700",
                     "#FF8C00", "#8B4513", "#6B8E23", "#BA55D3", "#A52A2A",
                     "#8B0000", "#4B0082", "#800080", "#2F4F4F", "#B22222",
                     "#696969")

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

color_mapping
# ---- Step 5: Drawing ----
group <- c("pre_post", "cell_type", "Primary_or_met")
color_palette <- c("#6A3D9A", "#33A02C")

p <- DimPlot(
  object = Mel,

    group.by = group[1],
  #cols = color_palette,     
  #cols = custom_ordered_colors,  # Choose to use the costumed colors
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
    nrow = 8,  
    title.position = "top"  
  ))

print(p_legend_bottom)
ggsave("UMAP_custom_colors.tiff", width = 10, height = 6, dpi = 300)

svg("UMAP_Mel_cell_type_Legend_bottom.svg", width = 10, height = 10)
print(p_legend_bottom)
dev.off()

# Dot Plot####
Mel$TEST_res_0.4 <- factor(
  Mel$TEST_res_0.4,
  levels = c(0,1,7,5,2,15,4,6,10,3,8,9,11,12,13,14) 
)

dot_p <- DotPlot(Mel,
                 features = marker_db$Macrophages,
                 group.by = "TEST_res_0.4") +
  RotatedAxis() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "italic"  
    )
  )

emf("Mel_Macro_dot.emf")
print(dot_p)
dev.off()

# GO
GO <- dotplot(Mel_GO_activeT$ego, 
        showCategory = 20,          # 显示前20个条目
        color = "p.adjust",         # 颜色映射调整p值
        label_format = 60) +        # 条目名称换行长度
  ggtitle("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # x轴标签倾斜45度

emf("GO_CD8cyc.emf", width = 10, height = 10)
print(GO)
dev.off()


# 数据准备：按细胞类型分组计算各组别占比
celltype_group_data <- analysis_groups %>% 
  count(cell_type, group) %>%               # 按细胞类型和组别计数
  group_by(cell_type) %>%                   # 以细胞类型为分组依据
  mutate(prop = n / sum(n)) %>%             # 计算各组别在细胞类型内的比例
  left_join(
    color_mapping, 
    by = "cell_type"
  ) %>%                                     # 合并颜色信息（可选）
  # 设置因子顺序以控制堆叠顺序
  mutate(
    group = factor(group, 
                   levels = c("Primary (Untreated)", 
                              "Metastatic (Untreated)", 
                              "Metastatic (Treated)"))
  )

# 方式1：堆积柱状图（所有细胞类型堆叠）
ggplot(celltype_group_data, 
       aes(x = cell_type, y = prop, fill = group)) +  # x轴为细胞类型，填充为组别
  geom_col(width = 0.7) +                            # 柱状图
  geom_text(
    aes(label = ifelse(prop >= 0.05, paste0(round(prop * 100, 1), "%"), "")),  # 仅显示≥5%的标签
    position = position_stack(vjust = 0.5),
    color = "white",
    size = 3
  ) +
  scale_fill_brewer(palette = "Set2") +               # 使用颜色方案（可替换为具体颜色值）
  labs(x = "Cell Type", y = "Proportion", fill = "Clinical Group") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 确保细胞类型按自定义顺序排序
celltype_group_data <- celltype_group_data %>% 
  mutate(
    cell_type = factor(cell_type, levels = ordered_celltypes)  # 强制排序
  )

# 生成水平堆积柱状图（翻转坐标系）
ggplot(celltype_group_data, 
       aes(y = cell_type, x = prop, fill = group)) +   # y轴为细胞类型（隐藏标签）
  geom_col(width = 0.7) +
  geom_text(
    aes(label = ifelse(prop >= 0.05, paste0(round(prop * 100, 1), "%"), "")),
    position = position_stack(vjust = 0.5),
    color = "white",
    size = 3
  ) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Proportion", y = NULL, fill = "Clinical Group") +   # 隐藏y轴标题
  theme_classic(base_size = 14) +
  theme(
    #axis.text.y = element_blank(),             # 隐藏y轴细胞类型名称
    #axis.ticks.y = element_blank(),            # 隐藏y轴刻度线
    legend.position = "right"                 # 图例置于底部
  )            


# 提取前6类中心坐标####
centroid_df <- FetchData(Mel, vars = c("umap_1", "umap_2", group[2])) |>
  group_by(across(all_of(group[2]))) |>
  summarise(
    x = median(umap_1)+ 2.5,
    y = median(umap_2) + 3,  
    .groups = "drop"
  ) |>
  slice_head(n = 6)

p <- DimPlot(
  object = Mel,
  
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
    nrow = 8,  
    title.position = "top"  
  ))
