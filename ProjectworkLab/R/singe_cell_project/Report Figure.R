# Draw graphs used in reports####
#Import Package
library(ggplot2)
library(reshape2)
library(ggrepel)
library(dplyr)
library(tidyverse)
library(patchwork)
library(Seurat)
library(harmony)
library(readr)
library(devEMF)
library(tidyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(speckle)
library(limma)
library(extrafont)
loadfonts(device = "win")

# Import Font
font_import()                   

# Load Data####
MBM <- readRDS("./data/MBM.RDS")
Mel <- readRDS("./data/Mel.RDS")

load("./Mel_plot.RData")
save.image("./Mel_plot.RData")

load("./Mel_meta_plot.RData")
save.image("./Mel_meta_plot.RData")

load("./MBM_plot.RData")
save.image("./MBM_plot.RData")

load("./Harmony.RData")
save.image("./Harmony.RData")

# Statistic----
# Extract metadata from Seurat objects and de-duplicate by patient[MBM/Mel]
patient_metadata <- MBM@meta.data %>% 
  as.data.frame() %>%             
  distinct(donor_id, .keep_all = TRUE) %>% 
  dplyr::select(Gender, pre_post) 

# Age statistical analysis
age_summary <- list(
  mean_age = mean(patient_metadata$age, na.rm = TRUE),
  median_age = median(patient_metadata$age, na.rm = TRUE)
)

# Sex statistical analysis
sex_distribution <- as.data.frame(table(patient_metadata$Gender)) %>%
  dplyr::rename(Gender = 1, Count = Freq)

# Treatment group statistics
pre_post_distribution <- as.data.frame(table(patient_metadata$pre_post)) %>%
  dplyr::rename(Gender = 1, Count = Freq)  

# P-value
# Chi-square test of the number of people in the treatment groups
pre_post_p <- chisq.test(pre_post_distribution$Count)$p.value

# Tests of differences between age groups
age_test <- tryCatch({
  t.test(age ~ pre_post, data = patient_metadata)
}, error = function(e) {
  wilcox.test(age ~ pre_post, data = patient_metadata)
})

# Fisher's exact test for sex distribution and treatment group
gender_group_table <- table(patient_metadata$sex, patient_metadata$pre_post)
gender_group_p <- if (all(dim(gender_group_table) >= 2)) {
  fisher.test(gender_group_table, simulate.p.value = TRUE)$p.value
} else NA

# Result Output
cat("============= Basic statistic =============\n")
cat("Age statistic (all patients)：\n")
print(age_summary)

cat("\n性别分布：\n")
print(sex_distribution)

cat("\n治疗分组分布：\n")
print(pre_post_distribution)

cat("\n============= Statistical tests =============\n")
cat("1. Treatment Group Proportion Test (χ²):\n",
    "P-value =", format.pval(pre_post_p, digits = 3), "\n\n")

cat("2. Test of differences between age groups：\n",
    "Method:", if (exists("statistic", age_test)) "Student's t-test" else "Wilcoxon rank sum test", "\n",
    "P-value =", format.pval(age_test$p.value, digits = 3), "\n\n")

cat("3. 性别-治疗分组独立性检验 (Fisher's exact):\n",
    "P-value =", if (!is.na(gender_group_p)) format.pval(gender_group_p, digits = 3) else "Not applicable", "\n")

# Figure3 Global Setting----
# Original colour vector
original_colors <- c("#4169E1", "#1E90FF", "#000080", "#32CD32", "#FFD700",
                     "#FF8C00", "#8B4513", "#6B8E23", "#BA55D3", "#A52A2A",
                     "#CD853F", "#4B0082", "#800080", "#2F4F4F", "#B22222",
                     "#696969")

# Create a new mapping data frame to ensure the correct colour correspondence
color_mapping <- data.frame(
  cell_type = ordered_celltypes,
  color = original_colors[1:length(target_order)] 
)
custom_ordered_colors <- as.character(color_mapping$color)

# Check colour
print(data.frame(
  Cluster = target_order,
  CellType = ordered_celltypes,
  AssignedColor = custom_ordered_colors
))

# Define the target clustering order
target_order <- c(0,1,7,5,2,15,4,6,10,3,8,9,11,12,13,14)


# Figure3A----
# Colour: Metastatic #6A3D9A; Primary #33A02C
figure3a_p <- DimPlot(
  object = Mel,
  group.by = "Primary_or_met",
  cols = c("#6A3D9A", "#33A02C"),
  label = FALSE,
  pt.size = 0.5,
  reduction = "umap"
) + 
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle("") +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_blank(),   
    plot.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  )

# Figure3B----
# Colour: Treated（Post）#CA0020 Untreated（Pre）#0571B0
figure3b_p <- DimPlot(
  object = Mel,
  group.by = "pre_post",
  cols = c("#CA0020", "#0571B0"),
  label = FALSE,
  pt.size = 0.5,
  reduction = "umap"
) + 
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle("") +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_blank(),   
    plot.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  )

# Figure3C----
figure3c_p <- DimPlot(
  object = Mel,
  group.by = "cell_type",
  cols = custom_ordered_colors, # For manually cell type
  label = FALSE,
  pt.size = 0.5,
  reduction = "umap"
) + 
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle("") +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_blank(),   
    plot.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  )

# Figure3D----
Mel$TEST_res_0.4 <- factor(
  Mel$TEST_res_0.4,
  levels = c(0,1,7,5,2,15,4,6,10,3,8,9,11,12,13,14) 
)

marker_gene = c("CD3D", "IL7R", "CD4", "CD8A",        # T cells
                "ISG15", "RSAD2", "IFIT1", "OAS1",    # Interferon-stimulated genes (ISGs)
                "CD19", "MS4A1", "CD79A",             # B cells
                "LYZ", "CD14", "FCGR3A",              # Macrophages
                "GNLY","NKG7",                        # NK cells
                "ACTA2","CTSK","CCND1","BCAN",        # Tumour-associated fibroblasts
                "TYR","MITF", "SOX10",                # Melanocytes
                "BRAF","CDKN2A","TP53")               # Malignant cells

figure3d_p <- DotPlot(Mel, 
                      features = marker_gene, 
                      group.by = "TEST_res_0.4") + 
  theme(
    axis.text.x = element_text(            # X-axis adjustment
      angle = 45,
      hjust = 1,
      face = "italic",
      family = "Arial",
      size = 8               # word size
    ),
    legend.position = "none"
  ) + 
  labs(x = NULL)

# Figure3E----
# Data preprocessing and factor order correction
meta_data <- Mel@meta.data[, c("Primary_or_met", "cell_type")]

# 计算各组的细胞类型比例（注意保留原始顺序）
cell_table <- table(meta_data$Primary_or_met, meta_data$cell_type)
cell_prop <- as.data.frame(prop.table(cell_table, margin = 1)) 
colnames(cell_prop) <- c("Primary_met", "cell_type", "Proportion")

# Add chi-square test code section
chisq_test_result <- chisq.test(cell_table)
fisher_test_result <- fisher.test(cell_table, simulate.p.value = TRUE) 

# Generate label column
cell_prop_labeled <- cell_prop %>%
  mutate(
    label = ifelse(Proportion >= 0.01, 
                   sprintf("%.1f%%", Proportion* 100), 
                   "")
  )

# Create a color vector
color_vec <- setNames(color_mapping$color, color_mapping$cell_type)

# Draw a stacked column chart
figure3e_p <- ggplot(cell_prop_labeled, 
                     aes(x = factor(Primary_met,levels = c("Primary", "Metastatic")),
                         y = Proportion, fill = cell_type)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  geom_text(aes(label = label), 
            position = position_stack(vjust = 0.5, reverse = TRUE),
            size = 3, color = "white") +
  scale_fill_manual(
    values = color_vec,
    guide = guide_legend(  
      keywidth = unit(3, "mm"), 
      keyheight = unit(3, "mm"), 
      override.aes = list(size = 0.5)  
    )
  ) +
  labs(x = "", y = "Proportion") +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 6, margin = margin(r = 6)), 
    legend.spacing.x = unit(0.5, "mm"),
    legend.key.spacing = unit(0.5, "mm"),
    legend.margin = margin(l = -5)
  )

# 添加显著性注释到图中
figure3e_p <- figure3e_p +
  annotate("text",
           x = 1.5,  # 两组中间位置
           y = 1.05, # y轴上方区域
           label = ifelse(chisq_test_result$p.value < 0.001,
                          "***",
                          sprintf("p = %.3f", chisq_test_result$p.value)),
           color = "black",
           size = 4)

# Figure3----
figure3_p <- 
  (figure3a_p | figure3b_p) /  
  (figure3c_p | figure3d_p) /
  (figure3e_p | plot_spacer())


# 设置A4纸张尺寸参数（竖版：21x29.7cm）
a4_width <- 21   # 厘米（短边）
a4_height <- 27.7 # 厘米（长边）
dpi_level <- 600  # 打印级分辨率

# 推荐矢量格式（PDF）
ggsave("./Report Plot/figure3_p.pdf", 
       plot = figure3_p,
       width = a4_width,
       height = a4_height,
       units = "cm",
       dpi = dpi_level,
       device = cairo_pdf, # 支持字体嵌入
       bg = "white")       # 透明背景可改为"transparent"

# 如需位图格式（TIFF无损压缩）
ggsave("./Report Plot/figure3_p.tiff",
       plot = figure3_p,
       width = a4_width, 
       height = a4_height,
       units = "cm",
       dpi = dpi_level,
       compression = "lzw", # 无损压缩
       bg = "white")

#Figure4A----
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

print(data.frame(
  Cluster = target_order,
  CellType = ordered_celltypes,
  AssignedColor = custom_ordered_colors
))

figure4a_p <- DimPlot(
  object = Mel_metastases,
  group.by = "cell_type",
  cols = custom_ordered_colors, # For manually cell type
  label = FALSE,
  pt.size = 0.5,
  reduction = "umap"
) + 
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle("") +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_blank(),   
    plot.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  )

#Figure4B----
marker_gene = c("CD3D", "IL7R", "CD4", "CD8A",        # T cells
                "ISG15", "RSAD2", "IFIT1", "OAS1",    # Interferon-stimulated genes (ISGs)
                "CD19", "MS4A1", "CD79A",             # B cells
                "LYZ", "CD14", "FCGR3A",              # Macrophages
                "GNLY","NKG7",                        # NK cells
                "ACTA2","CTSK","CCND1","BCAN",        # Tumour-associated fibroblasts
                "TYR","MITF", "SOX10",                # Melanocytes
                "BRAF","CDKN2A","TP53")               # Malignant cells

Mel_metastases$TEST_res_0.3 <- factor(
  Mel_metastases$TEST_res_0.3,
  levels = c("0","1","6","2","16","5","12","3","4","7","8","9","10","11","13","14","15") 
)

figure4b_p <- DotPlot(Mel_metastases, 
                      features = marker_gene, 
                      group.by = "TEST_res_0.3") + 
  theme(
    axis.text.x = element_text(           
      angle = 45,
      hjust = 1,
      face = "italic",
      family = "Arial",
      size = 8               
    ),
    legend.position = "none"  
  ) + 
  labs(x = NULL)

#Figure4C----
figure4c_p <- DimPlot(
  object = Mel_metastases,
  group.by = "pre_post",
  cols = c("#0571B0", "#CA0020"),
  label = FALSE,
  pt.size = 0.5,
  reduction = "umap"
) + 
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle("") +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_blank(),   
    plot.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  )

#Figure4D----
# 转换因子顺序（核心步骤）
Mel_metastases@meta.data <- Mel_metastases@meta.data %>% 
  mutate(
    pre_post = recode(pre_post, "Pre" = "Untreated", "Post" = "Treated"),
    pre_post = factor(pre_post, levels = c("Untreated", "Treated")),  # 确定顺序
    cell_type = factor(cell_type)  # 保持原顺序或按需排序
  )

# 计算分组占比
# 修改后的数据处理与绘图代码
# 修改后的数据处理与绘图代码
prop_data <- Mel_metastases@meta.data %>%
  group_by(cell_type, pre_post) %>%
  tally() %>%
  group_by(cell_type) %>%
  mutate(
    Proportion = n/sum(n)*100,
    # 动态生成标签（过滤小于25%的数值）
    label = if_else(Proportion >= 25, 
                    sprintf("%.1f%%", Proportion), 
                    "")
  ) %>%
  ungroup()

figure4d_p <- ggplot(prop_data, 
                     aes(y = cell_type, x = Proportion, fill = pre_post)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  # 添加智能比例标注
  geom_text(
    aes(label = label), 
    position = position_stack(vjust = 0.5, reverse = TRUE),
    size = 1.8, 
    color = "white",  # 根据背景色优化
    fontface = "bold"
  ) + 
  scale_fill_manual(values = c("Untreated" = "#0571B0", "Treated" = "#CA0020")) +
  labs(x = "Proportion (%)", y = "") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.y = element_blank(),
    legend.position = "none"
  )

#Figure4E----
# 直接统计每个治疗组的细胞数
Mel_meta_treatment_table <- as.data.frame(table(Mel_metastases@meta.data$ICB_exposed))

# 优化列名和格式
colnames(Mel_meta_treatment_table) <- c("Treatment", "Cell_Count")
Mel_meta_treatment_table$Proportion <- prop.table(Mel_meta_treatment_table$Cell_Count)

# 添加百分比显示列（保留两位小数）
Mel_meta_treatment_table$Percentage <- round(Mel_meta_treatment_table$Proportion * 100, 2)

# 按细胞数降序排列
Mel_meta_treatment_table <- Mel_meta_treatment_table[order(-Mel_meta_treatment_table$Cell_Count), ]

# 移除中间计算列
Mel_meta_treatment_table$Proportion <- NULL 

# 查看结果
print(Mel_meta_treatment_table)
figure4e <- Mel_meta_treatment_table
write.csv(figure4e,"./Report Plot/figure4e.csv")


#Figure4F----
# 数据预处理及因子顺序修正 
meta_data <- Mel_metastases@meta.data[, c("pre_post", "cell_type")]

# 计算各组的细胞类型比例（注意保留原始顺序）
cell_table <- table(meta_data$pre_post, meta_data$cell_type)
cell_prop <- as.data.frame(prop.table(cell_table, margin = 1)) 
colnames(cell_prop) <- c("pre_post", "cell_type", "Proportion")

# 增加卡方检验代码部分（基于原始频数数据cell_table）
chisq_test_result <- chisq.test(cell_table)
fisher_test_result <- fisher.test(cell_table, simulate.p.value = TRUE) 

# 生成标签列（忽略低于3%的标签）
cell_prop_labeled <- cell_prop %>%
  mutate(
    label = ifelse(Proportion >= 0.03, 
                   sprintf("%.1f%%", Proportion* 100), 
                   "")
  )

cell_prop_labeled <- cell_prop_labeled %>%
  mutate(pre_post = recode(pre_post, "Pre" = "Untreated", "Post" = "Treated"))

# 创建颜色向量（确保颜色映射正确）
color_vec <- setNames(color_mapping$color, color_mapping$cell_type)

# 绘制堆积柱状图
figure4f_p <- ggplot(cell_prop_labeled, 
                     aes(x = factor(pre_post, levels = c("Untreated", "Treated")),
                         y = Proportion, fill = cell_type)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  geom_text(aes(label = label), 
            position = position_stack(vjust = 0.5, reverse = TRUE),
            size = 3, color = "white") +
  scale_fill_manual(
    values = color_vec,
    guide = guide_legend(  
      keywidth = unit(3, "mm"), 
      keyheight = unit(3, "mm"), 
      override.aes = list(size = 0.5)  
    )
  ) +
  labs(x = "", y = "Proportion") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 7),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 6, margin = margin(r = 6)), 
    legend.spacing.x = unit(0.5, "mm"),
    legend.key.spacing = unit(0.5, "mm"),
    legend.margin = margin(l = -5)
  )

# 添加显著性注释到图中
figure4f_p <- figure4f_p +
  annotate("text",
           x = 1.5,  # 两组中间位置
           y = 1.05, # y轴上方区域
           label = ifelse(chisq_test_result$p.value < 0.001,
                          "***",
                          sprintf("p = %.3f", chisq_test_result$p.value)),
           color = "black",
           size = 4)

#Figure4----
figure4_p <- 
  (figure4a_p | figure4b_p) /  
  (figure4c_p | figure4d_p) /
  (plot_spacer() | figure4f_p)

# 设置A4纸张尺寸参数（竖版：21x29.7cm）
a4_width <- 21   # 厘米（短边）
a4_height <- 27.7 # 厘米（长边）
dpi_level <- 600  # 打印级分辨率

# 推荐矢量格式（PDF）
ggsave("./Report Plot/figure4_p.pdf", 
       plot = figure4_p,
       width = a4_width,
       height = a4_height,
       units = "cm",
       dpi = dpi_level,
       device = cairo_pdf, # 支持字体嵌入
       bg = "white")       # 透明背景可改为"transparent"

# 子图分别保存
ggsave("./Report Plot/figure4/figure4f_p.tiff",
       plot = figure4f_p,
       width = a4_width/2, 
       height = a4_height/3,
       units = "cm",
       dpi = dpi_level,
       compression = "lzw", # 无损压缩
       bg = "white")

# 如需位图格式（TIFF无损压缩）
ggsave("./Report Plot/figure4_p.tiff",
       plot = figure4_p,
       width = a4_width, 
       height = a4_height,
       units = "cm",
       dpi = dpi_level,
       compression = "lzw", # 无损压缩
       bg = "white")

#Figure5A----
# 提取目标亚群
t_cells <- subset(Mel_metastases, idents = c(0,1,6))

# 设置直观Cluster命名
t_cells$named_clusters <- factor(
  t_cells$seurat_clusters, 
  labels = c("CD8+ T cell", "CD4+ T cell", "Cycling CD8+ T cell")
)

group_colors <- c(
  "CD8+ T cell" = "#4169E1",  
  "CD4+ T cell" = "#1E90FF",  
  "Cycling CD8+ T cell" = "#000080" 
)

# activation
activation_genes <- c("CD69", "IL2RA", "MKI67", "IFNG", "GZMB", "PRF1")
figure5a_p <- DoHeatmap(t_cells,
                              features = activation_genes,
                              group.by = "named_clusters",
                              group.colors = group_colors,
                              slot = "scale.data",    # 使用标准化数据
                              disp.min = -2,          # 截断Z-score范围
                              disp.max = 2,
                              angle = 45,             # 调整标签角度
                              size = 4) + 
  scale_fill_gradientn(colors = c("#2E22EA", "#F9F8F6", "#FD1222"))

#Figure5B----
# exhaustion
exhaustion_genes <- c("PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT")
figure5b_p <- DoHeatmap(t_cells,
          features = exhaustion_genes,
          group.by = "named_clusters",
          group.colors = group_colors,
          cells = order(t_cells$named_clusters), # 按cluster排序
          size = 4,
          disp.min = -1.5,
          disp.max = 2.5) +
  scale_fill_gradientn(colors = c("#2E22EA", "#F9F8F6", "#FD1222"))

#Figure5C&D&E----
group_colors_treat <- c(
  "Untreated" = "#0571B0",  
  "Treated" = "#CA0020"
)

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
                  group.colors = group_colors_treat,
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
                  group.colors = group_colors_treat,
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

figure5c_p <- heatmap_list[[1]]
figure5d_p <- heatmap_list[[2]]
figure5e_p <- heatmap_list[[3]]

#Figure5F----
# 设置分组比较
Idents(t_cells) <- "pre_post" # 设置分组列
markers <- FindMarkers(t_cells,
                       ident.1 = "Treated",  
                       ident.2 = "Untreated",    
                       min.pct = 0.25,          
                       logfc.threshold = 0.25,   
                       test.use = "wilcox")      

# 添加基因列并计算转换值
volcano_data <- markers %>%
  tibble::rownames_to_column("gene") %>%
  mutate(
    log10p = -log10(p_val_adj),
    regulation = case_when(
      avg_log2FC > 0 & p_val_adj < 0.05 ~ "Up",
      avg_log2FC < 0 & p_val_adj < 0.05 ~ "Down",
      TRUE ~ "Not sig"
    )
  )

# 对应调整基因筛选部分
top_up <- volcano_data %>% 
  filter(regulation == "Up") %>% 
  arrange(p_val_adj) %>% 
  head(5)  

top_down <- volcano_data %>% 
  filter(regulation == "Down") %>% 
  arrange(p_val_adj) %>% 
  head(5)  

# 提取全部上调下调基因
all_up <- volcano_data %>% 
  filter(regulation == "Up") %>% 
  arrange(p_val_adj) 

all_down <- volcano_data %>% 
  filter(regulation == "Down") %>% 
  arrange(p_val_adj)

# 修改后的绘图代码
figure5f_p <- ggplot(volcano_data, aes(x = avg_log2FC, y = log10p)) +
  geom_point(aes(color = regulation), alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Down" = "#2E22EA", "Up" = "#FD1222", "Not sig" = "grey60")) +
  
  # 新增坐标轴优化
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +  # 确保0点在x轴中心对称
  geom_vline(xintercept = 0, color = "grey40", alpha = 0.8) +    # 高亮0基准线
  
  # 调整显著性阈值线
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey40") +
  
  # 减少标签数量至上下调各5个
  geom_text_repel(
    data = rbind(
      top_up %>% head(5),  # 仅显示上调Top5
      top_down %>% head(5) # 仅显示下调Top5
    ),
    aes(label = gene),
    size = 3.5,
    box.padding = 0.3,
    max.overlaps = 20,      # 减少标签重叠容忍度
    segment.color = "grey50", # 标签引导线颜色
    nudge_x = 0.1          # 横向微调标签位置
  ) +
  
  # 优化标签和主题
  labs(
    x = expression(Log[2]*" Fold Change"),
    y = expression(-Log[10]*"(Adj.P-value)"),
    color = "Regulation"
  ) +
  theme_classic(base_size = 14) +  # 增大基础字体
  theme(
    panel.grid.major = element_line(color = "grey93"), # 更浅的网格线
    legend.position = "top",
    axis.title = element_text(face = "bold"),          # 坐标轴标题加粗
    legend.title = element_text(face = "bold")
  )

# 保存结果
write.csv(top_up, "./data/Mel_T_cell/Top5_Upregulated_Genes.csv", row.names = FALSE)
write.csv(top_down, "./data/Mel_T_cell/Top5_Downregulated_Genes.csv", row.names = FALSE)

write.csv(all_up, "./data/Mel_T_cell/Upregulated_Genes.csv", row.names = FALSE)
write.csv(all_down, "./data/Mel_T_cell/Downregulated_Genes.csv", row.names = FALSE)

#Figure5G----
# 步骤1: 准备基因列表
gene_list <- all_up$gene  # 提取差异基因symbol
gene_list_down <- all_down$gene

# 步骤2: 将symbol转换为Entrez ID（enrichGO默认需要）
entrez_ids <- bitr(gene_list_down, 
                   fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = org.Hs.eg.db) %>% 
  pull(ENTREZID)

# 步骤3: GO富集分析
go_res <- enrichGO(
  gene          = entrez_ids,
  OrgDb         = org.Hs.eg.db,        # 指定人类数据库
  keyType       = "ENTREZID",          # 输入基因类型
  ont           = "BP",                # 选择生物学过程(Biological Process)
  pAdjustMethod = "BH",                # 多重检验校正方法
  pvalueCutoff  = 0.05,                # p值阈值
  qvalueCutoff  = 0.2,                 # q值阈值
  readable      = TRUE                 # 转换EntrezID为基因symbol
)

# 步骤4: 结果可视化
# 绘制前15显著GO项的点图
figure5g_1 <- ggplot(go_res@result[1:15, ], aes(x=GeneRatio, y=reorder(Description, GeneRatio))) + 
  geom_point(aes(size=Count, color=pvalue)) + 
  theme_bw() +
  theme(axis.title.y = element_blank())

figure5g_2 <- ggplot(go_res@result[1:15, ], aes(x=GeneRatio, y=reorder(Description, GeneRatio))) + 
  geom_point(aes(size=Count, color=pvalue)) + 
  theme_bw() +
  theme(axis.title.y = element_blank())

# 导出结果到CSV
write.csv(go_res@result, "./data/Mel_T_cell/Upregulated_GO_enrichment_results.csv", row.names = FALSE)
write.csv(go_res@result, "./data/Mel_T_cell/Downregulated_GO_enrichment_results.csv", row.names = FALSE)

#Figure5----
figure5_p_1 <- 
  (figure5a_p | figure5b_p) /  
  (figure5c_p) /
  (figure5d_p) /
  (figure5e_p)

# 设置A4纸张尺寸参数（竖版：21x29.7cm）
a4_width <- 21   # 厘米（短边）
a4_height <- 27.7 # 厘米（长边）
dpi_level <- 600  # 打印级分辨率

# 推荐矢量格式（PDF）
ggsave("./Report Plot/figure5/figure5_p_2.pdf", 
       plot = figure5_p_2,
       width = a4_width,
       height = a4_height,
       units = "cm",
       dpi = dpi_level,
       device = cairo_pdf, # 支持字体嵌入
       bg = "white")       # 透明背景可改为"transparent"

# 如需位图格式（TIFF无损压缩）
ggsave("./Report Plot/figure5/figure5_p_2.tiff",
       plot = figure5_p_2,
       width = a4_width, 
       height = a4_height,
       units = "cm",
       dpi = dpi_level,
       compression = "lzw", # 无损压缩
       bg = "white")

ggsave("./Report Plot/figure5/figure5g_2.tiff",
       plot = figure5g_2,
       width = a4_width, 
       height = a4_height/3,
       units = "cm",
       dpi = dpi_level,
       compression = "lzw", # 无损压缩
       bg = "white")

#Figure6A----
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

print(data.frame(
  Cluster = target_order,
  CellType = ordered_celltypes,
  AssignedColor = custom_ordered_colors
))

figure6a_p <- DimPlot(
  object = MBM,
  group.by = "cell_type",
  cols = custom_ordered_colors, # For manually cell type
  label = FALSE,
  pt.size = 0.5,
  reduction = "umap"
) + 
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle("") +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_blank(),   
    plot.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  )

#Figure6B----
marker_gene = c("CD3D", "IL7R", "CD4", "CD8A",        # T cells
                "ISG15", "RSAD2", "IFIT1", "OAS1",    # Interferon-stimulated genes (ISGs)
                "CD19", "MS4A1", "CD79A",             # B cells
                "LYZ", "CD14", "FCGR3A",              # Macrophages
                "GNLY","NKG7",                        # NK cells
                "ACTA2","CTSK","CCND1","BCAN",        # Tumour-associated fibroblasts
                "TYR","MITF", "SOX10",                # Melanocytes
                "BRAF","CDKN2A","TP53")               # Malignant cells

MBM$TEST_res_0.4 <- factor(
  MBM$TEST_res_0.4,
  levels = c(0,1,11,13,14,20,3,10,4,19,2,5,6,7,8,9,12,15,16,17,18) 
)

figure6b_p <- DotPlot(MBM, 
                      features = marker_gene, 
                      group.by = "TEST_res_0.4") + 
  theme(
    axis.text.x = element_text(           
      angle = 45,
      hjust = 1,
      face = "italic",
      family = "Arial",
      size = 8               
    ),
    legend.position = "none"  
  ) + 
  labs(x = NULL)

#Figure6C----
MBM@meta.data$pre_post <- factor(MBM@meta.data$pre_post, levels = c("Pre", "Post"))

figure6c_p <- DimPlot(
  object = MBM,
  group.by = "pre_post",
  cols = c("#0571B0", "#CA0020"),
  label = FALSE,
  pt.size = 0.5,
  reduction = "umap"
) + 
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle("") +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_blank(),   
    plot.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  )

#Figure6D----
# 转换因子顺序（核心步骤）
MBM@meta.data <- MBM@meta.data %>% 
  mutate(
    pre_post = recode(pre_post, "Pre" = "Untreated", "Post" = "Treated"),
    pre_post = factor(pre_post, levels = c("Untreated", "Treated")),  # 确定顺序
    cell_type = factor(cell_type)  # 保持原顺序或按需排序
  )

# 计算分组占比
# 修改后的数据处理与绘图代码
# 修改后的数据处理与绘图代码
prop_data <- MBM@meta.data %>%
  group_by(cell_type, pre_post) %>%
  tally() %>%
  group_by(cell_type) %>%
  mutate(
    Proportion = n/sum(n)*100,
    # 动态生成标签（过滤小于25%的数值）
    label = if_else(Proportion >= 25, 
                    sprintf("%.1f%%", Proportion), 
                    "")
  ) %>%
  ungroup()

figure6d_p <- ggplot(prop_data, 
                     aes(y = cell_type, x = Proportion, fill = pre_post)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  # 添加智能比例标注
  geom_text(
    aes(label = label), 
    position = position_stack(vjust = 0.5, reverse = TRUE),
    size = 1.8, 
    color = "white",  # 根据背景色优化
    fontface = "bold"
  ) + 
  scale_fill_manual(values = c("Untreated" = "#0571B0", "Treated" = "#CA0020")) +
  labs(x = "Proportion (%)", y = "") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.y = element_blank(),
    legend.position = "none"
  )

#Figure6E----
# 直接统计每个治疗组的细胞数
MBM_treatment_table <- as.data.frame(table(MBM@meta.data$`Immunotherapy #1`))

# 优化列名和格式
colnames(MBM_treatment_table) <- c("Treatment", "Cell_Count")
MBM_treatment_table$Proportion <- prop.table(MBM_treatment_table$Cell_Count)

# 添加百分比显示列（保留两位小数）
MBM_treatment_table$Percentage <- round(MBM_treatment_table$Proportion * 100, 2)

# 按细胞数降序排列
MBM_treatment_table <- MBM_treatment_table[order(-MBM_treatment_table$Cell_Count), ]

# 移除中间计算列
MBM_treatment_table$Proportion <- NULL 

# 查看结果
print(MBM_treatment_table)
figure6e <- MBM_treatment_table
write.csv(figure6e,"./Report Plot/figure6/figure6e.csv")

#Figure6F----
# 数据预处理及因子顺序修正 
meta_data <- MBM@meta.data[, c("pre_post", "cell_type")]

# 计算各组的细胞类型比例（注意保留原始顺序）
cell_table <- table(meta_data$pre_post, meta_data$cell_type)
cell_prop <- as.data.frame(prop.table(cell_table, margin = 1)) 
colnames(cell_prop) <- c("pre_post", "cell_type", "Proportion")

# 增加卡方检验代码部分（基于原始频数数据cell_table）
chisq_test_result <- chisq.test(cell_table)
fisher_test_result <- fisher.test(cell_table, simulate.p.value = TRUE) 

# 生成标签列（忽略低于3%的标签）
cell_prop_labeled <- cell_prop %>%
  mutate(
    label = ifelse(Proportion >= 0.03, 
                   sprintf("%.1f%%", Proportion* 100), 
                   "")
  )

cell_prop_labeled <- cell_prop_labeled %>%
  mutate(pre_post = recode(pre_post, "Pre" = "Untreated", "Post" = "Treated"))

# 创建颜色向量（确保颜色映射正确）
color_vec <- setNames(color_mapping$color, color_mapping$cell_type)

# 绘制堆积柱状图
figure6f_p <- ggplot(cell_prop_labeled, 
                     aes(x = factor(pre_post, levels = c("Untreated", "Treated")),
                         y = Proportion, fill = cell_type)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  geom_text(aes(label = label), 
            position = position_stack(vjust = 0.5, reverse = TRUE),
            size = 3, color = "white") +
  scale_fill_manual(
    values = color_vec,
    guide = guide_legend(  
      ncol = 1,  #  强制图例单列显示
      keywidth = unit(2, "mm"),  #  缩小色块宽度
      keyheight = unit(2, "mm"), #  缩小色块高度
      override.aes = list(size = 0.5),
      byrow = TRUE  #  按行填充条目（优化紧凑性）
    )
  ) +
  labs(x = "", y = "Proportion") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 7),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 5, margin = margin(r = 2)), # 缩小文本并减少右侧边距
    legend.spacing.x = unit(0.2, "mm"),  # 压缩条目水平间距
    legend.key.spacing = unit(0.2, "mm"), #  压缩色块与文本间距
    legend.margin = margin(l = -8, t = -2)  #  减少图例左侧和顶部外边距（负值收缩空间）
  )

# 添加显著性注释到图中
figure6f_p <- figure6f_p +
  annotate("text",
           x = 1.5,  # 两组中间位置
           y = 1.05, # y轴上方区域
           label = ifelse(chisq_test_result$p.value < 0.001,
                          "***",
                          sprintf("p = %.3f", chisq_test_result$p.value)),
           color = "black",
           size = 4)

#Figure6----
figure6_p <- 
  (figure6a_p | figure6b_p) /  
  (figure6c_p | figure6d_p ) /
  (plot_spacer()  | figure6f_p)

# 设置A4纸张尺寸参数（竖版：21x29.7cm）
a4_width <- 21   # 厘米（短边）
a4_height <- 27.7 # 厘米（长边）
dpi_level <- 600  # 打印级分辨率

# 推荐矢量格式（PDF）
ggsave("./Report Plot/figure6/figure6_p.pdf", 
       plot = figure6_p,
       width = a4_width,
       height = a4_height,
       units = "cm",
       dpi = dpi_level,
       device = cairo_pdf, # 支持字体嵌入
       bg = "white")       # 透明背景可改为"transparent"

# 如需位图格式（TIFF无损压缩）
ggsave("./Report Plot/figure6/figure6_p.tiff",
       plot = figure6_p,
       width = a4_width, 
       height = a4_height,
       units = "cm",
       dpi = dpi_level,
       compression = "lzw", # 无损压缩
       bg = "white")

ggsave("./Report Plot/figure6/figure6f_p.tiff",
       plot = figure6f_p,
       width = a4_width/2, 
       height = a4_height/3,
       units = "cm",
       dpi = dpi_level,
       compression = "lzw", # 无损压缩
       bg = "white")

#Figure7A----
# 提取目标亚群
t_cells <- subset(MBM, idents = c(0, 1, 11, 13))

# 设置直观Cluster命名
t_cells$named_clusters <- factor(
  t_cells$seurat_clusters, 
  labels = c("CD8+ T cell",
             "CD4+ T cell",
             "Activated CD4+ T cell",
             "Cycling CD8+ T cell")
)

group_colors <- c(
  "CD8+ T cell" = "#4169E1",  
  "CD4+ T cell" = "#1E90FF",  
  "Activated CD4+ T cell" = "#4682B4",
  "Cycling CD8+ T cell" = "#000080"
)

# activation
activation_genes <- c("CD69", "IL2RA", "MKI67", "IFNG", "GZMB", "PRF1")
figure7a_p <- DoHeatmap(t_cells,
                        features = activation_genes,
                        group.by = "named_clusters",
                        group.colors = group_colors,
                        slot = "scale.data",    # 使用标准化数据
                        disp.min = -2,          # 截断Z-score范围
                        disp.max = 2,
                        angle = 45,             # 调整标签角度
                        size = 4) + 
  scale_fill_gradientn(colors = c("#2E22EA", "#F9F8F6", "#FD1222"))


#Figure7B----
# exhaustion
exhaustion_genes <- c("PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT")
figure7b_p <- DoHeatmap(t_cells,
                        features = exhaustion_genes,
                        group.by = "named_clusters",
                        group.colors = group_colors,
                        cells = order(t_cells$named_clusters), # 按cluster排序
                        size = 4,
                        disp.min = -1.5,
                        disp.max = 2.5) +
  scale_fill_gradientn(colors = c("#2E22EA", "#F9F8F6", "#FD1222"))

#Figure7C&D&E&F----
group_colors_treat <- c(
  "Untreated" = "#0571B0",  
  "Treated" = "#CA0020"
)

# 按细胞类型拆分成四个子对象
split_objects <- SplitObject(t_cells, split.by = "named_clusters")

# 生成每个细胞类型的对比图
heatmap_list <- lapply(split_objects, function(sub_obj) { 
  # 获取当前细胞类型名称（确保split_objects的列表命名正确）
  cell_type <- unique(sub_obj@meta.data$named_clusters)
  
  # 排列顺序：确保Pre在前，Post在后（直接操作子对象）
  ordered_cells <- Cells(sub_obj)[order(sub_obj$pre_post, decreasing = FALSE)] 
  
  # 双面板图生成（直接使用子对象sub_obj）
  p1 <- DoHeatmap(sub_obj,
                  features = activation_genes,
                  cells = ordered_cells,
                  group.by = "pre_post",
                  group.colors = group_colors_treat,
                  slot = "scale.data", 
                  disp.min = -2, disp.max = 2,
                  size = 3.5, angle = 30) +
    scale_fill_gradientn(colors = c("#2E22EA", "#F9F8F6", "#FD1222"), guide = "none") +
    ggtitle(paste(cell_type, "Activation")) +
    theme(plot.title = element_text(size = 10),legend.position = "none")  # 统一标题字号
  
  p2 <- DoHeatmap(sub_obj,
                  features = exhaustion_genes,
                  cells = ordered_cells,
                  group.by = "pre_post",
                  group.colors = group_colors_treat,
                  slot = "scale.data",
                  disp.min = -1.5, disp.max = 2.5,
                  size = 3.5, angle = 30) +
    scale_fill_gradientn(colors = c("#2E22EA", "#F9F8F6", "#FD1222"), guide = "none") +
    ggtitle(paste(cell_type, "Exhaustion")) +
    theme(plot.title = element_text(size = 10),legend.position = "none") 
  
  # 组合双面板并调整比例
  p1 + p2 + plot_layout(widths = c(1, 0.8)) 
})





group_colors_treat <- c(
  "Untreated" = "#0571B0",  
  "Treated" = "#CA0020"
)

# 修改后的生成函数
heatmap_list <- lapply(split_objects, function(sub_obj) {
  # 步骤1：将pre_post中的"Pre/Post"替换为"Untreated/Treated"
  sub_obj@meta.data$pre_post <- ifelse(sub_obj@meta.data$pre_post == "Pre", "Untreated", "Treated")
  
  # 强制转换为因子以保持顺序
  sub_obj@meta.data$pre_post <- factor(sub_obj@meta.data$pre_post,
                                       levels = c("Untreated", "Treated"))
  
  cell_type <- unique(sub_obj@meta.data$named_clusters)
  
  # 步骤2：排序时基于新分组名称
  ordered_cells <- Cells(sub_obj)[order(sub_obj$pre_post, decreasing = FALSE)] 
  
  p1 <- DoHeatmap(sub_obj,
                  features = activation_genes,
                  cells = ordered_cells,
                  group.by = "pre_post",
                  group.colors = group_colors_treat,  # 直接使用新颜色
                  group.bar = TRUE,  # 确保显示分组色条
                  slot = "scale.data", 
                  disp.min = -2, disp.max = 2,
                  size = 3.5, angle = 30) +
    scale_fill_gradientn(colors = c("#2E22EA", "#F9F8F6", "#FD1222"), guide = "none") +
    ggtitle(paste(cell_type, "Activation")) +
    theme(plot.title = element_text(size = 10), legend.position = "none")
  
  p2 <- DoHeatmap(sub_obj,
                  features = exhaustion_genes,
                  cells = ordered_cells,
                  group.by = "pre_post",
                  group.colors = group_colors_treat,
                  group.bar = TRUE,  # 确保显示分组色条
                  slot = "scale.data",
                  disp.min = -1.5, disp.max = 2.5,
                  size = 3.5, angle = 30) +
    scale_fill_gradientn(colors = c("#2E22EA", "#F9F8F6", "#FD1222"), guide = "none") +
    ggtitle(paste(cell_type, "Exhaustion")) +
    theme(plot.title = element_text(size = 10), legend.position = "none")
  
  p1 + p2 + plot_layout(widths = c(1, 0.8))
})







figure7c_p <- heatmap_list[[1]]
figure7d_p <- heatmap_list[[2]]
figure7e_p <- heatmap_list[[3]]
figure7f_p <- heatmap_list[[4]]

#Figure7G----
# 设置分组比较
Idents(t_cells) <- "pre_post" # 设置分组列
markers <- FindMarkers(t_cells,
                       ident.1 = "Treated",  
                       ident.2 = "Untreated",    
                       min.pct = 0.25,          
                       logfc.threshold = 0.25,   
                       test.use = "wilcox")      

# 添加基因列并计算转换值
volcano_data <- markers %>%
  tibble::rownames_to_column("gene") %>%
  mutate(
    log10p = -log10(p_val_adj),
    regulation = case_when(
      avg_log2FC > 0 & p_val_adj < 0.05 ~ "Up",
      avg_log2FC < 0 & p_val_adj < 0.05 ~ "Down",
      TRUE ~ "Not sig"
    )
  )

# 对应调整基因筛选部分
top_up <- volcano_data %>% 
  filter(regulation == "Up") %>% 
  arrange(p_val_adj) %>% 
  head(5)  

top_down <- volcano_data %>% 
  filter(regulation == "Down") %>% 
  arrange(p_val_adj) %>% 
  head(5)  

# 提取全部上调下调基因
all_up <- volcano_data %>% 
  filter(regulation == "Up") %>% 
  arrange(p_val_adj) 

all_down <- volcano_data %>% 
  filter(regulation == "Down") %>% 
  arrange(p_val_adj)

# 修改后的绘图代码
figure7g_p <- ggplot(volcano_data, aes(x = avg_log2FC, y = log10p)) +
  geom_point(aes(color = regulation), alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Down" = "#2E22EA", "Up" = "#FD1222", "Not sig" = "grey60")) +
  
  # 新增坐标轴优化
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +  # 确保0点在x轴中心对称
  geom_vline(xintercept = 0, color = "grey40", alpha = 0.8) +    # 高亮0基准线
  
  # 调整显著性阈值线
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey40") +
  
  # 减少标签数量至上下调各5个
  geom_text_repel(
    data = rbind(
      top_up %>% head(5),  # 仅显示上调Top5
      top_down %>% head(5) # 仅显示下调Top5
    ),
    aes(label = gene),
    size = 3.5,
    box.padding = 0.3,
    max.overlaps = 20,      # 减少标签重叠容忍度
    segment.color = "grey50", # 标签引导线颜色
    nudge_x = 0.1          # 横向微调标签位置
  ) +
  
  # 优化标签和主题
  labs(
    x = expression(Log[2]*" Fold Change"),
    y = expression(-Log[10]*"(Adj.P-value)"),
    color = "Regulation"
  ) +
  theme_classic(base_size = 14) +  # 增大基础字体
  theme(
    panel.grid.major = element_line(color = "grey93"), # 更浅的网格线
    legend.position = "top",
    axis.title = element_text(face = "bold"),          # 坐标轴标题加粗
    legend.title = element_text(face = "bold")
  )

# 保存结果
write.csv(top_up, "./data/Mel_T_cell/Top5_Upregulated_Genes.csv", row.names = FALSE)
write.csv(top_down, "./data/Mel_T_cell/Top5_Downregulated_Genes.csv", row.names = FALSE)

write.csv(all_up, "./data/Mel_T_cell/Upregulated_Genes.csv", row.names = FALSE)
write.csv(all_down, "./data/Mel_T_cell/Downregulated_Genes.csv", row.names = FALSE)

#Figure7H----
# 步骤1: 准备基因列表
gene_list <- all_up$gene  # 提取差异基因symbol
gene_list_down <- all_down$gene

# 步骤2: 将symbol转换为Entrez ID（enrichGO默认需要）
entrez_ids <- bitr(gene_list, 
                   fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = org.Hs.eg.db) %>% 
  pull(ENTREZID)

# 步骤3: GO富集分析
go_res <- enrichGO(
  gene          = entrez_ids,
  OrgDb         = org.Hs.eg.db,        # 指定人类数据库
  keyType       = "ENTREZID",          # 输入基因类型
  ont           = "BP",                # 选择生物学过程(Biological Process)
  pAdjustMethod = "BH",                # 多重检验校正方法
  pvalueCutoff  = 0.05,                # p值阈值
  qvalueCutoff  = 0.2,                 # q值阈值
  readable      = TRUE                 # 转换EntrezID为基因symbol
)


go_res@result <- go_res@result %>%
  separate(GeneRatio, c("GeneInTerm", "GeneTotal"), sep = "/", convert = TRUE) %>%
  mutate(GeneRatio = GeneInTerm / GeneTotal)
go_res@result$Description <- gsub("\\(.*\\)", "", go_res@result$Description)

# 步骤4: 结果可视化
# 绘制前15显著GO项的点图
# 原始长名称： 
original_name <- "adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains"

# 建议缩短为（保留核心语义）：
short_name <- "adaptive immune response via IgSF somatic recombination"

# 定位并替换（避免误改其他名称）
go_res@result$Description[go_res@result$Description == original_name] <- short_name

figure7h_1 <- ggplot(go_res@result[1:15, ], aes(x=GeneRatio, y=reorder(Description, GeneRatio))) + 
  geom_point(aes(size=Count, color=pvalue)) + 
  theme_bw() +
  theme(axis.title.y = element_blank())

figure7h_2 <- ggplot(go_res@result[1:15, ], aes(x=GeneRatio, y=reorder(Description, GeneRatio))) + 
  geom_point(aes(size=Count, color=pvalue)) + 
  theme_bw() +
  theme(axis.title.y = element_blank())

# 导出结果到CSV
write.csv(go_res@result, "./data/Mel_T_cell/Upregulated_GO_enrichment_results.csv", row.names = FALSE)
write.csv(go_res@result, "./data/Mel_T_cell/Downregulated_GO_enrichment_results.csv", row.names = FALSE)

#Figure7----
figure7_p_1 <- 
  (figure7c_p) /
  (figure7d_p) /
  (figure7e_p) /
  (figure7f_p)

# 设置A4纸张尺寸参数（竖版：21x29.7cm）
a4_width <- 21   # 厘米（短边）
a4_height <- 29.7 # 厘米（长边）
dpi_level <- 600  # 打印级分辨率

# 推荐矢量格式（PDF）
ggsave("./Report Plot/figure7/figure7_p_1.pdf", 
       plot = figure7_p_1,
       width = a4_width,
       height = a4_height,
       units = "cm",
       dpi = dpi_level,
       device = cairo_pdf, # 支持字体嵌入
       bg = "white")       # 透明背景可改为"transparent"

# 如需位图格式（TIFF无损压缩）
ggsave("./Report Plot/figure7/figure7_p_1.tiff",
       plot = figure7_p_1,
       width = a4_width, 
       height = a4_height,
       units = "cm",
       dpi = dpi_level,
       compression = "lzw", # 无损压缩
       bg = "white")

ggsave("./Report Plot/figure7/figure7b_p.tiff",
       plot = figure7b_p,
       width = a4_width, 
       height = a4_height/2,
       units = "cm",
       dpi = dpi_level,
       compression = "lzw", # 无损压缩
       bg = "white")

# 提取post患者的元数据
post_meta <- MBM@meta.data[MBM@meta.data$pre_post == "Post", ]

# 提取所有免疫治疗列的非空值
therapy_values <- unlist(post_meta[, grep("Immunotherapy #", colnames(post_meta))])
therapy_non_na <- therapy_values[!is.na(therapy_values) & therapy_values != ""]

# 创建频次统计表
if(length(therapy_non_na) > 0) {
  therapy_counts <- sort(table(therapy_non_na), decreasing = TRUE)
  result_df <- data.frame(
    Treatment = names(therapy_counts),
    Frequency = as.numeric(therapy_counts),
    row.names = NULL
  )
} else {
  result_df <- data.frame(Treatment = character(), Frequency = integer())
}

# 打印表格型结果
print(knitr::kable(result_df, align = c("l","r")))

#Figure8:Pre-Processing----
Mel_meta_T <- subset(Mel_metastases, idents = c("0","1","6"))
MBM_T <- subset(MBM, idents = c("0","1", "11","13"))

unique(Mel_meta_T$cell_type)
unique(MBM_T$cell_type)

Mel_meta_T$group <- "Mel_metastases"
MBM_T$group <- "MBM"

merged_T <- merge(Mel_meta_T, MBM_T)

merged_T <- NormalizeData(merged_T) %>%
  FindVariableFeatures(nfeatures = 3000) %>%
  ScaleData(vars.to.regress = c("nCount_RNA")) %>% 
  RunPCA(npcs = 50, verbose = FALSE)

merged_T <- RunHarmony(
  merged_T,
  group.by.vars = "orig.ident",
  theta = 2, # Adjust batch correction strength
  lambda = 1,
  plot_convergence = TRUE
)

ElbowPlot(merged_T, ndims = 50)
merged_T <- FindNeighbors(object = merged_T, dims = 1:12) 
merged_T <- FindClusters(object = merged_T,
                               resolution = 0.3,
                               cluster.name = "TEST_res_0.3")

merged_T <- RunUMAP(merged_T,
                   dims = 1:12,
                   min.dist = 0.3,
                   spread = 1)
DimPlot(merged_T, reduction = "umap", group.by = "TEST_res_0.3")

marker_gene = c("CD3D", "IL7R", "CD4", "CD8A",        # T cells
                "ISG15", "RSAD2", "IFIT1", "OAS1",    # Interferon-stimulated genes (ISGs)
                "MKI67")          

DotPlot(merged_T, 
                      features = marker_gene, 
                      group.by = "TEST_res_0.3") + 
  theme(
    axis.text.x = element_text(            # X轴文本调整
      angle = 45,
      hjust = 1,
      face = "italic",
      family = "Arial",
      size = 8               # 新增字号设置（建议使用8-10pt）
    ),
    legend.position = "none"  # 禁用所有图例
  ) + 
  labs(x = NULL)              # 可选：移除默认X轴标题

merged_T <- subset(merged_T, subset = TEST_res_0.3 != 6)

#Figure8A----
figure8a_p <- DimPlot(
  object = merged_T,
  group.by = "TEST_res_0.3",
  cols = c("#003366","#000080","#3A6EA5","#0047AB","#87CEEB","#5D8AA8"),
  label = FALSE,
  pt.size = 0.5,
  reduction = "umap"
) + 
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle("") +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_blank(),   
    plot.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  )


#Figure8B----
figure8b_p <- DimPlot(
  object = merged_T,
  group.by = "group",
  cols = c("#6A3D9A", "#33A02C"),
  label = FALSE,
  pt.size = 0.5,
  reduction = "umap"
) + 
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle("") +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_blank(),   
    plot.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  )

#Figure8C----
figure8c_p <- DimPlot(
  object = merged_T,
  group.by = "pre_post",
  cols = c("#CA0020", "#0571B0"),
  label = FALSE,
  pt.size = 0.5,
  reduction = "umap"
) + 
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle("") +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_blank(),   
    plot.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  )

#Figure8D----
Idents(merged_T) <- "group"
print(table(Idents(merged_T)))

dge_results <- FindMarkers(
  object = merged_T,
  ident.1 = "MBM",   # 注意拼写一致性
  ident.2 = "Mel_metastases",
  test.use = "wilcox",              
  logfc.threshold = 0.25,           
  min.pct = 0.1,                    
  only.pos = FALSE                  
)

dge_results <- dge_results[dge_results$p_val > 0, ]

# 添加基因列并计算转换值
volcano_data <- dge_results %>%
  tibble::rownames_to_column("gene") %>%
  mutate(
    log10p = -log10(p_val_adj),
    regulation = case_when(
      avg_log2FC > 0 & p_val_adj < 0.05 ~ "Up",
      avg_log2FC < 0 & p_val_adj < 0.05 ~ "Down",
      TRUE ~ "Not sig"
    )
  )

# 对应调整基因筛选部分
top_up <- volcano_data %>% 
  filter(regulation == "Up") %>% 
  arrange(p_val_adj) %>% 
  head(5)  

top_down <- volcano_data %>% 
  filter(regulation == "Down") %>% 
  arrange(p_val_adj) %>% 
  head(5)  

# 提取全部上调下调基因
all_up <- volcano_data %>% 
  filter(regulation == "Up") %>% 
  arrange(p_val_adj) 

all_down <- volcano_data %>% 
  filter(regulation == "Down") %>% 
  arrange(p_val_adj)

# 修改后的绘图代码
figure8d_p <- ggplot(volcano_data, aes(x = avg_log2FC, y = log10p)) +
  geom_point(aes(color = regulation), alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Down" = "#2E22EA", "Up" = "#FD1222", "Not sig" = "grey60")) +
  
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +  # 确保0点在x轴中心对称
  geom_vline(xintercept = 0, color = "grey40", alpha = 0.8) +    # 高亮0基准线
  
  # 调整显著性阈值线
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey40") +
  
  # 选取四个目标基因
  geom_text_repel(
    data = subset(volcano_data, 
                  gene %in% c("MPP4", "CDC14C", "IKBKAP", "RAD23B")),  
    aes(label = gene),
    size = 3.5,
    box.padding = 0.3,
    max.overlaps = 20,
    segment.color = "grey50",
    nudge_x = 0.1,
    show.legend = FALSE  # 避免标注干扰图例
  )+
  
  # 优化标签和主题
  labs(
    x = expression(Log[2]*" Fold Change"),
    y = expression(-Log[10]*"(Adj.P-value)"),
    color = "Regulation"
  ) +
  theme_classic(base_size = 14) +  # 增大基础字体
  theme(
    panel.grid.major = element_line(color = "grey93"), # 更浅的网格线
    legend.position = "top",
    axis.title = element_text(face = "bold"),          # 坐标轴标题加粗
    legend.title = element_text(face = "bold")
  )

# 保存结果
write.csv(top_up, "./data/Merged_T_cell/Top5_Upregulated_Genes.csv", row.names = FALSE)
write.csv(top_down, "./data/Merged_T_cell/Top5_Downregulated_Genes.csv", row.names = FALSE)

write.csv(all_up, "./data/Merged_T_cell/Upregulated_Genes.csv", row.names = FALSE)
write.csv(all_down, "./data/Merged_T_cell/Downregulated_Genes.csv", row.names = FALSE)

#Figure8E----
# 步骤1: 准备基因列表
gene_list_up <- all_up$gene  # 提取差异基因symbol
gene_list_down <- all_down$gene

# 步骤2: 将symbol转换为Entrez ID（enrichGO默认需要）
entrez_ids <- bitr(gene_list_up,            #Choose up/down 
                   fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = org.Hs.eg.db) %>% 
  pull(ENTREZID)

# 步骤3: GO富集分析
go_res <- enrichGO(
  gene          = entrez_ids,
  OrgDb         = org.Hs.eg.db,        # 指定人类数据库
  keyType       = "ENTREZID",          # 输入基因类型
  ont           = "BP",                # 选择生物学过程(Biological Process)
  pAdjustMethod = "BH",                # 多重检验校正方法
  pvalueCutoff  = 0.05,                # p值阈值
  qvalueCutoff  = 0.2,                 # q值阈值
  readable      = TRUE                 # 转换EntrezID为基因symbol
)

go_res@result <- go_res@result %>%
  separate(GeneRatio, c("GeneInTerm", "GeneTotal"), sep = "/", convert = TRUE) %>%
  mutate(GeneRatio = GeneInTerm / GeneTotal)
go_res@result$Description <- gsub("\\(.*\\)", "", go_res@result$Description)

# 步骤4: 结果可视化
# 绘制前15显著GO项的点图
figure8e_1 <- ggplot(go_res@result[1:15, ], aes(x=GeneRatio, y=reorder(Description, GeneRatio))) + 
  geom_point(aes(size=Count, color=pvalue)) + 
  theme_bw() +
  theme(axis.title.y = element_blank())

figure8e_2 <- ggplot(go_res@result[1:15, ], aes(x=GeneRatio, y=reorder(Description, GeneRatio))) + 
  geom_point(aes(size=Count, color=pvalue)) + 
  theme_bw() +
  theme(axis.title.y = element_blank())

# 导出结果到CSV
write.csv(go_res@result, "./data/Merged_T_cell/Upregulated_GO_enrichment_results.csv", row.names = FALSE)
write.csv(go_res@result, "./data/Merged_T_cell/Downregulated_GO_enrichment_results.csv", row.names = FALSE)

#Figure8----
figure8_p_1 <- 
  (figure8a_p | figure8b_p) /  
  (figure8c_p | figure8d_p) 

figure8_p_2 <- 
  (figure8e_1) / 
  (figure8e_2)

# 设置A4纸张尺寸参数（竖版：21x29.7cm）
a4_width <- 21   # 厘米（短边）
a4_height <- 27.7 # 厘米（长边）
dpi_level <- 600  # 打印级分辨率

# 矢量格式（PDF）
ggsave("./Report Plot/figure8/figure8_p_1.pdf", 
       plot = figure8_p_1,
       width = a4_width,
       height = a4_height,
       units = "cm",
       dpi = dpi_level,
       device = cairo_pdf, # 支持字体嵌入
       bg = "white")       # "white"/"transparent"

# 位图格式（TIFF无损压缩）
ggsave("./Report Plot/figure8/figure8_p_1.tiff",
       plot = figure8_p_1,
       width = a4_width, 
       height = a4_height,
       units = "cm",
       dpi = dpi_level,
       compression = "lzw", # 无损压缩
       bg = "white")