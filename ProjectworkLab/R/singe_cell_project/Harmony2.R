#Import Package####
library(Seurat)
library(harmony)
library(dplyr)
library(readr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(reshape2)
source("./Script/My_Function.R")

# Import Data####
MBM <- readRDS("./data/MBM.RDS")
Mel_metastases <- readRDS("./data/Mel_metastases.RDS")

#Save/Load Data####
save.image("./Harmony.RData")
load("./Harmony.RData")

# 提取两个数据集的CD4+ T细胞
mcd4 <- subset(MBM, idents = c("1", "11"))  # MBM的CD4+/Activated CD4+
scd4 <- subset(Mel_metastases, idents = "1") # Mel_metastases的CD4+

# 添加分组信息
mcd4$group <- "MBM"
scd4$group <- "Mel_metastases"

# 验证修改后再合并（保留原始代码结构）
merged_cd4T <- merge(mcd4, scd4)

# Standardized Process
merged_cd4T <- NormalizeData(merged_cd4T) %>%
  FindVariableFeatures(nfeatures = 3000) %>%
  ScaleData(vars.to.regress = c("nCount_RNA")) %>% 
  RunPCA(npcs = 50, verbose = FALSE)

# Harmony integrates batch effects----
merged_cd4T <- RunHarmony(
  merged_cd4T,
  group.by.vars = "orig.ident",
  theta = 2, # Adjust batch correction strength
  lambda = 1,
  plot_convergence = TRUE
)

# UMAP


# 设定身份确认分类正确（检查是否仅包含目标两组）
Idents(merged_cd4T) <- "cell_type"
print(table(Idents(merged_cd4T)))  # 必须显示两组存在

# 执行差异分析（核心代码）
dge_results <- FindMarkers(
  object = merged_cd4T,
  ident.1 = "Activated CD4+ T cell",   # 注意拼写一致性
  ident.2 = "CD4+ T cell",
  test.use = "wilcox",              # 对单细胞数据更稳健的非参数检验
  logfc.threshold = 0.25,           # 降低阈值以捕捉更多ISG类基因
  min.pct = 0.1,                    # 确保低丰度但重要的活化基因不被过滤
  only.pos = FALSE                  # 保留双向差异基因
)

# 添加基因名称列方便后续操作
dge_results$gene <- rownames(dge_results)

# 筛选显著差异基因（根据文献调整阈值）
sig_genes <- subset(dge_results, p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)
print(paste("显著差异基因数量：", nrow(sig_genes)))

# 构建绘图数据框
plot_data <- dge_results %>%
  mutate(
    color = case_when(
      avg_log2FC > 0.5 & p_val_adj < 0.05 ~ "Upregulated",
      avg_log2FC < -0.5 & p_val_adj < 0.05 ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )

# 创建可交互的Enhanced Volcano
ggplot(plot_data, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = color), alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Upregulated" = "#E64B35", 
                                "Downregulated" = "#3182BD", 
                                "Not significant" = "#BDBDBD")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#666666") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "#666666") +
  ggtitle("Actived vs Other T Cells") +
  labs(x = "Log2(Fold Change)", y = "-Log10(Adjusted P-value)") +
  theme_classic(base_size = 14) +
  theme(legend.position = "right") +
  # 标记关键基因（需与用户数据匹配）
  ggrepel::geom_text_repel(
    data = subset(plot_data, gene %in% c("ISG15",
                                         "RSAD2",
                                         "IFI44L",
                                         "IFIT3",
                                         "IFIT2",
                                         "IFI6",
                                         "MX1")),
    aes(label = gene),
    box.padding = 0.5,
    max.overlaps = 100
  )





#----------------------For vlnplot-----------------------#
# 确保元数据列存在
if("cell_type" %in% colnames(mcd4@meta.data)) {
  # 修改MBM数据集中的标签
  mcd4$cell_type <- ifelse(
    mcd4$cell_type == "CD4+ T cell",  # 筛选条件
    "brain CD4+ T cell",              # 符合条件的替换值
    mcd4$cell_type                    # 不符合的保持原值
  )
  # 查看修改结果
  table(mcd4$cell_type, useNA = "always")
}

if("cell_type" %in% colnames(scd4@meta.data)) {
  # 修改Mel_metastases数据集中的标签
  scd4$cell_type <- ifelse(
    scd4$cell_type == "CD4+ T cell",
    "other CD4+ T cell",
    scd4$cell_type
  )
  # 查看修改结果
  table(scd4$cell_type, useNA = "always")
}






# Step 1：验证当前值分布
table(merged_cd4T$cell_type, useNA = "always")

# Step 2：精准替换目标标签
if("cell_type" %in% colnames(merged_cd4T@meta.data)) {
  # 创建新列防止原数据意外修改
  merged_cd4T@meta.data$new_cell_type <- as.character(merged_cd4T@meta.data$cell_type)
  
  # 限定严格匹配数字2（防止部分匹配错误）
  target_idx <- which(merged_cd4T@meta.data$new_cell_type == "2") 
  merged_cd4T@meta.data$new_cell_type[target_idx] <- "Actived brain CD4+ T cell"  # 注：生物学常用拼写为"Activated"
  
  # 核对修改后的分布
  print("修改后标签分布：")
  print(table(merged_cd4T@meta.data$new_cell_type, useNA = "always"))
  
  # 将修改列设置回原列名称（可选）
  merged_cd4T@meta.data$cell_type <- merged_cd4T@meta.data$new_cell_type
  merged_cd4T@meta.data$new_cell_type <- NULL
}

# Step 3：最终验证（确保没有覆盖其他有效数值）
table(merged_cd4T$cell_type)



# 绘制小提琴图
ISG15_vil_p <- VlnPlot(merged_cd4T, 
                       features = c("ISG15"), 
                       group.by = "cell_type",
                       pt.size = 0,    # 隐藏数据点
                       ncol = 3,       # 两列布局
                       split.by = NULL) + 
  theme_classic() +
  labs(title = "ISG15", 
       x = "", 
       y = "Expression Level") +
  theme(
    axis.text.x = element_blank()
  )


RSAD2_vil_p <- VlnPlot(merged_cd4T, 
                       features = c("RSAD2"), 
                       group.by = "cell_type",
                       pt.size = 0,    # 隐藏数据点
                       ncol = 2,       # 两列布局
                       split.by = NULL) + 
  theme_classic() +
  labs(title = "RSAD2", 
       x = "", 
       y = "Expression Level") +
  theme(
    axis.text.x = element_blank()
  )

IFIT1_vil_p <- VlnPlot(merged_cd4T, 
                       features = c("IFIT1"), 
                       group.by = "cell_type",
                       pt.size = 0,    # 隐藏数据点
                       ncol = 2,       # 两列布局
                       split.by = NULL) + 
  theme_classic() +
  labs(title = "IFIT1", 
       x = "", 
       y = "Expression Level")+
  theme(
    axis.text.x = element_blank()
  )

OAS1_vil_p <- VlnPlot(merged_cd4T, 
                      features = c("OAS1"), 
                      group.by = "cell_type",
                      pt.size = 0,    # 隐藏数据点
                      ncol = 2,       # 两列布局
                      split.by = NULL) + 
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