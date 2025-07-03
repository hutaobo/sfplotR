setwd("./example")

# 加载必要的 R 包
library(dplyr)
library(reshape2)
library(ggplot2)

# 使用正斜杠（推荐）
data <- readRDS("GSE250346_Seurat_GSE250346_CORRECTED_SEE_RDS_README_082024.rds")

metadata <- data@meta.data

# 定义样本列表及其对应名称（注意每个元素之间用逗号分隔）
sample_list <- c("VUHD090", "VUHD049", "VUHD038", "THD0008", "THD0011",
                 "VUHD069", "VUHD095", "VUHD113", "VUHD116A", "VUHD116B",
                 "TILD080LA", "TILD028LA", "TILD167LA", "TILD113LA", "TILD111LA",
                 "VUILD49LA", "TILD130LA", "VUILD110LA", "TILD117LA", "VUILD78LA",
                 "VUILD91LA", "VUILD48LA2", "VUILD102LA", "VUILD96LA", "VUILD48LA1",
                 "TILD315MA", "VUILD58MA", "TILD049MA", "TILD299MA", "TILD117MA2",
                 "VUILD141MA", "VUILD142MA", "VUILD106MA", "VUILD115MA", "TILD117MA1",
                 "TILD175MA", "VUILD78MA", "VUILD91MA", "VUILD104MA1", "VUILD105MA2",
                 "VUILD102MA", "VUILD107MA", "VUILD96MA", "VUILD104MA2", "VUILD105MA1")

# 初始化一个空列表，用于存放每个样本处理后的数据框
result_list <- list()

# 遍历每个样本
for (sample in sample_list) {
  # 根据样本名称筛选出对应的元数据（假设 metadata 数据框已定义）
  submeta <- metadata[metadata$sample == sample, ]

  # 调用自定义函数计算 cophenetic 距离
  # 注意：compute_cophenetic_distances_from_df 函数应已在环境中定义
  result <- compute_cophenetic_distances_from_df(submeta,
                                                 cluster_col = "final_CT",
                                                 x_col = "x_centroid",
                                                 y_col = "y_centroid")

  # 将结果中的 cophenetic 数据转换为矩阵格式，保留原有的行名和列名
  row_coph <- as.matrix(result$row_cophenetic_df)

  # 将矩阵转换为数据框，同时将矩阵的行名保存为新的一列 'row'
  df <- as.data.frame(row_coph)
  df$row <- rownames(row_coph)

  # 将数据框从宽格式转换为长格式，保留 'row' 列不变
  # 转换后 'column' 列为原矩阵的列名，'value' 列为对应的数值
  df_long <- melt(df, id.vars = "row", variable.name = "column", value.name = "value")

  # 添加当前样本的名称信息
  df_long$sample <- sample

  # 将当前样本的长格式数据存入结果列表
  result_list[[length(result_list) + 1]] <- df_long
}

# 合并所有样本的结果数据框（按行合并）
final_df <- do.call(rbind, result_list)

# 定义各组对应的样本列表
hd_samples <- c("VUHD090", "VUHD049", "VUHD038", "THD0008", "THD0011", "VUHD069", "VUHD095", "VUHD113", "VUHD116A", "VUHD116B")
la_samples <- c("TILD080LA", "TILD028LA", "TILD167LA", "TILD113LA", "TILD111LA", "VUILD49LA", "TILD130LA", "VUILD110LA", "TILD117LA", "VUILD78LA", "VUILD91LA", "VUILD48LA2", "VUILD102LA", "VUILD96LA", "VUILD48LA1")
ma_samples <- c("TILD315MA", "VUILD58MA", "TILD049MA", "TILD299MA", "TILD117MA2", "VUILD141MA", "VUILD142MA", "VUILD106MA", "VUILD115MA", "TILD117MA1", "TILD175MA", "VUILD78MA", "VUILD91MA", "VUILD104MA1", "VUILD105MA2", "VUILD102MA", "VUILD107MA", "VUILD96MA", "VUILD104MA2", "VUILD105MA1")

# 根据 sample 列创建一个新的分组列
final_df <- final_df %>%
  mutate(group = case_when(
    sample %in% hd_samples ~ "HD",
    sample %in% la_samples ~ "LA",
    sample %in% ma_samples ~ "MA",
    TRUE ~ "Other"
  ))

# 查看分组情况
print(table(final_df$group))

# 将最终合并的数据框导出为 CSV 文件，不保留行名
write.csv(final_df, file = "cophenetic_distances_searcher_D_score_in_all_samples.csv", row.names = FALSE)

# 读取 CSV 数据
final_df <- read.csv("cophenetic_distances_searcher_D_score_in_all_samples.csv", stringsAsFactors = FALSE)
View(final_df)

library(ggplot2)
library(dplyr)
library(grid)  # 用于 unit()

# 筛选数据，只保留 row 和 column 均为 "AT1"、"AT2" 或 "Capillary" 的记录
valid_levels <- c("AT1", "AT2", "Capillary")

df_subset <- final_df %>%
  filter(row %in% valid_levels, column %in% valid_levels)

# 绘制按 row 和 column 分面的 boxplot，同时缩小面板间距，并在 boxplot 上添加分散的点
p <- ggplot(df_subset, aes(x = group, y = value)) +
  geom_boxplot(outlier.shape = NA) +  # 隐藏 boxplot 自带的异常值点
  geom_jitter(aes(color = group), width = 0.2, size = 2, alpha = 0.7) +
  facet_grid(row ~ column) +
  theme_bw() +
  theme(panel.spacing = unit(0.5, "lines")) +  # 缩小各个面板之间的间隙
  labs(title = "Searcher's D Score", x = "Group", y = "Value")

# 根据 facet 的行数和列数自动计算 PDF 的宽度和高度
nrow_panels <- length(unique(df_subset$row))
ncol_panels <- length(unique(df_subset$column))

# 设定每个面板的宽度和高度（单位：英寸），这里示例为 3 英寸
panel_width <- 2
panel_height <- 2

pdf_width <- ncol_panels * panel_width
pdf_height <- nrow_panels * panel_height

# 保存图形到 PDF 文件，长宽根据计算结果自动调整
ggsave("Searcher's D Score plot RRS.pdf", plot = p, width = pdf_width, height = pdf_height, units = "in", limitsize = FALSE)


library(ggplot2)
library(dplyr)
library(grid)  # 用于 unit()

# 筛选数据，只保留 row 和 column 均为 "AT1"、"AT2" 或 "Capillary" 的记录
valid_levels <- c("Multiciliated", "Basal", "Goblet", "RASC", "Secretory", "Myofibroblasts", "PNEC")

df_subset <- final_df %>%
  filter(row %in% valid_levels, column %in% valid_levels)

# 绘制按 row 和 column 分面的 boxplot，同时缩小面板间距，并在 boxplot 上添加分散的点
p <- ggplot(df_subset, aes(x = group, y = value)) +
  geom_boxplot(outlier.shape = NA) +  # 隐藏 boxplot 自带的异常值点
  geom_jitter(aes(color = group), width = 0.2, size = 2, alpha = 0.7) +
  facet_grid(row ~ column) +
  theme_bw() +
  theme(panel.spacing = unit(0.5, "lines")) +  # 缩小各个面板之间的间隙
  labs(title = "Searcher's D Score", x = "Group", y = "Value")

# 根据 facet 的行数和列数自动计算 PDF 的宽度和高度
nrow_panels <- length(unique(df_subset$row))
ncol_panels <- length(unique(df_subset$column))

# 设定每个面板的宽度和高度（单位：英寸），这里示例为 3 英寸
panel_width <- 2
panel_height <- 2

pdf_width <- ncol_panels * panel_width
pdf_height <- nrow_panels * panel_height

# 保存图形到 PDF 文件，长宽根据计算结果自动调整
ggsave("Searcher's D Score plot ARS.pdf", plot = p, width = pdf_width, height = pdf_height, units = "in", limitsize = FALSE)


library(ggplot2)
library(dplyr)
library(gridExtra)  # 用于拼接图形
library(grid)       # 用于 unit()

# 1. 找出所有唯一的 facet 组合（根据 row 与 column 分组）
facet_ids <- final_df %>% distinct(row, column)

# 2. 对每个组合分别绘图，并在图中添加面板标题
plot_list <- lapply(seq_len(nrow(facet_ids)), function(i) {
  sub_data <- final_df %>%
    filter(row == facet_ids$row[i], column == facet_ids$column[i])

  # 构造面板标题，可以根据需要自定义
  panel_title <- paste("Row:", facet_ids$row[i], "\nColumn:", facet_ids$column[i])

  ggplot(sub_data, aes(x = group, y = value)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = group), width = 0.2, size = 2, alpha = 0.7) +
    labs(title = panel_title, x = "Group", y = "Value") +
    scale_y_continuous(limits = c(0, 1)) +  # 设置 y 轴范围
    theme_bw() +
    theme(panel.spacing = unit(0.5, "lines"),
          plot.title = element_text(size = 10, face = "bold"))
})

# 3. 计算面板的行数和列数
nrow_panels <- length(unique(final_df$row))
ncol_panels <- length(unique(final_df$column))

# 4. 设置每个面板的大小（例如：每个面板宽3英寸，高3英寸），并打开 PDF 设备
pdf_width <- ncol_panels * 3
pdf_height <- nrow_panels * 2.5

pdf("Searcher's D Score plot All.pdf", width = pdf_width, height = pdf_height)

# 5. 拼接所有图形到一个页面上
grid.arrange(grobs = plot_list, ncol = ncol_panels, nrow = nrow_panels)

# 6. 关闭设备
dev.off()
