library(reshape2)  # 加载 melt 所在包

cells_partitioned_by_annotation <- read.csv("Y:/long/publication_datasets/Vannan_2023_Lung_Fibrosis/Xenium/HE_annotations/cells_partitioned_by_annotation.csv")
metadata <- cells_partitioned_by_annotation[!duplicated(cells_partitioned_by_annotation$full_cell_id), ]

sample_list <- unique(metadata$sample)

result_list <- list()

for (sample in sample_list) {
  for (area in c("Minimally Remodeled Alveoli",
                 "Normal Alveoli",
                 "Advanced Remodeling",
                 "Remodeled Epithelium",
                 "Microscopic Honeycombing",
                 "Emphysema",
                 "Small Airway",
                 "Remnant Alveoli",
                 "Large Airway")) {
    # 筛选对应的元数据
    submeta <- metadata[metadata$sample == sample & metadata$Annotation_Type == area, ]

    # 如果数据为空或关键列缺失，则跳过
    if(nrow(submeta) != 0) {
      # 调用自定义函数计算 cophenetic 距离
      result <- compute_cophenetic_distances_from_df(submeta,
                                                     cluster_col = "final_CT",
                                                     x_col = "x_centroid",
                                                     y_col = "y_centroid")

      # 如果 result 或其中关键部分为空，也应跳过后续操作
      if(is.null(result) || is.null(result$row_cophenetic_df)) {
        next
      }

      row_coph <- as.matrix(result$row_cophenetic_df)
      df <- as.data.frame(row_coph)
      df$row <- rownames(row_coph)

      df_long <- melt(df, id.vars = "row", variable.name = "column", value.name = "value")
      df_long$sample <- sample
      df_long$area <- area

      result_list[[length(result_list) + 1]] <- df_long
    }
  }
}

final_df <- do.call(rbind, result_list)
write.csv(final_df, file = '分区域的Structure Score.csv', quote = FALSE, row.names = FALSE)


# 定义各组对应的样本列表
hd_samples <- c("VUHD090", "VUHD049", "VUHD038", "THD0008", "THD0011", "VUHD069", "VUHD095", "VUHD113", "VUHD116A", "VUHD116B")

la_samples <- c("VUILD102LF", "VUILD96LF", "VUILD104LF", "VUILD105LF", "VUILD48LF", "TILD117LF", "VUILD78LF", "VUILD91LF", "TILD049LF", "TILD111LF", "TILD167LF", "TILD299LF", "TILD315LF", "VUILD110", "VUILD49")

ma_samples <- c("VUILD102MF", "VUILD107MF", "VUILD96MF", "VUILD104MF", "VUILD105MF", "VUILD48MF", "TILD117MF", "VUILD78MF", "VUILD91MF", "TILD028MF", "TILD080MF", "TILD113MF", "TILD117MFB", "TILD130MF", "TILD175", "VUILD106", "VUILD115", "VUILD141", "VUILD58")


library(dplyr)
library(ggplot2)
library(patchwork)

final_df$group <- ifelse(final_df$sample %in% hd_samples, "Healthy",
                         ifelse(final_df$sample %in% la_samples, "Less",
                                ifelse(final_df$sample %in% ma_samples, "More", NA)))

# 如果还未设置因子水平，可以先设置：
final_df$group <- factor(final_df$group, levels = c("Healthy", "Less", "More"))

# 计算所有图统一的 y 轴范围
common_limits <- c(0, 1)

# 定义统一的颜色映射（也可以自己指定颜色）
color_scale <- scale_color_discrete(drop = FALSE, limits = c("Healthy", "Less", "More"))

# Normal Alveoli 箱线图（保留 y 轴）
plot_normal <- final_df %>%
  filter(area == "Normal Alveoli",
         row == "AT2",
         column == "Capillary") %>%
  ggplot(aes(x = "", y = value)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = group), width = 0.1, size = 2, alpha = 0.7) +
  labs(title = "Normal Alveoli", x = "", y = "Value") +
  scale_y_continuous(limits = common_limits) +
  color_scale +
  theme_minimal()

# Advanced Remodeling 箱线图（去除 y 轴）
plot_advanced <- final_df %>%
  filter(area == "Advanced Remodeling",
         row == "AT2",
         column == "Capillary") %>%
  ggplot(aes(x = "", y = value)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = group), width = 0.1, size = 2, alpha = 0.7) +
  labs(title = "Advanced Remodeling", x = "", y = "Value") +
  scale_y_continuous(limits = common_limits) +
  color_scale +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank())

# Remodeled Epithelium 箱线图（去除 y 轴）
plot_remodeled <- final_df %>%
  filter(area == "Remodeled Epithelium",
         row == "AT2",
         column == "Capillary") %>%
  ggplot(aes(x = "", y = value)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = group), width = 0.1, size = 2, alpha = 0.7) +
  labs(title = "Remodeled Epithelium", x = "", y = "Value") +
  scale_y_continuous(limits = common_limits) +
  color_scale +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank())

# Remnant Alveoli 箱线图（去除 y 轴）
plot_remnant <- final_df %>%
  filter(area == "Remnant Alveoli",
         row == "AT2",
         column == "Capillary") %>%
  ggplot(aes(x = "", y = value)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = group), width = 0.1, size = 2, alpha = 0.7) +
  labs(title = "Remnant Alveoli", x = "", y = "Value") +
  scale_y_continuous(limits = common_limits) +
  color_scale +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank())

# Minimally Remodeled Alveoli 箱线图（去除 y 轴）
plot_minimally <- final_df %>%
  filter(area == "Minimally Remodeled Alveoli",
         row == "AT2",
         column == "Capillary") %>%
  ggplot(aes(x = "", y = value)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = group), width = 0.1, size = 2, alpha = 0.7) +
  labs(title = "Minimally Remodeled Alveoli", x = "", y = "Value") +
  scale_y_continuous(limits = common_limits) +
  color_scale +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank())

# Microscopic Honeycombing 箱线图（去除 y 轴）
plot_microscopic <- final_df %>%
  filter(area == "Microscopic Honeycombing",
         row == "AT2",
         column == "Capillary") %>%
  ggplot(aes(x = "", y = value)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = group), width = 0.1, size = 2, alpha = 0.7) +
  labs(title = "Microscopic Honeycombing", x = "", y = "Value") +
  scale_y_continuous(limits = common_limits) +
  color_scale +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank())

# 将六个图排成一行（6列），并合并图例，统一放在最右边
combined_plot <- (plot_normal + plot_minimally + plot_remnant + plot_advanced + plot_microscopic + plot_remodeled) +
  plot_layout(ncol = 6, guides = "collect") &
  theme(legend.position = "right")

# 显示组合图
combined_plot
ggsave(filename = 'fig. 1h.pdf', plot = combined_plot, width = 7, height = 4, useDingbats = FALSE)
