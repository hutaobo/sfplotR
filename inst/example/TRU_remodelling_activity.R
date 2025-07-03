setwd("./example")

# 加载必要的 R 包
library(dplyr)
library(reshape2)
library(ggplot2)
library(cluster)

# 读取 Seurat 对象，建议使用正斜杠（/）
seurat_data <- readRDS("GSE250346_Seurat_GSE250346_CORRECTED_SEE_RDS_README_082024.rds")

# 提取 Seurat 对象的元数据
seurat_metadata <- seurat_data@meta.data

# 定义要筛选的样本名称
selected_sample <- c("VUILD110LA")

# 根据样本名称筛选对应的元数据
filtered_metadata <- seurat_metadata[seurat_metadata$sample == selected_sample, ]

source('../R/compute_cluster_nn_distance_df.R')

# 调用自定义函数计算细胞到各簇最近邻距离的数据框
nn_distance_df <- compute_cluster_nn_distance_df(filtered_metadata,
                                                 cluster_col = "final_CT",
                                                 x_col = "x_centroid",
                                                 y_col = "y_centroid")

# 将计算结果拆分为数值型距离数据和原始的细胞分组标签
numeric_distance_df <- nn_distance_df[, -which(names(nn_distance_df) == "cell_cluster")]
original_cluster_labels <- nn_distance_df[, "cell_cluster"]




library(ggplot2)
library(jpeg)
library(grid)

# 载入图片 (背景图)
img_path <- "C:/Users/taobo.hu/Downloads/morphology_focus.ome.tif - VUILD110.jpg"
img <- readJPEG(img_path)

# 获取图片尺寸
img_dim <- dim(img)  # [高度, 宽度, 色彩通道]

# 构建数据框
df <- data.frame(
  value = ifelse(original_cluster_labels %in% c('AT2', 'Capillary'),
                 apply(numeric_distance_df[, c('AT2', 'Capillary')], 1, max),
                 NA),
  x = filtered_metadata$x_centroid,
  y = filtered_metadata$y_centroid
)
df <- df[complete.cases(df),]

# 沿纵轴（Y轴）翻转180度
df$y <- img_dim[1] - df$y

# 绘制热力图叠加到图片上
p <- ggplot(df, aes(x = x, y = y)) +
  annotation_custom(
    rasterGrob(img,
               width = unit(1, "npc"),
               height = unit(1, "npc")),
    xmin = 0, xmax = img_dim[2], ymin = 0, ymax = img_dim[1]
  ) +
  stat_density2d(aes(fill = ..level.., alpha = ..level..),
                 geom = "polygon", color = NA) +
  scale_fill_gradientn(
    colors = rev(heat.colors(10)),
    oob = scales::squish  # 将超出范围的值压缩到范围内
  ) +
  scale_alpha(range = c(0, 0.7), guide = FALSE) +
  coord_fixed(ratio = 1, xlim = c(0, img_dim[2]), ylim = c(0, img_dim[1]), expand = FALSE) +
  theme_void() +
  labs(title = "TRU remodelling activity") +
  theme(
    plot.title = element_text(hjust = 0.5)  # 标题居中
  )

# 保存图片为PDF
output_path <- "C:/Users/taobo.hu/Downloads/TRU_remodelling_activity.pdf"
ggsave(output_path, plot = p, device = "pdf", width = 8, height = 6)
