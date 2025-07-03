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

# numeric_distance_df <- nn_distance_df[, c('AT1', 'AT2', 'Capillary')]

# 对 nn_distance_df 的每一行进行归一化
numeric_distance_df <- t(apply(numeric_distance_df, 1, function(row) {
  row_min <- min(row)
  row_max <- max(row)
  # 如果当前行的最大值和最小值相同，避免除0错误，这里返回全0
  if(row_max == row_min) {
    return(rep(0, length(row)))
  }
  # 使用公式归一化到 [-1, 1]
  2 * (row - row_min) / (row_max - row_min) - 1
}))

# 转换回 data frame
numeric_distance_df <- as.data.frame(numeric_distance_df)

# 查看结果的前几行
head(numeric_distance_df)




original_cluster_labels <- nn_distance_df[, "cell_cluster"]

# 如果尚未安装 cluster 包，请先安装
if (!require(cluster)) {
  install.packages("cluster")
  library(cluster)
}

# 设置随机种子以保证结果可重复
set.seed(123)

# 使用 CLARA 进行聚类（k 为你希望分成的簇数，请根据需要调整）
clara_result <- clara(numeric_distance_df, k = 20, metric = "euclidean")

# 假设 clara_result$clustering 是 CLARA 得到的簇标签
confusion_matrix <- table(Clara = clara_result$clustering, Original = original_cluster_labels)
# View(confusion_matrix)

# 构造包含 cell_id 和 group 两列的数据框
output_df <- data.frame(cell_id = filtered_metadata$cell_id,
                        group   = clara_result$clustering)

# 保存为 CSV 文件，文件名为 "clara_cluster_assignment.csv"，不保存行名
write.csv(output_df, file = "clara_cluster_assignment.csv", row.names = FALSE)



# 加载所需包
library(ggplot2)
library(patchwork)

# 假设 df 数据框中已经包含 x_centroid, y_centroid 以及两种聚类结果
# 添加聚类结果列
df$ClaraCluster <- as.factor(clara_result$clustering)
df$OriginalCluster <- as.factor(original_cluster_labels)

# 分别绘制两个图
p1 <- ggplot(df, aes(x = x_centroid, y = y_centroid, color = ClaraCluster)) +
  geom_point() +
  theme_bw() +
  labs(title = "CLARA 聚类结果",
       x = "x_centroid",
       y = "y_centroid")

p2 <- ggplot(df, aes(x = x_centroid, y = y_centroid, color = OriginalCluster)) +
  geom_point() +
  theme_bw() +
  labs(title = "原始聚类结果",
       x = "x_centroid",
       y = "y_centroid")

# 使用 patchwork 将两个图上下排列
combined_plot <- p1 / p2

# 保存为 PDF 文件，文件名为 Clustering_Comparison.pdf
ggsave(filename = "Clustering_Comparison.pdf", plot = p1, width = 12, height = 12)






# 1. 计算数值数据框的行间欧式距离，并进行层次聚类（采用 Ward.D2 方法）
distance_matrix <- dist(numeric_distance_df)
hierarchical_clustering <- hclust(distance_matrix, method = "ward.D2")

# 2. 自动确定最佳聚类数 k：使用 silhouette 方法（取 k 从 2 到 21，可根据数据情况调整范围）
silhouette_widths <- numeric(21)  # 初始化向量，假定 k 的最大值为 21
for (k in 2:21) {
  temp_clusters <- cutree(hierarchical_clustering, k = k)
  sil <- silhouette(temp_clusters, distance_matrix)
  silhouette_widths[k] <- mean(sil[, 3])
}
optimal_cluster_count <- which.max(silhouette_widths)
cat("自动确定的最佳聚类数 k =", optimal_cluster_count, "\n")

# 3. 根据最佳聚类数 k 对数据进行划分
hierarchical_clusters <- cutree(hierarchical_clustering, k = optimal_cluster_count)

# 4. 构造混淆矩阵，对比层次聚类结果与原始细胞分组标签
confusion_matrix <- table(Hierarchical_Clusters = hierarchical_clusters,
                          Original_Labels = original_cluster_labels)
print(confusion_matrix)

