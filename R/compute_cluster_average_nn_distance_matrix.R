compute_cluster_average_nn_distance_matrix <- function(df,
                                                       cluster_col = "Cluster",
                                                       x_col = "x",
                                                       y_col = "y") {
  # 检查必要的列是否存在
  required_cols <- c(cluster_col, x_col, y_col)
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(sprintf("数据框缺少必要的列: %s", paste(missing_cols, collapse = ", ")))
  }

  # 如果没有 cell_id 列，则使用行号作为 cell id
  cell_ids <- if ("cell_id" %in% colnames(df)) as.character(df$cell_id) else as.character(seq_len(nrow(df)))

  # 将聚类列转换为因子，并提取所有唯一分组
  clusters <- as.factor(df[[cluster_col]])
  unique_clusters <- levels(clusters)

  # 提取指定的坐标矩阵
  coords <- as.matrix(df[, c(x_col, y_col)])

  # 检查 RANN 包是否已安装
  if (!requireNamespace("RANN", quietly = TRUE)) {
    stop("需要安装 'RANN' 包，请先安装该包。")
  }

  # 对每个分组计算所有细胞到该分组中最近细胞的距离
  nn_dists_list <- lapply(unique_clusters, function(cluster_id) {
    cluster_mask <- clusters == cluster_id
    coords_cluster <- coords[cluster_mask, , drop = FALSE]
    if (nrow(coords_cluster) == 0) {
      rep(NA, nrow(df))
    } else {
      nn_result <- RANN::nn2(data = coords_cluster, query = coords, k = 1)
      nn_result$nn.dists[, 1]
    }
  })
  names(nn_dists_list) <- unique_clusters

  # 组合为矩阵，并设置行名为 cell_ids
  nn_dist_matrix <- do.call(cbind, nn_dists_list)
  rownames(nn_dist_matrix) <- cell_ids

  # 转换为 data frame 并添加每个细胞所属分组信息
  nn_dist_df <- as.data.frame(nn_dist_matrix)
  nn_dist_df$cell_cluster <- clusters

  # 计算每个细胞分组中，各分组距离的均值
  group_means <- aggregate(. ~ cell_cluster, data = nn_dist_df, FUN = mean, na.rm = TRUE)
  rownames(group_means) <- group_means$cell_cluster
  group_means$cell_cluster <- NULL

  # 删除整列均为 NA 的分组
  valid_cols <- colSums(is.na(group_means)) < nrow(group_means)
  group_means_clean <- group_means[, valid_cols, drop = FALSE]

  if (ncol(group_means_clean) == 0) {
    warning("所有分组的平均距离均为 NA，请检查数据。")
    return(data.frame())
  }

  return(group_means_clean)
}
