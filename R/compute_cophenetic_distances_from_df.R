compute_cophenetic_distances_from_df <- function(df,
                                                 cluster_col = "Cluster",
                                                 x_col = "x",
                                                 y_col = "y",
                                                 method = "average") {
  # 检查必要列是否存在
  if (!cluster_col %in% colnames(df)) {
    stop(sprintf("列 '%s' 不存在于数据框中，请检查输入数据。", cluster_col))
  }
  if (!all(c(x_col, y_col) %in% colnames(df))) {
    stop(sprintf("数据框必须包含 '%s' 和 '%s' 两列用于构成坐标。", x_col, y_col))
  }

  # 如果没有 cell_id 列，则使用行号作为 cell_id
  if ("cell_id" %in% colnames(df)) {
    cell_ids <- as.character(df$cell_id)
  } else {
    cell_ids <- as.character(1:nrow(df))
  }

  # 将 cluster 列转换为因子，并获取唯一分组
  clusters <- as.factor(df[[cluster_col]])
  unique_clusters <- levels(clusters)

  # 提取坐标矩阵，根据参数指定列名
  coords <- as.matrix(df[, c(x_col, y_col)])

  # 初始化存储最近邻距离的矩阵，行：cell_ids，列：unique_clusters
  df_nearest_cluster_dist <- matrix(NA, nrow = nrow(df), ncol = length(unique_clusters),
                                    dimnames = list(cell_ids, unique_clusters))

  # 加载 RANN 包用于最近邻搜索
  if (!requireNamespace("RANN", quietly = TRUE)) {
    stop("需要安装 'RANN' 包，请先安装该包。")
  }

  # 对每个分组计算所有细胞到该分组中最近细胞的距离
  for (c in unique_clusters) {
    mask_c <- clusters == c
    coords_c <- coords[mask_c, , drop = FALSE]

    if (nrow(coords_c) == 0) {
      df_nearest_cluster_dist[, c] <- NA
    } else {
      nn_result <- RANN::nn2(data = coords_c, query = coords, k = 1)
      df_nearest_cluster_dist[, c] <- nn_result$nn.dists[, 1]
    }
  }

  # 将最近邻距离矩阵转换为 data frame，并添加细胞所属分组信息
  df_nearest_cluster_dist_df <- as.data.frame(df_nearest_cluster_dist)
  df_nearest_cluster_dist_df$cell_cluster <- clusters

  # 对每个细胞分组，计算各分组距离的均值
  df_group_mean <- aggregate(. ~ cell_cluster,
                             data = df_nearest_cluster_dist_df,
                             FUN = mean, na.rm = TRUE)
  rownames(df_group_mean) <- df_group_mean$cell_cluster
  df_group_mean$cell_cluster <- NULL

  # 删除整列均为 NA 的分组
  df_group_mean_clean <- df_group_mean[, colSums(is.na(df_group_mean)) < nrow(df_group_mean), drop = FALSE]
  if (ncol(df_group_mean_clean) == 0) {
    warning("df_group_mean_clean 为空，请检查数据。")
    return(list(row_cophenetic_df = data.frame(),
                col_cophenetic_df = data.frame()))
  }

  # 对行和列分别进行层次聚类
  row_linkage <- hclust(dist(df_group_mean_clean), method = method)
  col_linkage <- hclust(dist(t(df_group_mean_clean)), method = method)

  # 计算 cophenetic 距离矩阵
  row_cophenetic <- cophenetic(row_linkage)
  col_cophenetic <- cophenetic(col_linkage)

  # 定义归一化函数
  normalize_matrix <- function(mat) {
    dmin <- min(mat, na.rm = TRUE)
    dmax <- max(mat, na.rm = TRUE)
    if (dmin == dmax) {
      return(mat)
    }
    (mat - dmin) / (dmax - dmin)
  }

  row_cophenetic_norm <- normalize_matrix(row_cophenetic)
  col_cophenetic_norm <- normalize_matrix(col_cophenetic)

  # 输出距离范围信息
  cat(sprintf("行方向 cophenetic 距离：原始范围 [%.4f, %.4f]，归一化后范围 [%.4f, %.4f]\n",
              min(row_cophenetic), max(row_cophenetic),
              min(row_cophenetic_norm), max(row_cophenetic_norm)))
  cat(sprintf("列方向 cophenetic 距离：原始范围 [%.4f, %.4f]，归一化后范围 [%.4f, %.4f]\n",
              min(col_cophenetic), max(col_cophenetic),
              min(col_cophenetic_norm), max(col_cophenetic_norm)))

  # 返回归一化后的 cophenetic 距离矩阵（行、列分别）
  return(list(row_cophenetic_df = row_cophenetic_norm,
              col_cophenetic_df = col_cophenetic_norm))
}
