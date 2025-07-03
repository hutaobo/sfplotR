compute_cluster_nn_distance_df <- function(df,
                                             cluster_col = "Cluster",
                                             x_col = "x",
                                             y_col = "y") {
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

    # 提取坐标矩阵
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

    return(df_nearest_cluster_dist_df)
  }