compute_cluster_average_nn_distance_matrix <- function(df,
                                                       cluster_col = "Cluster",
                                                       x_col = "x",
                                                       y_col = "y") {
  # Check that required columns exist.
  required_cols <- c(cluster_col, x_col, y_col)
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(sprintf("Data frame is missing required columns: %s",
                 paste(missing_cols, collapse = ", ")))
  }

  # Use row numbers as cell IDs when no cell_id column is present.
  cell_ids <- if ("cell_id" %in% colnames(df)) as.character(df$cell_id) else as.character(seq_len(nrow(df)))

  # Convert the cluster column to a factor and extract all groups.
  clusters <- as.factor(df[[cluster_col]])
  unique_clusters <- levels(clusters)

  # Extract the coordinate matrix.
  coords <- as.matrix(df[, c(x_col, y_col)])

  # Check that RANN is installed for nearest-neighbor search.
  if (!requireNamespace("RANN", quietly = TRUE)) {
    stop("The 'RANN' package is required but is not installed.")
  }

  # For each group, compute each cell's distance to the nearest cell in that group.
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

  # Combine distances into a matrix indexed by cell IDs.
  nn_dist_matrix <- do.call(cbind, nn_dists_list)
  rownames(nn_dist_matrix) <- cell_ids

  # Convert to a data frame and add each cell's own cluster.
  nn_dist_df <- as.data.frame(nn_dist_matrix)
  nn_dist_df$cell_cluster <- clusters

  # Compute mean nearest-neighbor distance from each cell cluster to every target cluster.
  group_means <- stats::aggregate(. ~ cell_cluster, data = nn_dist_df, FUN = mean, na.rm = TRUE)
  rownames(group_means) <- group_means$cell_cluster
  group_means$cell_cluster <- NULL

  # Remove target clusters whose distances are all NA.
  valid_cols <- colSums(is.na(group_means)) < nrow(group_means)
  group_means_clean <- group_means[, valid_cols, drop = FALSE]

  if (ncol(group_means_clean) == 0) {
    warning("All group mean distances are NA. Please check the input data.")
    return(data.frame())
  }

  return(group_means_clean)
}
