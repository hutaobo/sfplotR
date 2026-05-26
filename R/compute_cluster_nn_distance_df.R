compute_cluster_nn_distance_df <- function(df,
                                             cluster_col = "Cluster",
                                             x_col = "x",
                                             y_col = "y") {
    # Check that required columns exist.
    if (!cluster_col %in% colnames(df)) {
      stop(sprintf("Column '%s' does not exist in the data frame.", cluster_col))
    }
    if (!all(c(x_col, y_col) %in% colnames(df))) {
      stop(sprintf("Data frame must contain '%s' and '%s' coordinate columns.",
                   x_col, y_col))
    }

    # Use row numbers as cell IDs when no cell_id column is present.
    if ("cell_id" %in% colnames(df)) {
      cell_ids <- as.character(df$cell_id)
    } else {
      cell_ids <- as.character(1:nrow(df))
    }

    # Convert the cluster column to a factor and extract all groups.
    clusters <- as.factor(df[[cluster_col]])
    unique_clusters <- levels(clusters)

    # Extract the coordinate matrix.
    coords <- as.matrix(df[, c(x_col, y_col)])

    # Initialize the nearest-neighbor distance matrix.
    df_nearest_cluster_dist <- matrix(NA, nrow = nrow(df), ncol = length(unique_clusters),
                                      dimnames = list(cell_ids, unique_clusters))

    # Check that RANN is installed for nearest-neighbor search.
    if (!requireNamespace("RANN", quietly = TRUE)) {
      stop("The 'RANN' package is required but is not installed.")
    }

    # For each group, compute each cell's distance to the nearest cell in that group.
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

    # Convert to a data frame and add each cell's own cluster.
    df_nearest_cluster_dist_df <- as.data.frame(df_nearest_cluster_dist)
    df_nearest_cluster_dist_df$cell_cluster <- clusters

    return(df_nearest_cluster_dist_df)
  }
