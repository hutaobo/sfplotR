compute_cophenetic_distances_from_df <- function(df,
                                                 cluster_col = "Cluster",
                                                 x_col = "x",
                                                 y_col = "y",
                                                 method = "average",
                                                 print = FALSE) {
  # Check that the required columns exist
  if (!cluster_col %in% colnames(df)) {
    stop(sprintf("Column '%s' does not exist in the data frame. Please check your input.",
                 cluster_col))
  }
  if (!all(c(x_col, y_col) %in% colnames(df))) {
    stop(sprintf("Data frame must contain both '%s' and '%s' columns for coordinates.",
                 x_col, y_col))
  }

  # Use 'cell_id' column if present; otherwise, use row numbers
  if ("cell_id" %in% colnames(df)) {
    cell_ids <- as.character(df$cell_id)
  } else {
    cell_ids <- as.character(1:nrow(df))
  }

  # Convert the cluster column to factor and get unique levels
  clusters <- as.factor(df[[cluster_col]])
  unique_clusters <- levels(clusters)

  # Extract the coordinate matrix using the specified column names
  coords <- as.matrix(df[, c(x_col, y_col)])

  # Initialize a matrix to store nearest-neighbor distances:
  # rows = cell_ids, cols = unique_clusters
  df_nearest_cluster_dist <- matrix(
    NA,
    nrow = nrow(df),
    ncol = length(unique_clusters),
    dimnames = list(cell_ids, unique_clusters)
  )

  # Load RANN package for nearest-neighbor search
  if (!requireNamespace("RANN", quietly = TRUE)) {
    stop("The 'RANN' package is required but not installed. Please install it first.")
  }

  # For each cluster, compute distance from every cell to its nearest cell in that cluster
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

  # Convert the nearest-neighbor distance matrix to a data frame
  # and add a column for each cell's cluster
  df_nearest_cluster_dist_df <- as.data.frame(df_nearest_cluster_dist)
  df_nearest_cluster_dist_df$cell_cluster <- clusters

  # Compute the mean distance to each cluster for each cell cluster
  df_group_mean <- aggregate(
    . ~ cell_cluster,
    data = df_nearest_cluster_dist_df,
    FUN = mean,
    na.rm = TRUE
  )
  rownames(df_group_mean) <- df_group_mean$cell_cluster
  df_group_mean$cell_cluster <- NULL

  # Remove any cluster columns that are all NA
  df_group_mean_clean <- df_group_mean[,
    colSums(is.na(df_group_mean)) < nrow(df_group_mean),
    drop = FALSE
  ]
  if (ncol(df_group_mean_clean) == 0) {
    warning("df_group_mean_clean is empty. Please check your data.")
    return(list(
      row_cophenetic_df = data.frame(),
      col_cophenetic_df = data.frame()
    ))
  }

  # Perform hierarchical clustering on rows and columns
  row_linkage <- hclust(dist(df_group_mean_clean), method = method)
  col_linkage <- hclust(dist(t(df_group_mean_clean)), method = method)

  # Compute cophenetic distance matrices
  row_cophenetic <- cophenetic(row_linkage)
  col_cophenetic <- cophenetic(col_linkage)

  # Normalization helper function
  normalize_matrix <- function(mat) {
    dmin <- min(mat, na.rm = TRUE)
    dmax <- max(mat, na.rm = TRUE)
    if (dmin == dmax) return(mat)
    (mat - dmin) / (dmax - dmin)
  }

  row_cophenetic_norm <- normalize_matrix(row_cophenetic)
  col_cophenetic_norm <- normalize_matrix(col_cophenetic)

  if (print == TRUE) {
    # Print distance range information
    cat(sprintf(
      "Row cophenetic distances: original range [%.4f, %.4f], normalized range [%.4f, %.4f]\n",
      min(row_cophenetic), max(row_cophenetic),
      min(row_cophenetic_norm), max(row_cophenetic_norm)
    ))
    cat(sprintf(
      "Column cophenetic distances: original range [%.4f, %.4f], normalized range [%.4f, %.4f]\n",
      min(col_cophenetic), max(col_cophenetic),
      min(col_cophenetic_norm), max(col_cophenetic_norm)
    ))
  }

  # Return the normalized cophenetic distance matrices
  return(list(
    row_cophenetic_df = as.matrix(row_cophenetic_norm),
    col_cophenetic_df = as.matrix(col_cophenetic_norm)
  ))
}
