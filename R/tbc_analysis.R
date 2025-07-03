#!/usr/bin/env Rscript

# tbc_analysis.R — R version of the transcript-by-cell spatial analysis
#
# This script performs a spatial analysis by treating transcripts of each gene
# as a unique cell type and calculating spatial proximity scores
# (cophenetic distances) between every gene and every true cell type.
#
# It is designed to be memory‑efficient for large datasets by processing genes
# in parallel using the 'future' framework and now shows a live progress bar
# using the 'progressr' package.

# --- 1. SCRIPT SETUP ---

suppressPackageStartupMessages({
  library(future)
  library(future.apply)
  library(dplyr)
  library(data.table)
  library(tools)
  library(progressr)          # NEW: progress bar support
})

# Choose a handler for how the progress bar is rendered.
# For CLI scripts 'txt' is safe everywhere; use 'progress' in RStudio.
handlers("txt")

# --- 2. WORKER FUNCTION ---

#' Process a single gene to compute its cophenetic distances to cell types.
#'
#' This function is executed by each parallel worker. It calculates the spatial
#' relationship between the transcripts of a single gene and the existing cell
#' types.
#' @param gene The name of the gene (string) to process.
#' @param coords_df A data frame of transcript coordinates with columns 'x',
#'        'y', 'feature_name'.
#' @param cell_metadata_df A data frame of cell metadata with columns 'x',
#'        'y', 'celltype'.
#' @param celltype_coph_df The pre‑computed cophenetic distance matrix between
#'        cell types.
#' @param coph_method The hierarchical clustering method (e.g., "average").
#' @return A one‑row data frame of cophenetic distances for the gene, or NULL on
#'         error.
process_gene <- function(gene,
                         coords_df,
                         cell_metadata_df,
                         celltype_coph_df,
                         coph_method) {
  tryCatch({
    # Select transcripts for the current gene
    gene_coords <- coords_df[coords_df$feature_name == gene, ]

    # Case 1: No spatial coordinates for this gene
    if (nrow(gene_coords) == 0) {
      # Gene absent in original cell‑type matrix → return a row of NAs
      if (!gene %in% rownames(celltype_coph_df)) {
        empty_cols <- colnames(celltype_coph_df)
        res <- data.frame(matrix(NA, nrow = 1, ncol = length(empty_cols)))
        colnames(res) <- empty_cols
        rownames(res) <- gene
        return(res)
      }

      # Otherwise just return the pre‑computed row, minus self‑distance
      res <- celltype_coph_df[gene, , drop = FALSE]
      res <- res[, colnames(res) != gene, drop = FALSE]
      return(res)
    }

    # Case 2: Compute cophenetic distances for this gene's transcripts
    gene_coords_formatted <- data.frame(
      x        = gene_coords$x,
      y        = gene_coords$y,
      celltype = gene
    )

    combined_df <- rbind(cell_metadata_df, gene_coords_formatted)

    coph_results <- compute_cophenetic_distances_from_df(
      df          = combined_df,
      cluster_col = "celltype",
      x_col       = "x",
      y_col       = "y",
      method      = coph_method
    )

    # Convert 'dist' to a full matrix first
    row_coph_mat  <- as.matrix(coph_results$row_cophenetic_df)

    # Extract the row for this gene
    gene_coph_row <- row_coph_mat[gene, , drop = FALSE]

    # Remove the self-distance (always 0)
    gene_coph_row <- gene_coph_row[, colnames(gene_coph_row) != gene, drop = FALSE]

    return(gene_coph_row)
  }, error = function(e) {
    warning(sprintf("Worker failed on gene %s: %s", gene, e$message))
    return(NULL)
  })
}

# --- 3. MAIN ANALYSIS FUNCTION ---

#' Perform transcript‑by‑cell spatial analysis with a live progress bar.
transcript_by_cell_analysis <- function(cell_metadata,
                                        coords,
                                        sample_name   = NULL,
                                        output_folder = NULL,
                                        coph_method   = "average",
                                        n_jobs        = 8) {
  # --- Output directory & sample name ---
  if (is.null(sample_name)) sample_name <- basename(folder)
  if (is.null(output_folder))
    output_folder <- file.path(".", paste0("t_by_c_", sample_name))
  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

  message("--- Loading and preparing data ---")

  # Remove control probes
  coords <- coords[!grepl("NegControl|Unassigned", coords$feature_name,
                          ignore.case = TRUE), ]

  # Keep needed cols only
  cell_metadata_clean <- cell_metadata[, c("x", "y", "celltype")]

  # Gene universe after filtering
  genes <- sort(unique(coords$feature_name))

  # --- Cell‑type‑by‑cell‑type baseline map ---
  message("--- Computing initial celltype distance map (StructureMap) ---")
  celltype_coph <- compute_cophenetic_distances_from_df(
    df          = cell_metadata_clean,
    cluster_col = "celltype",
    x_col       = "x",
    y_col       = "y",
    method      = coph_method
  )$row_cophenetic_df

  plot_cophenetic_heatmap(matrix          = celltype_coph,
                          matrix_name     = "row_coph",
                          output_dir      = output_folder,
                          output_filename = paste0("StructureMap_of_",
                                                   sample_name, ".pdf"),
                          sample          = sample_name)

  write.csv(celltype_coph,
            file      = file.path(output_folder,
                                  paste0("StructureMap_table_",
                                         sample_name, ".csv")),
            row.names = TRUE)

  # --- Parallel gene processing with progress bar ---
  message(sprintf("--- Initializing %d workers to process %d genes ---",
                  n_jobs, length(genes)))
  plan(multisession, workers = n_jobs)

  results_list <- with_progress({
    p <- progressor(steps = length(genes))

    future_lapply(
      X   = genes,
      FUN = function(gene,
                     coords_df,
                     cell_metadata_df,
                     celltype_coph_df,
                     coph_method) {
        p(sprintf("Processed %s", gene))
        process_gene(gene,
                     coords_df        = coords_df,
                     cell_metadata_df = cell_metadata_df,
                     celltype_coph_df = celltype_coph_df,
                     coph_method      = coph_method)
      },
      coords_df        = coords,
      cell_metadata_df = cell_metadata_clean,
      celltype_coph_df = celltype_coph,
      coph_method      = coph_method,
      future.seed      = TRUE
    )
  })
  names(results_list) <- genes

  # --- Combine & save ---
  message("--- Combining and saving final results ---")

  results_list <- results_list[!sapply(results_list, is.null)]

  final_df <- data.table::rbindlist(results_list, idcol = "gene")
  setDF(final_df)
  rownames(final_df) <- final_df$gene
  final_df$gene <- NULL

  output_csv <- file.path(output_folder,
                          paste0("t_and_c_result_", sample_name, ".csv"))
  write.csv(final_df, file = output_csv, row.names = TRUE)

  message(paste("[DONE] Outputs saved to:", output_folder))
}

# --- End of script ---
