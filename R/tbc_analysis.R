#!/usr/bin/env Rscript

# tbc_analysis.R — R version of the transcript-by-cell spatial analysis
#
# This script performs a spatial analysis by treating transcripts of each gene
# as a unique cell type and calculating spatial proximity scores (cophenetic
# distances) between every gene and every true cell type.
#
# It is designed to be memory-efficient for large datasets by processing genes
# in parallel using the 'future' framework.

# --- 1. SCRIPT SETUP ---

# Install and load necessary packages
# Ensure you have these packages installed:
# install.packages(c("future", "future.apply", "dplyr", "data.table",
#                    "pheatmap", "RANN", "RColorBrewer", "latex2exp",
#                    "showtext"))

suppressPackageStartupMessages({
  library(future)
  library(future.apply)
  library(dplyr)
  library(data.table)
  library(tools)
})



# --- 2. WORKER FUNCTION ---

#' Process a single gene to compute its cophenetic distances to cell types.
#'
#' This function is executed by each parallel worker. It calculates the spatial
#' relationship between the transcripts of a single gene and the existing
#' cell types.
#'
#' @param gene The name of the gene (string) to process.
#' @param coords_df A data frame of transcript coordinates with columns 'x',
#'                  'y', 'feature_name'.
#' @param cell_metadata_df A data frame of cell metadata with columns 'x',
#'                         'y', 'celltype'.
#' @param celltype_coph_df The pre-computed cophenetic distance matrix
#'                         between cell types.
#' @param coph_method The hierarchical clustering method (e.g., "average").
#' @return A one-row data frame of cophenetic distances for the gene,
#'         or NULL on error.
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
      # If the gene is not in the original celltype matrix (e.g., a gene
      # with no cells), return a row of NAs.
      if (!gene %in% rownames(celltype_coph_df)) {
        empty_cols <- colnames(celltype_coph_df)
        res <- data.frame(matrix(NA, nrow = 1, ncol = length(empty_cols)))
        colnames(res) <- empty_cols
        rownames(res) <- gene
        return(res)
      }

      # Otherwise, return the pre-computed row from the celltype-by-celltype
      # matrix, removing the column corresponding to the gene itself.
      res <- celltype_coph_df[gene, , drop = FALSE]
      res <- res[, colnames(res) != gene, drop = FALSE]
      return(res)
    }

    # Case 2: Compute cophenetic distances for the gene's transcripts
    # Treat the gene's transcripts as a new temporary celltype
    gene_coords_formatted <- data.frame(
      x        = gene_coords$x,
      y        = gene_coords$y,
      celltype = gene
    )

    # Combine with the actual cell data
    combined_df <- rbind(
      cell_metadata_df,
      gene_coords_formatted
    )

    # Compute cophenetic distances on the combined data
    coph_results <- compute_cophenetic_distances_from_df(
      df          = combined_df,
      cluster_col = "celltype",
      x_col       = "x",
      y_col       = "y",
      method      = coph_method
    )

    # The result we need is the row for our gene from the 'row' matrix
    gene_coph_row <- coph_results$row_cophenetic_df[gene, , drop = FALSE]

    # Remove the column for the gene itself (distance to itself is 0)
    gene_coph_row <- gene_coph_row[
      , colnames(gene_coph_row) != gene,
      drop = FALSE
    ]

    return(gene_coph_row)

  }, error = function(e) {
    # Log any error without interrupting the parallel execution
    warning(sprintf("Worker failed on gene %s: %s", gene, e$message))
    return(NULL)
  })
}


# --- 3. MAIN ANALYSIS FUNCTION ---

#' Perform transcript-by-cell spatial analysis.
#'
#' This function loads spatial data, computes a base celltype-by-celltype
#' spatial distance matrix, and then iterates through all genes in parallel to
#' compute the spatial distance of each gene's transcripts to all cell types.
#'
#' @param sample_name  A name for the sample (defaults to the folder name).
#' @param output_folder Directory to save results.
#' @param coph_method  Clustering method for hclust ("average", "complete", etc.).
#' @param n_jobs       Number of parallel processes to use.
#' @param df           An optional data frame to map cell_id to a new "group"
#'                     annotation.
#' @return Creates output files in the specified directory.
transcript_by_cell_analysis <- function(cell_metadata,
                                        coords,
                                        sample_name   = NULL,
                                        output_folder = NULL,
                                        coph_method   = "average",
                                        n_jobs        = 4) {

  # --- Setup output directory and sample name ---
  if (is.null(sample_name)) {
    sample_name <- basename(folder)
  }
  if (is.null(output_folder)) {
    output_folder <- file.path(".", paste0("t_by_c_", sample_name))
  }
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }

  # --- Load Xenium Data ---
  # NOTE: R data loading for Xenium is different from Python.
  # The user needs to replace this section with their data loading code.
  # Common packages are Seurat, Voyager, or Giotto.
  # The goal is to produce two data frames:
  # 1. `cell_metadata`: Contains cell info with columns 'x', 'y',
  #    'cell_id', and a cluster/celltype column.
  # 2. `coords`: Contains transcript info with columns 'x', 'y',
  #    'feature_name', 'cell_id'.

  message("--- Loading and preparing data ---")
  message("NOTE: Please replace the placeholder data loading section with your code.")

  # #### START OF USER-DEFINED DATA LOADING ####
  # Example using placeholder data. Replace with your actual data loading.
  # We assume you have loaded your data into a Seurat object `so`.
  # so <- readRDS("path/to/your/seurat_object.rds")
  #
  # # 1. Prepare cell metadata
  # cell_metadata <- so@meta.data
  # spatial_coords <- as.data.frame(Embeddings(so, "spatial"))
  # cell_metadata$x <- spatial_coords$spatial_1
  # cell_metadata$y <- spatial_coords$spatial_2
  # cell_metadata$cell_id <- rownames(cell_metadata)
  #
  # # 2. Prepare transcript coordinates
  # # coords <- read.csv("path/to/transcripts.csv")
  # #### END OF USER-DEFINED DATA LOADING ####

  # For demonstration, let's create dummy data matching the required structure
  # genes <- unique(coords$feature_name)

  # Filter out control probes from transcripts
  coords <- coords[
    !grepl("NegControl|Unassigned", coords$feature_name,
           ignore.case = TRUE),
  ]

  # Select final columns
  cell_metadata_clean <- cell_metadata[, c("x", "y", "celltype")]

  # --- Compute and save initial celltype-by-celltype distances ---
  message("--- Computing initial celltype distance map (StructureMap) ---")
  celltype_coph <- compute_cophenetic_distances_from_df(
    df          = cell_metadata_clean,
    cluster_col = "celltype",
    x_col       = "x",
    y_col       = "y",
    method      = coph_method
  )$row_cophenetic_df

  # Plot and save the heatmap
  plot_cophenetic_heatmap(
    matrix          = celltype_coph,
    matrix_name     = "row_coph",
    output_dir      = output_folder,
    output_filename = paste0("StructureMap_of_", sample_name, ".pdf"),
    sample          = sample_name
  )

  # Save the table
  write.csv(
    celltype_coph,
    file      = file.path(output_folder,
                          paste0("StructureMap_table_", sample_name,
                                 ".csv")),
    row.names = TRUE
  )

  # --- Parallel processing of genes ---
  message(sprintf(
    "--- Initializing %d workers to process %d genes ---",
    n_jobs, length(genes)
  ))
  plan(multisession, workers = n_jobs)

  results_list <- future_lapply(
    X                = genes,
    FUN              = process_gene,
    coords_df        = coords,
    cell_metadata_df = cell_metadata_clean,
    celltype_coph_df = celltype_coph,
    coph_method      = coph_method,
    future.seed      = TRUE
  )
  names(results_list) <- genes

  # --- Combine results and save to file ---
  message("--- Combining and saving final results ---")

  # Filter out any NULL results from failed genes
  results_list <- results_list[!sapply(results_list, is.null)]

  # Efficiently combine list of data frames into one
  final_df <- data.table::rbindlist(results_list, idcol = "gene")
  setDF(final_df)
  rownames(final_df) <- final_df$gene
  final_df$gene    <- NULL

  # Write the final combined table
  output_csv <- file.path(
    output_folder,
    paste0("t_and_c_result_", sample_name, ".csv")
  )
  write.csv(final_df, file = output_csv, row.names = TRUE)

  message(paste("[DONE] Outputs saved to:", output_folder))
}
