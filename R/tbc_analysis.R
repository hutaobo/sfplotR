#!/usr/bin/env Rscript
# tbc_analysis.R — Transcript-by-cell spatial analysis (Linux-optimised)

# --- 1. SCRIPT SETUP ---------------------------------------------------------
suppressPackageStartupMessages({
  library(future)
  library(future.apply)
  library(dplyr)
  library(data.table)
  library(tools)
  library(progressr)  # live progress bar
})
# 将上限设置为 6 GiB (6 * 1024^3 bytes)
options(future.globals.maxSize = 100 * 1024^3)

# --- 2. WORKER FUNCTION ------------------------------------------------------

#' Compute cophenetic distances for a single gene vs. cell types
process_gene <- function(gene,
                         coords_df,
                         cell_metadata_df,
                         celltype_coph_df,
                         coph_method) {
  tryCatch({
    gene_coords <- coords_df[coords_df$feature_name == gene, ]

    # ── Case 1: gene has no coordinates ──────────────────────────────────
    if (nrow(gene_coords) == 0) {
      if (!gene %in% rownames(celltype_coph_df)) {
        res <- as.data.frame(matrix(NA_real_,
                                    nrow = 1,
                                    ncol = ncol(celltype_coph_df)))
        colnames(res) <- colnames(celltype_coph_df)
        rownames(res) <- gene
        return(res)
      }
      res <- celltype_coph_df[gene, , drop = FALSE]
      res <- res[, colnames(res) != gene, drop = FALSE]
      return(res)
    }

    # ── Case 2: Recalculate gene ↔︎ cell types distances ───────────────
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

    gene_coph_row <- as.matrix(coph_results$row_cophenetic_df)[gene, ,
                                                              drop = FALSE]
    gene_coph_row <- gene_coph_row[, colnames(gene_coph_row) != gene,
                                   drop = FALSE]
    gene_coph_row
  }, error = function(e) {
    warning(sprintf("Worker failed on gene %s: %s", gene, e$message))
    NULL
  })
}

# --- 3. MAIN ANALYSIS FUNCTION ----------------------------------------------

#' Transcript-by-cell spatial analysis with progress bar & forked parallelism
#'
#' @param n_jobs Number of parallel workers (<= `parallel::detectCores()`)
transcript_by_cell_analysis <- function(cell_metadata,
                                        coords,
                                        sample_name   = NULL,
                                        output_folder = NULL,
                                        coph_method   = "average",
                                        n_jobs        = 8) {

  if (is.null(sample_name)) sample_name <- "sample"
  if (is.null(output_folder))
    output_folder <- file.path(".", paste0("t_by_c_", sample_name))
  if (!dir.exists(output_folder))
    dir.create(output_folder, recursive = TRUE)

  message("--- Loading and preparing data ---")

  # ── Pre-processing ───────────────────────────────────────────────────
  coords <- coords[!grepl("NegControl|Unassigned",
                          coords$feature_name, ignore.case = TRUE), ]
  cell_metadata_clean <- cell_metadata[, c("x", "y", "celltype")]
  genes <- sort(unique(coords$feature_name))

  # ── Compute baseline celltype↔︎celltype distances ─────────────────────
  message("--- Computing initial cell-type StructureMap ---")
  celltype_coph <- compute_cophenetic_distances_from_df(
    df          = cell_metadata_clean,
    cluster_col = "celltype",
    x_col       = "x",
    y_col       = "y",
    method      = coph_method
  )$row_cophenetic_df

  plot_cophenetic_heatmap(
    matrix          = celltype_coph,
    matrix_name     = "row_coph",
    output_dir      = output_folder,
    output_filename = paste0("StructureMap_of_", sample_name, ".pdf"),
    sample          = sample_name
  )

  write.csv(celltype_coph,
            file.path(output_folder,
                      paste0("StructureMap_table_", sample_name, ".csv")),
            row.names = TRUE)

  # ── Process each gene in parallel ─────────────────────────────────────
  message(sprintf("--- Forking %d workers to process %d genes ---",
                  n_jobs, length(genes)))

  handlers("progress")  # CLI progress bar
  future::plan("multicore", workers = n_jobs)
  on.exit(future::plan("sequential"), add = TRUE)  # Restore sequential processing on exit

  results_list <- with_progress({
    p <- progressr::progressor(steps = length(genes))

    # --- START: MODIFIED SECTION ---
    # The key change is to use `future.globals` to explicitly pass large
    # objects and other necessary functions/variables to the workers.
    # This prevents `future` from bundling the large data inside the
    # function's environment, thus avoiding the size limit error.
    future.apply::future_lapply(
      X   = genes,
      FUN = function(gene) {
        # This function now only takes 'gene' as a direct argument.
        # All other variables (`p`, `process_gene`, `coords_df`, etc.)
        # are found in the global environment of the worker process.
        p(sprintf("Processed %s", gene))
        process_gene(gene,
                     coords_df,
                     cell_metadata_df,
                     celltype_coph_df,
                     coph_method)
      },
      future.seed      = TRUE,
      future.globals   = list(
          p = p,
          process_gene = process_gene,
          coords_df = coords,
          cell_metadata_df = cell_metadata_clean,
          celltype_coph_df = celltype_coph,
          coph_method = coph_method
      )
    )
    # --- END: MODIFIED SECTION ---
  })
  names(results_list) <- genes

  # ── Combine and save results ──────────────────────────────────────────
  message("--- Combining and saving final results ---")
  results_list <- results_list[!sapply(results_list, is.null)]

  final_df <- data.table::rbindlist(results_list, idcol = "gene")
  data.table::setDF(final_df)
  rownames(final_df) <- final_df$gene
  final_df$gene <- NULL

  write.csv(final_df,
            file = file.path(output_folder,
                             paste0("t_and_c_result_", sample_name, ".csv")),
            row.names = TRUE)
  message(paste("[DONE] Outputs saved to:", output_folder))
}

# --- End of script -----------------------------------------------------------
