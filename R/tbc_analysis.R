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

# --- 2. WORKER FUNCTION ------------------------------------------------------

#' Compute cophenetic distances for a single gene vs. cell types
process_gene <- function(gene,
                         coords_df,
                         cell_metadata_df,
                         celltype_coph_df,
                         coph_method) {
  tryCatch({
    gene_coords <- coords_df[coords_df$feature_name == gene, ]

    # ── Case 1: gene 没坐标 ────────────────────────────────────────────────
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

    # ── Case 2: 重新计算 gene ↔︎ cell types ───────────────────────────────
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

  # ── 预处理 ────────────────────────────────────────────────────────────────
  coords <- coords[!grepl("NegControl|Unassigned",
                          coords$feature_name, ignore.case = TRUE), ]
  cell_metadata_clean <- cell_metadata[, c("x", "y", "celltype")]
  genes <- sort(unique(coords$feature_name))

  # ── 计算 celltype↔︎celltype 基准距离 ──────────────────────────────────────
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

  # ── 并行处理各基因 ───────────────────────────────────────────────────────
  message(sprintf("--- Forking %d workers to process %d genes ---",
                  n_jobs, length(genes)))

  handlers("progress")  # CLI 进度条
  future::plan("multicore", workers = n_jobs)
  on.exit(future::plan("sequential"), add = TRUE)  # 运行结束后恢复单线程

  results_list <- with_progress({
    p <- progressr::progressor(steps = length(genes))

    future_lapply(
      X   = genes,
      FUN = function(gene,
                     coords_df,
                     cell_metadata_df,
                     celltype_coph_df,
                     coph_method) {
        p(sprintf("Processed %s", gene))
        process_gene(gene,
                     coords_df,
                     cell_metadata_df,
                     celltype_coph_df,
                     coph_method)
      },
      coords_df        = coords,
      cell_metadata_df = cell_metadata_clean,
      celltype_coph_df = celltype_coph,
      coph_method      = coph_method,
      future.seed      = TRUE   # reproducible RNG per worker
    )
  })
  names(results_list) <- genes

  # ── 汇总并保存 ────────────────────────────────────────────────────────────
  message("--- Combining and saving final results ---")
  results_list <- results_list[!sapply(results_list, is.null)]

  final_df <- data.table::rbindlist(results_list, idcol = "gene")
  setDF(final_df)
  rownames(final_df) <- final_df$gene
  final_df$gene <- NULL

  write.csv(final_df,
            file = file.path(output_folder,
                             paste0("t_and_c_result_", sample_name, ".csv")),
            row.names = TRUE)
  message(paste("[DONE] Outputs saved to:", output_folder))
}

# --- End of script -----------------------------------------------------------
