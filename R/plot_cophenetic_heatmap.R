plot_cophenetic_heatmap <- function(matrix,
                                    matrix_name = NULL,   # used to distinguish "row_coph" and "col_coph"
                                    output_dir = NULL,
                                    output_filename = NULL,
                                    figsize = c(9, 9),
                                    cellwidth = 11,
                                    cmap = "RdBu",
                                    linewidths = 0.5,     # pheatmap does not allow setting cell border width directly
                                    annot = TRUE,
                                    sample = "Sample",
                                    xlabel = NULL,
                                    ylabel = NULL,
                                    show_dendrogram = TRUE,
                                    return_xlabels = FALSE,
                                    parse_label = FALSE) {   # new parameter

  # 1) Set output directory
  if (is.null(output_dir)) {
    output_dir <- getwd()
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Set title, default filename, and axis labels based on matrix_name
  if (!is.null(matrix_name) && matrix_name == "row_coph") {
    title_str <- paste("Searcher's D score of", sample)
    default_filename <- paste0("Searcher's D score_of_", sample, ".pdf")
    xlabel <- "Searcher"
    ylabel <- "Searcher"
  } else if (!is.null(matrix_name) && matrix_name == "col_coph") {
    title_str <- paste("Findee's D score of", sample)
    default_filename <- paste0("Findee's D score_of_", sample, ".pdf")
    xlabel <- "Findee"
    ylabel <- "Findee"
  } else {
    title_str <- paste("SFplot of", sample)
    default_filename <- paste0("SFplot of ", sample, ".pdf")
    if (is.null(xlabel)) {
      xlabel <- "Findee"
    }
    if (is.null(ylabel)) {
      ylabel <- "Searcher"
    }
  }

  # If output_filename is not specified, use default
  if (is.null(output_filename)) {
    output_filename <- default_filename
  }
  output_file <- file.path(output_dir, output_filename)

  # Determine clustering based on show_dendrogram
  cluster_rows <- show_dendrogram
  cluster_cols <- show_dendrogram

  # Set color palette, default RdBu (requires RColorBrewer)
  if (cmap == "RdBu") {
    if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
      stop("Please install the RColorBrewer package first")
    }
    my_palette <- colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(100)
  } else {
    # If using other palettes, adjust manually
    my_palette <- colorRampPalette(colors = c("blue", "white", "red"))(100)
  }

  if(parse_label == TRUE) {
    label_parse <- function(breaks) {
      parse(text = latex2exp::TeX(breaks))
    }
    row_labels <- label_parse(rownames(matrix))
    col_labels <- label_parse(colnames(matrix))
  } else {
    row_labels <- rownames(matrix)
    col_labels <- colnames(matrix)
  }

  # Open PDF device; note pheatmap width/height are in inches
  library(showtext)
  font_add_google("Roboto Mono", "roboto mono")
  showtext_auto()

  cairo_pdf(file = output_file, width = figsize[1], height = figsize[2], family = "roboto mono")
  p <- pheatmap::pheatmap(matrix,
                          clustering_method = "average",
                          clustering_distance_rows = "euclidean",
                          clustering_distance_cols = "euclidean",
                          color = my_palette,
                          cluster_rows = cluster_rows,
                          cluster_cols = cluster_cols,
                          border_color = "white",
                          show_rownames = TRUE,
                          show_colnames = TRUE,
                          labels_row = row_labels,
                          labels_col = col_labels,
                          main = title_str,
                          cellwidth = cellwidth,    # adjust cell width as needed (pixels)
                          cellheight = cellwidth)   # adjust cell height as needed (pixels)
  print(p)
  dev.off()

  # If return_xlabels is TRUE, extract x-axis labels and return after closing device
  if(return_xlabels) {
    # Here assuming x-axis labels are stored in a grob named "col_names"
    x_labels_grob <- p$gtable$grobs[[which(p$gtable$layout$name == "col_names")]]
    x_labels <- x_labels_grob$label
    return(x_labels)
  }

  message("Heatmap saved to: ", output_file)
}
