plot_cophenetic_heatmap <- function(matrix,
                                    matrix_name = NULL,   # 用于区分 "row_coph" 和 "col_coph"
                                    output_dir = NULL,
                                    output_filename = NULL,
                                    figsize = c(9, 9),
                                    cellwidth = 11,
                                    cmap = "RdBu",
                                    linewidths = 0.5,     # pheatmap 中无法直接设置单元格边框宽度
                                    annot = TRUE,
                                    sample = "Sample",
                                    xlabel = NULL,
                                    ylabel = NULL,
                                    show_dendrogram = TRUE,
                                    return_xlabels = FALSE) {   # 新增参数

  # 1) 设置输出目录
  if (is.null(output_dir)) {
    output_dir <- getwd()
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # 根据 matrix_name 设置标题、默认文件名以及轴标签
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

  # 如果未指定输出文件名，则使用默认文件名
  if (is.null(output_filename)) {
    output_filename <- default_filename
  }
  output_file <- file.path(output_dir, output_filename)

  # 根据 show_dendrogram 参数确定是否进行聚类
  cluster_rows <- show_dendrogram
  cluster_cols <- show_dendrogram

  # 设置颜色调色板，默认为 RdBu（需依赖 RColorBrewer 包）
  if (cmap == "RdBu") {
    if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
      stop("请先安装 RColorBrewer 包")
    }
    my_palette <- colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(100)
  } else {
    # 若使用其它调色板，可自行调整
    my_palette <- colorRampPalette(colors = c("blue", "white", "red"))(100)
  }

  label_parse <- function(breaks) {
    parse(text = latex2exp::TeX(breaks))
  }
  row_labels <- label_parse(rownames(matrix))
  col_labels <- label_parse(colnames(matrix))

  # 打开 PDF 设备，注意 pheatmap 的 width/height 参数默认以英寸为单位
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
                          cellwidth = cellwidth,    # 根据需要调整单元格宽度（单位：像素）
                          cellheight = cellwidth)   # 根据需要调整单元格高度（单位：像素）
  print(p)
  dev.off()

  # 如果要求返回 x 轴标签，提取标签并关闭设备后返回
  if(return_xlabels) {
    # 这里假设 x 轴标签存储在名为 "col_labels" 的 grob 中
    x_labels_grob <- p$gtable$grobs[[which(p$gtable$layout$name == "col_names")]]
    x_labels <- x_labels_grob$label
    return(x_labels)
  }

  message("Heatmap saved to: ", output_file)
}
