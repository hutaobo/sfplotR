setwd("Y:/long/publication_datasets/Vannan_2023_Lung_Fibrosis/Rcode/sfplotR/example")

# 加载必要的 R 包
library(dplyr)
library(ggplot2)

# 使用正斜杠（推荐）
data <- readRDS("Y:/taobo/Downloads/TSOHP/GSE250346_Seurat_GSE250346_CORRECTED_SEE_RDS_README_082024.rds")

metadata <- data@meta.data
samplelist <- unique(metadata$sample)

# 定义各组对应的样本列表
hd_samples <- c("VUHD090", "VUHD049", "VUHD038", "THD0008", "THD0011", "VUHD069", "VUHD095", "VUHD113", "VUHD116A", "VUHD116B")
la_samples <- c("TILD080LA", "TILD028LA", "TILD167LA", "TILD113LA", "TILD111LA", "VUILD49LA", "TILD130LA", "VUILD110LA", "TILD117LA", "VUILD78LA", "VUILD91LA", "VUILD48LA2", "VUILD102LA", "VUILD96LA", "VUILD48LA1")
ma_samples <- c("TILD315MA", "VUILD58MA", "TILD049MA", "TILD299MA", "TILD117MA2", "VUILD141MA", "VUILD142MA", "VUILD106MA", "VUILD115MA", "TILD117MA1", "TILD175MA", "VUILD78MA", "VUILD91MA", "VUILD104MA1", "VUILD105MA2", "VUILD102MA", "VUILD107MA", "VUILD96MA", "VUILD104MA2", "VUILD105MA1")

# 定义样本列表及其对应名称
sample_lists <- list(hd_samples, la_samples, ma_samples)
sample_names <- c("Healthy", "Less", "More")

# 遍历每个样本列表及其名称
for (i in 1:1) {
  sample_list <- sample_lists[[i]]
  sample_name <- sample_names[i]

  # 用于存储当前样本组的计算结果
  row_cophenetic_list <- list()
  col_cophenetic_list <- list()

  # 遍历当前样本列表中的每个样本
  for (sample in sample_list) {
    submeta <- metadata[metadata$sample == sample, ]
    result <- compute_cophenetic_distances_from_df(submeta,
                                                   cluster_col = "final_CT",
                                                   x_col = "x_centroid",
                                                   y_col = "y_centroid")
    row_coph <- as.matrix(result$row_cophenetic_df)
    col_coph <- as.matrix(result$col_cophenetic_df)

    # 添加到列表中
    row_cophenetic_list[[length(row_cophenetic_list) + 1]] <- row_coph
    col_cophenetic_list[[length(col_cophenetic_list) + 1]] <- col_coph
  }

  # 对齐所有矩阵：确保它们拥有相同的行和列标签
  all_row_index <- sort(unique(unlist(lapply(row_cophenetic_list, rownames))))
  all_row_columns <- sort(unique(unlist(lapply(row_cophenetic_list, colnames))))
  aligned_row_cophenetic_list <- lapply(row_cophenetic_list, function(df) {
    matrix <- matrix(NA, nrow=length(all_row_index), ncol=length(all_row_columns),
                     dimnames=list(all_row_index, all_row_columns))
    matrix[rownames(df), colnames(df)] <- df
    return(matrix)
  })

  all_col_index <- sort(unique(unlist(lapply(col_cophenetic_list, rownames))))
  all_col_columns <- sort(unique(unlist(lapply(col_cophenetic_list, colnames))))
  aligned_col_cophenetic_list <- lapply(col_cophenetic_list, function(df) {
    matrix <- matrix(NA, nrow=length(all_col_index), ncol=length(all_col_columns),
                     dimnames=list(all_col_index, all_col_columns))
    matrix[rownames(df), colnames(df)] <- df
    return(matrix)
  })

  # 计算平均值（忽略 NA）
  average_row_values <- apply(simplify2array(aligned_row_cophenetic_list), c(1,2),
                              function(x) mean(x, na.rm=TRUE))
  average_row_cophenetic <- as.data.frame(average_row_values)
  rownames(average_row_cophenetic) <- all_row_index
  colnames(average_row_cophenetic) <- all_row_columns

  average_col_values <- apply(simplify2array(aligned_col_cophenetic_list), c(1,2),
                              function(x) mean(x, na.rm=TRUE))
  average_col_cophenetic <- as.data.frame(average_col_values)
  rownames(average_col_cophenetic) <- all_col_index
  colnames(average_col_cophenetic) <- all_col_columns

  # 将 NA 替换为 1
  average_row_cophenetic[is.na(average_row_cophenetic)] <- 1
  average_col_cophenetic[is.na(average_col_cophenetic)] <- 1

  # 绘制 row_cophenetic 热图
  ordered_names <- plot_cophenetic_heatmap(
    matrix = average_row_cophenetic,
    matrix_name = "row_coph",
    figsize = c(12, 12),
    output_dir = "./multi_plot",
    sample = sample_name,
    show_dendrogram = TRUE,
    return_xlabels = TRUE
  )
}

# 遍历每个样本列表及其名称
for (i in 2:3) {
  sample_list <- sample_lists[[i]]
  sample_name <- sample_names[i]

  # 用于存储当前样本组的计算结果
  row_cophenetic_list <- list()
  col_cophenetic_list <- list()

  # 遍历当前样本列表中的每个样本
  for (sample in sample_list) {
    submeta <- metadata[metadata$sample == sample, ]
    result <- compute_cophenetic_distances_from_df(submeta,
                                                   cluster_col = "final_CT",
                                                   x_col = "x_centroid",
                                                   y_col = "y_centroid")
    row_coph <- as.matrix(result$row_cophenetic_df)
    col_coph <- as.matrix(result$col_cophenetic_df)

    # 添加到列表中
    row_cophenetic_list[[length(row_cophenetic_list) + 1]] <- row_coph
    col_cophenetic_list[[length(col_cophenetic_list) + 1]] <- col_coph
  }

  # 对齐所有矩阵：确保它们拥有相同的行和列标签
  all_row_index <- sort(unique(unlist(lapply(row_cophenetic_list, rownames))))
  all_row_columns <- sort(unique(unlist(lapply(row_cophenetic_list, colnames))))
  aligned_row_cophenetic_list <- lapply(row_cophenetic_list, function(df) {
    matrix <- matrix(NA, nrow=length(all_row_index), ncol=length(all_row_columns),
                     dimnames=list(all_row_index, all_row_columns))
    matrix[rownames(df), colnames(df)] <- df
    return(matrix)
  })

  all_col_index <- sort(unique(unlist(lapply(col_cophenetic_list, rownames))))
  all_col_columns <- sort(unique(unlist(lapply(col_cophenetic_list, colnames))))
  aligned_col_cophenetic_list <- lapply(col_cophenetic_list, function(df) {
    matrix <- matrix(NA, nrow=length(all_col_index), ncol=length(all_col_columns),
                     dimnames=list(all_col_index, all_col_columns))
    matrix[rownames(df), colnames(df)] <- df
    return(matrix)
  })

  # 计算平均值（忽略 NA）
  average_row_values <- apply(simplify2array(aligned_row_cophenetic_list), c(1,2),
                              function(x) mean(x, na.rm=TRUE))
  average_row_cophenetic <- as.data.frame(average_row_values)
  rownames(average_row_cophenetic) <- all_row_index
  colnames(average_row_cophenetic) <- all_row_columns

  average_col_values <- apply(simplify2array(aligned_col_cophenetic_list), c(1,2),
                              function(x) mean(x, na.rm=TRUE))
  average_col_cophenetic <- as.data.frame(average_col_values)
  rownames(average_col_cophenetic) <- all_col_index
  colnames(average_col_cophenetic) <- all_col_columns

  # 输出结果
  cat("Average Row Cophenetic Distance for current sample group:\n")
  print(average_row_cophenetic)
  cat("\nAverage Column Cophenetic Distance for current sample group:\n")
  print(average_col_cophenetic)

  # 将 NA 替换为 1
  average_row_cophenetic[is.na(average_row_cophenetic)] <- 1
  average_col_cophenetic[is.na(average_col_cophenetic)] <- 1

  average_row_cophenetic_ordered <- average_row_cophenetic[ordered_names, ordered_names]

  # 绘制 row_cophenetic 热图
  plot_cophenetic_heatmap(
    matrix = average_row_cophenetic_ordered,
    matrix_name = "row_coph",
    figsize = c(12, 12),
    output_dir = "./multi_plot",
    sample = sample_name,
    show_dendrogram = FALSE
  )
}
