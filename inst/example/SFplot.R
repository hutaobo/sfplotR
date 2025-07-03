setwd("Y:/long/publication_datasets/Vannan_2023_Lung_Fibrosis/Rcode/sfplotR/example")

# 加载必要的 R 包
library(dplyr)
library(ggplot2)

# 使用正斜杠（推荐）
data <- readRDS("Y:/taobo/Downloads/TSOHP/GSE250346_Seurat_GSE250346_CORRECTED_SEE_RDS_README_082024.rds")

metadata <- data@meta.data
samplelist <- unique(metadata$sample)

for (sample in samplelist) {
  submeta <- metadata[metadata$sample == sample, ]
  result <- compute_cluster_average_nn_distance_matrix(submeta,
                                                       cluster_col = "final_CT",
                                                       x_col = "x_centroid",
                                                       y_col = "y_centroid")
  plot_cophenetic_heatmap(
    matrix = result,
    matrix_name = "",
    figsize = c(12, 12),
    output_dir = "./output",
    sample = sample
  )
}
