library(readr)
X874141_transcripts <- read_csv("Y:/long/projects/Pakagedevelopment/SpatialMap/notebook/874141_transcripts.csv")


# 合并cell和transcript的StructureMap ------------------------------------------

# 读取 CSV 文件
transcript = subset(X874141_transcripts, grepl("^CXCL5", feature_name))
transcript$feature_name <- paste0("$\\textit{", transcript$feature_name, "}$")

cells = read.csv("Y:/long/projects/NPC/Ximum5000 group1 and group 2 raw data/874141_celltype.csv")
cells$feature_name <- paste0("$\\textbf{", gsub(" ", "~", cells$feature_name), "}$")

df = rbind(transcript, cells)

source("Y:/long/projects/NPC/rcode/compute_cophenetic_distances_from_df.R")
source("Y:/long/projects/Pakagedevelopment/SpatialMap/sfplotR/R/plot_cophenetic_heatmap.R")

# 计算 cophenetic 距离
result <- compute_cophenetic_distances_from_df(
  df,
  cluster_col = "feature_name",
  x_col = "x",
  y_col = "y",
  method = "average"
)

Sys.setlocale("LC_CTYPE", "en_US.UTF-8")

# 绘制 cophenetic 热图
plot_cophenetic_heatmap(
  result[['row_cophenetic_df']],
  figsize = c(15, 15),
  matrix_name = "row_coph",
  output_dir = './',
  sample = "NPC_CXCL5"
)

