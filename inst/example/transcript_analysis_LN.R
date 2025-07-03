library(readr)
library(dplyr)

cells = read_csv('/Volumes/mn-moldia/long/10X_datasets/Xenium/Xenium_5K/Xenium_Prime_Human_Lymph_Node_Reactive_FFPE_outs/cells.csv.gz')

celltype <- read_csv('/Volumes/mn-moldia/long/10X_datasets/Xenium/Xenium_5K/Xenium_Prime_Human_Lymph_Node_Reactive_FFPE_outs/Xenium_Prime_Human_Lymph_Node_Reactive_FFPE_cell_types.csv')

cells <- cells %>%
  left_join(celltype %>% select(cell_id, group), by = "cell_id")
cells$group <- paste0("$\\textbf{", gsub(" ", "~", cells$group), "}$")

source("/Volumes/mn-moldia/long/projects/NPC/rcode/compute_cophenetic_distances_from_df.R")
source("/Volumes/mn-moldia/long/projects/Pakagedevelopment/CellGPS/sfplotR/R/plot_cophenetic_heatmap.R")

# 计算 cophenetic 距离
result <- compute_cophenetic_distances_from_df(
  cells,
  cluster_col = "group",
  x_col = "x_centroid",
  y_col = "y_centroid",
  method = "average"
)

# 绘制 cophenetic 热图
plot_cophenetic_heatmap(
  result[['row_cophenetic_df']],
  figsize = c(15, 15),
  cellwidth = 13,
  matrix_name = "row_coph",
  output_dir = '~/Downloads/',
  sample = "LN"
)


# 读取transcript ------------------------------------------------------------

library(arrow)
library(dplyr)

# 打开 Parquet 数据集（不会立刻加载所有数据）
ds <- open_dataset("/Volumes/mn-moldia/long/10X_datasets/Xenium/Xenium_5K/Xenium_Prime_Human_Lymph_Node_Reactive_FFPE_outs/transcripts.parquet", format = "parquet")

# 通过 dplyr 操作进行筛选或选择部分列，例如只读取部分列和符合条件的行
freq_table <- ds %>%
  count(feature_name) %>%   # 计算每个 feature_name 的出现次数
  collect()                 # 收集结果到内存

freq_table = subset(freq_table, freq_table$feature_name start_with 'CCR  / CCL/ CXCR/ CXCL')

subset_ds <- ds %>%
  select(x_location, y_location, z_location, cell_id, feature_name) %>%
  semi_join(freq_table, by = "feature_name")

df <- collect(subset_ds)
df$feature_name <- paste0("$\\textit{", df$feature_name, "}$")


# 合并cell和transcript的StructureMap ------------------------------------------

matrix <- rbind(setNames(cells[, c("x_centroid","y_centroid","cell_id","group")], c("x","y","cell_id","group")),
                setNames(df[, c("x_location","y_location","cell_id","feature_name")], c("x","y","cell_id","group")))

source("/Volumes/mn-moldia/long/projects/NPC/rcode/compute_cophenetic_distances_from_df.R")
source("/Volumes/mn-moldia/long/projects/Pakagedevelopment/SpatialMap/sfplotR/R/plot_cophenetic_heatmap.R")

# 计算 cophenetic 距离
result <- compute_cophenetic_distances_from_df(
  matrix,
  cluster_col = "group",
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
  output_dir = '~/Downloads/',
  sample = "5K_LN"
)
