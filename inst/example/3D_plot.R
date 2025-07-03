setwd("Y:/long/publication_datasets/Vannan_2023_Lung_Fibrosis/Rcode/sfplotR/example")

# 加载必要的 R 包
library(dplyr)
library(reshape2)
library(ggplot2)

# 使用正斜杠（推荐）
data <- readRDS("Y:/taobo/Downloads/TSOHP/GSE250346_Seurat_GSE250346_CORRECTED_SEE_RDS_README_082024.rds")

metadata <- data@meta.data

# 定义样本列表及其对应名称（注意每个元素之间用逗号分隔）
sample_list <- c("VUHD090", "VUHD049", "VUHD038", "THD0008", "THD0011",
                 "VUHD069", "VUHD095", "VUHD113", "VUHD116A", "VUHD116B",
                 "TILD080LA", "TILD028LA", "TILD167LA", "TILD113LA", "TILD111LA",
                 "VUILD49LA", "TILD130LA", "VUILD110LA", "TILD117LA", "VUILD78LA",
                 "VUILD91LA", "VUILD48LA2", "VUILD102LA", "VUILD96LA", "VUILD48LA1",
                 "TILD315MA", "VUILD58MA", "TILD049MA", "TILD299MA", "TILD117MA2",
                 "VUILD141MA", "VUILD142MA", "VUILD106MA", "VUILD115MA", "TILD117MA1",
                 "TILD175MA", "VUILD78MA", "VUILD91MA", "VUILD104MA1", "VUILD105MA2",
                 "VUILD102MA", "VUILD107MA", "VUILD96MA", "VUILD104MA2", "VUILD105MA1")

# 初始化一个空列表，用于存放每个样本处理后的数据框
result_list <- list()

# 遍历每个样本
for (sample in sample_list) {
  # 根据样本名称筛选出对应的元数据（假设 metadata 数据框已定义）
  submeta <- metadata[metadata$sample == sample, ]

  # 调用自定义函数计算 cophenetic 距离
  # 注意：compute_cophenetic_distances_from_df 函数应已在环境中定义
  result <- compute_cophenetic_distances_from_df(submeta,
                                                 cluster_col = "final_CT",
                                                 x_col = "x_centroid",
                                                 y_col = "y_centroid")

  # 将结果中的 cophenetic 数据转换为矩阵格式，保留原有的行名和列名
  row_coph <- as.matrix(result$row_cophenetic_df)

  # 将矩阵转换为数据框，同时将矩阵的行名保存为新的一列 'row'
  df <- as.data.frame(row_coph)
  df$row <- rownames(row_coph)

  # 将数据框从宽格式转换为长格式，保留 'row' 列不变
  # 转换后 'column' 列为原矩阵的列名，'value' 列为对应的数值
  df_long <- melt(df, id.vars = "row", variable.name = "column", value.name = "value")

  # 添加当前样本的名称信息
  df_long$sample <- sample

  # 将当前样本的长格式数据存入结果列表
  result_list[[length(result_list) + 1]] <- df_long
}

# 合并所有样本的结果数据框（按行合并）
final_df <- do.call(rbind, result_list)

# 将最终合并的数据框导出为 CSV 文件，不保留行名
write.csv(final_df, file = "cophenetic_distances_searcher_D_score_in_all_samples.csv", row.names = FALSE)

# 读取 CSV 数据
final_df <- read.csv("cophenetic_distances_searcher_D_score_in_all_samples.csv", stringsAsFactors = FALSE)
View(final_df)

# 定义各组对应的样本列表
hd_samples <- c("VUHD090", "VUHD049", "VUHD038", "THD0008", "THD0011", "VUHD069", "VUHD095", "VUHD113", "VUHD116A", "VUHD116B")
la_samples <- c("TILD080LA", "TILD028LA", "TILD167LA", "TILD113LA", "TILD111LA", "VUILD49LA", "TILD130LA", "VUILD110LA", "TILD117LA", "VUILD78LA", "VUILD91LA", "VUILD48LA2", "VUILD102LA", "VUILD96LA", "VUILD48LA1")
ma_samples <- c("TILD315MA", "VUILD58MA", "TILD049MA", "TILD299MA", "TILD117MA2", "VUILD141MA", "VUILD142MA", "VUILD106MA", "VUILD115MA", "TILD117MA1", "TILD175MA", "VUILD78MA", "VUILD91MA", "VUILD104MA1", "VUILD105MA2", "VUILD102MA", "VUILD107MA", "VUILD96MA", "VUILD104MA2", "VUILD105MA1")

# 根据 sample 列创建一个新的分组列
final_df <- final_df %>%
  mutate(group = case_when(
    sample %in% hd_samples ~ "HD",
    sample %in% la_samples ~ "LA",
    sample %in% ma_samples ~ "MA",
    TRUE ~ "Other"
  ))

# 查看分组情况
print(table(final_df$group))

# 1. 安装并加载必要的包
if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")
if (!require("caret")) install.packages("caret")
if (!require("randomForest")) install.packages("randomForest")

library(dplyr)
library(tidyr)
library(caret)
library(randomForest)

# 2. 数据预处理：生成统一特征并删除重复行，只保留 value, sample, group, feature 四列
final_df_clean <- final_df %>%
  mutate(feature = ifelse(row < column,
                          paste(row, column, sep = "_"),
                          paste(column, row, sep = "_"))) %>%
  distinct(sample, feature, .keep_all = TRUE) %>%
  select(value, sample, group, feature)

# 3. 转换为宽格式，每个样本一行，每个 feature 一列；对于重复取均值
df_wide <- final_df_clean %>%
  pivot_wider(names_from = feature, values_from = value, values_fn = mean)

# 4. 删除包含 NA 值的列
df_wide <- df_wide %>% select_if(~ !any(is.na(.)))

# 5. 由于部分变量名中可能有特殊字符或空格，利用 make.names() 将列名转换为合法变量名
names(df_wide) <- make.names(names(df_wide))

# 检查数据结构
str(df_wide)

# 6. 将 group 转换为因子
df_wide$group <- as.factor(df_wide$group)

# 7. 分割数据集为训练集和测试集（70% 训练，30% 测试）
set.seed(123)
trainIndex <- createDataPartition(df_wide$group, p = 0.7, list = FALSE)
trainData <- df_wide[trainIndex, ]
testData  <- df_wide[-trainIndex, ]

# 8. 利用 RFE 进行特征选择
# 设定 RFE 控制参数，使用随机森林作为评估函数
control_rfe <- rfeControl(functions = rfFuncs, method = "cv", number = 5)

# 训练时排除 group 和 sample 列（sample 只是标识符）
train_predictors <- trainData %>% select(-group, -sample)
train_response <- trainData$group

set.seed(123)
rfe_results <- rfe(x = train_predictors,
                   y = train_response,
                   sizes = seq(1, ncol(train_predictors)),
                   rfeControl = control_rfe)

# 查看 RFE 结果及被选中的特征
print(rfe_results)
selected_features <- predictors(rfe_results)
print(selected_features)

# 9. 构造最终模型的公式
# 由于已经用 make.names() 将变量名转换为合法形式，这里直接使用
formula_str <- paste("group ~", paste(selected_features, collapse = " + "))
final_formula <- as.formula(formula_str)
print(final_formula)

# 10. 利用 caret 建立随机森林模型
control <- trainControl(method = "cv", number = 5)
set.seed(123)
model_rf <- train(final_formula,
                  data = trainData %>% select(-sample),  # 排除 sample 列
                  method = "rf",
                  trControl = control)

# 输出模型信息及调参结果
print(model_rf)

# 11. 模型评估：在测试集上进行预测并输出混淆矩阵
predictions <- predict(model_rf, newdata = testData)
confMat <- confusionMatrix(predictions, testData$group)
print(confMat)
