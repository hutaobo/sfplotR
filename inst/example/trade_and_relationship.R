install.packages("comtradr")
library(comtradr)

# 获取美国自德国、法国、日本、墨西哥的年度进口额
df <- ct_get_data(
  reporter = "USA", partner = c("DEU","FRA","JPN","MEX"),
  commodity_code = "TOTAL", start_date = 2018, end_date = 2023,
  flow_direction = "import"
)

install.packages("tradestatistics")
library(tradestatistics)




install.packages(c("haven", "dplyr", "tidyr"))
library(haven)   # 读取 .sav
library(dplyr)   # 数据操作
library(tidyr)   # 数据重塑

raw <- read_sav(
  "C:/Users/taobo.hu/Downloads/Pew-Research-Center-Global-Attitudes-Spring-2022-Survey-Data/Pew Research Center Global Attitudes Spring 2022 Dataset.sav"
)
