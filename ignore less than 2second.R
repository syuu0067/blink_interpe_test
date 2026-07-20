#除外用コード
rm(list=ls())
library(dplyr)

#データ取込
subject <- "s2"
data_num <- 1
eject_num <- 3
data_num <- sprintf("%02d",data_num)
file_name <- paste(subject,data_num,"_mw_dat.csv", sep="")
data <- read.csv(file_name)

#実行列
ejected_dat <- data %>%
  dplyr::filter(trial !=eject_num)

##保存
write.csv(ejected_dat,file_name,row.names = FALSE)

