# 解析開始 --------------------------------------------------------------------
rm(list=ls())
library(ggplot2)
library(tidyverse)
library(patchwork)
library(dplyr)
library(zoo)
##解析
subject <- "s2"
mw_filename <- paste(subject,"_mw_all.csv",sep="")
rf_filename <- paste(subject,"_rf_all.csv",sep="")
mw_dat <- read.csv(mw_filename)
rf_dat <- read.csv(rf_filename)
##図形描写
mw_median <- mw_dat %>%
  group_by(time) %>%
  summarise(
    median_pupil = median(pupil, na.rm = TRUE),
    sd_pupil = sd(pupil, na.rm = TRUE),
    n_trial = n_distinct(trial),
    .groups = "drop"
  )

rf_median <- rf_dat %>%
  group_by(time) %>%
  summarise(
    median_pupil = median(pupil, na.rm = TRUE),
    sd_pupil = sd(pupil, na.rm = TRUE),
    n_trial = n_distinct(trial),
    .groups = "drop"
  )

model_mw <- lm(mw_median$median_pupil ~ mw_median$time, data = mw_median)
model_rf <- lm(rf_median$median_pupil ~ rf_median$time, data = rf_median)

y_range <- range(
  c(mw_median$median_pupil,
    rf_median$median_pupil),
  na.rm = TRUE
)

{MW_median <- ggplot(mw_median, aes(x = time, y = median_pupil)) +
    geom_line() +
    coord_cartesian(ylim = y_range) +
    labs(
      x = "Time (s)",
      y = "median pupil",
      title = "MW_medain"
    ) +
    geom_smooth(method = "lm", se = TRUE, color = "red") +
    theme_bw()
  
  RF_median <- ggplot(rf_median, aes(x = time, y = median_pupil,)) +
    geom_line() +
    coord_cartesian(ylim = y_range) +
    labs(
      x = "Time (s)",
      y = "median pupil",
      title = "RF_medain"
    ) +
    geom_smooth(method = "lm", se = TRUE, color = "blue") +
    theme_bw()
  
  RF_median + MW_median}
##平均を使ったノンパラ
mw_mean_trial <- mw_dat %>%
  group_by(trial) %>%
  summarise(
    mean_pupil_trial = mean(pupil, na.rm = TRUE),
    sd_pupil_trial = sd(pupil, na.rm = TRUE),
    n_trial = n_distinct(trial),
    .groups = "drop"
  )

rf_mean_trial <- rf_dat %>%
  group_by(trial) %>%
  summarise(
    mean_pupil_trial = mean(pupil, na.rm = TRUE),
    sd_pupil_trial = sd(pupil, na.rm = TRUE),
    n_trial = n_distinct(trial),
    .groups = "drop"
  )

##これいらんかも。AIはデフォルト推奨、これは吉史先生推奨
library(exactRankTests)
wilcox.exact(mw_mean_trial$mean_pupil_trial,rf_mean_trial$mean_pupil_trial)
wilcox.exact(mw_mean_trial$sd_pupil_trial,rf_mean_trial$sd_pupil_trial)
t.test(mw_mean_trial$mean_pupil_trial,rf_mean_trial$mean_pupil_trial)

#瞬目持続時間を出す。下処理内のデータフレームに手を出すからエラーの原因かも
rf_dat <- rf_dat %>%
  arrange(trial, time) %>%
  group_by(trial) %>%
  mutate(
    blink = !is.na(for_blink_duration),
    blink_count = cumsum(
      blink != lag(blink, default = first(blink))
    )
  ) %>%
  ungroup()

mw_dat <- mw_dat %>%
  arrange(trial, time) %>%
  group_by(trial) %>%
  mutate(
    blink = !is.na(for_blink_duration),
    blink_count = cumsum(
      blink != lag(blink, default = first(blink))
    )
  ) %>%
  ungroup()


rf_blink <- rf_dat %>%
  arrange(trial,time_rel) %>%
  group_by(trial) %>%
  summarise(blink_sample_counts=(sum(is.na(for_blink_duration))-3)*4,
            nmber_of_blink = max(blink_count)/2) %>%
  mutate(blink_duration_rate = (blink_sample_counts/9000)*100)

mw_blink <- mw_dat %>%
  arrange(trial,time_rel) %>%
  group_by(trial) %>%
  summarise(blink_sample_counts=(sum(is.na(for_blink_duration))-3)*4,
            number_of_blink = max(blink_count)/2) %>%
  mutate(blink_duration_rate = (blink_sample_counts/9000)*100)

#engbertの適応範囲を変えた場合
rf_dat <- rf_dat %>%
  arrange(trial, time) %>%
  group_by(trial) %>%
  mutate(
    blink_count = cumsum(
      blink_duration != lag(blink_duration, default = first(blink_duration))
    )
  ) %>%
  ungroup()

mw_dat <- mw_dat %>%
  arrange(trial, time) %>%
  group_by(trial) %>%
  mutate(
    blink_count = cumsum(
      blink_duration != lag(blink_duration, default = first(blink_duration))
    )
  ) %>%
  ungroup()


rf_blink <- rf_dat %>%
  arrange(trial,time_rel) %>%
  group_by(trial) %>%
  summarise(blink_sample_counts=(sum(is.na(for_blink_duration))-3)*4,
            nmber_of_blink = max(blink_count)/2) %>%
  mutate(blink_duration_rate = (blink_sample_counts/9000)*100)

mw_blink <- mw_dat %>%
  arrange(trial,time_rel) %>%
  group_by(trial) %>%
  summarise(blink_sample_counts=(sum(is.na(for_blink_duration))-3)*4,
            number_of_blink = max(blink_count)/2) %>%
  mutate(blink_duration_rate = (blink_sample_counts/9000)*100)

#瞬目持続時間の平均、SD
data.frame(
  RF = c(
    Mean = mean(rf_blink$blink_sample_counts),
    SD = sd(rf_blink$blink_sample_counts)
  ),
  MW = c(
    Mean = mean(mw_blink$blink_sample_counts) ,
    SD = sd(mw_blink$blink_sample_counts)
  )
)

wilcox.exact(mw_blink$blink_sample_counts,rf_blink$blink_sample_counts)

##blink rate
data.frame(
  RF = c(
    minutes_blink_rate = mean(rf_blink$nmber_of_blink)*6.67,
    SD = sd(rf_blink$nmber_of_blink)
  ),
  MW = c(
    minutes_blink_rate = mean(mw_blink$number_of_blink)*6.67,
    SD = sd(mw_blink$number_of_blink)
  )
)

wilcox.exact(mw_blink$number_of_blink,rf_blink$nmber_of_blink)

##水平座標　平均
mw_gaze_mean <- mw_dat %>%
  group_by(trial) %>%
  summarise(
    mean_gaze_x = mean(x, na.rm = TRUE),
    mean_gaze_y = mean(y, na.rm = TRUE),
    sd_gaze_x = sd(x, na.rm = TRUE),
    sd_gaze_y = sd(y, na.rm = TRUE),
    n_trial = n_distinct(trial),
    .groups = "drop"
  )

rf_gaze_mean <- rf_dat %>%
  group_by(trial) %>%
  summarise(
    mean_gaze_x = mean(x, na.rm = TRUE),
    mean_gaze_y = mean(y, na.rm = TRUE),
    sd_gaze_x = sd(x, na.rm = TRUE),
    sd_gaze_y = sd(y, na.rm = TRUE),
    n_trial = n_distinct(trial),
    .groups = "drop"
  )

data.frame(
  RF = c(
    Gaze_X_Pos.Mean = mean(rf_gaze_mean$mean_gaze_x),
    Gaze_Y_Pos.Mean = mean(rf_gaze_mean$mean_gaze_y),
    Gaze_X_Pos.std.Dev = sd(rf_gaze_mean$mean_gaze_x),
    Gaze_Y_Pos.std.Dev = sd(rf_gaze_mean$mean_gaze_y)
  ),
  MW = c(
    Gaze_X_Pos.Mean = mean(mw_gaze_mean$mean_gaze_x),
    Gaze_Y_Pos.Mean = mean(mw_gaze_mean$mean_gaze_y),
    Gaze_X_Pos.std.Dev = sd(mw_gaze_mean$mean_gaze_x),
    Gaze_Y_Pos.std.Dev = sd(mw_gaze_mean$mean_gaze_y)
  )
)

wilcox.exact(mw_gaze_mean$mean_gaze_x,rf_gaze_mean$mean_gaze_x)
wilcox.exact(mw_gaze_mean$sd_gaze_x,rf_gaze_mean$sd_gaze_x)

wilcox.exact(mw_gaze_mean$mean_gaze_y,rf_gaze_mean$mean_gaze_y)
wilcox.exact(mw_gaze_mean$sd_gaze_y,rf_gaze_mean$sd_gaze_y)

