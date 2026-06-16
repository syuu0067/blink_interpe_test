rm(list=ls())
Sys.setenv(EDFAPI = "C:/Program Files (x86)/SR Research/EyeLink/EDF_Access_API")
library(eyelinkReader)

blink_raw <- read_edf(
  "mwtest2.edf",
  import_samples = TRUE,
  start_marker = "REFOCUS",
  end_marker = "MW_REPORT"
)

library(dplyr)
library(zoo)
library(dplyr)

samples <- blink_raw$samples %>%
  select(
    trial,
    eye,
    time,
    time_rel,
    x = gxR,
    y = gyR,
    pupil = paR
  )

sr=250
screen_center_x <- 1920 / 2  
screen_center_y <- 1080 / 2  
samples$x <- samples$x - screen_center_x
samples$y <- samples$y - screen_center_y

# ---- monitor settings ----
diag_inch <- 23
aspect_w <- 16
aspect_h <- 9
view_dist_cm <- 50

res_x <- 1920
res_y <- 1080



samples_interp <- samples %>%
  arrange(trial, time_rel) %>%
  group_by(trial) %>%
  mutate(
    is_zero_gaze = is.na(x),
    is_zero_pupil = !is.na(pupil) & pupil == 0,
    
    na_group = cumsum(
      is_zero_gaze != lag(is_zero_gaze, default = first(is_zero_gaze))
    ),
    zero_group = cumsum(
      is_zero_pupil != lag(is_zero_pupil, default = first(is_zero_pupil))
    )
  ) %>%
  ungroup()

zero_segments <- samples_interp %>%
  filter(is_zero_pupil) %>%
  group_by(trial, zero_group) %>%
  summarise(
    zero_start = min(time_rel, na.rm = TRUE),
    zero_end = max(time_rel, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    bad_start = zero_start - 100,
    bad_end = zero_end + 300
  )

samples_flagged <- samples_interp %>%
  group_by(trial) %>%
  group_modify(~ {
    s <- .x
    z <- zero_segments %>% filter(trial == .y$trial)
    
    if (nrow(z) == 0) {
      s$blink_window <- FALSE
    } else {
      s$blink_window <- sapply(
        s$time_rel,
        function(t) any(t >= z$bad_start & t <= z$bad_end)
      )
    }
    
    s
  }) %>%
  ungroup()

samples_na <- samples_flagged %>%
  mutate(
    pupil_for_interp = if_else(blink_window, NA_real_, pupil)
  )

samples_interp_final <- samples_na %>%
  arrange(trial, time_rel) %>%
  group_by(trial) %>%
  mutate(
    pupil_interp = na.approx(
      pupil_for_interp,
      x = time_rel,
      na.rm = FALSE
    )
  ) %>%
  mutate(
    gaze_x = na.approx(
      x,
      x = time_rel,
      na.rm = FALSE
    )
  ) %>%
  mutate(
    gaze_y = na.approx(
      y,
      x = time_rel,
      na.rm = FALSE
    )
  ) %>%
  ungroup()

samples_interp_final <- samples_interp_final%>%
filter(!is.na(pupil_interp)) 

analysis_data <- samples_interp_final %>%
  arrange(trial, time_rel) %>%
  group_by(trial) %>%
  mutate(
    pupil_diff = pupil_interp - lag(pupil_interp)
  ) %>%
  ungroup()
  
pupil_velocity_average = mean(analysis_data$pupil_diff,na.rm=TRUE)

SD = sd(analysis_data$pupil_diff,na.rm = TRUE)
σ = 1.5

low_threshold = pupil_velocity_average - SD*σ
high_threshold = pupil_velocity_average + SD*σ


raw_artifacts <- which(analysis_data$pupil_diff<low_threshold | analysis_data$pupil_diff>high_threshold)
a <- diff(raw_artifacts)
artifacts_point <- sort(unique(c(raw_artifacts,raw_artifacts)))


artifacts=artifacts_point
for (i in length(a):1) {
  if(a[i]<=4 & a[i]>1)
    for (d in (a[i]-1):1) {
      artifacts <- append(artifacts,raw_artifacts[i]+d,after = i)
    }
}

artifacts <- sort(unique(artifacts))

analysis_artifact_interp <- analysis_data %>%
  mutate(
    artifact_velocity = row_number() %in% artifacts,
    pupil_for_artifact_interp = if_else(
      artifact_velocity,
      NA_real_,
      pupil_interp
    )
  ) %>%
  arrange(trial, time_rel) %>%
  group_by(trial) %>%
  mutate(
    pupil_interp2 = na.approx(
      pupil_for_artifact_interp,
      x = time_rel,
      na.rm = FALSE
    )
  ) %>%
  ungroup()


analysis_check <- analysis_artifact_interp %>%
  mutate(
    artifact_velocity = row_number() %in% artifacts
  )

blink_times <- analysis_check %>%
  filter(blink_window) %>%
  pull(time_rel)


gaze_velocity <- function(x, sampling_rate) {
  dt <- 1 / sampling_rate
  n <- length(x)
  v <- rep(NA_real_, n)
  
  for (i in 3:(n - 2)) {
    v[i] <- (x[i + 2] + x[i + 1] - x[i - 1] - x[i - 2]) / (6 * dt)
  }
  
  v
}

median_based_sd <- function(v) {
  v <- v[is.finite(v)]
  sqrt(median((v - median(v))^2))
}

analysis_check <- analysis_check %>%
arrange(trial, time_rel) %>%
group_by(trial) %>%
  mutate(
vx = gaze_velocity(gaze_x,sr),
vy = gaze_velocity(gaze_y,sr)
)%>%
ungroup()

gaze_interp <- function(v,x,samplig_rate){
  dt <- 1/samplig_rate
  v2 <- v
  v2 <- v2[-c(1, 2, length(v2), (length(v2)-1))]
  x0 <- x[2]
  v3 = 0
  result <- c()

  for (i in 1:length(v2)) {
    v3<- v3+v2[i]
    result <- c(result,v3*dt+x0)
  }
  x_rec <- c(x[1],x[2],result,x[(length(x)-1)],x[length(x)])
  x_rec
}


detect_candidates <- function(x, y, sampling_rate = 250, lambda = 5, min_samples = 3) {

  vx <- gaze_velocity(x,sr)
  vy <- gaze_velocity(y,sr)
  
  sigma_x <- median_based_sd(vx)
  sigma_y <- median_based_sd(vy)
  
  eta_x <- lambda * sigma_x
  eta_y <- lambda * sigma_y
  
  cand <- (vx / eta_x)^2 + (vy / eta_y)^2 > 1
  cand[is.na(cand)] <- FALSE
  
  starts <- which(diff(c(FALSE, cand)) == 1)
  ends   <- which(diff(c(cand, FALSE)) == -1)
  
  result <- lapply(seq_along(starts), function(k) {
    s <- starts[k]
    e <- ends[k]
    
    if ((e - s + 1) < min_samples) return(NULL)
    
    amp <- sqrt((x[e] - x[s])^2 + (y[e] - y[s])^2)
    pv  <- max(sqrt(vx[s:e]^2 + vy[s:e]^2), na.rm = TRUE)
    
    data.frame(
      onset = s,
      offset = e,
      duration_samples = e - s + 1,
      duration_ms = (e - s + 1) * (1000 / sampling_rate),
      amplitude = amp,
      peak_vel = pv
    )
  })
  
  events <- do.call(rbind, Filter(Negate(is.null), result))
  
}
events <- analysis_check %>%
  arrange(trial, time_rel) %>%
  group_by(trial) %>%
  group_modify(~ {
    ev <- detect_candidates(
      x = .x$gaze_x,
      y = .x$gaze_y,
      sampling_rate = sr
    )
    
    if (is.null(ev)) return(tibble())
    
    ev
      
  }) %>%
  ungroup()


second_dat <- analysis_check %>%
  group_by(trial) %>%
  mutate(
    pupil_for_interp_second = pupil_interp2
  ) %>%
  group_modify(~ {
    s <- .x
    ev <- events %>% filter(trial == .y$trial)
    
    if (nrow(ev) > 0) {
      for (i in 1:nrow(ev)) {
        start <- max(1, ev$onset[i])
        end <- min(nrow(s), ev$offset[i])
        s$pupil_for_interp_second[start:end] <- NA
      }
    }
    
    s
  }) %>%
  ungroup()

final_dat <- second_dat %>%
  arrange(trial, time_rel) %>%
  group_by(trial) %>%
  mutate(
    Pupil_final = na.approx(
       pupil_for_interp_second,
      x = time_rel,
      na.rm = FALSE
    )
  ) %>%
  ungroup()



library(ggplot2)
library(tidyverse)
library(patchwork)

dat<-second_dat

gorgx<-ggplot(data=dat)+geom_line( aes(x=time_rel, y=xv), color='red')
gorgx<-gorgx+ylim(-(blink_raw$display_coords[3]/2), (blink_raw$display_coords[3]/2))
gorgy<-ggplot(data=dat)+geom_line( aes(x=time_rel, y=yv), color='blue')
gorgy<-gorgy+ylim(-(blink_raw$display_coords[4]/2), (blink_raw$display_coords[4]/2))
gppl<-ggplot(data=dat)+geom_line(aes(x=time_rel, y=pupil_interp2), color='green')

gorgx / gorgy / gppl

tmpwid<-100
tblink<-1
prd<-(blink_raw$blinks$sttime_rel[tblink]/(1000/blink_raw$headers$rec_sample_rate)-tmpwid):(blink_raw$blinks$entime_rel[tblink]/(1000/blink_raw$headers$rec_sample_rate)+tmpwid)
plot(prd,dat$pupil_second[prd])

onsets<-events$onset*(1000/sr)
offsets<-events$offset*(1000/sr)
outdir <- "plot_dat"
for (i in 1:nrow(events)){
  nms<-i
  prd1<-(events$onset[nms]-tmpwid):(events$offset[nms]+tmpwid)
  prd2<-(events$onset[nms]):(events$offset[nms])
  xx<-c(events$onset[nms]-tmpwid, events$offset[nms]+tmpwid)
  yy<-c(min(dat$pupil_second[prd1],na.rm = T), max(dat$pupil_second[prd1],na.rm=T))
  png(
    filename = file.path(outdir, paste0("microsaccade_", i, ".png")),
    width = 800,
    height = 600
  )
  plot(prd1, dat$pupil_second[prd1], xlim=xx, ylim=yy)
  par (new =T)
  plot(prd2, dat$pupil_second[prd2], col='red', xlim=xx, ylim=yy)
  res[nms,]
  dev.off()
}

