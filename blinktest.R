Sys.setenv(EDFAPI = "C:/Program Files (x86)/SR Research/EyeLink/EDF_Access_API")
library(eyelinkReader)

blink_raw <- read_edf(
  "blinktest.edf",
  import_samples = TRUE
)

library(dplyr)
library(zoo)

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

screen_center_x <- 1920 / 2  
screen_center_y <- 1080 / 2  
analysis_check$x <- analysis_check$x - screen_center_x
analysis_check$y <- analysis_check$y - screen_center_y

# ---- monitor settings ----
diag_inch <- 23
aspect_w <- 16
aspect_h <- 9
view_dist_cm <- 50

res_x <- 1920
res_y <- 1080

samples <- samples %>%
  bind_cols(deg_pos)

samples_interp <- samples %>%
  arrange(trial, time_rel) %>%
  group_by(trial) %>%
  mutate(
    is_zero_pupil = !is.na(pupil) & pupil == 0,
    
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
  ungroup()

analysis_data <- samples_interp_final %>%
  filter(!is.na(pupil_interp))


x <- diff(analysis_data$pupil_interp)
pupil_velocity_average = mean(x)

SD = sd(x,na.rm = TRUE)
σ = 1.5

low_threshold = pupil_velocity_average - SD*σ
high_threshold = pupil_velocity_average + SD*σ


raw_artifacts <- which(x<low_threshold | x>high_threshold)
a <- diff(raw_artifacts)
artifacts_point <- sort(unique(c(raw_artifacts,raw_artifacts+1)))

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


analysis_check <- analysis_data %>%
  mutate(
    artifact_velocity = row_number() %in% artifacts
  )

blink_times <- analysis_check %>%
  filter(blink_window) %>%
  pull(time_rel)

analysis_check <- analysis_check %>%
  mutate(
    within_50ms_from_blink = sapply(
      time_rel,
      function(t) any(abs(t - blink_times) <= 50)
    ),
    
    blink_related_artifact = blink_window |
      (artifact_velocity & within_50ms_from_blink),
    
    pupil_restored = if_else(
      blink_related_artifact,
      pupil,
      pupil_interp
    )
  )

gaze_velocity <- function(x, sampling_rate) {
  x <- analysis_check$x_deg
  dt <- 1 / 250
  n <- length(x)
  v <- rep(NA_real_, n)
  
  if (n < 5) return(v)
  
  for (i in 3:(n - 2)) {
    v[i] <- (x[i + 2] + x[i + 1] - x[i - 1] - x[i - 2]) / (6 * dt)
  }
  
  v
}

interped_gaze <- function(v) {
  v <- v[is.finite(v)]
  sqrt(median(v^2) - median(v)^2)
}

detect_candidates <- function(x, y, sampling_rate = 250, lambda = 5) {
  vx <- gaze_velocity(x, sampling_rate)
  vy <- gaze_velocity(y, sampling_rate)
  
  sigma_x <- interped_gaze(vx)
  sigma_y <- interped_gaze(vy)
  
  eta_x <- lambda * sigma_x
  eta_y <- lambda * sigma_y
  
  cand <- (vx / eta_x)^2 + (vy / eta_y)^2 > 1
  cand[is.na(cand)] <- FALSE
  
  saccade_like_raw <- (vx / eta_x)^2 + (vy / eta_y)^2 > 1
  
  analysis_check$saccade_like_raw <- c(FALSE, saccade_like_raw)
  list(
    candidate = cand,
    vx = vx,
    vy = vy,
    eta_x = eta_x,
    eta_y = eta_y
  )
}




analysis_check <- analysis_check %>%
  arrange(trial, time_rel) %>%
  group_by(trial) %>%
  mutate(
    saccade_group = cumsum(
      saccade_like_raw != lag(saccade_like_raw, default = first(saccade_like_raw))
    )
  ) %>%
  group_by(trial, saccade_group) %>%
  mutate(
    saccade_like = saccade_like_raw & n() >= 2
  ) %>%
  ungroup()

imrl<- analysis_check[analysis_check$artifact_velocity==TRUE,]

blink_temp<- analysis_check[analysis_check$blink_window==TRUE,]
saccade_temp<- blink_temp[blink_temp$saccade_like==TRUE,]
saccade_temp2<-analysis_check[analysis_check$saccade_like==TRUE,]
plot(analysis_check$x,analysis_check$y)


plot(analysis_check$x_deg, analysis_check$y_deg, col = ifelse(!is.na(analysis_check$blink_window) & analysis_check$blink_window, ifelse(analysis_check$saccade_like,"green","red"), ifelse(!is.na(analysis_check$saccade_like)& analysis_check$saccade_like ,"blue","black")))


t <- seq(0, 2 * pi, length = 200)
ellipse_x <- threshold_x * cos(t)
ellipse_y <- threshold_y * sin(t)
lines(ellipse_x, ellipse_y, col = "blue", lwd = 2, lty = 2)

