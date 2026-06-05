# Microsaccade detection based on Engbert & Mergenthaler (2006)
# Paper: "Microsaccades are triggered by low retinal image slip"
# PNAS 103(18): 7192-7197

# ---- パラメータ ----
SR      <- 500   # サンプリングレート (Hz)
DT      <- 1/SR  # 時間ステップ (s) = 0.002s
LAMBDA  <- 5     # 閾値の倍率（論文推奨値）
MIN_DUR <- 3     # 最小持続サンプル数（3サンプル = 6ms）

# ---- 式[1]: 位置列 → 速度列 ----
calc_velocity <- function(x, dt) {
  n <- length(x)
  v <- rep(NA, n)
  # 両端2サンプルはNAのまま
  for (i in 3:(n - 2)) {
    v[i] <- (x[i+2] + x[i+1] - x[i-1] - x[i-2]) / (6 * dt)
  }
  return(v)
}

# ---- 中央値ベースの標準偏差 (σ) ----
# σ² = median((v - median(v))²)
calc_sigma <- function(v) {
  v <- v[!is.na(v)]
  sqrt(median((v - median(v))^2))
}

# ---- マイクロサッカード検出 ----
detect_microsaccades <- function(x, y, sr = SR, lambda = LAMBDA, min_dur = MIN_DUR) {
  dt <- 1 / sr

  # 速度計算
  vx <- calc_velocity(x, dt)
  vy <- calc_velocity(y, dt)

  # 各軸の閾値 η = λ·σ
  eta_x <- lambda * calc_sigma(vx)
  eta_y <- lambda * calc_sigma(vy)

  # 楕円閾値: (vx/ηx)² + (vy/ηy)² > 1
  above <- (vx / eta_x)^2 + (vy / eta_y)^2 > 1
  above[is.na(above)] <- FALSE

  # 閾値超過区間を抽出
  starts <- which(diff(c(FALSE, above)) ==  1)
  ends   <- which(diff(c(above, FALSE)) == -1)

  result <- lapply(seq_along(starts), function(k) {
    s <- starts[k]
    e <- ends[k]
    if ((e - s + 1) < min_dur) return(NULL)  # 最小持続時間フィルタ

    amp <- sqrt((x[e] - x[s])^2 + (y[e] - y[s])^2)
    pv  <- sqrt(max(vx[s:e]^2 + vy[s:e]^2, na.rm = TRUE))

    data.frame(onset = s, offset = e,
               duration_ms = (e - s + 1) * (1000/sr),
               amplitude_deg = amp,
               peak_vel = pv)
  })
# 
  do.call(rbind, Filter(Negate(is.null), result))
}

# ---- 使用例（ダミーデータ） ----
set.seed(42)
n   <- 1500  # 3秒分 @ 500Hz
t   <- seq(0, by = DT, length.out = n)
# ランダムウォーク（ドリフト）を模擬
x   <- cumsum(rnorm(n, 0, 0.01))
y   <- cumsum(rnorm(n, 0, 0.01))
# サッカードを1個埋め込む（500サンプル目付近）
x[500:510] <- x[500:510] + seq(0, 0.5, length.out = 11)
y[500:510] <- y[500:510] + seq(0, 0.3, length.out = 11)

ms <- detect_microsaccades(x, y)
print(ms)
