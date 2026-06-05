Sys.setenv(EDFAPI = "C:/Program Files (x86)/SR Research/EyeLink/EDF_Access_API")
library(eyelinkReader)
# フォルダへのパス
DIR <- r"(C:\Users\shuma\OneDrive - 学校法人　金沢工業大学\ドキュメント\pupil)"

# ファイルパス（作業ディレクトリにある場合）
edf_file <- read_edf("mwtest.edf")

edf <- read_edf(
  "mwtest.edf",
  import_samples = TRUE
)

event <- edf$events


edSys.which("make")

eyelinkReader::compiled_library_status()

