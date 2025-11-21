################################################################################
# FAST Weather Preprocessing 
################################################################################

################################################################################


load_if_needed <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}

load_if_needed(c(
  "zoo",
  "RcppRoll",
  "decisionSupport",
  "compiler",
  "data.table"
))

######################################################################

process_weather_data <- function(
    file_path = "weather_koeln-bonn_final_compressed.rds",
    scenarios = c("historical", "ssp126", "ssp245", "ssp370", "ssp585"),
    base_temp = 1) {
  
  
  weather_combined <- readRDS(file_path)
  setDT(weather_combined)
  
  # Sort once by id_season, date (assumes day ordering exists)
  setkey(weather_combined, id_season, yday)
  
  
  # Compute rolling means by reference using fast grouping
  weather_combined[, Ts_5cm_smooth := frollmean(Ts_5cm, n = 7, align = "right", fill = NA_real_), by = id_season]
  weather_combined[, PAR := Ra * 0.45]
  weather_combined[, GDD_daily := pmax(0, Tavg - base_temp)]
  

  # Pre-extract scenario name from id_season only once
  weather_combined[, scenario := tstrsplit(id_season, "--", fixed = TRUE, keep = 1L)]
  
  # Split in one vectorized call; no regex inside lapply
  weather_precomputed <- split(weather_combined, by = "scenario", keep.by = FALSE)
  
  # Inside each scenario, further split by id_season but keep as shallow copies
  weather_precomputed <- lapply(weather_precomputed, function(dt) {
    split(dt, by = "id_season", keep.by = FALSE)
  })a
  
  
  return(weather_precomputed)
}


  weather_precomputed <- process_weather_data(
    file_path = "weather_koeln-bonn_final_compressed.rds",
    base_temp = 1
  )
## Please explain better what the script actually does 
# Give titles, structure it also give comments what every section does 