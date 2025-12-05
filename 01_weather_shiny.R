# FAST Weather Preprocessing for Onion Climate Model ----
#
# This script:
#   - Loads required packages (installing them if missing).
#   - Reads compressed weather data (historical + SSP scenarios).
#   - Assumes derived variables (e.g. Ts_5cm_smooth, PAR, GDD_daily)
#     have already been computed and stored in the input file.
#   - Splits the data into a nested list:
#       result[[scenario]][[id_season]] -> data.table of daily weather

# Load required packages ----

load_if_needed <- function(pkgs) {
  for (pkg in pkgs) { 
    
    # Install package if it is not available
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
    }
    
    # Attach package, suppressing startup messages
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE)
    )
  }
}

load_if_needed(c(
  "zoo",             # time series tools (often used with weather data)
  "RcppRoll",        # fast rolling operations (alternate to data.table tools)
  "decisionSupport", # Monte Carlo / decision tools (used elsewhere in project)
  "compiler",        # byte-code compilation for speed (used in other scripts)
  "data.table"       # fast data manipulation and rolling functions
))


# Process weather data function ----

process_weather_data <- function(
    file_path = "weather_koeln-bonn_final_compressed.rds",
    scenarios = c("historical", "ssp126", "ssp245", "ssp370", "ssp585")
) {
  # 1) Read and convert to data.table ----------------------------------------
  weather_combined <- readRDS(file_path)
  data.table::setDT(weather_combined)  # in-place conversion
  
  # Sort once by season and day-of-year (yday) for correct temporal order
  data.table::setkey(weather_combined, id_season, yday)
  
  # NOTE:
  # 2) No longer creating derived variables here.
  #    We assume columns such as Ts_5cm_smooth, PAR, GDD_daily, ET0_mm, etc.
  #    are already present in 'weather_combined' from upstream processing.
  
  # 3) Extract scenario name from id_season -----------------------------------
  
  weather_combined[
    ,
    scenario := data.table::tstrsplit(id_season, "--", fixed = TRUE, keep = 1L)
  ]
  
  # 4) Split into a nested list: scenario -> id_season ------------------------
  
  # First split by scenario:
  #   result[[scenario]] -> all seasons for that scenario
  weather_precomputed <- split(
    weather_combined,
    by      = "scenario",
    keep.by = FALSE
  )
  
  # Then, within each scenario, split by id_season:
  #   result[[scenario]][[id_season]] -> one season's daily weather
  weather_precomputed <- lapply(weather_precomputed, function(dt) {
    split(dt, by = "id_season", keep.by = FALSE)
  })
  
  return(weather_precomputed)
}


# Run preprocessing ----

# This call reads the default weather file and prepares `weather_precomputed`.
weather_precomputed <- process_weather_data(
  file_path = "weather_koeln-bonn_final_compressed.rds"
)
