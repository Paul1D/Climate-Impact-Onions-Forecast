# Helper Functions for Onion Climate Model ----
#
# This script defines:
#   - A package loader helper: load_if_needed()
#   - A factory function: helper_function()
#     which compiles and returns (or attaches) a set of helper functions:
#
#   Generic utilities:
#     - sat()              : smooth saturating transform
#     - consec_counter()   : count consecutive TRUE values (e.g. dry days)
#     - fast_chance_event(): simple stochastic gate (legacy)
#     - expected_loss()    : risk × loss_fraction, in [0,1]
#
#   Sowing rule:
#     - compute_onion_sowing_yday_fast(): discover sowing DOY from weather
#
#   Stress modules:
#     - get_drought_stress()
#     - get_extreme_rain_stress()
#     - get_hail_stress()
#     - get_botrytis_stress() 
#     - get_downy_mildew_stress()
#     - get_fusarium_stress()
#     - get_thrips_stress()
#
#   Biomass:
#     - calc_bio_vectorized(): daily biomass increment from PAR, LAI, T, water
#
# At the end, helper_function(attach_to_global = TRUE) is called so that all
# helpers are available by name in the global environment.
# All functions are based on inputs from the weather data file. Threshholds, impact factors and optimums are defined in the input table.



# Load required packages ----

load_if_needed <- function(pkgs) {
  for (pkg in pkgs) {
    
    # Install package if it is not installed yet
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
    }
    
    # Attach package quietly
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE)
    )
  }
}

load_if_needed(c(
  "zoo",
  "RcppRoll",
  "decisionSupport",
  "compiler",
  "data.table"
))


# Helper Function: Stress, ----

helper_function <- function(attach_to_global = FALSE) {
  
  # Generic utility functions ----
  
  # Saturating transform: x >= 0 -> [0,1) with slope controlled by `impact`
  sat <- compiler::cmpfun(function(x, impact) {
    pmin(1 - exp(-impact * pmax(x, 0)), 1)
  })
  
  # Count consecutive TRUEs in a logical vector (e.g. dry days in a row)
  consec_counter <- compiler::cmpfun(function(x_is_true) {
    x <- as.logical(x_is_true)
    n <- length(x)
    out <- integer(n)
    if (n == 0L) return(out)
    k <- 0L
    for (i in seq_len(n)) {
      if (isTRUE(x[i])) {
        k <- k + 1L
      } else {
        k <- 0L
      }
      out[i] <- k
    }
    out
  })
  
  # Stochastic gate: returns value_if with probability p, otherwise value_if_not
  fast_chance_event <- compiler::cmpfun(function(p, value_if, value_if_not = 1) {
    if (is.na(p)) return(value_if_not)
    if (runif(1) < p) value_if else value_if_not
  })
  
  # Expected loss: risk × loss_fraction, NA-safe, bounded to [0,1]
  expected_loss <- compiler::cmpfun(function(risk, loss_fraction) {
    if (is.na(risk) || is.na(loss_fraction)) return(0)
    pmin(pmax(risk, 0), 1) * pmin(pmax(loss_fraction, 0), 1)
  })
  
  
  # Onion sowing date (vectorized) ----
  
  compute_onion_sowing_yday_fast <- compiler::cmpfun(function(df,
                                                              T_soil_min_sowing_p,
                                                              warm_days_needed_sowing_p,
                                                              frost_buffer_days_sowing_p,
                                                              min_yday_sowing_c,
                                                              max_yday_sowing_c,
                                                              rain_window_days_sowing_p,
                                                              rain_sum_max_mm_sowing_p,
                                                              min_dry_days_sowing_c) {
    n <- nrow(df)
    if (n == 0L) return(max_yday_sowing_c)
    
    # Window lengths (converted to integers)
    n_warm  <- max(1L, as.integer(round(warm_days_needed_sowing_p)))
    n_frost <- max(1L, as.integer(round(frost_buffer_days_sowing_p)))
    n_rain  <- max(1L, as.integer(round(rain_window_days_sowing_p)))
    n_dry   <- max(1L, as.integer(round(min_dry_days_sowing_c)))
    
    # Basic flags
    frost_flag <- df$Tmin < -2                     # frost occurrence
    is_dry     <- is.na(df$Prec) | df$Prec == 0    # no rain recorded
    
    # 1) Warm spell: soil temperature at 5 cm above threshold
    warm_ok <- data.table::frollsum(
      df$Ts_5cm_smooth >= T_soil_min_sowing_p,
      n      = n_warm,
      align  = "left",
      fill   = 0L
    ) == n_warm
    
    # 2) Frost buffer after warm spell
    frost_sum <- data.table::frollsum(
      frost_flag,
      n      = n_frost,
      align  = "left",
      fill   = 0L
    )
    frost_ok  <- data.table::shift(frost_sum, n = n_warm, fill = 0L) == 0L
    
    # 3) Preceding rain (before current day)
    Prec0 <- ifelse(is.na(df$Prec), 0, df$Prec)
    rain_prev <- data.table::frollsum(
      Prec0,
      n      = n_rain,
      align  = "right",
      fill   = 0
    )
    rain_prev_excl_today <- data.table::shift(rain_prev, n = 1L, fill = 0)
    
    # 4) Dry streak (consecutive dry days up to current day)
    dry_streak <- consec_counter(is_dry)
    
    # 5) Calendar window
    in_window <- df$yday >= min_yday_sowing_c & df$yday <= max_yday_sowing_c
    
    # Candidate days: satisfy all conditions
    cand <- warm_ok & frost_ok &
      (rain_prev_excl_today <= rain_sum_max_mm_sowing_p) &
      (dry_streak >= n_dry) & in_window
    
    idx <- which(cand)
    if (length(idx) == 0L) return(max_yday_sowing_c)
    
    # Place sowing at the end of the warm spell
    end_idx <- idx[1] + n_warm - 1L
    end_idx <- min(end_idx, n)
    df$yday[end_idx]
  })
  
  
  # Stress functions (abiotic & biotic) ----
  
  # Drought stress: combines RH, high temperature and dry spell length
  get_drought_stress <- compiler::cmpfun(function(Tavg,
                                                  RH_mean,
                                                  days_consec_dry,
                                                  rh_drought_threshold_p,
                                                  Tavg_drought_threshold_p,
                                                  impact_rh_drought_t,
                                                  impact_days_dry_drought_t,
                                                  impact_temp_drought_t) {
    # RH component
    R_RH <- ifelse(
      RH_mean < rh_drought_threshold_p,
      1 - exp(- (impact_rh_drought_t / 2) * (rh_drought_threshold_p - RH_mean)),
      0
    )
    R_RH <- pmin(R_RH, 1)
    
    # Temperature component
    R_T  <- ifelse(
      Tavg > Tavg_drought_threshold_p,
      1 - exp(- (impact_temp_drought_t / 2) * (Tavg - Tavg_drought_threshold_p)),
      0
    )
    R_T <- pmin(R_T, 1)
    
    # Dry streak component
    R_D <- 1 - exp(- (impact_days_dry_drought_t / 2) * days_consec_dry)
    R_D <- pmin(R_D, 1)
    
    pmin(R_RH * R_T * R_D, 1)
  })
  
  # Extreme rain stress: weighted count of extreme rainfall days
  get_extreme_rain_stress <- compiler::cmpfun(function(Prec,
                                                       prec_extreme_rain_medium_p,
                                                       prec_extreme_rain_high_p,
                                                       impact_days_extreme_rain_t) {
    extreme_days_medium <- sum(
      Prec >= prec_extreme_rain_medium_p & Prec < prec_extreme_rain_high_p,
      na.rm = TRUE
    )
    extreme_days_high   <- sum(
      Prec >= prec_extreme_rain_high_p,
      na.rm = TRUE
    )
    weighted_extreme_days <- extreme_days_medium + 2 * extreme_days_high
    risk_extreme_rain <- 1 - exp(- (impact_days_extreme_rain_t / 2) * weighted_extreme_days)
    pmin(risk_extreme_rain, 1)
  })
  
  # Hail stress: proxy using warm + rainy days
  get_hail_stress <- compiler::cmpfun(function(Tavg,
                                               Prec,
                                               prec_hail_threshold_p,
                                               Tavg_hail_threshold_p,
                                               impact_days_hail_t) {
    hail_days <- sum(
      Prec >= prec_hail_threshold_p & Tavg >= Tavg_hail_threshold_p,
      na.rm = TRUE
    )
    risk_hail <- 1 - exp(- (impact_days_hail_t) * hail_days)
    pmin(risk_hail, 1)
  })
  
  # Botrytis stress: high RH, many wet days, suitable T range
  
  get_botrytis_stress <- compiler::cmpfun(function(Tmin,      
                                                   RH_max,   
                                                   days_consec_wet,
                                                   rh_botrytis_threshold_p,
                                                   impact_rh_botrytis_t,
                                                   impact_days_wet_botrytis_t,
                                                   Topt_botrytis_p,
                                                   Twidth_botrytis_p) {

    
    # 1) Luftfeuchte-Gate: ohne sehr hohe RH kein Risiko
    R_H <- ifelse(
      RH_max >= rh_botrytis_threshold_p,  # z.B. 90%
      1 - exp(- (impact_rh_botrytis_t / 2) * (RH_max - rh_botrytis_threshold_p)),
      0
    )
    R_H <- pmin(R_H, 1)
    
    # 2) length of wet period 
    R_D <- 1 - exp(- (impact_days_wet_botrytis_t / 2) * days_consec_wet)
    R_D <- pmin(R_D, 1)
    
    # 3) Temperature during the night
    R_T <- exp(- (Tmin - Topt_botrytis_p)^2 / (2 * Twidth_botrytis_p^2))
    R_T <- pmin(R_T, 1)
    
    # 4) Gesamt-Risiko
    pmin(R_H * R_D * R_T, 1)
  })
  
  
  # Downy mildew stress: similar structure, different parameters
  get_downy_mildew_stress <- compiler::cmpfun(function(Tmin,
                                                       RH_max,   
                                                       days_consec_wet,
                                                       rh_mildew_threshold_p,
                                                       impact_rh_mildew_t,
                                                       impact_days_wet_mildew_t,
                                                       Topt_mildew_p,
                                                       Twidth_mildew_p) {
  
    
    # 1) High RH required 
    R_H <- ifelse(
      RH_max >= rh_mildew_threshold_p,
      1 - exp(- (impact_rh_mildew_t / 2) * (RH_max - rh_mildew_threshold_p)),
      0
    )
    R_H <- pmin(R_H, 1)
    
    # 2) amount of wet days--> therefore wet nights
    R_D <- 1 - exp(- (impact_days_wet_mildew_t / 2) * days_consec_wet)
    R_D <- pmin(R_D, 1)
    
    # 3) Temperature optimum during the night 
    R_T <- exp(- (Tmin - Topt_mildew_p)^2 / (2 * Twidth_mildew_p^2))
    R_T <- pmin(R_T, 1)
    
    pmin(R_H * R_D * R_T, 1)
  })
  
  
  # Fusarium stress: soil T + RH + wet days
  get_fusarium_stress <- compiler::cmpfun(function(Ts_5cm,
                                                   RH_mean,
                                                   days_consec_wet,
                                                   rh_fusarium_threshold_p,
                                                   impact_rh_fusarium_t,
                                                   impact_days_wet_fusarium_t,
                                                   Tsopt_fusarium_p,
                                                   Tswidth_fusarium_p) {
    R_H <- ifelse(
      RH_mean >= rh_fusarium_threshold_p,
      1 - exp(- (impact_rh_fusarium_t / 2) * (RH_mean - rh_fusarium_threshold_p)),
      0
    )
    R_H <- pmin(R_H, 1)
    
    R_D <- 1 - exp(- (impact_days_wet_fusarium_t / 2) * days_consec_wet)
    R_D <- pmin(R_D, 1)
    
    R_T <- exp(- (Ts_5cm - Tsopt_fusarium_p)^2 / (2 * Tswidth_fusarium_p^2))
    R_T <- pmin(R_T, 1)
    
    pmin(R_H * R_D * R_T, 1)
  })
  
  # Thrips stress: dry air, suitable T, dry spells
  get_thrips_stress <- compiler::cmpfun(function(Tavg,
                                                 RH_mean,
                                                 days_consec_dry,
                                                 rh_thrips_threshold_p,
                                                 Topt_thrips_p,
                                                 Twidth_thrips_p,
                                                 impact_rh_thrips_t,
                                                 impact_days_dry_thrips_t) {
    R_H <- ifelse(
      RH_mean < rh_thrips_threshold_p,
      1 - exp(- (impact_rh_thrips_t / 2) * (rh_thrips_threshold_p - RH_mean)),
      0
    )
    R_H <- pmin(R_H, 1)
    
    R_T <- exp(- (Tavg - Topt_thrips_p)^2 / (2 * Twidth_thrips_p^2))
    R_T <- pmin(R_T, 1)
    
    R_D <- 1 - exp(- (impact_days_dry_thrips_t / 2) * days_consec_dry)
    R_D <- pmin(R_D, 1)
    
    pmin(R_H * R_T * R_D, 1)
  })
  
  
  # Biomass calculation (vectorized) ----
  
  calc_bio_vectorized <- compiler::cmpfun(function(PAR, LAI, Tavg, Prec,
                                                   f_T_1_lower, f_T_1_upper,
                                                   f_T_0_lower, f_T_0_upper,
                                                   f_W_1_lower, f_W_1_upper,
                                                   f_W_0_5, LUE_onion, lec_k) {
    # Temperature response
    f_T <- ifelse(Tavg >= f_T_1_lower & Tavg <= f_T_1_upper, 1,
                  ifelse(Tavg <  f_T_0_lower | Tavg >  f_T_0_upper, 0, 0.5))
    
    # Water response
    f_W <- ifelse(Prec >= f_W_1_lower & Prec <= f_W_1_upper, 1,
                  ifelse(Prec < f_W_0_5, 0.5, 0.7))
    
    # Daily biomass increment via Beer-Lambert + LUE
    LUE_onion * PAR * (1 - exp(-lec_k * LAI)) * f_T * f_W
  })
  

# total_biomass       -> g DM / m² (Summe)
# PAR: MJ m^-2 d^-1
# LUE_onion_p: g DM MJ^-1
# calc_bio_vectorized -> g DM m^-2 d^-1
  
  
  
  # Collect helpers and optionally attach to global environment ----
  
  helpers <- list(
    sat                         = sat,
    expected_loss               = expected_loss,
    fast_chance_event           = fast_chance_event,
    consec_counter              = consec_counter,
    compute_onion_sowing_yday_fast = compute_onion_sowing_yday_fast,
    get_drought_stress          = get_drought_stress,
    get_extreme_rain_stress     = get_extreme_rain_stress,
    get_hail_stress             = get_hail_stress,
    get_botrytis_stress         = get_botrytis_stress,
    get_downy_mildew_stress     = get_downy_mildew_stress,
    get_thrips_stress           = get_thrips_stress,
    get_fusarium_stress         = get_fusarium_stress,
    calc_bio_vectorized         = calc_bio_vectorized
  )
  
  if (attach_to_global) {
    list2env(helpers, envir = .GlobalEnv)
    message("✅ Helper functions attached to global environment.")
  }
  
  helpers
}


# Attach helpers for interactive use ----

helper_function(attach_to_global = TRUE)
