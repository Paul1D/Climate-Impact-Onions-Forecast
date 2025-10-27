################################################################################
# Helper Function load 
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

################################################################################
# Optimized Helper Functions 
################################################################################

helper_function <- function(attach_to_global = FALSE) {
  
  # Fast monolithic counter (faster than rle for long vectors)
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
  
  # Simple stochastic gate
  fast_chance_event <- compiler::cmpfun(function(p, value_if, value_if_not = 1) {
    if (is.na(p)) return(value_if_not)
    if (runif(1) < p) value_if else value_if_not
  })
  
  # Saturation 
  sat <- compiler::cmpfun(function(x, impact) pmin(1 - exp(-impact * pmax(x, 0)), 1))
  
  ################################################################################
  # Onion sowing — fully vectorized 
  ################################################################################
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
    
    # --- Sanitize rolling window lengths
    n_warm  <- max(1L, as.integer(round(warm_days_needed_sowing_p)))
    n_frost <- max(1L, as.integer(round(frost_buffer_days_sowing_p)))
    n_rain  <- max(1L, as.integer(round(rain_window_days_sowing_p)))
    n_dry   <- max(1L, as.integer(round(min_dry_days_sowing_c)))
    
    # Flags
    frost_flag <- df$Tmin < -2
    is_dry     <- is.na(df$Prec) | df$Prec == 0
    
    # Rolling warm spell: consecutive warm days with Ts >= threshold
    warm_ok <- data.table::frollsum(df$Ts_5cm_smooth >= T_soil_min_sowing_p,
                                    n = n_warm,
                                    align = "left", fill = 0L) == n_warm
    
    # No frost in the buffer after warm spell (shift by warm window)
    frost_sum <- data.table::frollsum(frost_flag,
                                      n = n_frost,
                                      align = "left", fill = 0L)
    frost_ok <- data.table::shift(frost_sum, n = n_warm, fill = 0L) == 0L
    
    # Rain sum in the *previous* n_rain days before current day (exclude today)
    Prec0 <- ifelse(is.na(df$Prec), 0, df$Prec)
    rain_prev <- data.table::frollsum(Prec0, n = n_rain,
                                      align = "right", fill = 0)
    rain_prev_excl_today <- data.table::shift(rain_prev, n = 1L, fill = 0)
    
    # Dry streak length up to i
    dry_streak <- consec_counter(is_dry)
    
    # Candidate sowing window filter
    in_window <- df$yday >= min_yday_sowing_c & df$yday <= max_yday_sowing_c
    cand <- warm_ok & frost_ok &
      (rain_prev_excl_today <= rain_sum_max_mm_sowing_p) &
      (dry_streak >= n_dry) & in_window
    
    idx <- which(cand)
    if (length(idx) == 0L) return(max_yday_sowing_c)
    
    # Return yday at the end of the warm spell
    end_idx <- idx[1] + n_warm - 1L
    end_idx <- min(end_idx, n)
    df$yday[end_idx]
  })
  
  
  ################################################################################
  # Stress functions — lean branching, vector-friendly internals
  ################################################################################
  
  get_seedbed_stress <- compiler::cmpfun(function(Prec, Ts_5cm, RH_mean, day_consec_wet,
                                                  prec_seedbed_medium_p,
                                                  prec_seedbed_high_p,
                                                  Ts_min_seedbed_p,
                                                  RH_high_seedbed_p,
                                                  risk_seedbed_medium_t,
                                                  risk_seedbed_high_t,
                                                  impact_Ts_seedbed_t,
                                                  impact_RH_seedbed_t,
                                                  risk_additional_day = 0.05) {
    R_P <- if (Prec > prec_seedbed_high_p) risk_seedbed_high_t
    else if (Prec > prec_seedbed_medium_p) risk_seedbed_medium_t else 0
    R_Ts <- ifelse(Ts_5cm < Ts_min_seedbed_p,
                   sat(Ts_min_seedbed_p - Ts_5cm, impact_Ts_seedbed_t/2), 0)
    R_RH <- ifelse(RH_mean > RH_high_seedbed_p,
                   sat(RH_mean - RH_high_seedbed_p, impact_RH_seedbed_t/2), 0)
    pmin((R_P + R_Ts + R_RH) + (day_consec_wet - 1) * risk_additional_day, 1)
  })
  
  get_drought_stress <- compiler::cmpfun(function(Prec, Tavg, RH_mean, days_consec_dry,
                                                  prec_drought_threshold_p,
                                                  Tavg_drought_threshold_p,
                                                  RH_drought_threshold_p,
                                                  impact_prec_drought_t,
                                                  impact_days_dry_drought_t,
                                                  impact_temp_drought_t,
                                                  impact_RH_drought_t) {
    R_P  <- ifelse(Prec < prec_drought_threshold_p,
                   sat(prec_drought_threshold_p - Prec, impact_prec_drought_t/2), 0)
    R_T  <- ifelse(Tavg > Tavg_drought_threshold_p,
                   sat(Tavg - Tavg_drought_threshold_p, impact_temp_drought_t/2), 0)
    R_D  <- sat(days_consec_dry, impact_days_dry_drought_t/2)
    R_RH <- ifelse(RH_mean < RH_drought_threshold_p,
                   sat(RH_drought_threshold_p - RH_mean, impact_RH_drought_t/2), 0)
    pmin(R_P * R_T * R_D * R_RH, 1)
  })
  
  get_extreme_rain_stress <- compiler::cmpfun(function(Prec,
                                                       prec_extreme_rain_medium_p,
                                                       prec_extreme_rain_high_p,
                                                       impact_days_extreme_rain_t) {
    extreme_days_medium <- sum(Prec >= prec_extreme_rain_medium_p & 
                                 Prec <  prec_extreme_rain_high_p, na.rm = TRUE)
    extreme_days_high   <- sum(Prec >= prec_extreme_rain_high_p, na.rm = TRUE)
    sat(extreme_days_medium + 2 * extreme_days_high, impact_days_extreme_rain_t/2)
  })
  
  get_hail_stress <- compiler::cmpfun(function(Tavg, Prec,
                                               prec_hail_threshold_p,
                                               Tavg_hail_threshold_p,
                                               impact_days_hail_t) {
    hail_days <- sum(Prec >= prec_hail_threshold_p & 
                       Tavg >= Tavg_hail_threshold_p, na.rm = TRUE)
    sat(hail_days, impact_days_hail_t/2)
  })
  
  get_botrytis_stress <- compiler::cmpfun(function(Tavg, RH_mean, Prec, days_consec_wet,
                                                   Topt_botrytis_p, Twidth_botrytis_p,
                                                   RH_botrytis_threshold_p,
                                                   prec_botrytis_threshold_p,
                                                   wet_days_needed_p,
                                                   impact_RH_botrytis_t,
                                                   impact_prec_botrytis_t,
                                                   impact_days_wet_botrytis_t) {
    R_T <- exp(-((Tavg - Topt_botrytis_p)^2) / (2 * Twidth_botrytis_p^2))
    R_RH <- ifelse(RH_mean > RH_botrytis_threshold_p,
                   sat(RH_mean - RH_botrytis_threshold_p, impact_RH_botrytis_t/2), 0)
    R_P  <- ifelse(Prec >= prec_botrytis_threshold_p,
                   sat(Prec - prec_botrytis_threshold_p, impact_prec_botrytis_t/2), 0)
    R_D  <- sat(days_consec_wet - wet_days_needed_p, impact_days_wet_botrytis_t/2)
    pmin(R_T * R_RH * R_P * R_D, 1)
  })
  
  get_downy_mildew_stress <- compiler::cmpfun(function(Tavg, RH_mean, days_consec_wet,
                                                       Topt_mildew_p, Twidth_mildew_p,
                                                       RH_mildew_threshold_p,
                                                       wet_days_needed_p,
                                                       impact_RH_mildew_t,
                                                       impact_days_wet_mildew_t) {
    R_T  <- exp(-((Tavg - Topt_mildew_p)^2) / (2 * Twidth_mildew_p^2))
    R_RH <- ifelse(RH_mean >= RH_mildew_threshold_p,
                   sat(RH_mean - RH_mildew_threshold_p, impact_RH_mildew_t/2), 0)
    R_D  <- sat(days_consec_wet - wet_days_needed_p, impact_days_wet_mildew_t/2)
    pmin(R_T * R_RH * R_D, 1)
  })
  
  get_thrips_stress <- compiler::cmpfun(function(Tavg, Prec, days_consec_dry,
                                                 prec_thrips_threshold_p,
                                                 Topt_thrips_p, Twidth_thrips_p,
                                                 impact_prec_thrips_t,
                                                 impact_days_dry_thrips_t) {
    R_P <- ifelse(Prec < prec_thrips_threshold_p,
                  sat(prec_thrips_threshold_p - Prec, impact_prec_thrips_t/2), 0)
    R_T <- exp(-((Tavg - Topt_thrips_p)^2) / (2 * Twidth_thrips_p^2))
    R_D <- sat(days_consec_dry, impact_days_dry_thrips_t/2)
    pmin(R_P * R_T * R_D, 1)
  })
  
  get_onion_fly_stress <- compiler::cmpfun(function(Tavg, Prec, days_consec_wet,
                                                    prec_onion_fly_threshold_p,
                                                    Topt_onion_fly_p,
                                                    Twidth_onion_fly_p,
                                                    impact_prec_onion_fly_t,
                                                    impact_days_wet_onion_fly_t) {
    R_P <- ifelse(Prec >= prec_onion_fly_threshold_p,
                  sat(Prec - prec_onion_fly_threshold_p, impact_prec_onion_fly_t/2), 0)
    R_T <- exp(-((Tavg - Topt_onion_fly_p)^2) / (2 * Twidth_onion_fly_p^2))
    R_D <- sat(days_consec_wet, impact_days_wet_onion_fly_t/2)
    pmin(R_P * R_T * R_D, 1)
  })
  
  get_wireworm_stress <- compiler::cmpfun(function(Tavg, Ts_5cm, Prec, days_consec_wet,
                                                   prec_wireworm_threshold_p,
                                                   Topt_wireworm_p, Twidth_wireworm_p,
                                                   Ts_opt_wireworm_p,
                                                   Ts_width_wireworm_p,
                                                   impact_prec_wireworm_t,
                                                   impact_days_wet_wireworm_t) {
    R_P  <- ifelse(Prec >= prec_wireworm_threshold_p,
                   sat(Prec - prec_wireworm_threshold_p, impact_prec_wireworm_t/2), 0)
    R_Ta <- exp(-((Tavg - Topt_wireworm_p)^2) / (2 * Twidth_wireworm_p^2))
    R_Ts <- exp(-((Ts_5cm - Ts_opt_wireworm_p)^2) / (2 * Ts_width_wireworm_p^2))
    R_D  <- sat(days_consec_wet, impact_days_wet_wireworm_t/2)
    pmin(R_P * R_Ta * R_Ts * R_D, 1)
  })
  
  get_fusarium_stress <- compiler::cmpfun(function(Tavg, Prec, days_consec_wet,
                                                   prec_fusarium_threshold_p,
                                                   impact_prec_fusarium_t,
                                                   impact_days_wet_fusarium_t,
                                                   Topt_fusarium_p,
                                                   Twidth_fusarium_p) {
    R_P <- ifelse(Prec >= prec_fusarium_threshold_p,
                  sat(Prec - prec_fusarium_threshold_p, impact_prec_fusarium_t/2), 0)
    R_D <- sat(days_consec_wet, impact_days_wet_fusarium_t/2)
    R_T <- exp(-((Tavg - Topt_fusarium_p)^2) / (2 * Twidth_fusarium_p^2))
    pmin(R_P * R_D * R_T, 1)
  })
  
  ################################################################################
  # Biomass — vectorized
  ################################################################################
  calc_bio_vectorized <- compiler::cmpfun(function(PAR, LAI, Tavg, Prec,
                                                   f_T_1_lower, f_T_1_upper,
                                                   f_T_0_lower, f_T_0_upper,
                                                   f_W_1_lower, f_W_1_upper,
                                                   f_W_0_5, LUE_onion, lec_k) {
    f_T <- ifelse(Tavg >= f_T_1_lower & Tavg <= f_T_1_upper, 1,
                  ifelse(Tavg <  f_T_0_lower | Tavg >  f_T_0_upper, 0, 0.5))
    f_W <- ifelse(Prec >= f_W_1_lower & Prec <= f_W_1_upper, 1,
                  ifelse(Prec < f_W_0_5, 0.5, 0.7))
    LUE_onion * PAR * (1 - exp(-lec_k * LAI)) * f_T * f_W
  })
  
  helpers <- list(
    sat = sat,
    fast_chance_event = fast_chance_event,
    consec_counter = consec_counter,
    compute_onion_sowing_yday_fast = compute_onion_sowing_yday_fast,
    get_seedbed_stress = get_seedbed_stress,
    get_drought_stress = get_drought_stress,
    get_extreme_rain_stress = get_extreme_rain_stress,
    get_hail_stress = get_hail_stress,
    get_botrytis_stress = get_botrytis_stress,
    get_downy_mildew_stress = get_downy_mildew_stress,
    get_thrips_stress = get_thrips_stress,
    get_onion_fly_stress = get_onion_fly_stress,
    get_wireworm_stress = get_wireworm_stress,
    get_fusarium_stress = get_fusarium_stress,
    calc_bio_vectorized = calc_bio_vectorized
  )
  
  if (attach_to_global) {
    list2env(helpers, envir = .GlobalEnv)
    message("✅ Helper functions attached to global environment.")
  }
  helpers
}

# Attach helpers to global for plainNames model
helper_function(attach_to_global = TRUE)
  