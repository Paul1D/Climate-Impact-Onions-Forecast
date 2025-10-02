 library(tidyverse)
 library(decisionSupport)
 library(dplyr)



##### Input Table Onion model


input_variables <- read.csv("input_table_onion.csv", header = TRUE, sep = ";")


onion_climate_impact <- function(){ # Start of onion_climate_impact function
  
  # 1st step randomly selecet a season
  # Loop through each scenario and assign to List
  
  weather_scenario_list <- list()
  
  for (scenario in scenarios) {
    # Find all matching names
    matching_names <- names(weather_list)[grepl(paste0("^", scenario, "--"), names(weather_list))]
    
    # Randomly select one
    selected_name <- sample(matching_names, 1)
    
    # Store in the list with name like "weather_historical"
    weather_scenario_list[[paste0("weather_", scenario)]] <- weather_list[[selected_name]]
    
    # Optional: print the selected name
    cat("Selected for", scenario, ":", selected_name, "\n")
  }
  
  ## just for name storage
  
  weather_scenario_names <- weather_scenario_list
  
  
  ### Several steps are needed in order to calculate PAR, these are explained in the following
  
  # Calculate Extraterrestrial Radiation (Ra)
  #
  # This function calculates Ra (in MJ/m²/day) based on latitude and day of the year,
  # following FAO-56 methodology.
  # This is to calcualte extra terrestrial radiation, which is required in order to calculate PAR
  
  # Functions
  calc_Ra <- function(latitude_koeln_bonn_c, yday) {
    Gsc <- 0.0820
    phi <- latitude_koeln_bonn_c * pi / 180
    dr <- 1 + 0.033 * cos(2 * pi * yday / 365)
    delta <- 0.409 * sin((2 * pi * yday / 365) - 1.39)
    omega_s <- acos(-tan(phi) * tan(delta))
    Ra <- (24 * 60 / pi) * Gsc * dr *
      (omega_s * sin(phi) * sin(delta) + cos(phi) * cos(delta) * sin(omega_s))
    return(Ra)
  }
  
  calc_Rs <- function(Ra, deltaT) {
    return(Krs_c * Ra * deltaT)
  }
  
  
  # Process each scenario
  weather_scenario_list <- lapply(names(weather_scenario_list), function(scenario) {
    df <- weather_scenario_list[[scenario]]
    
    # --- Calculate PAR ---
    Ra <- calc_Ra(latitude_koeln_bonn_c, df$yday)
    deltaT <- df$Tmax - df$Tmin
    Rs <- calc_Rs(Ra, deltaT)
    df$PAR <- Rs * 0.45
    
    # --- Calculate daily GDD ---
    df$GDD_daily <- pmax(0, ((df$Tmax + df$Tmin) / 2) - base_temp_p)
    
    # Sort and filter days after planting
    df_sorted <- df[order(df$yday), ]
    df_after_planting <- df_sorted[df_sorted$yday >= planting_yday_p, ]
    
    
    # PHASE 1: Emergence
    
    gdd_emergence_cumsum <- cumsum(df_after_planting$GDD_daily)
    emergence_index <- which(gdd_emergence_cumsum >= GDD_field_emergence_required_p)[1]
    
    if (!is.na(emergence_index)) {
      emergence_yday <- df_after_planting$yday[emergence_index]
      df$in_emergence_phase <- df$yday %in% df_after_planting$yday[1:emergence_index]
    } else {
      emergence_yday <- NA
      df$in_emergence_phase <- FALSE
    }
    
    
    # PHASE 2: Vegetative
    
    if (!is.na(emergence_yday)) {
      df_after_emergence <- df[df$yday >= emergence_yday, ]
      gdd_veg_cumsum <- cumsum(df_after_emergence$GDD_daily)
      bulbing_index <- which(gdd_veg_cumsum >= GDD_vegetative_required_p)[1]
      
      if (!is.na(bulbing_index)) {
        bulbing_yday <- df_after_emergence$yday[bulbing_index]
        df$in_vegetative_phase <- df$yday >= emergence_yday & df$yday < bulbing_yday
      } else {
        bulbing_yday <- NA
        df$in_vegetative_phase <- FALSE
      }
    } else {
      bulbing_yday <- NA
      df$in_vegetative_phase <- FALSE
    }
    
    
    # PHASE 3: Bulbing
    
    if (!is.na(bulbing_yday)) {
      df_after_bulbing <- df[df$yday >= bulbing_yday, ]
      gdd_bulbing_cumsum <- cumsum(df_after_bulbing$GDD_daily)
      maturation_index <- which(gdd_bulbing_cumsum >= GDD_bulbing_required_p)[1]
      
      if (!is.na(maturation_index)) {
        maturation_yday <- df_after_bulbing$yday[maturation_index]
        df$in_bulbing_phase <- df$yday >= bulbing_yday & df$yday < maturation_yday
      } else {
        maturation_yday <- NA
        df$in_bulbing_phase <- FALSE
      }
    } else {
      maturation_yday <- NA
      df$in_bulbing_phase <- FALSE
    }
    
    
    # PHASE 4: Maturation + Harvest Window
    
    if (!is.na(maturation_yday)) {
      df_after_maturation <- df[df$yday >= maturation_yday, ]
      gdd_mat_cumsum <- cumsum(df_after_maturation$GDD_daily)
      harvest_index <- which(gdd_mat_cumsum >= GDD_maturation_required_p)[1]
      
      if (!is.na(harvest_index)) {
        harvest_yday <- df_after_maturation$yday[harvest_index]
        df$in_maturation_phase <- df$yday >= maturation_yday & df$yday < harvest_yday
        df$in_harvest_window <- df$yday >= harvest_yday
      } else {
        harvest_yday <- NA
        df$in_maturation_phase <- FALSE
        df$in_harvest_window <- FALSE
      }
    } else {
      harvest_yday <- NA
      df$in_maturation_phase <- FALSE
      df$in_harvest_window <- FALSE
    }
    
    
    # Final annotation
    
    df$emergence_yday <- emergence_yday
    df$bulbing_yday <- bulbing_yday
    df$maturation_yday <- maturation_yday
    df$harvest_yday <- harvest_yday
    
    return(df)
  })
  
  # Reapply names (optional but safe)
  names(weather_scenario_list) <- names(weather_scenario_names)
  
  
  
  ###### Functions for yield reducing factors #########
  
  ## Abiotic stress factors
  
  # Seedbed preparation
  #Prec in mm/day 

  get_seedbed_stress <- function(Prec,
                               day_consec_wet,
                               prec_seedbed_medium_p,
                               prec_seedbed_high_p,
                               risk_seedbed_medium_t,
                               risk_seedbed_high_t,
                               risk_additional_day = 0.05) {
  
  # base risk depending on precipitation thresholds
  if (Prec > prec_seedbed_high_p) {
    risk_seedbed <- risk_seedbed_high_t
  } else if (Prec > prec_seedbed_medium_p) {
    risk_seedbed <- risk_seedbed_medium_t
  } else {
    risk_seedbed <- 0
  }
 
  # scale risk by consecutive wet days
  risk_seedbed <- risk_seedbed + (day_consec_wet - 1) * risk_additional_day ###Rethink depending on how its calculated
  
  # cap at 1 (100%)
  risk_seedbed <- pmin(risk_seedbed, 1)
  
  return(risk_seedbed)
}

  
#Drought Stress
  
  get_drought_stress <- function(Prec,
                                 Tavg,
                                 days_consec_dry,
                                 prec_drought_threshold_p,
                                 Tavg_drought_threshold_p,         # threshold temp where drought stress starts
                                 impact_prec_drought_t,      # stakeholder input 0–100
                                 impact_days_dry_drought_t,          # stakeholder input 0.1–1
                                 impact_temp_drought_t     # stakeholder input 0.1–1
                                ) {
    
    # --- 1. Precipitation component ---
    # Active if Prec below threshold → the drier, the higher the risk
    R_P <- ifelse(Prec < prec_drought_threshold_p,
                  1 - exp(- (impact_prec_drought_t / 2) * (prec_drought_threshold_p - Prec)),
                  0)
    R_P <- pmin(R_P, 1)
    
    # --- 2. Temperature component ---
    # Higher temps = higher risk, saturating to 1
    R_T <- ifelse(Tavg > Tavg_drought_threshold_p,
                  1 - exp(- (impact_temp_drought_t / 2) * (Tavg - Tavg_drought_threshold_p)),
                  0)
    R_T <- pmin(R_T, 1)
    
    # --- 3. Consecutive dry days component ---
    R_D <- 1 - exp(- (impact_days_dry_drought_t / 2) * days_consec_dry)
    
    # --- Combined drought risk (multiplicative) ---
    risk <- R_P * R_T * R_D
    risk_drought <- pmin(risk, 1)   # saf25ety cap
    
    return(risk_drought)
  }
  
  # Hail Stress
  
  get_hail_stress <- function(Tavg,
                              Prec,
                              prec_hail_threshold_p,
                              Tavg_hail_threshold_p,
                              impact_days_hail_t) {
    
    # hail-trigger days = days with both high rainfall AND high temperature
    hail_days <- sum(Prec >= prec_hail_threshold_p &
                       Tavg >= Tavg_hail_threshold_p, na.rm = TRUE)
    
    # convert to risk with exponential saturation
    risk_hail <- 1 - exp(- (impact_days_hail_t / 2) * hail_days)
    
    return(pmin(risk_hail, 1))
  }
  
  
  # Extreme Rainfall
  get_extreme_rain_stress <- function(Prec,
                                      prec_extreme_rain_medium_p,
                                      prec_extreme_rain_high_p,
                                      impact_days_extreme_rain_t) {impact_days_extreme_rain_t
    
    # count days exceeding thresholds
    extreme_days_medium <- sum(Prec >= prec_extreme_rain_medium_p &
                                 Prec < prec_extreme_rain_high_p, na.rm = TRUE)
    extreme_days_high   <- sum(Prec >= prec_extreme_rain_high_p, na.rm = TRUE)
    
    # combine into weighted "extreme days"
    weighted_extreme_days <- extreme_days_medium + 2 * extreme_days_high
    
    # convert to risk with exponential saturation
    risk_extreme_rain <- 1 - exp(- (impact_days_extreme_rain_t / 2) * weighted_extreme_days)
    
    return(pmin(risk_extreme_rain, 1))

  }
  ## Harvest Risk
  
  
  
  ### Biotic Stress Factors
  
  ## Fungal Pathogens
  
  # Risk for Botrytis

    get_botrytis_stress <- function(Tavg,
                                    Prec,
                                    days_consec_wet,
                                    prec_botrytis_threshold_p      = 1,   # mm threshold for "rainy day"
                                    impact_prec_botrytis_t      = 50,  # 0–1; calibrates effect of daily rain amount
                                    impact_days_wet_botrytis_t      = 40,  # 0–1; calibrates effect of consecutive wet days
                                    Topt_botrytis_p        = 18,  # °C; optimum temperature for Botrytis
                                    Twidth_botrytis_p      = 4   # °C; width around optimum
                                    ) {
      
      # --- 1. Precipitation component ---
      R_P <- ifelse(Prec >= prec_botrytis_threshold_p,
                    1 - exp(- (impact_prec_botrytis_t / 2) * (Prec - prec_botrytis_threshold_p)),
                    0)
      R_P <- pmin(R_P, 1)
      
      # --- 2. Consecutive wet days component ---
      R_D <- 1 - exp(- (impact_days_wet_botrytis_t / 2) * days_consec_wet)
      R_D <- pmin(R_D, 1)
      
      # --- 3. Temperature component (optimum around 15–20 °C) ---
      R_T <- exp(- (Tavg - Topt_botrytis_p)^2 / (2 * Twidth_botrytis_p^2))
      R_T <- pmin(R_T, 1)
      
      # --- Combined risk (multiplicative) ---
      risk_botrytis <- pmin(R_P * R_D * R_T, 1)
      
      # --- Chance event ---
      return(risk_botrytis)
    }
    
  
  
    get_fusarium_stress <- function(Tavg,
                                    Prec,
                                    days_consec_wet,
                                    prec_fusarium_threshold_p,       # mm threshold for "wet" day
                                    impact_prec_fusarium_t,          # 0–1; calibrates effect of daily rain
                                    impact_days_wet_fusarium_t,      # 0–1; calibrates effect of consecutive wet days
                                    Topt_fusarium_p,                 # °C; optimum for Fusarium infection
                                    Twidth_fusarium_p               # °C; width around optimum
                                    ) {
      
      # --- 1. Precipitation component (proxy for soil/root wetness) ---
      R_P <- ifelse(Prec >= prec_fusarium_threshold_p,
                    1 - exp(- (impact_prec_fusarium_t / 2) * (Prec - prec_fusarium_threshold_p)),
                    0)
      R_P <- pmin(R_P, 1)
      
      # --- 2. Consecutive wet days component ---
      R_D <- 1 - exp(- (impact_days_wet_fusarium_t / 2) * days_consec_wet)
      R_D <- pmin(R_D, 1)
      
      # --- 3. Temperature component (optimum ~27 °C for onion basal rot) ---
      R_T <- exp(- (Tavg - Topt_fusarium_p)^2 / (2 * Twidth_fusarium_p^2))
      R_T <- pmin(R_T, 1)
      
      # --- Combined risk (multiplicative) ---
      risk_fusarium <- pmin(R_P * R_D * R_T, 1)
      
      # --- Chance event ---
      return(risk_fusarium)
    }
  
  ## Downy Mildew 
  
  get_downy_mildew_stress <- function(Tavg,
                                      Prec,
                                      days_consec_wet,
                                      prec_mildew_threshold_p,      # mm threshold for "wet" day
                                      impact_prec_mildew_t,       # 0–1; calibrates effect of daily rain
                                      impact_days_wet_mildew_t,   # 0–1; calibrates effect of consecutive wet days
                                      Topt_mildew_p,              # °C; optimum for downy mildew infection
                                      Twidth_mildew_p            # °C; width around optimum
                                      ) {
    
    # --- 1. Precipitation component (proxy for leaf wetness) ---
    R_P <- ifelse(Prec >= prec_mildew_threshold_p,
                  1 - exp(- (impact_prec_mildew_t / 2) * (Prec - prec_mildew_threshold_p)),
                  0)
    R_P <- pmin(R_P, 1)
    
    # --- 2. Consecutive wet days component ---
    R_D <- 1 - exp(- (impact_days_wet_mildew_t / 2) * days_consec_wet)
    R_D <- pmin(R_D, 1)
    
    # --- 3. Temperature component (optimum ~12 °C) ---
    R_T <- exp(- (Tavg - Topt_mildew_p)^2 / (2 * Twidth_mildew_p^2))
    R_T <- pmin(R_T, 1)
    
    # --- Combined risk (multiplicative) ---
    risk_mildew <- pmin(R_P * R_D * R_T, 1)
    
    # --- Chance event ---
    return(risk_mildew)
  }
  
  
  ## Animal Pressure
  get_thrips_stress <- function(Tavg,
                                Prec,
                                days_consec_dry,
                                prec_thrips_threshold_p,        # mm/day; rainfall threshold below which dryness contributes
                                Topt_thrips_p,                  # °C; optimum for thrips development
                                Twidth_thrips_p,                # °C; width around optimum
                                impact_prec_thrips_t,           # 0–1; effect of dryness
                                impact_days_dry_thrips_t       # 0–1; effect of consecutive dry days
                                ) {     # proportional yield loss
    
    # --- 1. Precipitation component (dryness effect) ---
    R_P <- ifelse(Prec < prec_thrips_threshold_p,
                  1 - exp(- (impact_prec_thrips_t / 2) * (prec_thrips_threshold_p - Prec)),
                  0)
    R_P <- pmin(R_P, 1)
    
    # --- 2. Temperature component (optimum for thrips ~28 °C) ---
    R_T <- exp(- (Tavg - Topt_thrips_p)^2 / (2 * Twidth_thrips_p^2))
    R_T <- pmin(R_T, 1)
    
    # --- 3. Consecutive dry days component ---
    R_D <- 1 - exp(- (impact_days_dry_thrips_t / 2) * days_consec_dry)
    R_D <- pmin(R_D, 1)
  
  # --- Combined risk ---
  risk_thrips <- pmin(R_P * R_T * R_D, 1)
  
  return(risk_thrips)
}

  
  get_onion_fly_stress <- function(Tavg,
                                   Prec,
                                   days_consec_wet,
                                   prec_onion_fly_threshold_p,   # mm/day; rainfall threshold for wetness
                                   Topt_onion_fly_p,             # °C; optimum for activity
                                   Twidth_onion_fly_p,           # °C; width around optimum
                                   impact_prec_onion_fly_t,      # 0–1; effect of wetness
                                   impact_days_wet_onion_fly_t  # 0–1; effect of consecutive wet days
                                   ) {
    
    # --- 1. Precipitation component (wetness effect) ---
    R_P <- ifelse(Prec >= prec_onion_fly_threshold_p,
                  1 - exp(- (impact_prec_onion_fly_t / 2) * (Prec - prec_onion_fly_threshold_p)),
                  0)
    R_P <- pmin(R_P, 1)
    
    # --- 2. Temperature component (optimum ~23 °C) ---
    R_T <- exp(- (Tavg - Topt_onion_fly_p)^2 / (2 * Twidth_onion_fly_p^2))
    R_T <- pmin(R_T, 1)
    
    # --- 3. Consecutive wet days ---
    R_D <- 1 - exp(- (impact_days_wet_onion_fly_t / 2) * days_consec_wet)
    R_D <- pmin(R_D, 1)
    
    # --- Combined risk ---
    risk_onion_fly <- pmin(R_P * R_T * R_D, 1)
    
    return(risk_onion_fly * yield_reduction_onion_fly_t)
  }
  
  
  get_wireworm_stress <- function(Tavg,
                                  Prec,
                                  days_consec_wet,
                                  prec_wireworm_threshold_p,     # mm/day; rainfall threshold for wetness
                                  Topt_wireworm_p,               # °C; optimum for activity
                                  Twidth_wireworm_p,             # °C; width around optimum
                                  impact_prec_wireworm_t,        # 0–1; effect of wetness
                                  impact_days_wet_wireworm_t    # 0–1; effect of consecutive wet days
                                  ) {
    
    # --- 1. Precipitation component (wetness effect) ---
    R_P <- ifelse(Prec >= prec_wireworm_threshold_p,
                  1 - exp(- (impact_prec_wireworm_t / 2) * (Prec - prec_wireworm_threshold_p)),
                  0)
    R_P <- pmin(R_P, 1)
    
    # --- 2. Temperature component (optimum ~18 °C) ---
    R_T <- exp(- (Tavg - Topt_wireworm_p)^2 / (2 * Twidth_wireworm_p^2))
    R_T <- pmin(R_T, 1)
    
    # --- 3. Consecutive wet days ---
    R_D <- 1 - exp(- (impact_days_wet_wireworm_t / 2) * days_consec_wet)
    R_D <- pmin(R_D, 1)
    
    # --- Combined risk ---
    risk_wireworm <- pmin(R_P * R_T * R_D, 1)
    
    return(risk_wireworm * yield_reduction_wireworm_t)
  }
  
  
  #####  Apply stresss functions
  #####  Apply stress functions
  season_risks <- lapply(weather_scenario_list, function(df) {
    
    safe_mean <- function(x, condition) {
      if (any(condition, na.rm = TRUE)) mean(x[condition], na.rm = TRUE) else 0
    }
    safe_max  <- function(x, condition) {
      if (any(condition, na.rm = TRUE)) max(x[condition], na.rm = TRUE) else 0
    }
    
    ##### --- EMERGENCE PHASE ---
    emergence_filter <- df$in_emergence_phase == TRUE
    emergence_Tavg <- safe_mean(df$Tavg, emergence_filter)
    emergence_Prec <- safe_mean(df$Prec, emergence_filter)
    emergence_consec_wet <- safe_max(df$day_consec_wet, emergence_filter)
    
    # risks
    emergence_seed_bed_risk <- get_seedbed_stress(emergence_Prec, emergence_consec_wet,
                                                  prec_seedbed_medium_p,
                                                  prec_seedbed_high_p,
                                                  risk_seedbed_medium_t,
                                                  risk_seedbed_high_t)
    emergence_drought_risk <- get_drought_stress(emergence_Prec, emergence_Tavg, emergence_consec_wet,
                                                 prec_drought_threshold_p,
                                                 Tavg_drought_threshold_p,
                                                 impact_prec_drought_t,
                                                 impact_days_dry_drought_t,
                                                 impact_temp_drought_t)
    emergence_extreme_rain_risk <- get_extreme_rain_stress(df$Prec[emergence_filter],
                                                           prec_extreme_rain_medium_p,
                                                           prec_extreme_rain_high_p,
                                                           impact_days_extreme_rain_t)
    emergence_hail_risk <- get_hail_stress(df$Tavg[emergence_filter],
                                           df$Prec[emergence_filter],
                                           prec_hail_threshold_p,
                                           Tavg_hail_threshold_p,
                                           impact_days_hail_t)
    emergence_fusarium_risk <- get_fusarium_stress(emergence_Tavg, emergence_Prec, emergence_consec_wet,
                                                   prec_fusarium_threshold_p,
                                                   impact_prec_fusarium_t,
                                                   impact_days_wet_fusarium_t,
                                                   Topt_fusarium_p,
                                                   Twidth_fusarium_p)
    emergence_onion_fly_risk <- get_onion_fly_stress(emergence_Tavg, emergence_Prec, emergence_consec_wet,
                                                     prec_onion_fly_threshold_p,
                                                     Topt_onion_fly_p,
                                                     Twidth_onion_fly_p,
                                                     impact_prec_onion_fly_t,
                                                     impact_days_wet_onion_fly_t)
    emergence_wireworm_risk <- get_wireworm_stress(emergence_Tavg, emergence_Prec, emergence_consec_wet,
                                                   prec_wireworm_threshold_p,
                                                   Topt_wireworm_p,
                                                   Twidth_wireworm_p,
                                                   impact_prec_wireworm_t,
                                                   impact_days_wet_wireworm_t)
    
    # yield reductions
    yield_reduction_seedbed_emergence <- chance_event(chance = emergence_seed_bed_risk, value_if = yield_reduction_seedbed_t, value_if_not = 1)
    yield_reduction_drought_emergence <- chance_event(chance = emergence_drought_risk, value_if = yield_reduction_drought_t, value_if_not = 1)
    yield_reduction_extreme_rain_emergence <- chance_event(chance = emergence_extreme_rain_risk, value_if = yield_reduction_extreme_rain_t, value_if_not = 1)
    yield_reduction_hail_emergence <- chance_event(chance = emergence_hail_risk, value_if = yield_reduction_hail_p, value_if_not = 1)
    yield_reduction_fusarium_emergence <- chance_event(chance = emergence_fusarium_risk, value_if = yield_reduction_fusarium_t, value_if_not = 1)
    yield_reduction_onion_fly_emergence <- chance_event(chance = emergence_onion_fly_risk, value_if = yield_reduction_onion_fly_t, value_if_not = 1)
    yield_reduction_wireworm_emergence <- chance_event(chance = emergence_wireworm_risk, value_if = yield_reduction_wireworm_t, value_if_not = 1)
    
    ##### --- VEGETATIVE PHASE ---
    vegetative_filter <- df$in_vegetative_phase == TRUE
    vegetative_Tavg <- safe_mean(df$Tavg, vegetative_filter)
    vegetative_Prec <- safe_mean(df$Prec, vegetative_filter)
    vegetative_consec_wet <- safe_max(df$day_consec_wet, vegetative_filter)
    
    vegetative_drought_risk <- get_drought_stress(vegetative_Prec, vegetative_Tavg, vegetative_consec_wet,
                                                  prec_drought_threshold_p,
                                                  Tavg_drought_threshold_p,
                                                  impact_prec_drought_t,
                                                  impact_days_dry_drought_t,
                                                  impact_temp_drought_t)
    vegetative_extreme_rain_risk <- get_extreme_rain_stress(df$Prec[vegetative_filter],
                                                            prec_extreme_rain_medium_p,
                                                            prec_extreme_rain_high_p,
                                                            impact_days_extreme_rain_t)
    vegetative_hail_risk <- get_hail_stress(df$Tavg[vegetative_filter],
                                            df$Prec[vegetative_filter],
                                            prec_hail_threshold_p,
                                            Tavg_hail_threshold_p,
                                            impact_days_hail_t)
    vegetative_fusarium_risk <- get_fusarium_stress(vegetative_Tavg, vegetative_Prec, vegetative_consec_wet,
                                                    prec_fusarium_threshold_p,
                                                    impact_prec_fusarium_t,
                                                    impact_days_wet_fusarium_t,
                                                    Topt_fusarium_p,
                                                    Twidth_fusarium_p)
    vegetative_downy_mildew_risk <- get_downy_mildew_stress(vegetative_Tavg, vegetative_Prec, vegetative_consec_wet,
                                                            prec_mildew_threshold_p,
                                                            impact_prec_mildew_t,
                                                            impact_days_wet_mildew_t,
                                                            Topt_mildew_p,
                                                            Twidth_mildew_p)
    vegetative_thrips_risk <- get_thrips_stress(vegetative_Tavg, vegetative_Prec, vegetative_consec_wet,
                                                prec_thrips_threshold_p,
                                                Topt_thrips_p,
                                                Twidth_thrips_p,
                                                impact_prec_thrips_t,
                                                impact_days_dry_thrips_t)
    vegetative_onion_fly_risk <- get_onion_fly_stress(vegetative_Tavg, vegetative_Prec, vegetative_consec_wet,
                                                      prec_onion_fly_threshold_p,
                                                      Topt_onion_fly_p,
                                                      Twidth_onion_fly_p,
                                                      impact_prec_onion_fly_t,
                                                      impact_days_wet_onion_fly_t)
    
    # yield reductions
    yield_reduction_drought_vegetative <- chance_event(chance = vegetative_drought_risk, value_if = yield_reduction_drought_t, value_if_not = 1)
    yield_reduction_extreme_rain_vegetative <- chance_event(chance = vegetative_extreme_rain_risk, value_if = yield_reduction_extreme_rain_t, value_if_not = 1)
    yield_reduction_hail_vegetative <- chance_event(chance = vegetative_hail_risk, value_if = yield_reduction_hail_p, value_if_not = 1)
    yield_reduction_fusarium_vegetative <- chance_event(chance = vegetative_fusarium_risk, value_if = yield_reduction_fusarium_t, value_if_not = 1)
    yield_reduction_downy_mildew_vegetative <- chance_event(chance = vegetative_downy_mildew_risk, value_if = yield_reduction_downy_mildew_t, value_if_not = 1)
    yield_reduction_thrips_vegetative <- chance_event(chance = vegetative_thrips_risk, value_if = yield_reduction_thrips_t, value_if_not = 1)
    yield_reduction_onion_fly_vegetative <- chance_event(chance = vegetative_onion_fly_risk, value_if = yield_reduction_onion_fly_t, value_if_not = 1)
    
    ##### --- BULBING PHASE ---
    bulbing_filter <- df$in_bulbing_phase == TRUE
    bulbing_Tavg <- safe_mean(df$Tavg, bulbing_filter)
    bulbing_Prec <- safe_mean(df$Prec, bulbing_filter)
    bulbing_consec_wet <- safe_max(df$day_consec_wet, bulbing_filter)
    
    bulbing_drought_risk <- get_drought_stress(bulbing_Prec, bulbing_Tavg, bulbing_consec_wet,
                                               prec_drought_threshold_p,
                                               Tavg_drought_threshold_p,
                                               impact_prec_drought_t,
                                               impact_days_dry_drought_t,
                                               impact_temp_drought_t)
    bulbing_extreme_rain_risk <- get_extreme_rain_stress(df$Prec[bulbing_filter],
                                                         prec_extreme_rain_medium_p,
                                                         prec_extreme_rain_high_p,
                                                         impact_days_extreme_rain_t)
    bulbing_hail_risk <- get_hail_stress(df$Tavg[bulbing_filter],
                                         df$Prec[bulbing_filter],
                                         prec_hail_threshold_p,
                                         Tavg_hail_threshold_p,
                                         impact_days_hail_t)
    bulbing_botrytis_risk <- get_botrytis_stress(bulbing_Tavg, bulbing_Prec, bulbing_consec_wet,
                                                 prec_botrytis_threshold_p,
                                                 impact_prec_botrytis_t,
                                                 impact_days_wet_botrytis_t,
                                                 Topt_botrytis_p,
                                                 Twidth_botrytis_p)
    bulbing_fusarium_risk <- get_fusarium_stress(bulbing_Tavg, bulbing_Prec, bulbing_consec_wet,
                                                 prec_fusarium_threshold_p,
                                                 impact_prec_fusarium_t,
                                                 impact_days_wet_fusarium_t,
                                                 Topt_fusarium_p,
                                                 Twidth_fusarium_p)
    bulbing_downy_mildew_risk <- get_downy_mildew_stress(bulbing_Tavg, bulbing_Prec, bulbing_consec_wet,
                                                         prec_mildew_threshold_p,
                                                         impact_prec_mildew_t,
                                                         impact_days_wet_mildew_t,
                                                         Topt_mildew_p,
                                                         Twidth_mildew_p)
    bulbing_onion_fly_risk <- get_onion_fly_stress(bulbing_Tavg, bulbing_Prec, bulbing_consec_wet,
                                                   prec_onion_fly_threshold_p,
                                                   Topt_onion_fly_p,
                                                   Twidth_onion_fly_p,
                                                   impact_prec_onion_fly_t,
                                                   impact_days_wet_onion_fly_t)
    
    # yield reductions
    yield_reduction_drought_bulbing <- chance_event(chance = bulbing_drought_risk, value_if = yield_reduction_drought_t, value_if_not = 1)
    yield_reduction_extreme_rain_bulbing <- chance_event(chance = bulbing_extreme_rain_risk, value_if = yield_reduction_extreme_rain_t, value_if_not = 1)
    yield_reduction_hail_bulbing <- chance_event(chance = bulbing_hail_risk, value_if = yield_reduction_hail_p, value_if_not = 1)
    yield_reduction_botrytis_bulbing <- chance_event(chance = bulbing_botrytis_risk, value_if = yield_reduction_botrytis_t, value_if_not = 1)
    yield_reduction_fusarium_bulbing <- chance_event(chance = bulbing_fusarium_risk, value_if = yield_reduction_fusarium_t, value_if_not = 1)
    yield_reduction_downy_mildew_bulbing <- chance_event(chance = bulbing_downy_mildew_risk, value_if = yield_reduction_downy_mildew_t, value_if_not = 1)
    yield_reduction_onion_fly_bulbing <- chance_event(chance = bulbing_onion_fly_risk, value_if = yield_reduction_onion_fly_t, value_if_not = 1)
    
    ##### --- MATURATION PHASE ---
    maturation_filter <- df$in_maturation_phase == TRUE
    maturation_Tavg <- safe_mean(df$Tavg, maturation_filter)
    maturation_Prec <- safe_mean(df$Prec, maturation_filter)
    maturation_consec_wet <- safe_max(df$day_consec_wet, maturation_filter)
    
    maturation_drought_risk <- get_drought_stress(maturation_Prec, maturation_Tavg, maturation_consec_wet,
                                                  prec_drought_threshold_p,
                                                  Tavg_drought_threshold_p,
                                                  impact_prec_drought_t,
                                                  impact_days_dry_drought_t,
                                                  impact_temp_drought_t)
    maturation_extreme_rain_risk <- get_extreme_rain_stress(df$Prec[maturation_filter],
                                                            prec_extreme_rain_medium_p,
                                                            prec_extreme_rain_high_p,
                                                            impact_days_extreme_rain_t)
    maturation_hail_risk <- get_hail_stress(df$Tavg[maturation_filter],
                                            df$Prec[maturation_filter],
                                            prec_hail_threshold_p,
                                            Tavg_hail_threshold_p,
                                            impact_days_hail_t)
    maturation_botrytis_risk <- get_botrytis_stress(maturation_Tavg, maturation_Prec, maturation_consec_wet,
                                                    prec_botrytis_threshold_p,
                                                    impact_prec_botrytis_t,
                                                    impact_days_wet_botrytis_t,
                                                    Topt_botrytis_p,
                                                    Twidth_botrytis_p)
    maturation_fusarium_risk <- get_fusarium_stress(maturation_Tavg, maturation_Prec, maturation_consec_wet,
                                                    prec_fusarium_threshold_p,
                                                    impact_prec_fusarium_t,
                                                    impact_days_wet_fusarium_t,
                                                    Topt_fusarium_p,
                                                    Twidth_fusarium_p)
    maturation_downy_mildew_risk <- get_downy_mildew_stress(maturation_Tavg, maturation_Prec, maturation_consec_wet,
                                                            prec_mildew_threshold_p,
                                                            impact_prec_mildew_t,
                                                            impact_days_wet_mildew_t,
                                                            Topt_mildew_p,
                                                            Twidth_mildew_p)
    maturation_onion_fly_risk <- get_onion_fly_stress(maturation_Tavg, maturation_Prec, maturation_consec_wet,
                                                      prec_onion_fly_threshold_p,
                                                      Topt_onion_fly_p,
                                                      Twidth_onion_fly_p,
                                                      impact_prec_onion_fly_t,
                                                      impact_days_wet_onion_fly_t)
    
    # yield reductions
    yield_reduction_drought_maturation <- chance_event(chance = maturation_drought_risk, value_if = yield_reduction_drought_t, value_if_not = 1)
    yield_reduction_extreme_rain_maturation <- chance_event(chance = maturation_extreme_rain_risk, value_if = yield_reduction_extreme_rain_t, value_if_not = 1)
    yield_reduction_hail_maturation <- chance_event(chance = maturation_hail_risk, value_if = yield_reduction_hail_p, value_if_not = 1)
    yield_reduction_botrytis_maturation <- chance_event(chance = maturation_botrytis_risk, value_if = yield_reduction_botrytis_t, value_if_not = 1)
    yield_reduction_fusarium_maturation <- chance_event(chance = maturation_fusarium_risk, value_if = yield_reduction_fusarium_t, value_if_not = 1)
    yield_reduction_downy_mildew_maturation <- chance_event(chance = maturation_downy_mildew_risk, value_if = yield_reduction_downy_mildew_t, value_if_not = 1)
    yield_reduction_onion_fly_maturation <- chance_event(chance = maturation_onion_fly_risk, value_if = yield_reduction_onion_fly_t, value_if_not = 1)
    
    ##### --- RETURN LIST ---
    return(list(
      # Emergence
      emergence_seed_bed_risk = emergence_seed_bed_risk,
      emergence_drought_risk = emergence_drought_risk,
      emergence_extreme_rain_risk = emergence_extreme_rain_risk,
      emergence_hail_risk = emergence_hail_risk,
      emergence_fusarium_risk = emergence_fusarium_risk,
      emergence_onion_fly_risk = emergence_onion_fly_risk,
      emergence_wireworm_risk = emergence_wireworm_risk,
      yield_reduction_seedbed_emergence = yield_reduction_seedbed_emergence,
      yield_reduction_drought_emergence = yield_reduction_drought_emergence,
      yield_reduction_extreme_rain_emergence = yield_reduction_extreme_rain_emergence,
      yield_reduction_hail_emergence = yield_reduction_hail_emergence,
      yield_reduction_fusarium_emergence = yield_reduction_fusarium_emergence,
      yield_reduction_onion_fly_emergence = yield_reduction_onion_fly_emergence,
      yield_reduction_wireworm_emergence = yield_reduction_wireworm_emergence,
      
      # Vegetative
      vegetative_drought_risk = vegetative_drought_risk,
      vegetative_extreme_rain_risk = vegetative_extreme_rain_risk,
      vegetative_hail_risk = vegetative_hail_risk,
      vegetative_fusarium_risk = vegetative_fusarium_risk,
      vegetative_downy_mildew_risk = vegetative_downy_mildew_risk,
      vegetative_thrips_risk = vegetative_thrips_risk,
      vegetative_onion_fly_risk = vegetative_onion_fly_risk,
      yield_reduction_drought_vegetative = yield_reduction_drought_vegetative,
      yield_reduction_extreme_rain_vegetative = yield_reduction_extreme_rain_vegetative,
      yield_reduction_hail_vegetative = yield_reduction_hail_vegetative,
      yield_reduction_fusarium_vegetative = yield_reduction_fusarium_vegetative,
      yield_reduction_downy_mildew_vegetative = yield_reduction_downy_mildew_vegetative,
      yield_reduction_thrips_vegetative = yield_reduction_thrips_vegetative,
      yield_reduction_onion_fly_vegetative = yield_reduction_onion_fly_vegetative,
      
      # Bulbing
      bulbing_drought_risk = bulbing_drought_risk,
      bulbing_extreme_rain_risk = bulbing_extreme_rain_risk,
      bulbing_hail_risk = bulbing_hail_risk,
      bulbing_botrytis_risk = bulbing_botrytis_risk,
      bulbing_fusarium_risk = bulbing_fusarium_risk,
      bulbing_downy_mildew_risk = bulbing_downy_mildew_risk,
      bulbing_onion_fly_risk = bulbing_onion_fly_risk,
      yield_reduction_drought_bulbing = yield_reduction_drought_bulbing,
      yield_reduction_extreme_rain_bulbing = yield_reduction_extreme_rain_bulbing,
      yield_reduction_hail_bulbing = yield_reduction_hail_bulbing,
      yield_reduction_botrytis_bulbing = yield_reduction_botrytis_bulbing,
      yield_reduction_fusarium_bulbing = yield_reduction_fusarium_bulbing,
      yield_reduction_downy_mildew_bulbing = yield_reduction_downy_mildew_bulbing,
      yield_reduction_onion_fly_bulbing = yield_reduction_onion_fly_bulbing,
      
      # Maturation
      maturation_drought_risk = maturation_drought_risk,
      maturation_extreme_rain_risk = maturation_extreme_rain_risk,
      maturation_hail_risk = maturation_hail_risk,
      maturation_botrytis_risk = maturation_botrytis_risk,
      maturation_fusarium_risk = maturation_fusarium_risk,
      maturation_downy_mildew_risk = maturation_downy_mildew_risk,
      maturation_onion_fly_risk = maturation_onion_fly_risk,
      yield_reduction_drought_maturation = yield_reduction_drought_maturation,
      yield_reduction_extreme_rain_maturation = yield_reduction_extreme_rain_maturation,
      yield_reduction_hail_maturation = yield_reduction_hail_maturation,
      yield_reduction_botrytis_maturation = yield_reduction_botrytis_maturation,
      yield_reduction_fusarium_maturation = yield_reduction_fusarium_maturation,
      yield_reduction_downy_mildew_maturation = yield_reduction_downy_mildew_maturation,
      yield_reduction_onion_fly_maturation = yield_reduction_onion_fly_maturation
    ))
  })
  
  
  
  
  
  # Combine weather data and seasonal risk info
  weather_scenario_list <- Map(function(weather_df, risks) {
    # Add the risk variables as new columns (repeated to match number of rows)
    for (name in names(risks)) {
      weather_df[[name]] <- risks[[name]]
    }
    return(weather_df)
  }, weather_scenario_list, season_risks)
  
  
  
  ## Calculate Biomass Growth per Phase
  
  calculate_biomass_daily_emergence <- function(PAR, LAI_emergence_p, Tavg, Prec, emergence_phase_stress, LUE_onion_p) {
    
    
    # Temperature effect
    f_T <- ifelse(Tavg >= f_T_1_lower_p & Tavg <= f_T_1_upper_p, 1,
                  ifelse(Tavg < f_T_0_lower_p | Tavg > f_T_0_upper_p, 0, 0.5))
    
    # Water availability
    f_W <- ifelse(Prec >= f_W_1_lower_p & Prec <= f_W_1_upper_p, 1,
                  ifelse(Prec < f_W_0.5_p, 0.5, 0.7))
    
  
    
    
    # Vectorized light interception & biomass increment
    delta_B <- LUE_onion_p * PAR * (1 - exp(-lec_k_c * LAI_emergence_p)) * f_T * f_W 
    return(delta_B)
  }
  
  calculate_biomass_daily_vegetative <- function(PAR, LAI_veg_p, Tavg,
                                                 Prec, vegetative_phase_stress, LUE_onion_p) {
    
    # Temperature effect
    f_T <- ifelse(Tavg >= f_T_1_lower_p & Tavg <= f_T_1_upper_p, 1,
                  ifelse(Tavg < f_T_0_lower_p | Tavg > f_T_0_upper_p, 0, 0.5))
    
    # Water availability
    f_W <- ifelse(Prec >= f_W_1_lower_p & Prec <= f_W_1_upper_p, 1,
                  ifelse(Prec < f_W_0.5_p, 0.5, 0.7))
    
    
    delta_B <- LUE_onion_p * PAR * (1 - exp(-lec_k_c * LAI_veg_p)) * f_T * f_W 
    return(delta_B)
  }
  
  calculate_biomass_daily_bulbing <- function(PAR, LAI_bulbing_p, Tavg, Prec, bulbing_phase_stress, LUE_onion_p) {
    
    
    # Temperature effect
    f_T <- ifelse(Tavg >= f_T_1_lower_p & Tavg <= f_T_1_upper_p, 1,
                  ifelse(Tavg < f_T_0_lower_p | Tavg > f_T_0_upper_p, 0, 0.5))
    
    # Water availability
    f_W <- ifelse(Prec >= f_W_1_lower_p & Prec <= f_W_1_upper_p, 1,
                  ifelse(Prec < f_W_0.5_p, 0.5, 0.7))
   
    
    delta_B <- LUE_onion_p * PAR * (1 - exp(-lec_k_c * LAI_bulbing_p)) * f_T * f_W 
    return(delta_B)
  }
  
  calculate_biomass_daily_maturation <- function(PAR, LAI_maturation_p, Tavg, Prec, maturation_phase_stress, LUE_onion_p) {
    
    
    # Temperature effect
    f_T <- ifelse(Tavg >= f_T_1_lower_p & Tavg <= f_T_1_upper_p, 1,
                  ifelse(Tavg < f_T_0_lower_p | Tavg > f_T_0_upper_p, 0, 0.5))
    
    # Water availability
    f_W <- ifelse(Prec >= f_W_1_lower_p & Prec <= f_W_1_upper_p, 1,
                  ifelse(Prec < f_W_0.5_p, 0.5, 0.7))
    
  
    
    delta_B <- LUE_onion_p * PAR * (1 - exp(-lec_k_c * LAI_maturation_p)) * f_T * f_W 
    return(delta_B)
  }
  
  get_scalar <- function(df, column, filter) {
    if (column %in% names(df) && any(filter, na.rm = TRUE)) {
      value <- unique(df[[column]][filter])
      # Handles cases where the filter might not find any data, or if there are multiple unique values
      if (length(value) == 1) {
        return(value)
      } else {
        return(0)
      }
    } else {
      return(0)
    }
  }
  
  
  #### Apply functions for all phases
  
  #### Apply functions for all phases
  
  biomass_all_scenarios <- lapply(weather_scenario_list, function(df) {
    
    # Emergence phase
    df_emergence <- df[df$in_emergence_phase, ]
    biomass_emergence <- if (nrow(df_emergence) > 0) {
      sum(calculate_biomass_daily_emergence(
        PAR = df_emergence$PAR,
        LAI_emergence_p = LAI_emergence_p,
        Tavg = df_emergence$Tavg,
        Prec = df_emergence$Prec,
        LUE_onion_p = LUE_onion_p
      ), na.rm = TRUE)
    } else 0
    
    # Vegetative phase
    df_veg <- df[df$in_vegetative_phase, ]
    biomass_veg <- if (nrow(df_veg) > 0) {
      sum(calculate_biomass_daily_vegetative(
        PAR = df_veg$PAR,
        LAI_veg_p = LAI_veg_p,
        Tavg = df_veg$Tavg,
        Prec = df_veg$Prec,
        LUE_onion_p = LUE_onion_p
      ), na.rm = TRUE)
    } else 0
    
    # Bulbing phase
    df_bulbing <- df[df$in_bulbing_phase, ]
    biomass_bulbing <- if (nrow(df_bulbing) > 0) {
      sum(calculate_biomass_daily_bulbing(
        PAR = df_bulbing$PAR,
        LAI_bulbing_p = LAI_bulbing_p,
        Tavg = df_bulbing$Tavg,
        Prec = df_bulbing$Prec,
        LUE_onion_p = LUE_onion_p
      ), na.rm = TRUE)
    } else 0
    
    # Maturation phase
    df_maturation <- df[df$in_maturation_phase, ]
    biomass_maturation <- if (nrow(df_maturation) > 0) {
      sum(calculate_biomass_daily_maturation(
        PAR = df_maturation$PAR,
        LAI_maturation_p = LAI_maturation_p,
        Tavg = df_maturation$Tavg,
        Prec = df_maturation$Prec,
        LUE_onion_p = LUE_onion_p
      ), na.rm = TRUE)
    } else 0
    
    # Combine biomass across phases
    total_biomass_current_scenario <- biomass_emergence + biomass_veg + biomass_bulbing + biomass_maturation
    
    # Raw yield before stress
    raw_yield_per_ha <- (total_biomass_current_scenario * onions_per_ha_p * (1 - dry_onion_weight_p)) * HI_onions_p / 1000 / 1000
    
    # Collect all yield reduction multipliers (== value_if / value_if_not from chance_event)
    yield_reductions <- df[1, grepl("^yield_reduction_", names(df))]
    
    # Combined multiplier: product of all
    combined_yield_multiplier <- prod(unlist(yield_reductions), na.rm = TRUE)
    
    # Final yield after applying stress
    final_yield_per_ha <- raw_yield_per_ha * combined_yield_multiplier
    
    # Return everything
    final_output <- c(
      list(
        emergence_biomass = biomass_emergence,
        vegetative_biomass = biomass_veg,
        bulbing_biomass = biomass_bulbing,
        maturation_biomass = biomass_maturation,
        total_biomass = total_biomass_current_scenario,
        raw_yield_per_ha = raw_yield_per_ha,
        final_yield_per_ha = final_yield_per_ha,
        combined_yield_multiplier = combined_yield_multiplier
      ),
      as.list(yield_reductions)  # keep all individual multipliers
    )
    
    return(final_output)
  })
}

# Run the Monte Carlo simulation using the model function
onion_mc_simulation <- mcSimulation(estimate = as.estimate(input_variables),
                                    model_function = onion_climate_impact,
                                    numberOfModelRuns = num_simulations_c,
                                    functionSyntax = "plainNames")

