# Climate Impact Onion Model
library(chillR)
library(tidyverse)
library(decisionSupport)


flist <- list.files('future_weather/', full.names = TRUE)

future_weather <- purrr::map(flist, function(file_path) {
  
  # Extract filename
  filename <- basename(file_path) # just the file name, no directory
  
  # Split by '.' and extract parts
  name_parts <- str_split(filename, '\\.')[[1]]
  
  # Check if name_parts is long enough
  if (length(name_parts) < 5) {
    stop(paste("Filename", filename, "does not have enough parts."))
  }
  
  # Read and process the file
  read.csv(file_path) %>%
    select(-X) %>%
    mutate(
      ssp = name_parts[3],
      gcm = name_parts[4],
      scenario_year = name_parts[5],
      yday = lubridate::yday(DATE),
      season = ifelse(yday >= 0, Year + 1, Year)
    )
})

future_weather <- do.call('rbind', future_weather) 
future_weather$id = paste(future_weather$ssp, future_weather$gcm, future_weather$scenario_year, sep = '--')


#weed out incomplete seasons
drop_list <- future_weather %>% 
  group_by(id, season) %>% 
  summarise(n = n()) %>% 
  filter(n < 365)

future_weather <- future_weather %>% 
  filter(!(paste(id, season) %in% paste(drop_list$id, drop_list$season)))

hist_weather <- read.csv('weather_2020_koeln-bonn.csv') %>% 
  mutate(scenario_year = 2020, 
         ssp = 'historical',
         gcm = 'historical',
         yday = lubridate::yday(DATE),
         season = ifelse(yday >= 0,
                         yes = Year  +1,
                         no = Year),
         id = paste(ssp, gcm, scenario_year, sep = '--'))
drop_list <- hist_weather %>% 
  group_by(id, season) %>% 
  summarise(n = n()) %>% 
  filter(n < 365)
hist_weather <- hist_weather %>% 
  filter(!(paste(id, season) %in% paste(drop_list$id, drop_list$season)))

#combine hist and future weather
weather_combined <- future_weather %>% 
  rbind(hist_weather)

#add unique id to each seaoson
weather_combined$id_seaon <- paste(weather_combined$id, weather_combined$season, sep = '--')

#write.csv(weather_combined, 'weather_onion_koeln-bonn.csv')

weather_combined <- read.csv("weather_onion_koeln-bonn.csv")

weather_combined$Tavg <- ( weather_combined$Tmax + weather_combined$Tmin ) / 2

# Initialize the new column with zeros
weather_combined$day_consec_dry <- 0

dry_days <- weather_combined$Prec == 0
rle_dry <- rle(dry_days)
consec_dry <- unlist(lapply(seq_along(rle_dry$lengths), function(i) {
  if (rle_dry$values[i]) {
    seq_len(rle_dry$lengths[i])
  } else {
    rep(0, rle_dry$lengths[i])
  }
}))
weather_combined$day_consec_dry <- consec_dry


weather_combined$day_consec_wet <- 0

# Logical vector: TRUE wenn wet (regen), FALSE sonst
wet_days <- weather_combined$Prec > 0

# run length encoding auf wet_days
rle_wet <- rle(wet_days)


consec_wet <- inverse.rle(list(
  lengths = rle_wet$lengths,
  values = ifelse(rle_wet$values, seq_len(max(rle_wet$lengths)), 0)
))


consec_wet <- unlist(lapply(seq_along(rle_wet$lengths), function(i) {
  if (rle_wet$values[i]) {
    seq_len(rle_wet$lengths[i])
  } else {
    rep(0, rle_wet$lengths[i])
  }
}))

weather_combined$day_consec_wet <- consec_wet



#split by each season
weather_list <- split(x = weather_combined, f = weather_combined$id_seaon) 

# Define the scenario prefixes and corresponding variable names
scenarios <- c("historical", "ssp126", "ssp245", "ssp370", "ssp585")



# Weather Lists justs needs to be calculated once, model picks different random season each run.



##### Input Table Onion model


input_variables <- read.csv("input_table_onion.csv", header = TRUE, sep = ";")

make_variables <- function(est,n=1)
{ x<-random(rho=est, n=n)
for(i in colnames(x)) assign(i,as.numeric(x[1,i]),envir=.GlobalEnv)
} 

make_variables(as.estimate(input_variables))

onion_climate_impact <- function(){

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
  calc_Ra <- function(latitude_koeln_bonn, yday) {
    Gsc <- 0.0820
    phi <- latitude_koeln_bonn * pi / 180
    dr <- 1 + 0.033 * cos(2 * pi * yday / 365)
    delta <- 0.409 * sin((2 * pi * yday / 365) - 1.39)
    omega_s <- acos(-tan(phi) * tan(delta))
    Ra <- (24 * 60 / pi) * Gsc * dr *
      (omega_s * sin(phi) * sin(delta) + cos(phi) * cos(delta) * sin(omega_s))
    return(Ra)
  }
  
  calc_Rs <- function(Ra, deltaT) {
    return(Krs * Ra * deltaT)
  }
  
  
  # Process each scenario
  weather_scenario_list <- lapply(names(weather_scenario_list), function(scenario) {
    df <- weather_scenario_list[[scenario]]
    
    # --- Calculate PAR ---
    Ra <- calc_Ra(latitude_koeln_bonn, df$yday)
    deltaT <- df$Tmax - df$Tmin
    Rs <- calc_Rs(Ra, deltaT)
    df$PAR <- Rs * 0.45
    
    # --- Calculate daily GDD ---
    df$GDD_daily <- pmax(0, ((df$Tmax + df$Tmin) / 2) - base_temp)
    
    # Sort and filter days after planting
    df_sorted <- df[order(df$yday), ]
    df_after_planting <- df_sorted[df_sorted$yday >= planting_yday, ]
    
 
    # PHASE 1: Emergence
    
    gdd_emergence_cumsum <- cumsum(df_after_planting$GDD_daily)
    emergence_index <- which(gdd_emergence_cumsum >= GDD_field_emergence_required)[1]
    
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
      bulbing_index <- which(gdd_veg_cumsum >= GDD_vegetative_required)[1]
      
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
      maturation_index <- which(gdd_bulbing_cumsum >= GDD_bulbing_required)[1]
      
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
      harvest_index <- which(gdd_mat_cumsum >= GDD_maturation_required)[1]
      
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
  
  get_seedbed_stress <- function(Prec, Tavg, day_consec_wet) {
    # Trockenstress: zu wenig Regen + heißen Tagen
    if (day_consec_wet < 1 & Prec < 0.5 & Tavg > 25) {
      return(chance_event(
        chance = vv(var_mean = risk_seedbed_high, var_CV = var_CV, n = 1),
        value_if = yield_reduction_seedbed,
        value_if_not = 0, n = 1))
    } else if (day_consec_wet < 1 & Prec < 1) {
      return(chance_event(
        chance = vv(var_mean = risk_seedbed_medium, var_CV = var_CV, n = 1),
        value_if = yield_reduction_seedbed,
        value_if_not = 0, n = 1))
    }
    return(0)
  }
  
#Drought Stress 
  
  
  get_drought_stress <- function(Prec, Tavg, day_consec_wet) {
    if (Prec < 0.5 & Tavg > 25) {
      return(chance_event(
        chance = vv(var_mean = risk_drought_high, var_CV = var_CV, n = 1),
        value_if = yield_reduction_drought,
        value_if_not = 0, n = 1))
    } else if (Prec < 1 & day_consec_wet < 2) {
      return(chance_event(
        chance = vv(var_mean = risk_drought_medium, var_CV = var_CV, n = 1),
        value_if = yield_reduction_drought,
        value_if_not = 0, n = 1))
    }
    return(0)
  }
  

# Extreme Rainfall

  get_extreme_rain_stress <- function(Prec) {
    if (Prec > 40) {
      return(chance_event(
        chance = vv(var_mean = risk_extreme_rain_high, var_CV = var_CV, n = 1),
        value_if = yield_reduction_extreme_rain,
        value_if_not = 0, n = 1))
    } else if (Prec > 20) {
      return(chance_event(
        chance = vv(var_mean = risk_extreme_rain_medium, var_CV = var_CV, n = 1),
        value_if = yield_reduction_extreme_rain,
        value_if_not = 0, n = 1))
    }
    return(0)
  }
  

# Wind stress, not sure if really required 

# Harvest complications 

get_onion_harvest_risk <- function(Prec,
                                   DOY_harvest_window,
                                   prec_threshold_rain = 5,
                                   drought_threshold = 2) {
  
  # Subset Niederschlag auf Erntezeitfenster
  prec_window <- Prec[DOY_harvest_window]
  
  # Anzahl Tage mit starkem Regen (z. B. >5 mm)
  heavy_rain_days <- sum(prec_window > prec_threshold_rain, na.rm = TRUE)
  
  # Anzahl sehr trockener Tage (z. B. <2 mm)
  dry_days <- sum(prec_window < drought_threshold, na.rm = TRUE)
  
  # Risikoskalen definieren (0 = kein Risiko, 3 = hohes Risiko)
  rain_risk <- cut(heavy_rain_days, breaks = c(-1, 0, 2, 4, Inf), labels = 0:3)
  drought_risk <- cut(dry_days, breaks = c(-1, 1, 3, 5, Inf), labels = 0:3)
  
  # Ergebnisliste
  list(
    rain_risk = as.numeric(as.character(rain_risk)),
    drought_risk = as.numeric(as.character(drought_risk)),
    total_harvest_risk_score = sum(as.numeric(as.character(c(rain_risk, drought_risk))))
  )
}

    
### Biotic Stress Factors   

## Weed Pressure
# Low PAR , little bit of rain helps emergence, this needs to be updated any maybe different conditions for different phases 

get_weed_pressure_stress <- function(PAR, Prec, day_consec_wet) {
  if (PAR < 300) {
    if (Prec >= 1 & day_consec_wet >= 2) {
      return(chance_event(
        chance = vv(var_mean = risk_weed_high, var_CV = var_CV, n = 1),
        value_if = yield_reduction_weed,
        value_if_not = 0,
        n = 1
      ))
    } else if (Prec >= 0.5) {
      return(chance_event(
        chance = vv(var_mean = risk_weed_medium, var_CV = var_CV, n = 1),
        value_if = yield_reduction_weed,
        value_if_not = 0,
        n = 1
      ))
    }
  }
  return(0)
}

## Fungal Pathogens 

# Risk for Botrytis

get_botrytis_stress <- function(Tavg, Prec, day_consec_wet) {
  if (Prec >= 1 & day_consec_wet >= 2) {
    return(chance_event(
      chance = vv(var_mean = risk_neck_rot_high, var_CV = var_CV, n = 1),
      value_if = yield_reduction_neck_rot,
      value_if_not = 0,
      n = 1
    ))
    
  } else if (Prec >= 1) {
    return(chance_event(
      chance = vv(var_mean = risk_neck_rot_medium, var_CV = var_CV, n = 1),
      value_if = yield_reduction_neck_rot,
      value_if_not = 0,
      n = 1
    ))
  } else {
    return(0)
  }
}



# Fusarium (Fungi, soil, pathogen, wet & warm increases risk of infection)
get_fusarium_stress <- function(Tavg, Prec, day_consec_wet) {
  if (Tavg >= 15 & Tavg <= 28) {
    if (Prec >= 1 & day_consec_wet >= 2) {
      return(chance_event(
        chance = vv(var_mean = risk_fusarium_high, var_CV = var_CV, n = 1),
        value_if = yield_reduction_fusarium,
        value_if_not = 0,
        n = 1
      ))
    } else if (Prec >= 1) {
      return(chance_event(
        chance = vv(var_mean = risk_fusarium_medium, var_CV = var_CV, n = 1),
        value_if = yield_reduction_fusarium,
        value_if_not = 0,
        n = 1
      ))
    }
  }
  return(0)
}

# Downy Mildew (Peronospora destructor) – cold nights , wet periods 

get_downy_mildew_stress <- function(Tavg, Prec, day_consec_wet) {
  if (Tavg >= 10 & Tavg <= 22) {
    if (Prec >= 1 & day_consec_wet >= 2) {
      return(chance_event(
        chance = vv(var_mean = risk_downy_mildew_high, var_CV = var_CV, n = 1),
        value_if = yield_reduction_downy_mildew,
        value_if_not = 0,
        n = 1
      ))
    } else if (Prec >= 1) {
      return(chance_event(
        chance = vv(var_mean = risk_downy_mildew_medium, var_CV = var_CV, n = 1),
        value_if = yield_reduction_downy_mildew,
        value_if_not = 0,
        n = 1
      ))
    }
  }
  return(0)
}

## Animal Pressure

# Thripse – heiß & trocken
get_thrips_stress <- function(Tavg, Prec) {
  if (Tavg > 20 & Prec < 1) {
    return(chance_event(
      chance = vv(var_mean = risk_thrips_high, var_CV = var_CV, n = 1),
      value_if = yield_reduction_thrips,
      value_if_not = 0,
      n = 1
    ))
  }
  return(0)
}

# Zwiebel­fliege – warm, leicht feucht fördert Larvenaktivität
get_onion_fly_stress <- function(Tavg, Prec) {
  if (Tavg >= 12) {
    if (Prec >= 1) {
      return(chance_event(
        chance = vv(var_mean = risk_onion_fly_high, var_CV = var_CV, n = 1),
        value_if = yield_reduction_onion_fly,
        value_if_not = 0,
        n = 1
      ))
    } else {
      return(chance_event(
        chance = vv(var_mean = risk_onion_fly_medium, var_CV = var_CV, n = 1),
        value_if = yield_reduction_onion_fly,
        value_if_not = 0,
        n = 1
      ))
    }
  }
  return(0)
}

# Glasflügelzikade – trocken & warm bevorzugt
get_leafhopper_stress <- function(Tavg, Prec) {
  if (Tavg >= 15 & Prec < 1) {
    return(chance_event(
      chance = vv(var_mean = risk_leafhopper_medium, var_CV = var_CV, n = 1),
      value_if = yield_reduction_leafhopper,
      value_if_not = 0,
      n = 1
    ))
  }
  return(0)
}

# Drahtwurm – moderate Temperaturen, feuchte Böden (Regen als Proxy)
get_wireworm_stress <- function(Tavg, Prec) {
  if (Tavg >= 10 & Tavg <= 25 & Prec >= 1) {
    return(chance_event(
      chance = vv(var_mean = risk_wireworm_medium, var_CV = var_CV, n = 1),
      value_if = yield_reduction_wireworm,
      value_if_not = 0,
      n = 1
    ))
  }
  return(0)
}

#####  Apply stresss functions

season_risks <- lapply(weather_scenario_list, function(df) {
  
  # Helper to calculate mean values safely
  safe_mean <- function(x, condition) mean(x[condition], na.rm = TRUE)
  safe_max  <- function(x, condition) max(x[condition], na.rm = TRUE)
  
  # --- EMERGENCE PHASE STRESSORS
  emergence_filter <- df$in_emergence_phase == TRUE
  
  emergence_Tavg <- safe_mean(df$Tavg, emergence_filter)
  emergence_Prec <- safe_mean(df$Prec, emergence_filter)
  emergence_PAR  <- safe_mean(df$PAR, emergence_filter)
  emergence_consec_wet <- safe_max(df$day_consec_wet, emergence_filter)
  
  emergence_stress <- 0
  if (any(emergence_filter, na.rm = TRUE)) {
    emergence_stress <- emergence_stress +
      get_weed_pressure_stress(PAR = emergence_PAR, Prec = emergence_Prec, day_consec_wet = emergence_consec_wet) +
      get_fusarium_stress(Prec = emergence_Prec, Tavg = emergence_Tavg, day_consec_wet = emergence_consec_wet) +
      get_onion_fly_stress(Tavg = emergence_Tavg, Prec = emergence_Prec) +
      get_wireworm_stress(Tavg = emergence_Tavg, Prec = emergence_Prec)
  }
  
  # --- VEGETATIVE PHASE STRESSORS
  vegetative_filter <- df$in_vegetative_phase == TRUE
  
  vegetative_Tavg <- safe_mean(df$Tavg, vegetative_filter)
  vegetative_Prec <- safe_mean(df$Prec, vegetative_filter)
  vegetative_PAR  <- safe_mean(df$PAR, vegetative_filter)
  vegetative_consec_wet <- safe_max(df$day_consec_wet, vegetative_filter)
  
  vegetative_stress <- 0
  if (any(vegetative_filter, na.rm = TRUE)) {
    vegetative_stress <- vegetative_stress +
      get_fusarium_stress(Prec = vegetative_Prec, Tavg = vegetative_Tavg, day_consec_wet = vegetative_consec_wet) +
      get_downy_mildew_stress(Prec = vegetative_Prec, Tavg = vegetative_Tavg, day_consec_wet = vegetative_consec_wet) +
      get_thrips_stress(Tavg = vegetative_Tavg, Prec = vegetative_Prec) +
      get_leafhopper_stress(Tavg = vegetative_Tavg, Prec = vegetative_Prec) +
      get_onion_fly_stress(Tavg = vegetative_Tavg, Prec = vegetative_Prec)
  }
  
  # --- BULBING PHASE STRESSORS
  bulbing_filter <- df$in_bulbing_phase == TRUE
  
  bulbing_Tavg <- safe_mean(df$Tavg, bulbing_filter)
  bulbing_Prec <- safe_mean(df$Prec, bulbing_filter)
  bulbing_PAR  <- safe_mean(df$PAR, bulbing_filter)
  bulbing_consec_wet <- safe_max(df$day_consec_wet, bulbing_filter)
  
  bulbing_stress <- 0
  if (any(bulbing_filter, na.rm = TRUE)) {
    bulbing_stress <- bulbing_stress +
      get_botrytis_stress(Prec = bulbing_Prec, Tavg = bulbing_Tavg, day_consec_wet = bulbing_consec_wet) +
      get_fusarium_stress(Prec = bulbing_Prec, Tavg = bulbing_Tavg, day_consec_wet = bulbing_consec_wet) +
      get_downy_mildew_stress(Prec = bulbing_Prec, Tavg = bulbing_Tavg, day_consec_wet = bulbing_consec_wet) +
      get_onion_fly_stress(Tavg = bulbing_Tavg, Prec = bulbing_Prec)
  }
  
  # --- MATURATION PHASE STRESSORS
  maturation_filter <- df$in_maturation_phase == TRUE
  
  maturation_Tavg <- safe_mean(df$Tavg, maturation_filter)
  maturation_Prec <- safe_mean(df$Prec, maturation_filter)
  maturation_PAR  <- safe_mean(df$PAR, maturation_filter)
  maturation_consec_wet <- safe_max(df$day_consec_wet, maturation_filter)
  
  maturation_stress <- 0
  if (any(maturation_filter, na.rm = TRUE)) {
    maturation_stress <- maturation_stress +
      get_botrytis_stress(Prec = maturation_Prec, Tavg = maturation_Tavg, day_consec_wet = maturation_consec_wet) +
      get_fusarium_stress(Prec = maturation_Prec, Tavg = maturation_Tavg, day_consec_wet = maturation_consec_wet) +
      get_downy_mildew_stress(Prec = maturation_Prec, Tavg = maturation_Tavg, day_consec_wet = maturation_consec_wet) +
      get_onion_fly_stress(Tavg = maturation_Tavg, Prec = maturation_Prec)
  }
  
  # --- RETURN: Named list per phase
  return(list(
    emergence_phase_stress = emergence_stress,
    vegetative_phase_stress = vegetative_stress,
    bulbing_phase_stress = bulbing_stress,
    maturation_phase_stress = maturation_stress
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
calculate_biomass_daily <- function(PAR, LAI, Tavg, Prec, stress_factor, LUE_onion) {
  k <- 0.5  # Light extinction coefficient
  
  # Temperature effect
  f_T <- ifelse(Tavg >= 10 & Tavg <= 25, 1,
                ifelse(Tavg < 5 | Tavg > 35, 0, 0.5))
  
  # Water availability
  f_W <- ifelse(Prec >= 2 & Prec <= 10, 1,
                ifelse(Prec < 1, 0.5, 0.7))
  
  # Stress (same scalar applied to all days in a phase)
  f_S <- 1 - stress_factor
  
  # Vectorized light interception & biomass increment
  delta_B <- LUE_onion * PAR * (1 - exp(-k * LAI)) * f_T * f_W * f_S
  return(delta_B)
}

# Loop over each weather scenario
# Finaler Code für Biomasseberechnung über alle Szenarien hinweg

biomass_all_scenarios <- lapply(weather_scenario_list, function(df) {
  
  calculate_phase_biomass <- function(df_phase, LAI_value, stress_value) {
    if (nrow(df_phase) == 0) return(0)
    
    delta_B <- calculate_biomass_daily(
      PAR = df_phase$PAR,
      LAI = LAI_value,
      Tavg = df_phase$Tavg,
      Prec = df_phase$Prec,
      stress_factor = stress_value,
      LUE_onion = LUE_onion
    )
    
    return(sum(delta_B, na.rm = TRUE))
  }
  
  # Datensätze für jede Phase extrahieren
  df_emergence   <- df[df$in_emergence_phase   == TRUE, ]
  df_veg         <- df[df$in_vegetative_phase  == TRUE, ]
  df_bulbing     <- df[df$in_bulbing_phase     == TRUE, ]
  df_maturation  <- df[df$in_maturation_phase  == TRUE, ]
  
  # Stresswerte aus dem DataFrame extrahieren (einheitlich über die Phase)
  stress_emergence   <- unique(na.omit(df$emergence_phase_stress))
  stress_veg         <- unique(na.omit(df$vegetative_phase_stress))
  stress_bulbing     <- unique(na.omit(df$bulbing_phase_stress))
  stress_maturation  <- unique(na.omit(df$maturation_phase_stress))
  
  # Nur einen Wert pro Stress erlauben (sonst Warnung)
  stress_emergence   <- if (length(stress_emergence) == 1) stress_emergence else mean(stress_emergence, na.rm = TRUE)
  stress_veg         <- if (length(stress_veg) == 1) stress_veg else mean(stress_veg, na.rm = TRUE)
  stress_bulbing     <- if (length(stress_bulbing) == 1) stress_bulbing else mean(stress_bulbing, na.rm = TRUE)
  stress_maturation  <- if (length(stress_maturation) == 1) stress_maturation else mean(stress_maturation, na.rm = TRUE)
  
  # Biomasse pro Phase berechnen
  biomass_emergence   <- calculate_phase_biomass(df_emergence,   LAI_emergence,   stress_emergence)
  biomass_veg         <- calculate_phase_biomass(df_veg,         LAI_veg,         stress_veg)
  biomass_bulbing     <- calculate_phase_biomass(df_bulbing,     LAI_bulbing,     stress_bulbing)
  biomass_maturation  <- calculate_phase_biomass(df_maturation,  LAI_maturation,  stress_maturation)
  
  # Zusammenfassung
  list(
    emergence   = biomass_emergence,
    vegetative  = biomass_veg,
    bulbing     = biomass_bulbing,
    maturation  = biomass_maturation,
    total       = biomass_emergence + biomass_veg + biomass_bulbing + biomass_maturation
  )
})


print(unique(LAI_value))
