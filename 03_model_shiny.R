
################################################################################
# Onion Climate Impact Model
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

enableJIT(3)

################################################################################
# Pre-compute phase indices
################################################################################
precompute_phase_data <- function(df, sow) {
  last_y <- max(df$yday, na.rm = TRUE)
  list(
    em = list(idx = df$yday >  sow & df$yday <= (sow + 15),
              range = (sow + 1):min(sow + 15, last_y)),
    vg = list(idx = df$yday > (sow + 15) & df$yday <= (sow + 50),
              range = (sow + 16):min(sow + 50, last_y)),
    bl = list(idx = df$yday > (sow + 50) & df$yday <= (sow + 90),
              range = (sow + 51):min(sow + 90, last_y)),
    mt = list(idx = df$yday > (sow + 90),
              range = (sow + 91):last_y)
  )
}

################################################################################
# Batch risk computation (means & maxima computed once per phase)
################################################################################
compute_all_risks <- compiler::cmpfun(function(phase_data, df, params) {
  # Helper to compute stats per phase
  stats_for <- function(idx) list(
    T   = mean(df$Tavg[idx],    na.rm = TRUE),
    P   = mean(df$Prec[idx],    na.rm = TRUE),
    RH  = mean(df$RH_mean[idx], na.rm = TRUE),
    Ts  = mean(df$Ts_5cm[idx],  na.rm = TRUE),
    wet = max(df$day_consec_wet[idx], na.rm = TRUE),
    dry = max(df$day_consec_dry[idx], na.rm = TRUE)
  )
  
  ems <- stats_for(phase_data$em$idx)
  vgs <- stats_for(phase_data$vg$idx)
  bls <- stats_for(phase_data$bl$idx)
  mts <- stats_for(phase_data$mt$idx)
  
  # Vectors needed for some stresses
  em_Prec <- df$Prec[phase_data$em$idx]
  em_Tavg <- df$Tavg[phase_data$em$idx]
  vg_Prec <- df$Prec[phase_data$vg$idx]
  vg_Tavg <- df$Tavg[phase_data$vg$idx]
  bl_Prec <- df$Prec[phase_data$bl$idx]
  bl_Tavg <- df$Tavg[phase_data$bl$idx]
  mt_Prec <- df$Prec[phase_data$mt$idx]
  mt_Tavg <- df$Tavg[phase_data$mt$idx]
  
  # Emergence risks
  risks_em <- c(
    seedbed = get_seedbed_stress(ems$P, ems$Ts, ems$RH, ems$wet,
                                 params$prec_seedbed_medium_p, params$prec_seedbed_high_p,
                                 params$Ts_min_seedbed_p, params$RH_high_seedbed_p,
                                 params$risk_seedbed_medium_t, params$risk_seedbed_high_t,
                                 params$impact_Ts_seedbed_t, params$impact_RH_seedbed_t),
    drought = get_drought_stress(ems$P, ems$T, ems$RH, ems$dry,
                                 params$prec_drought_threshold_p, params$Tavg_drought_threshold_p,
                                 params$RH_drought_threshold_p, params$impact_prec_drought_t,
                                 params$impact_days_dry_drought_t, params$impact_temp_drought_t,
                                 params$impact_RH_drought_t),
    exrain  = get_extreme_rain_stress(em_Prec, params$prec_extreme_rain_medium_p,
                                      params$prec_extreme_rain_high_p, params$impact_days_extreme_rain_t),
    hail    = get_hail_stress(em_Tavg, em_Prec, params$prec_hail_threshold_p,
                              params$Tavg_hail_threshold_p, params$impact_days_hail_t),
    fusarium= get_fusarium_stress(ems$T, ems$P, ems$wet,
                                  params$prec_fusarium_threshold_p, params$impact_prec_fusarium_t,
                                  params$impact_days_wet_fusarium_t, params$Topt_fusarium_p,
                                  params$Twidth_fusarium_p),
    onfly   = get_onion_fly_stress(ems$T, ems$P, ems$wet,
                                   params$prec_onion_fly_threshold_p, params$Topt_onion_fly_p,
                                   params$Twidth_onion_fly_p, params$impact_prec_onion_fly_t,
                                   params$impact_days_wet_onion_fly_t),
    wire    = get_wireworm_stress(ems$T, ems$Ts, ems$P, ems$wet,
                                  params$prec_wireworm_threshold_p, params$Topt_wireworm_p,
                                  params$Twidth_wireworm_p, params$Ts_opt_wireworm_p,
                                  params$Ts_width_wireworm_p, params$impact_prec_wireworm_t,
                                  params$impact_days_wet_wireworm_t)
  )
  
  # Vegetative risks
  risks_vg <- c(
    drought = get_drought_stress(vgs$P, vgs$T, vgs$RH, vgs$dry,
                                 params$prec_drought_threshold_p, params$Tavg_drought_threshold_p,
                                 params$RH_drought_threshold_p, params$impact_prec_drought_t,
                                 params$impact_days_dry_drought_t, params$impact_temp_drought_t,
                                 params$impact_RH_drought_t),
    exrain  = get_extreme_rain_stress(vg_Prec, params$prec_extreme_rain_medium_p,
                                      params$prec_extreme_rain_high_p, params$impact_days_extreme_rain_t),
    hail    = get_hail_stress(vg_Tavg, vg_Prec, params$prec_hail_threshold_p,
                              params$Tavg_hail_threshold_p, params$impact_days_hail_t),
    fusarium= get_fusarium_stress(vgs$T, vgs$P, vgs$wet,
                                  params$prec_fusarium_threshold_p, params$impact_prec_fusarium_t,
                                  params$impact_days_wet_fusarium_t, params$Topt_fusarium_p,
                                  params$Twidth_fusarium_p),
    mildew  = get_downy_mildew_stress(vgs$T, vgs$RH, vgs$wet,
                                      params$Topt_mildew_p, params$Twidth_mildew_p,
                                      params$RH_mildew_threshold_p, params$wet_days_needed_mildew_p,
                                      params$impact_RH_mildew_t, params$impact_days_wet_mildew_t),
    thrips  = get_thrips_stress(vgs$T, vgs$P, vgs$dry,
                                params$prec_thrips_threshold_p, params$Topt_thrips_p,
                                params$Twidth_thrips_p, params$impact_prec_thrips_t,
                                params$impact_days_dry_thrips_t),
    onfly   = get_onion_fly_stress(vgs$T, vgs$P, vgs$wet,
                                   params$prec_onion_fly_threshold_p, params$Topt_onion_fly_p,
                                   params$Twidth_onion_fly_p, params$impact_prec_onion_fly_t,
                                   params$impact_days_wet_onion_fly_t)
  )
  
  # Bulbing risks
  risks_bl <- c(
    drought = get_drought_stress(bls$P, bls$T, bls$RH, bls$dry,
                                 params$prec_drought_threshold_p, params$Tavg_drought_threshold_p,
                                 params$RH_drought_threshold_p, params$impact_prec_drought_t,
                                 params$impact_days_dry_drought_t, params$impact_temp_drought_t,
                                 params$impact_RH_drought_t),
    exrain  = get_extreme_rain_stress(bl_Prec, params$prec_extreme_rain_medium_p,
                                      params$prec_extreme_rain_high_p, params$impact_days_extreme_rain_t),
    hail    = get_hail_stress(bl_Tavg, bl_Prec, params$prec_hail_threshold_p,
                              params$Tavg_hail_threshold_p, params$impact_days_hail_t),
    botrytis= get_botrytis_stress(bls$T, bls$RH, bls$P, bls$wet,
                                  params$Topt_botrytis_p, params$Twidth_botrytis_p,
                                  params$RH_botrytis_threshold_p, params$prec_botrytis_threshold_p,
                                  params$wet_days_needed_botrytis_p, params$impact_RH_botrytis_t,
                                  params$impact_prec_botrytis_t, params$impact_days_wet_botrytis_t),
    fusarium= get_fusarium_stress(bls$T, bls$P, bls$wet,
                                  params$prec_fusarium_threshold_p, params$impact_prec_fusarium_t,
                                  params$impact_days_wet_fusarium_t, params$Topt_fusarium_p,
                                  params$Twidth_fusarium_p),
    mildew  = get_downy_mildew_stress(bls$T, bls$RH, bls$wet,
                                      params$Topt_mildew_p, params$Twidth_mildew_p,
                                      params$RH_mildew_threshold_p, params$wet_days_needed_mildew_p,
                                      params$impact_RH_mildew_t, params$impact_days_wet_mildew_t),
    onfly   = get_onion_fly_stress(bls$T, bls$P, bls$wet,
                                   params$prec_onion_fly_threshold_p, params$Topt_onion_fly_p,
                                   params$Twidth_onion_fly_p, params$impact_prec_onion_fly_t,
                                   params$impact_days_wet_onion_fly_t)
  )
  
  # Maturation risks
  risks_mt <- c(
    drought = get_drought_stress(mts$P, mts$T, mts$RH, mts$dry,
                                 params$prec_drought_threshold_p, params$Tavg_drought_threshold_p,
                                 params$RH_drought_threshold_p, params$impact_prec_drought_t,
                                 params$impact_days_dry_drought_t, params$impact_temp_drought_t,
                                 params$impact_RH_drought_t),
    exrain  = get_extreme_rain_stress(mt_Prec, params$prec_extreme_rain_medium_p,
                                      params$prec_extreme_rain_high_p, params$impact_days_extreme_rain_t),
    hail    = get_hail_stress(mt_Tavg, mt_Prec, params$prec_hail_threshold_p,
                              params$Tavg_hail_threshold_p, params$impact_days_hail_t),
    botrytis= get_botrytis_stress(mts$T, mts$RH, mts$P, mts$wet,
                                  params$Topt_botrytis_p, params$Twidth_botrytis_p,
                                  params$RH_botrytis_threshold_p, params$prec_botrytis_threshold_p,
                                  params$wet_days_needed_botrytis_p, params$impact_RH_botrytis_t,
                                  params$impact_prec_botrytis_t, params$impact_days_wet_botrytis_t),
    fusarium= get_fusarium_stress(mts$T, mts$P, mts$wet,
                                  params$prec_fusarium_threshold_p, params$impact_prec_fusarium_t,
                                  params$impact_days_wet_fusarium_t, params$Topt_fusarium_p,
                                  params$Twidth_fusarium_p),
    mildew  = get_downy_mildew_stress(mts$T, mts$RH, mts$wet,
                                      params$Topt_mildew_p, params$Twidth_mildew_p,
                                      params$RH_mildew_threshold_p, params$wet_days_needed_mildew_p,
                                      params$impact_RH_mildew_t, params$impact_days_wet_mildew_t),
    onfly   = get_onion_fly_stress(mts$T, mts$P, mts$wet,
                                   params$prec_onion_fly_threshold_p, params$Topt_onion_fly_p,
                                   params$Twidth_onion_fly_p, params$impact_prec_onion_fly_t,
                                   params$impact_days_wet_onion_fly_t)
  )
  
  list(em = risks_em, vg = risks_vg, bl = risks_bl, mt = risks_mt)
})

helper_function(attach_to_global = TRUE)


################################################################################
# Main model 
################################################################################
onion_climate_impact <- compiler::cmpfun(function() {
  
  # Sample exactly one season data.table per scenario
  weather_scenario_list <- lapply(weather_precomputed, function(list_per_scenario) {
    sample(list_per_scenario, 1)[[1]]
  })
  
  # Collect params 
  params <- list(
    # Sowing
    T_soil_min_sowing_p = T_soil_min_sowing_p,
    warm_days_needed_sowing_p = warm_days_needed_sowing_p,
    frost_buffer_days_sowing_p = frost_buffer_days_sowing_p,
    min_yday_sowing_c = min_yday_sowing_c,
    max_yday_sowing_c = max_yday_sowing_c,
    rain_window_days_sowing_p = rain_window_days_sowing_p,
    rain_sum_max_mm_sowing_p = rain_sum_max_mm_sowing_p,
    min_dry_days_sowing_c = min_dry_days_sowing_c,
    
    # Seedbed
    prec_seedbed_medium_p = prec_seedbed_medium_p,
    prec_seedbed_high_p = prec_seedbed_high_p,
    Ts_min_seedbed_p = Ts_min_seedbed_p,
    RH_high_seedbed_p = RH_high_seedbed_p,
    risk_seedbed_medium_t = risk_seedbed_medium_t,
    risk_seedbed_high_t = risk_seedbed_high_t,
    impact_Ts_seedbed_t = impact_Ts_seedbed_t,
    impact_RH_seedbed_t = impact_RH_seedbed_t,
    
    # Drought
    prec_drought_threshold_p = prec_drought_threshold_p,
    Tavg_drought_threshold_p = Tavg_drought_threshold_p,
    RH_drought_threshold_p = RH_drought_threshold_p,
    impact_prec_drought_t = impact_prec_drought_t,
    impact_days_dry_drought_t = impact_days_dry_drought_t,
    impact_temp_drought_t = impact_temp_drought_t,
    impact_RH_drought_t = impact_RH_drought_t,
    
    # Extreme rain
    prec_extreme_rain_medium_p = prec_extreme_rain_medium_p,
    prec_extreme_rain_high_p = prec_extreme_rain_high_p,
    impact_days_extreme_rain_t = impact_days_extreme_rain_t,
    
    # Hail
    prec_hail_threshold_p = prec_hail_threshold_p,
    Tavg_hail_threshold_p = Tavg_hail_threshold_p,
    impact_days_hail_t = impact_days_hail_t,
    
    # Botrytis
    Topt_botrytis_p = Topt_botrytis_p,
    Twidth_botrytis_p = Twidth_botrytis_p,
    RH_botrytis_threshold_p = RH_botrytis_threshold_p,
    prec_botrytis_threshold_p = prec_botrytis_threshold_p,
    wet_days_needed_botrytis_p = wet_days_needed_botrytis_p,
    impact_RH_botrytis_t = impact_RH_botrytis_t,
    impact_prec_botrytis_t = impact_prec_botrytis_t,
    impact_days_wet_botrytis_t = impact_days_wet_botrytis_t,
    
    # Mildew
    Topt_mildew_p = Topt_mildew_p,
    Twidth_mildew_p = Twidth_mildew_p,
    RH_mildew_threshold_p = RH_mildew_threshold_p,
    wet_days_needed_mildew_p = wet_days_needed_mildew_p,
    impact_RH_mildew_t = impact_RH_mildew_t,
    impact_days_wet_mildew_t = impact_days_wet_mildew_t,
    
    # Thrips
    prec_thrips_threshold_p = prec_thrips_threshold_p,
    Topt_thrips_p = Topt_thrips_p,
    Twidth_thrips_p = Twidth_thrips_p,
    impact_prec_thrips_t = impact_prec_thrips_t,
    impact_days_dry_thrips_t = impact_days_dry_thrips_t,
    
    # Onion fly
    prec_onion_fly_threshold_p = prec_onion_fly_threshold_p,
    Topt_onion_fly_p = Topt_onion_fly_p,
    Twidth_onion_fly_p = Twidth_onion_fly_p,
    impact_prec_onion_fly_t = impact_prec_onion_fly_t,
    impact_days_wet_onion_fly_t = impact_days_wet_onion_fly_t,
    
    # Wireworm
    prec_wireworm_threshold_p = prec_wireworm_threshold_p,
    Topt_wireworm_p = Topt_wireworm_p,
    Twidth_wireworm_p = Twidth_wireworm_p,
    Ts_opt_wireworm_p = Ts_opt_wireworm_p,
    Ts_width_wireworm_p = Ts_width_wireworm_p,
    impact_prec_wireworm_t = impact_prec_wireworm_t,
    impact_days_wet_wireworm_t = impact_days_wet_wireworm_t,
    
    # Fusarium
    prec_fusarium_threshold_p = prec_fusarium_threshold_p,
    impact_prec_fusarium_t = impact_prec_fusarium_t,
    impact_days_wet_fusarium_t = impact_days_wet_fusarium_t,
    Topt_fusarium_p = Topt_fusarium_p,
    Twidth_fusarium_p = Twidth_fusarium_p,
    
    # Yield reductions
    yield_reduction_seedbed_t = yield_reduction_seedbed_t,
    yield_reduction_drought_t = yield_reduction_drought_t,
    yield_reduction_extreme_rain_t = yield_reduction_extreme_rain_t,
    yield_reduction_hail_t = yield_reduction_hail_t,
    yield_reduction_fusarium_t = yield_reduction_fusarium_t,
    yield_reduction_onion_fly_t = yield_reduction_onion_fly_t,
    yield_reduction_wireworm_t = yield_reduction_wireworm_t,
    yield_reduction_botrytis_t = yield_reduction_botrytis_t,
    yield_reduction_downy_mildew_t = yield_reduction_downy_mildew_t,
    yield_reduction_thrips_t = yield_reduction_thrips_t,
    
    # Biomass
    LAI_emergence_p = LAI_emergence_p,
    LAI_veg_p = LAI_veg_p,
    LAI_bulbing_p = LAI_bulbing_p,
    LAI_maturation_p = LAI_maturation_p,
    f_T_1_lower_p = f_T_1_lower_p,
    f_T_1_upper_p = f_T_1_upper_p,
    f_T_0_lower_p = f_T_0_lower_p,
    f_T_0_upper_p = f_T_0_upper_p,
    f_W_1_lower_p = f_W_1_lower_p,
    f_W_1_upper_p = f_W_1_upper_p,
    f_W_0.5_p = f_W_0.5_p,
    LUE_onion_p = LUE_onion_p,
    lec_k_c = lec_k_c,
    onions_per_ha_p = onions_per_ha_p,
    dry_onion_weight_p = dry_onion_weight_p,
    HI_onions_p = HI_onions_p
  )
  
  results <- vector("list", length(weather_scenario_list))
  names(results) <- names(weather_scenario_list)
  
  for (sc in names(weather_scenario_list)) {
    df <- weather_scenario_list[[sc]]
    setDT(df)  # ensure data.table ops
    
    # Consecutive counters (fast)
    df$day_consec_wet <- consec_counter(!is.na(df$Prec) & df$Prec > 0)
    df$day_consec_dry <- consec_counter(is.na(df$Prec) | df$Prec == 0)
    
    
    # Sowing day (vectorized)
    sow <- compute_onion_sowing_yday_fast(
      df = df,
      T_soil_min_sowing_p = params$T_soil_min_sowing_p,
      warm_days_needed_sowing_p = params$warm_days_needed_sowing_p,
      frost_buffer_days_sowing_p = params$frost_buffer_days_sowing_p,
      min_yday_sowing_c = params$min_yday_sowing_c,
      max_yday_sowing_c = params$max_yday_sowing_c,
      rain_window_days_sowing_p = params$rain_window_days_sowing_p,
      rain_sum_max_mm_sowing_p = params$rain_sum_max_mm_sowing_p,
      min_dry_days_sowing_c = params$min_dry_days_sowing_c
    )
    
    # Phase masks
    phase_data <- precompute_phase_data(df, sow)
    
    # Risks in one call
    all_risks <- compute_all_risks(phase_data, df, params)
    
    # Yield reduction multipliers (single draws)
    multipliers <- c(
      fast_chance_event(all_risks$em["seedbed"],  params$yield_reduction_seedbed_t),
      fast_chance_event(all_risks$em["drought"],  params$yield_reduction_drought_t),
      fast_chance_event(all_risks$em["exrain"],   params$yield_reduction_extreme_rain_t),
      fast_chance_event(all_risks$em["hail"],     params$yield_reduction_hail_t),
      fast_chance_event(all_risks$em["fusarium"], params$yield_reduction_fusarium_t),
      fast_chance_event(all_risks$em["onfly"],    params$yield_reduction_onion_fly_t),
      fast_chance_event(all_risks$em["wire"],     params$yield_reduction_wireworm_t),
      
      fast_chance_event(all_risks$vg["drought"],  params$yield_reduction_drought_t),
      fast_chance_event(all_risks$vg["exrain"],   params$yield_reduction_extreme_rain_t),
      fast_chance_event(all_risks$vg["hail"],     params$yield_reduction_hail_t),
      fast_chance_event(all_risks$vg["fusarium"], params$yield_reduction_fusarium_t),
      fast_chance_event(all_risks$vg["mildew"],   params$yield_reduction_downy_mildew_t),
      fast_chance_event(all_risks$vg["thrips"],   params$yield_reduction_thrips_t),
      fast_chance_event(all_risks$vg["onfly"],    params$yield_reduction_onion_fly_t),
      
      fast_chance_event(all_risks$bl["drought"],  params$yield_reduction_drought_t),
      fast_chance_event(all_risks$bl["exrain"],   params$yield_reduction_extreme_rain_t),
      fast_chance_event(all_risks$bl["hail"],     params$yield_reduction_hail_t),
      fast_chance_event(all_risks$bl["botrytis"], params$yield_reduction_botrytis_t),
      fast_chance_event(all_risks$bl["fusarium"], params$yield_reduction_fusarium_t),
      fast_chance_event(all_risks$bl["mildew"],   params$yield_reduction_downy_mildew_t),
      fast_chance_event(all_risks$bl["onfly"],    params$yield_reduction_onion_fly_t),
      
      fast_chance_event(all_risks$mt["drought"],  params$yield_reduction_drought_t),
      fast_chance_event(all_risks$mt["exrain"],   params$yield_reduction_extreme_rain_t),
      fast_chance_event(all_risks$mt["hail"],     params$yield_reduction_hail_t),
      fast_chance_event(all_risks$mt["botrytis"], params$yield_reduction_botrytis_t),
      fast_chance_event(all_risks$mt["fusarium"], params$yield_reduction_fusarium_t),
      fast_chance_event(all_risks$mt["mildew"],   params$yield_reduction_downy_mildew_t),
      fast_chance_event(all_risks$mt["onfly"],    params$yield_reduction_onion_fly_t)
    )
    
    combined_yield_multiplier <- prod(multipliers, na.rm = TRUE)
    
    # Biomass (vectorized per phase)
    biomass_em <- sum(calc_bio_vectorized(
      df$PAR[phase_data$em$idx], params$LAI_emergence_p,
      df$Tavg[phase_data$em$idx], df$Prec[phase_data$em$idx],
      params$f_T_1_lower_p, params$f_T_1_upper_p, params$f_T_0_lower_p, params$f_T_0_upper_p,
      params$f_W_1_lower_p, params$f_W_1_upper_p, params$f_W_0.5_p,
      params$LUE_onion_p, params$lec_k_c), na.rm = TRUE)
    
    biomass_vg <- sum(calc_bio_vectorized(
      df$PAR[phase_data$vg$idx], params$LAI_veg_p,
      df$Tavg[phase_data$vg$idx], df$Prec[phase_data$vg$idx],
      params$f_T_1_lower_p, params$f_T_1_upper_p, params$f_T_0_lower_p, params$f_T_0_upper_p,
      params$f_W_1_lower_p, params$f_W_1_upper_p, params$f_W_0.5_p,
      params$LUE_onion_p, params$lec_k_c), na.rm = TRUE)
    
    biomass_bl <- sum(calc_bio_vectorized(
      df$PAR[phase_data$bl$idx], params$LAI_bulbing_p,
      df$Tavg[phase_data$bl$idx], df$Prec[phase_data$bl$idx],
      params$f_T_1_lower_p, params$f_T_1_upper_p, params$f_T_0_lower_p, params$f_T_0_upper_p,
      params$f_W_1_lower_p, params$f_W_1_upper_p, params$f_W_0.5_p,
      params$LUE_onion_p, params$lec_k_c), na.rm = TRUE)
    
    biomass_mt <- sum(calc_bio_vectorized(
      df$PAR[phase_data$mt$idx], params$LAI_maturation_p,
      df$Tavg[phase_data$mt$idx], df$Prec[phase_data$mt$idx],
      params$f_T_1_lower_p, params$f_T_1_upper_p, params$f_T_0_lower_p, params$f_T_0_upper_p,
      params$f_W_1_lower_p, params$f_W_1_upper_p, params$f_W_0.5_p,
      params$LUE_onion_p, params$lec_k_c), na.rm = TRUE)
    
    total_biomass <- biomass_em + biomass_vg + biomass_bl + biomass_mt
    
    raw_yield <- (total_biomass * params$onions_per_ha_p *
                    (1 - params$dry_onion_weight_p)) * params$HI_onions_p / 1e6
    
    final_yield <- raw_yield * combined_yield_multiplier
    
    results[[sc]] <- list(
      sowing_yday = sow,
      harvest_yday = max(df$yday, na.rm = TRUE),
      total_biomass = total_biomass,
      raw_yield_per_ha = raw_yield,
      final_yield_per_ha = final_yield,
      combined_yield_multiplier = combined_yield_multiplier,
      
      # expose multipliers for diagnostics
      m_seedbed_em = multipliers[1],
      m_drought_em = multipliers[2],
      m_exrain_em  = multipliers[3],
      m_hail_em    = multipliers[4],
      m_fus_em     = multipliers[5],
      m_onfly_em   = multipliers[6],
      m_wire_em    = multipliers[7],
      
      m_drought_vg = multipliers[8],
      m_exrain_vg  = multipliers[9],
      m_hail_vg    = multipliers[10],
      m_fus_vg     = multipliers[11],
      m_mildew_vg  = multipliers[12],
      m_thrips_vg  = multipliers[13],
      m_onfly_vg   = multipliers[14],
      
      m_drought_bl = multipliers[15],
      m_exrain_bl  = multipliers[16],
      m_hail_bl    = multipliers[17],
      m_botrytis_bl= multipliers[18],
      m_fus_bl     = multipliers[19],
      m_mildew_bl  = multipliers[20],
      m_onfly_bl   = multipliers[21],
      
      m_drought_mt = multipliers[22],
      m_exrain_mt  = multipliers[23],
      m_hail_mt    = multipliers[24],
      m_botrytis_mt= multipliers[25],
      m_fus_mt     = multipliers[26],
      m_mildew_mt  = multipliers[27],
      m_onfly_mt   = multipliers[28]
    )
  }
  
  results
})

################################################################################
# Monte Carlo Simulation
################################################################################

input_variables <- read.csv("input_table_onion_final.csv", header = TRUE, sep = ";")

runtime <- system.time({
  onion_mc_simulation <- mcSimulation(
    estimate = as.estimate(input_variables),
    model_function = onion_climate_impact,
    numberOfModelRuns = 1000,
    functionSyntax = "plainNames"
  )
})

cat("\n⏱️ Optimized simulation finished in",
    round(runtime["elapsed"], 2), "seconds (", 
    round(runtime["elapsed"]/60, 2), "minutes)\n")

