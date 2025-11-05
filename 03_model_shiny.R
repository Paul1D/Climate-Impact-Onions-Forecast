################################################################################
# Onion Climate Impact Model (fixed + consistent version)
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
# Pre-compute phase indices based on cumulative GDD
################################################################################
precompute_phase_data <- compiler::cmpfun(function(df, sow,
                                                   GDD_field_emergence_required_p,
                                                   GDD_vegetative_required_p,
                                                   GDD_bulbing_required_p,
                                                   GDD_maturation_required_p) {
  n <- nrow(df)
  if (n == 0L) return(NULL)
  
  df[, GDD_cum := 0]
  after_sow <- df$yday >= sow
  df$GDD_cum[after_sow] <- cumsum(df$GDD_daily[after_sow])
  
  gdd_em <- GDD_field_emergence_required_p
  gdd_vg <- GDD_vegetative_required_p
  gdd_bl <- GDD_bulbing_required_p
  gdd_mt <- GDD_maturation_required_p
  
  em_idx <- after_sow & df$GDD_cum <= gdd_em
  vg_idx <- after_sow & df$GDD_cum > gdd_em & df$GDD_cum <= gdd_vg
  bl_idx <- after_sow & df$GDD_cum > gdd_vg & df$GDD_cum <= gdd_bl
  mt_idx <- after_sow & df$GDD_cum > gdd_bl & df$GDD_cum <= gdd_mt
  
  list(
    em = list(idx = em_idx, range = which(em_idx)),
    vg = list(idx = vg_idx, range = which(vg_idx)),
    bl = list(idx = bl_idx, range = which(bl_idx)),
    mt = list(idx = mt_idx, range = which(mt_idx))
  )
})

################################################################################
# Batch risk computation (means & maxima computed once per phase)
################################################################################
compute_all_risks <- compiler::cmpfun(function(phase_data, df, params) {
  
  stats_for <- function(idx) list(
    T   = mean(df$Tavg[idx],    na.rm = TRUE),
    P   = mean(df$Prec[idx],    na.rm = TRUE),
    RH  = mean(df$RH_mean[idx], na.rm = TRUE),
    Ts  = mean(df$Ts_5cm[idx],  na.rm = TRUE),
    wet = if (any(idx, na.rm = TRUE)) max(df$day_consec_wet[idx], na.rm = TRUE) else 0,
    dry = if (any(idx, na.rm = TRUE)) max(df$day_consec_dry[idx], na.rm = TRUE) else 0
  )
  
  ems <- stats_for(phase_data$em$idx)
  vgs <- stats_for(phase_data$vg$idx)
  bls <- stats_for(phase_data$bl$idx)
  mts <- stats_for(phase_data$mt$idx)

  
  # === EMERGENCE ===
  risks_em <- c(
    drought = get_drought_stress(ems$T, ems$RH, ems$dry,
                                 params$rh_drought_threshold_p, params$Tavg_drought_threshold_p,
                                 params$impact_rh_drought_t, params$impact_days_dry_drought_t, params$impact_temp_drought_t),
    exrain = get_extreme_rain_stress(df$Prec[phase_data$em$idx],
                                     params$prec_extreme_rain_medium_p, params$prec_extreme_rain_high_p,
                                     params$impact_days_extreme_rain_t),
    hail = get_hail_stress(df$Tavg[phase_data$em$idx], df$Prec[phase_data$em$idx],
                           params$prec_hail_threshold_p, params$Tavg_hail_threshold_p, params$impact_days_hail_t),
    fusarium = get_fusarium_stress(ems$Ts, ems$RH, ems$wet,
                                   params$rh_fusarium_threshold_p, params$impact_rh_fusarium_t,
                                   params$impact_days_wet_fusarium_t, params$Tsopt_fusarium_p, params$Tswidth_fusarium_p)
  )
  
  # === VEGETATIVE ===
  risks_vg <- c(
    drought = get_drought_stress(vgs$T, vgs$RH, vgs$dry,
                                 params$rh_drought_threshold_p, params$Tavg_drought_threshold_p,
                                 params$impact_rh_drought_t, params$impact_days_dry_drought_t, params$impact_temp_drought_t),
    exrain = get_extreme_rain_stress(df$Prec[phase_data$vg$idx],
                                     params$prec_extreme_rain_medium_p, params$prec_extreme_rain_high_p,
                                     params$impact_days_extreme_rain_t),
    hail = get_hail_stress(df$Tavg[phase_data$vg$idx], df$Prec[phase_data$vg$idx],
                           params$prec_hail_threshold_p, params$Tavg_hail_threshold_p, params$impact_days_hail_t),
    mildew = get_downy_mildew_stress(vgs$T, vgs$RH, vgs$wet,
                                     params$rh_mildew_threshold_p, params$impact_rh_mildew_t,
                                     params$impact_days_wet_mildew_t, params$Topt_mildew_p, params$Twidth_mildew_p),
    thrips = get_thrips_stress(vgs$T, vgs$RH, vgs$dry,
                               params$rh_thrips_threshold_p, params$Topt_thrips_p, params$Twidth_thrips_p,
                               params$impact_rh_thrips_t, params$impact_days_dry_thrips_t),
    fusarium = get_fusarium_stress(vgs$Ts, vgs$RH, vgs$wet,
                                   params$rh_fusarium_threshold_p, params$impact_rh_fusarium_t,
                                   params$impact_days_wet_fusarium_t, params$Tsopt_fusarium_p, params$Tswidth_fusarium_p)
  )
  
  # === BULBING ===
  risks_bl <- c(
    drought = get_drought_stress(bls$T, bls$RH, bls$dry,
                                 params$rh_drought_threshold_p, params$Tavg_drought_threshold_p,
                                 params$impact_rh_drought_t, params$impact_days_dry_drought_t, params$impact_temp_drought_t),
    botrytis = get_botrytis_stress(bls$T, bls$RH, bls$wet,
                                   params$rh_botrytis_threshold_p, params$impact_rh_botrytis_t,
                                   params$impact_days_wet_botrytis_t, params$Topt_botrytis_p, params$Twidth_botrytis_p),
    mildew = get_downy_mildew_stress(bls$T, bls$RH, bls$wet,
                                     params$rh_mildew_threshold_p, params$impact_rh_mildew_t,
                                     params$impact_days_wet_mildew_t, params$Topt_mildew_p, params$Twidth_mildew_p),
    fusarium = get_fusarium_stress(bls$Ts, bls$RH, bls$wet,
                                   params$rh_fusarium_threshold_p, params$impact_rh_fusarium_t,
                                   params$impact_days_wet_fusarium_t, params$Tsopt_fusarium_p, params$Tswidth_fusarium_p)
  )
  
  # === MATURATION ===
 risks_mt <- c(
    drought = get_drought_stress(mts$T, mts$RH, mts$dry,
                                 params$rh_drought_threshold_p, params$Tavg_drought_threshold_p,
                                 params$impact_rh_drought_t, params$impact_days_dry_drought_t, params$impact_temp_drought_t),
    exrain = get_extreme_rain_stress(df$Prec[phase_data$mt$idx],
                                     params$prec_extreme_rain_medium_p, params$prec_extreme_rain_high_p,
                                     params$impact_days_extreme_rain_t),
    hail = get_hail_stress(df$Tavg[phase_data$mt$idx], df$Prec[phase_data$mt$idx],
                           params$prec_hail_threshold_p, params$Tavg_hail_threshold_p, params$impact_days_hail_t),
    botrytis = get_botrytis_stress(mts$T, mts$RH, mts$wet,
                                   params$rh_botrytis_threshold_p, params$impact_rh_botrytis_t,
                                   params$impact_days_wet_botrytis_t, params$Topt_botrytis_p, params$Twidth_botrytis_p),
    fusarium = get_fusarium_stress(mts$Ts, mts$RH, mts$wet,
                                   params$rh_fusarium_threshold_p, params$impact_rh_fusarium_t,
                                   params$impact_days_wet_fusarium_t, params$Tsopt_fusarium_p, params$Tswidth_fusarium_p),
    mildew = get_downy_mildew_stress(mts$T, mts$RH, mts$wet,
                                     params$rh_mildew_threshold_p, params$impact_rh_mildew_t,
                                     params$impact_days_wet_mildew_t, params$Topt_mildew_p, params$Twidth_mildew_p)
    
  )
  list(em = risks_em, vg = risks_vg, bl = risks_bl, mt = risks_mt)
})

helper_function(attach_to_global = TRUE)

################################################################################
# Main model
################################################################################
onion_climate_impact <- compiler::cmpfun(function() {
  
  weather_scenario_list <- lapply(weather_precomputed, function(list_per_scenario) {
    sample(list_per_scenario, 1)[[1]]
  })
  
  params <- list(
    # --- Sowing
    T_soil_min_sowing_p = T_soil_min_sowing_p,
    warm_days_needed_sowing_p = warm_days_needed_sowing_p,
    frost_buffer_days_sowing_p = frost_buffer_days_sowing_p,
    min_yday_sowing_c = min_yday_sowing_c,
    max_yday_sowing_c = max_yday_sowing_c,
    rain_window_days_sowing_p = rain_window_days_sowing_p,
    rain_sum_max_mm_sowing_p = rain_sum_max_mm_sowing_p,
    min_dry_days_sowing_c = min_dry_days_sowing_c,
    
    # --- Drought
    rh_drought_threshold_p = rh_drought_threshold_p,
    Tavg_drought_threshold_p = Tavg_drought_threshold_p,
    impact_rh_drought_t = impact_rh_drought_t,
    impact_days_dry_drought_t = impact_days_dry_drought_t,
    impact_temp_drought_t = impact_temp_drought_t,
    
    # --- Extreme rainfall
    prec_extreme_rain_medium_p = prec_extreme_rain_medium_p,
    prec_extreme_rain_high_p = prec_extreme_rain_high_p,
    impact_days_extreme_rain_t = impact_days_extreme_rain_t,
    
    # --- Hail
    prec_hail_threshold_p = prec_hail_threshold_p,
    Tavg_hail_threshold_p = Tavg_hail_threshold_p,
    impact_days_hail_t = impact_days_hail_t,
    
    # --- Botrytis
    rh_botrytis_threshold_p = rh_botrytis_threshold_p,
    impact_rh_botrytis_t = impact_rh_botrytis_t,
    impact_days_wet_botrytis_t = impact_days_wet_botrytis_t,
    Topt_botrytis_p = Topt_botrytis_p,
    Twidth_botrytis_p = Twidth_botrytis_p,
    
    # --- Mildew
    rh_mildew_threshold_p = rh_mildew_threshold_p,
    impact_rh_mildew_t = impact_rh_mildew_t,
    impact_days_wet_mildew_t = impact_days_wet_mildew_t,
    Topt_mildew_p = Topt_mildew_p,
    Twidth_mildew_p = Twidth_mildew_p,
    
    # --- Thrips
    rh_thrips_threshold_p = rh_thrips_threshold_p,
    Topt_thrips_p = Topt_thrips_p,
    Twidth_thrips_p = Twidth_thrips_p,
    impact_rh_thrips_t = impact_rh_thrips_t,
    impact_days_dry_thrips_t = impact_days_dry_thrips_t,
    
    # --- Fusarium
    rh_fusarium_threshold_p = rh_fusarium_threshold_p,
    impact_rh_fusarium_t = impact_rh_fusarium_t,
    impact_days_wet_fusarium_t = impact_days_wet_fusarium_t,
    Tsopt_fusarium_p = Tsopt_fusarium_p,
    Tswidth_fusarium_p = Tswidth_fusarium_p,
    
    # --- Yield reductions
    yield_reduction_drought_t = yield_reduction_drought_t,
    yield_reduction_extreme_rain_t = yield_reduction_extreme_rain_t,
    yield_reduction_hail_t = yield_reduction_hail_t,
    yield_reduction_fusarium_t = yield_reduction_fusarium_t,
    yield_reduction_botrytis_t = yield_reduction_botrytis_t,
    yield_reduction_downy_mildew_t = yield_reduction_downy_mildew_t,
    yield_reduction_thrips_t = yield_reduction_thrips_t,

    
    # --- Biomass / yield
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
    dry_onion_weight_t = dry_onion_weight_t,
    HI_onions_t = HI_onions_t
  )
  
  results <- vector("list", length(weather_scenario_list))
  names(results) <- names(weather_scenario_list)
  
  for (sc in names(weather_scenario_list)) {
    df <- weather_scenario_list[[sc]]
    setDT(df)
    
    df$day_consec_wet <- consec_counter(!is.na(df$Prec) & df$Prec > 0)
    df$day_consec_dry <- consec_counter(is.na(df$Prec) | df$Prec == 0)
    
    sow <- compute_onion_sowing_yday_fast(df,
                                          params$T_soil_min_sowing_p, params$warm_days_needed_sowing_p,
                                          params$frost_buffer_days_sowing_p, params$min_yday_sowing_c,
                                          params$max_yday_sowing_c, params$rain_window_days_sowing_p,
                                          params$rain_sum_max_mm_sowing_p, params$min_dry_days_sowing_c)
    
    phase_data <- precompute_phase_data(df, sow,
                                        GDD_field_emergence_required_p, GDD_vegetative_required_p,
                                        GDD_bulbing_required_p, GDD_maturation_required_p)
    
    all_risks <- compute_all_risks(phase_data, df, params)
    
    multipliers <- c(
      fast_chance_event(all_risks$em["drought"],  params$yield_reduction_drought_t),
      fast_chance_event(all_risks$em["exrain"],   params$yield_reduction_extreme_rain_t),
      fast_chance_event(all_risks$em["hail"],     params$yield_reduction_hail_t),
      fast_chance_event(all_risks$em["fusarium"], params$yield_reduction_fusarium_t),
      
      fast_chance_event(all_risks$vg["drought"],  params$yield_reduction_drought_t),
      fast_chance_event(all_risks$vg["exrain"],   params$yield_reduction_extreme_rain_t),
      fast_chance_event(all_risks$vg["hail"],     params$yield_reduction_hail_t),
      fast_chance_event(all_risks$vg["fusarium"], params$yield_reduction_fusarium_t),
      fast_chance_event(all_risks$vg["mildew"],   params$yield_reduction_downy_mildew_t),
      fast_chance_event(all_risks$vg["thrips"],   params$yield_reduction_thrips_t),
      
      fast_chance_event(all_risks$bl["drought"],  params$yield_reduction_drought_t),
      fast_chance_event(all_risks$bl["exrain"],   params$yield_reduction_extreme_rain_t),
      fast_chance_event(all_risks$bl["hail"],     params$yield_reduction_hail_t),
      fast_chance_event(all_risks$bl["botrytis"], params$yield_reduction_botrytis_t),
      fast_chance_event(all_risks$bl["fusarium"], params$yield_reduction_fusarium_t),
      fast_chance_event(all_risks$bl["mildew"],   params$yield_reduction_downy_mildew_t),
      
      fast_chance_event(all_risks$mt["drought"],  params$yield_reduction_drought_t),
      fast_chance_event(all_risks$mt["exrain"],   params$yield_reduction_extreme_rain_t),
      fast_chance_event(all_risks$mt["hail"],     params$yield_reduction_hail_t),
      fast_chance_event(all_risks$mt["botrytis"], params$yield_reduction_botrytis_t),
      fast_chance_event(all_risks$mt["fusarium"], params$yield_reduction_fusarium_t),
      fast_chance_event(all_risks$mt["mildew"],   params$yield_reduction_downy_mildew_t)
    )
    
    combined_yield_multiplier <- prod(multipliers, na.rm = TRUE)
    
    biomass_em <- sum(calc_bio_vectorized(
      df$PAR[phase_data$em$idx], params$LAI_emergence_p,
      df$Tavg[phase_data$em$idx], df$Prec[phase_data$em$idx],
      params$f_T_1_lower_p, params$f_T_1_upper_p,
      params$f_T_0_lower_p, params$f_T_0_upper_p,
      params$f_W_1_lower_p, params$f_W_1_upper_p,
      params$f_W_0.5_p, params$LUE_onion_p, params$lec_k_c), na.rm = TRUE)
    
    biomass_vg <- sum(calc_bio_vectorized(
      df$PAR[phase_data$vg$idx], params$LAI_veg_p,
      df$Tavg[phase_data$vg$idx], df$Prec[phase_data$vg$idx],
      params$f_T_1_lower_p, params$f_T_1_upper_p,
      params$f_T_0_lower_p, params$f_T_0_upper_p,
      params$f_W_1_lower_p, params$f_W_1_upper_p,
      params$f_W_0.5_p, params$LUE_onion_p, params$lec_k_c), na.rm = TRUE)
    
    biomass_bl <- sum(calc_bio_vectorized(
      df$PAR[phase_data$bl$idx], params$LAI_bulbing_p,
      df$Tavg[phase_data$bl$idx], df$Prec[phase_data$bl$idx],
      params$f_T_1_lower_p, params$f_T_1_upper_p,
      params$f_T_0_lower_p, params$f_T_0_upper_p,
      params$f_W_1_lower_p, params$f_W_1_upper_p,
      params$f_W_0.5_p, params$LUE_onion_p, params$lec_k_c), na.rm = TRUE)
    
    biomass_mt <- sum(calc_bio_vectorized(
      df$PAR[phase_data$mt$idx], params$LAI_maturation_p,
      df$Tavg[phase_data$mt$idx], df$Prec[phase_data$mt$idx],
      params$f_T_1_lower_p, params$f_T_1_upper_p,
      params$f_T_0_lower_p, params$f_T_0_upper_p,
      params$f_W_1_lower_p, params$f_W_1_upper_p,
      params$f_W_0.5_p, params$LUE_onion_p, params$lec_k_c), na.rm = TRUE)
    
    total_biomass <- biomass_em + biomass_vg + biomass_bl + biomass_mt
    
    raw_yield <- (total_biomass * params$onions_per_ha_p *
                    (1 - params$dry_onion_weight_t)) * params$HI_onions_t / 1e6
    
    final_yield <- raw_yield * combined_yield_multiplier
    
    results[[sc]] <- list(
      sowing_yday = sow,
      harvest_yday <- if (any(phase_data$mt$idx)) {
        max(df$yday[phase_data$mt$idx])
      } else {
        max(df$yday[after_sow])  # fallback: last sowed day
      }
      ,
      total_biomass = total_biomass,
      raw_yield_per_ha = raw_yield,
      final_yield_per_ha = final_yield,
      combined_yield_multiplier = combined_yield_multiplier,
      
      # expose multipliers for diagnostics
      m_drought_em = multipliers[1],
      m_exrain_em  = multipliers[2],
      m_hail_em    = multipliers[3],
      m_fus_em     = multipliers[4],
    
      
      m_drought_vg = multipliers[5],
      m_exrain_vg  = multipliers[6],
      m_hail_vg    = multipliers[7],
      m_fus_vg     = multipliers[8],
      m_mildew_vg  = multipliers[9],
      m_thrips_vg  = multipliers[10],
   
      
      m_drought_bl = multipliers[11],
      m_exrain_bl  = multipliers[12],
      m_hail_bl    = multipliers[13],
      m_botrytis_bl= multipliers[14],
      m_fus_bl     = multipliers[15],
      m_mildew_bl  = multipliers[16],
     
      
      m_drought_mt = multipliers[17],
      m_exrain_mt  = multipliers[18],
      m_hail_mt    = multipliers[19],
      m_botrytis_mt= multipliers[20],
      m_fus_mt     = multipliers[21],
      m_mildew_mt  = multipliers[22]
      
    )
  }
  
 return(results)
})
################################################################################
# Monte Carlo Simulation
################################################################################

input_variables <- read.csv("input_table_onion.csv", header = TRUE, sep = ";")

  onion_mc_simulation <- mcSimulation(
    estimate = as.estimate(input_variables),
    model_function = onion_climate_impact,
    numberOfModelRuns = num_simulations_c,
    functionSyntax = "plainNames"
  )

