# Onion Climate Impact Model ----
#
# This script implements the core onion climate impact model used in the
# Monte Carlo simulations. It does the following:
#
# 1. Loads required packages and enables JIT compilation for speed.
# 2. Defines a function `precompute_phase_data()` that:
#      - uses cumulative GDD to split the season into four phases:
#        emergence (em), vegetative (vg), bulbing (bl), maturation (mt).
# 3. Defines a function `compute_all_risks()` that:
#      - aggregates weather (T, RH, rain, soil T, wet/dry streaks) by phase
#      - computes phase-level risk scores for drought, extreme rain, hail,
#        Botrytis, downy mildew, Fusarium, and thrips using functions.
# 4. Calls `helper_function(attach_to_global = TRUE)` to load all stress and
#    biomass helpers from the helper script.
# 5. Defines the main model function `onion_climate_impact()`:
#      - samples one season per scenario from `weather_precomputed`
#      - computes sowing date
#      - computes phase-wise risks and resulting biomass multipliers
#      - computes potential vs stressed biomass and yield
#      - returns results per scenario (historical, ssp126, ssp245, ssp370, ssp585)
# 6. Runs a Monte Carlo simulation via `mcSimulation()` using decisionSupport.

 

# Load required packages ----

load_if_needed <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
    }
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

# Enable just-in-time compilation for extra speed (level 3 = aggressive)
compiler::enableJIT(3)


# Phase definition via cumulative GDD ----

# Pre-compute which days belong to which growth phase (em, vg, bl, mt)
# using cumulative Growing Degree Days (GDD) after sowing.
precompute_phase_data <- compiler::cmpfun(function(df, sow,
                                                   GDD_field_emergence_required_p,
                                                   GDD_vegetative_required_p,
                                                   GDD_bulbing_required_p,
                                                   GDD_maturation_required_p) {
  n <- nrow(df)
  if (n == 0L) return(NULL)
  
  # Initialize cumulative GDD
  df[, GDD_cum := 0]
  
  # Only accumulate after sowing
  after_sow <- df$yday >= sow
  df$GDD_cum[after_sow] <- cumsum(df$GDD_daily[after_sow])
  
  # Phase thresholds (GDD-based)
  gdd_em <- GDD_field_emergence_required_p
  gdd_vg <- GDD_vegetative_required_p
  gdd_bl <- GDD_bulbing_required_p
  gdd_mt <- GDD_maturation_required_p
  
  # Logical indices for each phase
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


# Batch risk computation per phase ----

# Compute all stress risks per phase (em, vg, bl, mt) in one go.
# Uses phase means (T, RH, soil T, rain) and maximum wet/dry streaks.
compute_all_risks <- compiler::cmpfun(function(phase_data, df, params) {
  
  # Helper to compute phase-level summary statistics
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
  
  # --- EMERGENCE phase risks ---
  risks_em <- c(
    drought = get_drought_stress(
      ems$T, ems$RH, ems$dry,
      params$rh_drought_threshold_p, params$Tavg_drought_threshold_p,
      params$impact_rh_drought_t, params$impact_days_dry_drought_t, params$impact_temp_drought_t
    ),
    exrain = get_extreme_rain_stress(
      df$Prec[phase_data$em$idx],
      params$prec_extreme_rain_medium_p, params$prec_extreme_rain_high_p,
      params$impact_days_extreme_rain_t
    ),
    hail = get_hail_stress(
      df$Tavg[phase_data$em$idx], df$Prec[phase_data$em$idx],
      params$prec_hail_threshold_p, params$Tavg_hail_threshold_p, params$impact_days_hail_t
    ),
    fusarium = get_fusarium_stress(
      ems$Ts, ems$RH, ems$wet,
      params$rh_fusarium_threshold_p, params$impact_rh_fusarium_t,
      params$impact_days_wet_fusarium_t, params$Tsopt_fusarium_p, params$Tswidth_fusarium_p
    )
  )
  
  # --- VEGETATIVE phase risks ---
  risks_vg <- c(
    drought = get_drought_stress(
      vgs$T, vgs$RH, vgs$dry,
      params$rh_drought_threshold_p, params$Tavg_drought_threshold_p,
      params$impact_rh_drought_t, params$impact_days_dry_drought_t, params$impact_temp_drought_t
    ),
    exrain = get_extreme_rain_stress(
      df$Prec[phase_data$vg$idx],
      params$prec_extreme_rain_medium_p, params$prec_extreme_rain_high_p,
      params$impact_days_extreme_rain_t
    ),
    hail = get_hail_stress(
      df$Tavg[phase_data$vg$idx], df$Prec[phase_data$vg$idx],
      params$prec_hail_threshold_p, params$Tavg_hail_threshold_p, params$impact_days_hail_t
    ),
    mildew = get_downy_mildew_stress(
      vgs$T, vgs$RH, vgs$wet,
      params$rh_mildew_threshold_p, params$impact_rh_mildew_t,
      params$impact_days_wet_mildew_t, params$Topt_mildew_p, params$Twidth_mildew_p
    ),
    thrips = get_thrips_stress(
      vgs$T, vgs$RH, vgs$dry,
      params$rh_thrips_threshold_p, params$Topt_thrips_p, params$Twidth_thrips_p,
      params$impact_rh_thrips_t, params$impact_days_dry_thrips_t
    ),
    fusarium = get_fusarium_stress(
      vgs$Ts, vgs$RH, vgs$wet,
      params$rh_fusarium_threshold_p, params$impact_rh_fusarium_t,
      params$impact_days_wet_fusarium_t, params$Tsopt_fusarium_p, params$Tswidth_fusarium_p
    )
  )
  
  # --- BULBING phase risks ---
  risks_bl <- c(
    drought = get_drought_stress(
      bls$T, bls$RH, bls$dry,
      params$rh_drought_threshold_p, params$Tavg_drought_threshold_p,
      params$impact_rh_drought_t, params$impact_days_dry_drought_t, params$impact_temp_drought_t
    ),
    botrytis = get_botrytis_stress(
      bls$T, bls$RH, bls$wet,
      params$rh_botrytis_threshold_p, params$impact_rh_botrytis_t,
      params$impact_days_wet_botrytis_t, params$Topt_botrytis_p, params$Twidth_botrytis_p
    ),
    mildew = get_downy_mildew_stress(
      bls$T, bls$RH, bls$wet,
      params$rh_mildew_threshold_p, params$impact_rh_mildew_t,
      params$impact_days_wet_mildew_t, params$Topt_mildew_p, params$Twidth_mildew_p
    ),
    fusarium = get_fusarium_stress(
      bls$Ts, bls$RH, bls$wet,
      params$rh_fusarium_threshold_p, params$impact_rh_fusarium_t,
      params$impact_days_wet_fusarium_t, params$Tsopt_fusarium_p, params$Tswidth_fusarium_p
    )
  )
  
  # --- MATURATION phase risks ---
  risks_mt <- c(
    drought = get_drought_stress(
      mts$T, mts$RH, mts$dry,
      params$rh_drought_threshold_p, params$Tavg_drought_threshold_p,
      params$impact_rh_drought_t, params$impact_days_dry_drought_t, params$impact_temp_drought_t
    ),
    exrain = get_extreme_rain_stress(
      df$Prec[phase_data$mt$idx],
      params$prec_extreme_rain_medium_p, params$prec_extreme_rain_high_p,
      params$impact_days_extreme_rain_t
    ),
    hail = get_hail_stress(
      df$Tavg[phase_data$mt$idx], df$Prec[phase_data$mt$idx],
      params$prec_hail_threshold_p, params$Tavg_hail_threshold_p, params$impact_days_hail_t
    ),
    botrytis = get_botrytis_stress(
      mts$T, mts$RH, mts$wet,
      params$rh_botrytis_threshold_p, params$impact_rh_botrytis_t,
      params$impact_days_wet_botrytis_t, params$Topt_botrytis_p, params$Twidth_botrytis_p
    ),
    fusarium = get_fusarium_stress(
      mts$Ts, mts$RH, mts$wet,
      params$rh_fusarium_threshold_p, params$impact_rh_fusarium_t,
      params$impact_days_wet_fusarium_t, params$Tsopt_fusarium_p, params$Tswidth_fusarium_p
    ),
    mildew = get_downy_mildew_stress(
      mts$T, mts$RH, mts$wet,
      params$rh_mildew_threshold_p, params$impact_rh_mildew_t,
      params$impact_days_wet_mildew_t, params$Topt_mildew_p, params$Twidth_mildew_p
    )
  )
  
  list(em = risks_em, vg = risks_vg, bl = risks_bl, mt = risks_mt)
})


# Load compiled helper functions (stress + biomass + sowing) ----

# This assumes you've sourced the helper script where helper_function() is defined.
helper_function(attach_to_global = TRUE)


# Main onion climate model function ----

onion_climate_impact <- compiler::cmpfun(function() {
  
  # 1) For each scenario, randomly sample one season (one id_season) ----
  weather_scenario_list <- lapply(weather_precomputed, function(list_per_scenario) {
    sample(list_per_scenario, 1)[[1]]
  })
  
  # 2) Bundle parameters into a single list for passing around ----
  params <- list(
    # --- Sowing
    T_soil_min_sowing_p        = T_soil_min_sowing_p,
    warm_days_needed_sowing_p  = warm_days_needed_sowing_p,
    frost_buffer_days_sowing_p = frost_buffer_days_sowing_p,
    min_yday_sowing_c          = min_yday_sowing_c,
    max_yday_sowing_c          = max_yday_sowing_c,
    rain_window_days_sowing_p  = rain_window_days_sowing_p,
    rain_sum_max_mm_sowing_p   = rain_sum_max_mm_sowing_p,
    min_dry_days_sowing_c      = min_dry_days_sowing_c,
    
    # --- Drought
    rh_drought_threshold_p     = rh_drought_threshold_p,
    Tavg_drought_threshold_p   = Tavg_drought_threshold_p,
    impact_rh_drought_t        = impact_rh_drought_t,
    impact_days_dry_drought_t  = impact_days_dry_drought_t,
    impact_temp_drought_t      = impact_temp_drought_t,
    
    # --- Extreme rainfall
    prec_extreme_rain_medium_p = prec_extreme_rain_medium_p,
    prec_extreme_rain_high_p   = prec_extreme_rain_high_p,
    impact_days_extreme_rain_t = impact_days_extreme_rain_t,
    
    # --- Hail
    prec_hail_threshold_p      = prec_hail_threshold_p,
    Tavg_hail_threshold_p      = Tavg_hail_threshold_p,
    impact_days_hail_t         = impact_days_hail_t,
    
    # --- Botrytis
    rh_botrytis_threshold_p    = rh_botrytis_threshold_p,
    impact_rh_botrytis_t       = impact_rh_botrytis_t,
    impact_days_wet_botrytis_t = impact_days_wet_botrytis_t,
    Topt_botrytis_p            = Topt_botrytis_p,
    Twidth_botrytis_p          = Twidth_botrytis_p,
    
    # --- Mildew
    rh_mildew_threshold_p      = rh_mildew_threshold_p,
    impact_rh_mildew_t         = impact_rh_mildew_t,
    impact_days_wet_mildew_t   = impact_days_wet_mildew_t,
    Topt_mildew_p              = Topt_mildew_p,
    Twidth_mildew_p            = Twidth_mildew_p,
    
    # --- Thrips
    rh_thrips_threshold_p      = rh_thrips_threshold_p,
    Topt_thrips_p              = Topt_thrips_p,
    Twidth_thrips_p            = Twidth_thrips_p,
    impact_rh_thrips_t         = impact_rh_thrips_t,
    impact_days_dry_thrips_t   = impact_days_dry_thrips_t,
    
    # --- Fusarium
    rh_fusarium_threshold_p    = rh_fusarium_threshold_p,
    impact_rh_fusarium_t       = impact_rh_fusarium_t,
    impact_days_wet_fusarium_t = impact_days_wet_fusarium_t,
    Tsopt_fusarium_p           = Tsopt_fusarium_p,
    Tswidth_fusarium_p         = Tswidth_fusarium_p,
    
    # --- Yield reductions (fractional impact if event happens)
    yield_reduction_drought_t       = yield_reduction_drought_t,
    yield_reduction_extreme_rain_t  = yield_reduction_extreme_rain_t,
    yield_reduction_hail_t          = yield_reduction_hail_t,
    yield_reduction_fusarium_t      = yield_reduction_fusarium_t,
    yield_reduction_botrytis_t      = yield_reduction_botrytis_t,
    yield_reduction_downy_mildew_t  = yield_reduction_downy_mildew_t,
    yield_reduction_thrips_t        = yield_reduction_thrips_t,
    
    # --- Biomass / yield / canopy
    LAI_emergence_p = LAI_emergence_p,
    LAI_veg_p       = LAI_veg_p,
    LAI_bulbing_p   = LAI_bulbing_p,
    LAI_maturation_p= LAI_maturation_p,
    f_T_1_lower_p   = f_T_1_lower_p,
    f_T_1_upper_p   = f_T_1_upper_p,
    f_T_0_lower_p   = f_T_0_lower_p,
    f_T_0_upper_p   = f_T_0_upper_p,
    f_W_1_lower_p   = f_W_1_lower_p,
    f_W_1_upper_p   = f_W_1_upper_p,
    f_W_0.5_p       = f_W_0.5_p,
    LUE_onion_p     = LUE_onion_p,
    lec_k_c         = lec_k_c,
    onions_per_ha_p = onions_per_ha_p,
    dry_onion_weight_t = dry_onion_weight_t,
    HI_onions_t        = HI_onions_t
  )
  
  # 3) Prepare result container per scenario ----
  results <- vector("list", length(weather_scenario_list))
  names(results) <- names(weather_scenario_list)
  
  # 4) Loop over scenarios ----
  for (sc in names(weather_scenario_list)) {
    df <- weather_scenario_list[[sc]]
    data.table::setDT(df)
    
    # Derived streak variables (wet / dry days)
    df$day_consec_wet <- consec_counter(!is.na(df$Prec) & df$Prec > 0)
    df$day_consec_dry <- consec_counter(is.na(df$Prec) | df$Prec == 0)
    
    # Sowing date for this season and scenario
    sow <- compute_onion_sowing_yday_fast(
      df,
      params$T_soil_min_sowing_p, params$warm_days_needed_sowing_p,
      params$frost_buffer_days_sowing_p, params$min_yday_sowing_c,
      params$max_yday_sowing_c, params$rain_window_days_sowing_p,
      params$rain_sum_max_mm_sowing_p, params$min_dry_days_sowing_c
    )
    
    # Phase data (em, vg, bl, mt) based on cumulative GDD
    phase_data <- precompute_phase_data(
      df, sow,
      GDD_field_emergence_required_p, GDD_vegetative_required_p,
      GDD_bulbing_required_p, GDD_maturation_required_p
    )
    
    # Phase-wise risk scores
    all_risks <- compute_all_risks(phase_data, df, params)
    
    # 5) Convert risks into random multipliers (event may or may not occur) ----
    
    # Emergence phase multipliers
    em_multipliers <- c(
      fast_chance_event(all_risks$em["drought"],  params$yield_reduction_drought_t),
      fast_chance_event(all_risks$em["exrain"],   params$yield_reduction_extreme_rain_t),
      fast_chance_event(all_risks$em["hail"],     params$yield_reduction_hail_t),
      fast_chance_event(all_risks$em["fusarium"], params$yield_reduction_fusarium_t)
    )
    
    # Vegetative phase multipliers
    vg_multipliers <- c(
      fast_chance_event(all_risks$vg["drought"],  params$yield_reduction_drought_t),
      fast_chance_event(all_risks$vg["exrain"],   params$yield_reduction_extreme_rain_t),
      fast_chance_event(all_risks$vg["hail"],     params$yield_reduction_hail_t),
      fast_chance_event(all_risks$vg["fusarium"], params$yield_reduction_fusarium_t),
      fast_chance_event(all_risks$vg["mildew"],   params$yield_reduction_downy_mildew_t),
      fast_chance_event(all_risks$vg["thrips"],   params$yield_reduction_thrips_t)
    )
    
    # Bulbing phase multipliers
    bl_multipliers <- c(
      fast_chance_event(all_risks$bl["drought"],  params$yield_reduction_drought_t),
      fast_chance_event(all_risks$bl["exrain"],   params$yield_reduction_extreme_rain_t),
      fast_chance_event(all_risks$bl["hail"],     params$yield_reduction_hail_t),
      fast_chance_event(all_risks$bl["botrytis"], params$yield_reduction_botrytis_t),
      fast_chance_event(all_risks$bl["fusarium"], params$yield_reduction_fusarium_t),
      fast_chance_event(all_risks$bl["mildew"],   params$yield_reduction_downy_mildew_t)
    )
    
    # Maturation phase multipliers
    mt_multipliers <- c(
      fast_chance_event(all_risks$mt["drought"],  params$yield_reduction_drought_t),
      fast_chance_event(all_risks$mt["exrain"],   params$yield_reduction_extreme_rain_t),
      fast_chance_event(all_risks$mt["hail"],     params$yield_reduction_hail_t),
      fast_chance_event(all_risks$mt["botrytis"], params$yield_reduction_botrytis_t),
      fast_chance_event(all_risks$mt["fusarium"], params$yield_reduction_fusarium_t),
      fast_chance_event(all_risks$mt["mildew"],   params$yield_reduction_downy_mildew_t)
    )
    
    # Phase-wise biomass multipliers (product of individual stress effects)
    em_bio_multiplier <- prod(em_multipliers, na.rm = TRUE)
    vg_bio_multiplier <- prod(vg_multipliers, na.rm = TRUE)
    bl_bio_multiplier <- prod(bl_multipliers, na.rm = TRUE)
    mt_bio_multiplier <- prod(mt_multipliers, na.rm = TRUE)
    
    # 6) Compute potential biomass (without stress) per phase ----
    
    biomass_em_pot <- sum(calc_bio_vectorized(
      df$PAR[phase_data$em$idx], params$LAI_emergence_p,
      df$Tavg[phase_data$em$idx], df$Prec[phase_data$em$idx],
      params$f_T_1_lower_p, params$f_T_1_upper_p,
      params$f_T_0_lower_p, params$f_T_0_upper_p,
      params$f_W_1_lower_p, params$f_W_1_upper_p,
      params$f_W_0.5_p, params$LUE_onion_p, params$lec_k_c
    ), na.rm = TRUE)
    
    biomass_vg_pot <- sum(calc_bio_vectorized(
      df$PAR[phase_data$vg$idx], params$LAI_veg_p,
      df$Tavg[phase_data$vg$idx], df$Prec[phase_data$vg$idx],
      params$f_T_1_lower_p, params$f_T_1_upper_p,
      params$f_T_0_lower_p, params$f_T_0_upper_p,
      params$f_W_1_lower_p, params$f_W_1_upper_p,
      params$f_W_0.5_p, params$LUE_onion_p, params$lec_k_c
    ), na.rm = TRUE)
    
    biomass_bl_pot <- sum(calc_bio_vectorized(
      df$PAR[phase_data$bl$idx], params$LAI_bulbing_p,
      df$Tavg[phase_data$bl$idx], df$Prec[phase_data$bl$idx],
      params$f_T_1_lower_p, params$f_T_1_upper_p,
      params$f_T_0_lower_p, params$f_T_0_upper_p,
      params$f_W_1_lower_p, params$f_W_1_upper_p,
      params$f_W_0.5_p, params$LUE_onion_p, params$lec_k_c
    ), na.rm = TRUE)
    
    biomass_mt_pot <- sum(calc_bio_vectorized(
      df$PAR[phase_data$mt$idx], params$LAI_maturation_p,
      df$Tavg[phase_data$mt$idx], df$Prec[phase_data$mt$idx],
      params$f_T_1_lower_p, params$f_T_1_upper_p,
      params$f_T_0_lower_p, params$f_T_0_upper_p,
      params$f_W_1_lower_p, params$f_W_1_upper_p,
      params$f_W_0.5_p, params$LUE_onion_p, params$lec_k_c
    ), na.rm = TRUE)
    
    # 7) Apply stress multipliers to potential biomass ----
    
    biomass_em <- biomass_em_pot * em_bio_multiplier
    biomass_vg <- biomass_vg_pot * vg_bio_multiplier
    biomass_bl <- biomass_bl_pot * bl_bio_multiplier
    biomass_mt <- biomass_mt_pot * mt_bio_multiplier
    
    # Potential yield (no stress): based on sum of potential biomass
    potential_biomass <- biomass_em_pot + biomass_vg_pot + biomass_bl_pot + biomass_mt_pot
    potential_yield <- potential_biomass * params$onions_per_ha_p *
      (1 - params$dry_onion_weight_t) * params$HI_onions_t / 1e6  # kg/ha -> t/ha
    
    # Realized total biomass with stress
    total_biomass <- biomass_em + biomass_vg + biomass_bl + biomass_mt
    
    # Realized final yield (t/ha)
    final_yield <- (total_biomass * params$onions_per_ha_p *
                      (1 - params$dry_onion_weight_t)) * params$HI_onions_t / 1e6
    
    # Harvest yday: last day in maturation phase, or last day after sowing
    after_sow_idx <- df$yday >= sow
    harvest_yday <- if (any(phase_data$mt$idx)) {
      max(df$yday[phase_data$mt$idx])
    } else {
      max(df$yday[after_sow_idx], na.rm = TRUE)
    }
    
    # 8) Store scenario results ----
    results[[sc]] <- list(
      sowing_yday        = sow,
      harvest_yday       = harvest_yday,
      total_biomass      = total_biomass,
      raw_yield_per_ha   = potential_yield,
      final_yield_per_ha = final_yield,
      
      em_bio_multiplier  = em_bio_multiplier,
      vg_bio_multiplier  = vg_bio_multiplier,
      bl_bio_multiplier  = bl_bio_multiplier,
      mt_bio_multiplier  = mt_bio_multiplier,
      
      # Expose individual event multipliers for diagnostics (per phase)
      m_drought_em = em_multipliers[1],
      m_exrain_em  = em_multipliers[2],
      m_hail_em    = em_multipliers[3],
      m_fus_em     = em_multipliers[4],
      
      m_drought_vg = vg_multipliers[1],
      m_exrain_vg  = vg_multipliers[2],
      m_hail_vg    = vg_multipliers[3],
      m_fus_vg     = vg_multipliers[4],
      m_mildew_vg  = vg_multipliers[5],
      m_thrips_vg  = vg_multipliers[6],
      
      m_drought_bl = bl_multipliers[1],
      m_exrain_bl  = bl_multipliers[2],
      m_hail_bl    = bl_multipliers[3],
      m_botrytis_bl= bl_multipliers[4],
      m_fus_bl     = bl_multipliers[5],
      m_mildew_bl  = bl_multipliers[6],
      
      m_drought_mt = mt_multipliers[1],
      m_exrain_mt  = mt_multipliers[2],
      m_hail_mt    = mt_multipliers[3],
      m_botrytis_mt= mt_multipliers[4],
      m_fus_mt     = mt_multipliers[5],
      m_mildew_mt  = mt_multipliers[6]
    )
  }
  
  return(results)
})


# Monte Carlo simulation with decisionSupport ----

input_variables <- read.csv("input_table_onion.csv", header = TRUE, sep = ";")

onion_mc_simulation <- mcSimulation(
  estimate        = as.estimate(input_variables),
  model_function  = onion_climate_impact,
  numberOfModelRuns = 100,
  functionSyntax  = "plainNames"
)
