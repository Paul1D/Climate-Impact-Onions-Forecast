# Onion Climate Impact Model ----
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
#        Botrytis, downy mildew, Fusarium, and thrips.
# 4. Loads all helper functions (stress + biomass + sowing).
# 5. Defines the main model function `onion_climate_impact()`:
#      - samples one season per scenario
#      - computes sowing date
#      - computes stress risks and biomass multipliers
#      - computes potential vs realized yield
# 6. Runs Monte Carlo simulation.

#
# Extended with:
# - Geisenheim-style irrigation scheduling (soil water balance with TAW/RAW)
# - Rainfed vs irrigated yield outputs for direct comparison.
#
# Scientific basis:
# - FAO-56 Kc–ET0 approach for ETc and soil water balance.
# - nFK (nutzbare Feldkapazität) as plant-available water in the root zone.
# - Geisenheim irrigation scheduling & onion-specific kc values.
# - Onion: shallow rooting, high sensitivity to water stress, frequent light irrigations recommended.


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
  em_idx <- after_sow & df$GDD_cum <= GDD_field_emergence_required_p
  vg_idx <- after_sow & df$GDD_cum > GDD_field_emergence_required_p &
    df$GDD_cum <= GDD_vegetative_required_p
  bl_idx <- after_sow & df$GDD_cum > GDD_vegetative_required_p &
    df$GDD_cum <= GDD_bulbing_required_p
  mt_idx <- after_sow & df$GDD_cum > GDD_bulbing_required_p &
    df$GDD_cum <= GDD_maturation_required_p
  
  list(
    em = list(idx = em_idx, range = which(em_idx)),
    vg = list(idx = vg_idx, range = which(vg_idx)),
    bl = list(idx = bl_idx, range = which(bl_idx)),
    mt = list(idx = mt_idx, range = which(mt_idx))
  )
})


# Batch risk computation per phase ----
# Compute all stress risks per phase (em, vg, bl, mt) in one go.
# Uses phase means (T, RH, soil T, rain) and maximum wet/dry/hot streaks.
compute_all_risks <- compiler::cmpfun(function(phase_data, df, params) {
  
  stats_for <- function(idx) list(
    T        = mean(df$Tavg[idx],    na.rm = TRUE),
    Tmin     = mean(df$Tmin[idx],    na.rm = TRUE),   # nighttime temp proxy
    RHmean   = mean(df$RH_mean[idx], na.rm = TRUE),
    RHmax    = mean(df$RH_max[idx],  na.rm = TRUE),   # nighttime humidity proxy
    Ts       = mean(df$Ts_5cm[idx],  na.rm = TRUE),
    wet      = if (any(idx, na.rm = TRUE)) max(df$day_consec_wet[idx]) else 0,
    dry      = if (any(idx, na.rm = TRUE)) max(df$day_consec_dry[idx]) else 0,
    hot      = if (any(idx, na.rm = TRUE)) max(df$day_consec_hot[idx]) else 0
  )
  
  ems <- stats_for(phase_data$em$idx)
  vgs <- stats_for(phase_data$vg$idx)
  bls <- stats_for(phase_data$bl$idx)
  mts <- stats_for(phase_data$mt$idx)
  
  # --- EMERGENCE phase risks ---
  risks_em <- c(
    exrain = get_extreme_rain_stress(
      df$Prec[phase_data$em$idx],
      params$prec_extreme_rain_medium_p,
      params$prec_extreme_rain_high_p,
      params$impact_days_extreme_rain_t
    ),
    hail = get_hail_stress(
      df$Tavg[phase_data$em$idx], df$Prec[phase_data$em$idx],
      params$prec_hail_threshold_p, params$Tavg_hail_threshold_p,
      params$impact_days_hail_t
    ),
    fusarium = get_fusarium_stress(
      ems$Ts, ems$RHmean, ems$wet,
      params$rh_fusarium_threshold_p,
      params$impact_rh_fusarium_t,
      params$impact_days_wet_fusarium_t,
      params$Tsopt_fusarium_p,
      params$Tswidth_fusarium_p
    )
  )
  
  # --- VEGETATIVE phase risks ---
  risks_vg <- c(
    exrain = get_extreme_rain_stress(
      df$Prec[phase_data$vg$idx],
      params$prec_extreme_rain_medium_p,
      params$prec_extreme_rain_high_p,
      params$impact_days_extreme_rain_t
    ),
    hail = get_hail_stress(
      df$Tavg[phase_data$vg$idx], df$Prec[phase_data$vg$idx],
      params$prec_hail_threshold_p, params$Tavg_hail_threshold_p,
      params$impact_days_hail_t
    ),
    mildew = get_downy_mildew_stress(
      vgs$Tmin, vgs$RHmax, vgs$wet,
      params$rh_mildew_threshold_p,
      params$impact_rh_mildew_t,
      params$impact_days_wet_mildew_t,
      params$Topt_mildew_p,
      params$Twidth_mildew_p
    ),
    thrips = get_thrips_stress(
      vgs$T, vgs$RHmean, vgs$dry,
      params$rh_thrips_threshold_p,
      params$Topt_thrips_p,
      params$Twidth_thrips_p,
      params$impact_rh_thrips_t,
      params$impact_days_dry_thrips_t
    ),
    fusarium = get_fusarium_stress(
      vgs$Ts, vgs$RHmean, vgs$wet,
      params$rh_fusarium_threshold_p,
      params$impact_rh_fusarium_t,
      params$impact_days_wet_fusarium_t,
      params$Tsopt_fusarium_p,
      params$Tswidth_fusarium_p
    ),
    # heat risk (vg)
    heat = get_heat_stress(
      mean(df$Tmax[phase_data$vg$idx], na.rm = TRUE),
      vgs$hot,
      params$Tmax_heat_threshold_p,
      params$impact_temp_heat_t,
      params$impact_days_heat_t
    )
  )
  
  # --- BULBING phase risks ---
  risks_bl <- c(
    exrain = get_extreme_rain_stress(
      df$Prec[phase_data$bl$idx],
      params$prec_extreme_rain_medium_p,
      params$prec_extreme_rain_high_p,
      params$impact_days_extreme_rain_t
    ),
    hail = get_hail_stress(
      df$Tavg[phase_data$bl$idx], df$Prec[phase_data$bl$idx],
      params$prec_hail_threshold_p, params$Tavg_hail_threshold_p,
      params$impact_days_hail_t
    ),
    botrytis = get_botrytis_stress(
      bls$Tmin, bls$RHmax, bls$wet,
      params$rh_botrytis_threshold_p,
      params$impact_rh_botrytis_t,
      params$impact_days_wet_botrytis_t,
      params$Topt_botrytis_p,
      params$Twidth_botrytis_p
    ),
    mildew = get_downy_mildew_stress(
      bls$Tmin, bls$RHmax, bls$wet,
      params$rh_mildew_threshold_p,
      params$impact_rh_mildew_t,
      params$impact_days_wet_mildew_t,
      params$Topt_mildew_p,
      params$Twidth_mildew_p
    ),
    fusarium = get_fusarium_stress(
      bls$Ts, bls$RHmean, bls$wet,
      params$rh_fusarium_threshold_p,
      params$impact_rh_fusarium_t,
      params$impact_days_wet_fusarium_t,
      params$Tsopt_fusarium_p,
      params$Tswidth_fusarium_p
    ),
    heat = get_heat_stress(
      mean(df$Tmax[phase_data$bl$idx], na.rm = TRUE),
      bls$hot,
      params$Tmax_heat_threshold_p,
      params$impact_temp_heat_t,
      params$impact_days_heat_t
    ),
    thrips = get_thrips_stress(
      bls$T, bls$RHmean, bls$dry,
      params$rh_thrips_threshold_p,
      params$Topt_thrips_p,
      params$Twidth_thrips_p,
      params$impact_rh_thrips_t,
      params$impact_days_dry_thrips_t
    )
  )
  
  # --- MATURATION phase risks ---
  risks_mt <- c(
    exrain = get_extreme_rain_stress(
      df$Prec[phase_data$mt$idx],
      params$prec_extreme_rain_medium_p,
      params$prec_extreme_rain_high_p,
      params$impact_days_extreme_rain_t
    ),
    hail = get_hail_stress(
      df$Tavg[phase_data$mt$idx], df$Prec[phase_data$mt$idx],
      params$prec_hail_threshold_p, params$Tavg_hail_threshold_p,
      params$impact_days_hail_t
    ),
    botrytis = get_botrytis_stress(
      mts$Tmin, mts$RHmax, mts$wet,
      params$rh_botrytis_threshold_p,
      params$impact_rh_botrytis_t,
      params$impact_days_wet_botrytis_t,
      params$Topt_botrytis_p,
      params$Twidth_botrytis_p
    ),
    fusarium = get_fusarium_stress(
      mts$Ts, mts$RHmean, mts$wet,
      params$rh_fusarium_threshold_p,
      params$impact_rh_fusarium_t,
      params$impact_days_wet_fusarium_t,
      params$Tsopt_fusarium_p,
      params$Tswidth_fusarium_p
    ),
    mildew = get_downy_mildew_stress(
      mts$Tmin, mts$RHmax, mts$wet,
      params$rh_mildew_threshold_p,
      params$impact_rh_mildew_t,
      params$impact_days_wet_mildew_t,
      params$Topt_mildew_p,
      params$Twidth_mildew_p
    ),
    heat = get_heat_stress(
      mean(df$Tmax[phase_data$mt$idx], na.rm = TRUE),
      mts$hot,
      params$Tmax_heat_threshold_p,
      params$impact_temp_heat_t,
      params$impact_days_heat_t
    ),
    thrips = get_thrips_stress(
      mts$T, mts$RHmean, mts$dry,
      params$rh_thrips_threshold_p,
      params$Topt_thrips_p,
      params$Twidth_thrips_p,
      params$impact_rh_thrips_t,
      params$impact_days_dry_thrips_t
    )
  )
  
  list(em = risks_em, vg = risks_vg, bl = risks_bl, mt = risks_mt)
})


# Load compiled helper functions (stress + biomass + sowing) ----
helper_function(attach_to_global = TRUE)


# Irrigation scheduling (Geisenheim-style, phase-based) ----
# Uses a simple root-zone water balance:
#   Dr(t) = Dr(t-1) + ETc - P_eff - Irrig
# with ETc = kc * ET0 and Dr constrained between 0 and TAW.
compute_irrigation_geisenheim <- compiler::cmpfun(function(df, sow_yday, params, phase_data) {
  n <- nrow(df)
  irrig <- numeric(n)
  
  # map each day to a numeric phase: 1=em, 2=vg, 3=bl, 4=mt
  phase_id <- integer(n)
  phase_id[phase_data$em$idx] <- 1L
  phase_id[phase_data$vg$idx] <- 2L
  phase_id[phase_data$bl$idx] <- 3L
  phase_id[phase_data$mt$idx] <- 4L
  
  # per-phase kc and rooting depth (cm)
  kc_phase   <- c(params$kc_em, params$kc_vg, params$kc_bl, params$kc_mt)
  root_phase <- c(params$root_em_cm, params$root_vg_cm, params$root_bl_cm, params$root_mt_cm)
  
  # soil parameters (canonical local names)
  nFK_mm_per_cm      <- params$nFK_mm_per_cm
  fc_init_frac       <- params$fc_init_frac
  fc_target_frac     <- params$fc_target_frac
  deplete_trigger_fr <- params$deplete_trigger_fr
  max_irrig_mm       <- params$max_irrig_mm
  
  # state
  Dr       <- NA_real_     # soil water depletion (mm)
  prev_TAW <- NA_real_     # previous total available water (mm)
  
  for (i in seq_len(n)) {
    # before sowing: skip
    if (df$yday[i] < sow_yday) next
    
    pid <- phase_id[i]
    if (pid == 0L) next  # outside em/vg/bl/mt
    
    kc   <- kc_phase[pid]
    root <- root_phase[pid]
    
    if (is.na(kc) || is.na(root) || kc <= 0 || root <= 0) next
    
    # total available water (TAW) in current root zone
    TAW <- root * nFK_mm_per_cm
    
    # initialize depletion at (1 - fc_init_frac)*TAW (e.g. 10% depleted)
    if (is.na(Dr)) {
      Dr <- (1 - fc_init_frac) * TAW
    } else if (!is.na(prev_TAW) && prev_TAW > 0 && TAW > 0 && !isTRUE(all.equal(TAW, prev_TAW))) {
      # rescale depletion if rooting depth changed
      Dr <- Dr * (TAW / prev_TAW)
    }
    prev_TAW <- TAW
    
    # crop evapotranspiration and effective precipitation
    ETc  <- df$ET0_mm[i] * kc
    Peff <- df$Prec[i]   # you can later add an efficiency factor here
    
    # update depletion
    Dr <- Dr + ETc - Peff
    if (Dr < 0)   Dr <- 0
    if (Dr > TAW) Dr <- TAW
    
    # irrigation trigger (RAW threshold)
    if (Dr >= deplete_trigger_fr * TAW) {
      target_D <- (1 - fc_target_frac) * TAW
      gift     <- Dr - target_D
      if (gift < 0) gift <- 0
      if (gift > max_irrig_mm) gift <- max_irrig_mm
      
      if (gift > 0) {
        irrig[i] <- gift
        Dr <- Dr - gift
      }
    }
  }
  
  irrig
})


# MAIN MODEL FUNCTION ----------------------------------------------------

onion_climate_impact <- compiler::cmpfun(function() {
  
  # 1) For each scenario, randomly sample one season (one id_season) ----
  weather_scenario_list <- lapply(weather_precomputed, function(list_per_scenario) {
    sample(list_per_scenario, 1)[[1]]
  })
  
  # 2) Bundle parameters into a single list for passing around ----
  #    We map from your input variable names (*_p, *_c, *_t) to
  #    internal canonical names used by the functions.
  params <- list(
    # irrigation (soil water & kc / root by phase)
    nFK_mm_per_cm      = nfk_mm_per_cm_t,   # from input_table
    fc_init_frac       = fc_init_frac_c,
    fc_target_frac     = fc_target_frac_c,
    deplete_trigger_fr = deplete_trigger_fr_t,
    max_irrig_mm       = max_irrig_mm_p,
    
    root_em_cm         = root_em_cm_p,
    root_vg_cm         = root_vg_cm_p,
    root_bl_cm         = root_bl_cm_p,
    root_mt_cm         = root_mt_cm_p,
    
    kc_em              = kc_em_p,
    kc_vg              = kc_vg_p,
    kc_bl              = kc_bl_p,
    kc_mt              = kc_mt_p,
    
    # extreme rain
    prec_extreme_rain_medium_p = prec_extreme_rain_medium_p,
    prec_extreme_rain_high_p   = prec_extreme_rain_high_p,
    impact_days_extreme_rain_t = impact_days_extreme_rain_t,
    
    # hail
    prec_hail_threshold_p      = prec_hail_threshold_p,
    Tavg_hail_threshold_p      = Tavg_hail_threshold_p,
    impact_days_hail_t         = impact_days_hail_t,
    
    # heat
    Tmax_heat_threshold_p      = Tmax_heat_threshold_p,
    impact_temp_heat_t         = impact_temp_heat_t,
    impact_days_heat_t         = impact_days_heat_t,
    
    # botrytis
    rh_botrytis_threshold_p    = rh_botrytis_threshold_p,
    impact_rh_botrytis_t       = impact_rh_botrytis_t,
    impact_days_wet_botrytis_t = impact_days_wet_botrytis_t,
    Topt_botrytis_p            = Topt_botrytis_p,
    Twidth_botrytis_p          = Twidth_botrytis_p,
    
    # mildew
    rh_mildew_threshold_p      = rh_mildew_threshold_p,
    impact_rh_mildew_t         = impact_rh_mildew_t,
    impact_days_wet_mildew_t   = impact_days_wet_mildew_t,
    Topt_mildew_p              = Topt_mildew_p,
    Twidth_mildew_p            = Twidth_mildew_p,
    
    # thrips
    rh_thrips_threshold_p      = rh_thrips_threshold_p,
    Topt_thrips_p              = Topt_thrips_p,
    Twidth_thrips_p            = Twidth_thrips_p,
    impact_rh_thrips_t         = impact_rh_thrips_t,
    impact_days_dry_thrips_t   = impact_days_dry_thrips_t,
    
    # fusarium
    rh_fusarium_threshold_p    = rh_fusarium_threshold_p,
    impact_rh_fusarium_t       = impact_rh_fusarium_t,
    impact_days_wet_fusarium_t = impact_days_wet_fusarium_t,
    Tsopt_fusarium_p           = Tsopt_fusarium_p,
    Tswidth_fusarium_p         = Tswidth_fusarium_p,
    
    # yield losses
    yield_reduction_extreme_rain_t = yield_reduction_extreme_rain_t,
    yield_reduction_hail_t         = yield_reduction_hail_t,
    yield_reduction_fusarium_t     = yield_reduction_fusarium_t,
    yield_reduction_botrytis_t     = yield_reduction_botrytis_t,
    yield_reduction_downy_mildew_t = yield_reduction_downy_mildew_t,
    yield_reduction_thrips_t       = yield_reduction_thrips_t,
    yield_reduction_heat_t         = yield_reduction_heat_t,
    
    # biomass / radiation use
    LAI_emergence_p  = (leaf_area_em_per_plant_p / 10000) * (onions_per_ha_p / 10000),
    LAI_veg_p        = (leaf_area_veg_per_plant_p / 10000) * (onions_per_ha_p / 10000),
    LAI_bulbing_p    = (leaf_area_bl_per_plant_p / 10000) * (onions_per_ha_p / 10000),
    LAI_maturation_p = (leaf_area_mt_per_plant_p / 10000) * (onions_per_ha_p / 10000),
    
    f_T_1_lower_p   = f_T_1_lower_p,
    f_T_1_upper_p   = f_T_1_upper_p,
    f_T_0_lower_p   = f_T_0_lower_p,
    f_T_0_upper_p   = f_T_0_upper_p,
    f_W_1_lower_p   = f_W_1_lower_p,
    f_W_1_upper_p   = f_W_1_upper_p,
    f_W_0.5_p       = f_W_0.5_p,
    
    LUE_onion_p        = LUE_onion_p,
    lec_k_c            = lec_k_c,
    onions_per_ha_p    = onions_per_ha_p,
    HI_onions_t        = HI_onions_t
  )
  
  results <- vector("list", length(weather_scenario_list))
  names(results) <- names(weather_scenario_list)
  
  # 4) Loop over scenarios ----
  for (sc in names(weather_scenario_list)) {
    
    df <- weather_scenario_list[[sc]]
    setDT(df)
    
    # Derived streak variables (wet / dry days) – based on precipitation only
    df$day_consec_wet <- consec_counter(!is.na(df$Prec) & df$Prec > 0)
    df$day_consec_dry <- consec_counter(is.na(df$Prec) | df$Prec == 0)
    
    # HEAT STREAKS ---------------------------------------------------
    df$hot_day        <- df$Tmax >= params$Tmax_heat_threshold_p
    df$day_consec_hot <- consec_counter(df$hot_day)
    
    # Sowing date (fixed parameter)
    sow <- round(sowing_day_p)
    
    # Phase data (em, vg, bl, mt) based on cumulative GDD
    phase_data <- precompute_phase_data(
      df, sow,
      GDD_field_emergence_required_p,
      GDD_vegetative_required_p,
      GDD_bulbing_required_p,
      GDD_maturation_required_p
    )
    
    # Irrigation schedule (mm/day)
    df$Irrig_mm <- compute_irrigation_geisenheim(df, sow, params, phase_data)
    
    # Phase-wise risk scores (stresses do NOT depend on irrigation here)
    all_risks <- compute_all_risks(phase_data, df, params)
    
    # 5) Convert risks into random multipliers (event may or may not occur) ----
    
    # Emergence phase multipliers
    em_multipliers <- c(
      fast_chance_event(all_risks$em["exrain"],   params$yield_reduction_extreme_rain_t),
      fast_chance_event(all_risks$em["hail"],     params$yield_reduction_hail_t),
      fast_chance_event(all_risks$em["fusarium"], params$yield_reduction_fusarium_t)
    )
    
    # Vegetative phase multipliers
    vg_multipliers <- c(
      fast_chance_event(all_risks$vg["heat"],     params$yield_reduction_heat_t),
      fast_chance_event(all_risks$vg["exrain"],   params$yield_reduction_extreme_rain_t),
      fast_chance_event(all_risks$vg["hail"],     params$yield_reduction_hail_t),
      fast_chance_event(all_risks$vg["fusarium"], params$yield_reduction_fusarium_t),
      fast_chance_event(all_risks$vg["mildew"],   params$yield_reduction_downy_mildew_t),
      fast_chance_event(all_risks$vg["thrips"],   params$yield_reduction_thrips_t)
 
    )
    
    # Bulbing phase multipliers
    bl_multipliers <- c(
      fast_chance_event(all_risks$bl["heat"],     params$yield_reduction_heat_t),
      fast_chance_event(all_risks$bl["exrain"],   params$yield_reduction_extreme_rain_t),
      fast_chance_event(all_risks$bl["hail"],     params$yield_reduction_hail_t),
      fast_chance_event(all_risks$bl["fusarium"], params$yield_reduction_fusarium_t),
      fast_chance_event(all_risks$bl["mildew"],   params$yield_reduction_downy_mildew_t),
      fast_chance_event(all_risks$bl["botrytis"], params$yield_reduction_botrytis_t),
      fast_chance_event(all_risks$bl["thrips"],   params$yield_reduction_thrips_t)
    )
    
    # Maturation phase multipliers
    mt_multipliers <- c(
      fast_chance_event(all_risks$mt["heat"],     params$yield_reduction_heat_t),
      fast_chance_event(all_risks$mt["exrain"],   params$yield_reduction_extreme_rain_t),
      fast_chance_event(all_risks$mt["hail"],     params$yield_reduction_hail_t),
      fast_chance_event(all_risks$mt["fusarium"], params$yield_reduction_fusarium_t),
      fast_chance_event(all_risks$mt["mildew"],   params$yield_reduction_downy_mildew_t),
      fast_chance_event(all_risks$mt["botrytis"], params$yield_reduction_botrytis_t),
      fast_chance_event(all_risks$mt["thrips"],   params$yield_reduction_thrips_t)
    )
    
    # Phase-wise biomass multipliers (product of individual stress effects)
    em_bio_multiplier <- prod(em_multipliers)
    vg_bio_multiplier <- prod(vg_multipliers)
    bl_bio_multiplier <- prod(bl_multipliers)
    mt_bio_multiplier <- prod(mt_multipliers)
    
    
    # -----------------------------------------------------------------
    # 6A) RAINFED biomass and yield (no irrigation) -------------------
    # -----------------------------------------------------------------
    
    water_rf <- df$Prec  # only rainfall
    
    biomass_em_pot_rf <- sum(calc_bio_vectorized(
      df$PAR[phase_data$em$idx], params$LAI_emergence_p,
      df$Tavg[phase_data$em$idx], water_rf[phase_data$em$idx],
      params$f_T_1_lower_p, params$f_T_1_upper_p,
      params$f_T_0_lower_p, params$f_T_0_upper_p,
      params$f_W_1_lower_p, params$f_W_1_upper_p,
      params$f_W_0.5_p, params$LUE_onion_p, params$lec_k_c
    ))
    
    biomass_vg_pot_rf <- sum(calc_bio_vectorized(
      df$PAR[phase_data$vg$idx], params$LAI_veg_p,
      df$Tavg[phase_data$vg$idx], water_rf[phase_data$vg$idx],
      params$f_T_1_lower_p, params$f_T_1_upper_p,
      params$f_T_0_lower_p, params$f_T_0_upper_p,
      params$f_W_1_lower_p, params$f_W_1_upper_p,
      params$f_W_0.5_p, params$LUE_onion_p, params$lec_k_c
    ))
    
    biomass_bl_pot_rf <- sum(calc_bio_vectorized(
      df$PAR[phase_data$bl$idx], params$LAI_bulbing_p,
      df$Tavg[phase_data$bl$idx], water_rf[phase_data$bl$idx],
      params$f_T_1_lower_p, params$f_T_1_upper_p,
      params$f_T_0_lower_p, params$f_T_0_upper_p,
      params$f_W_1_lower_p, params$f_W_1_upper_p,
      params$f_W_0.5_p, params$LUE_onion_p, params$lec_k_c
    ))
    
    biomass_mt_pot_rf <- sum(calc_bio_vectorized(
      df$PAR[phase_data$mt$idx], params$LAI_maturation_p,
      df$Tavg[phase_data$mt$idx], water_rf[phase_data$mt$idx],
      params$f_T_1_lower_p, params$f_T_1_upper_p,
      params$f_T_0_lower_p, params$f_T_0_upper_p,
      params$f_W_1_lower_p, params$f_W_1_upper_p,
      params$f_W_0.5_p, params$LUE_onion_p, params$lec_k_c
    ))
    
    biomass_em_rf <- biomass_em_pot_rf * em_bio_multiplier
    biomass_vg_rf <- biomass_vg_pot_rf * vg_bio_multiplier
    biomass_bl_rf <- biomass_bl_pot_rf * bl_bio_multiplier
    biomass_mt_rf <- biomass_mt_pot_rf * mt_bio_multiplier
    
    potential_biomass_g_m2_rf <- biomass_em_pot_rf + biomass_vg_pot_rf +
      biomass_bl_pot_rf + biomass_mt_pot_rf
    realized_biomass_g_m2_rf  <- biomass_em_rf + biomass_vg_rf +
      biomass_bl_rf + biomass_mt_rf
    
    potential_biomass_t_ha_rf <- potential_biomass_g_m2_rf * 0.01
    realized_biomass_t_ha_rf  <- realized_biomass_g_m2_rf  * 0.01
    
    potential_yield_DM_t_ha_rf <- potential_biomass_t_ha_rf * params$HI_onions_t
    final_yield_DM_t_ha_rf     <- realized_biomass_t_ha_rf  * params$HI_onions_t
    
    
    # -----------------------------------------------------------------
    # 6B) IRRIGATED biomass and yield (Prec + Irrig_mm) ---------------
    # -----------------------------------------------------------------
    
    water_ir <- df$Prec + df$Irrig_mm  # rainfall + irrigation
    
    biomass_em_pot_ir <- sum(calc_bio_vectorized(
      df$PAR[phase_data$em$idx], params$LAI_emergence_p,
      df$Tavg[phase_data$em$idx], water_ir[phase_data$em$idx],
      params$f_T_1_lower_p, params$f_T_1_upper_p,
      params$f_T_0_lower_p, params$f_T_0_upper_p,
      params$f_W_1_lower_p, params$f_W_1_upper_p,
      params$f_W_0.5_p, params$LUE_onion_p, params$lec_k_c
    ))
    
    biomass_vg_pot_ir <- sum(calc_bio_vectorized(
      df$PAR[phase_data$vg$idx], params$LAI_veg_p,
      df$Tavg[phase_data$vg$idx], water_ir[phase_data$vg$idx],
      params$f_T_1_lower_p, params$f_T_1_upper_p,
      params$f_T_0_lower_p, params$f_T_0_upper_p,
      params$f_W_1_lower_p, params$f_W_1_upper_p,
      params$f_W_0.5_p, params$LUE_onion_p, params$lec_k_c
    ))
    
    biomass_bl_pot_ir <- sum(calc_bio_vectorized(
      df$PAR[phase_data$bl$idx], params$LAI_bulbing_p,
      df$Tavg[phase_data$bl$idx], water_ir[phase_data$bl$idx],
      params$f_T_1_lower_p, params$f_T_1_upper_p,
      params$f_T_0_lower_p, params$f_T_0_upper_p,
      params$f_W_1_lower_p, params$f_W_1_upper_p,
      params$f_W_0.5_p, params$LUE_onion_p, params$lec_k_c
    ))
    
    biomass_mt_pot_ir <- sum(calc_bio_vectorized(
      df$PAR[phase_data$mt$idx], params$LAI_maturation_p,
      df$Tavg[phase_data$mt$idx], water_ir[phase_data$mt$idx],
      params$f_T_1_lower_p, params$f_T_1_upper_p,
      params$f_T_0_lower_p, params$f_T_0_upper_p,
      params$f_W_1_lower_p, params$f_W_1_upper_p,
      params$f_W_0.5_p, params$LUE_onion_p, params$lec_k_c
    ))
    
    biomass_em_ir <- biomass_em_pot_ir * em_bio_multiplier
    biomass_vg_ir <- biomass_vg_pot_ir * vg_bio_multiplier
    biomass_bl_ir <- biomass_bl_pot_ir * bl_bio_multiplier
    biomass_mt_ir <- biomass_mt_pot_ir * mt_bio_multiplier
    
    potential_biomass_g_m2_ir <- biomass_em_pot_ir + biomass_vg_pot_ir +
      biomass_bl_pot_ir + biomass_mt_pot_ir
    realized_biomass_g_m2_ir  <- biomass_em_ir + biomass_vg_ir +
      biomass_bl_ir + biomass_mt_ir
    
    potential_biomass_t_ha_ir <- potential_biomass_g_m2_ir * 0.01
    realized_biomass_t_ha_ir  <- realized_biomass_g_m2_ir  * 0.01
    
    potential_yield_DM_t_ha_ir <- potential_biomass_t_ha_ir * params$HI_onions_t
    final_yield_DM_t_ha_ir     <- realized_biomass_t_ha_ir  * params$HI_onions_t
    
    
    # HARVEST DAY ------------------------------------------------------
    
    after_sow_idx <- df$yday >= sow
    harvest_yday <- if (any(phase_data$mt$idx)) {
      max(df$yday[phase_data$mt$idx])
    } else {
      max(df$yday[after_sow_idx], na.rm = TRUE)
    }
    
    # 8) Store scenario results ----
    results[[sc]] <- list(
      sowing_yday  = sow,
      harvest_yday = harvest_yday,
      
      # rainfed yields
      raw_yield_per_ha_rainfed   = potential_yield_DM_t_ha_rf,
      final_yield_per_ha_rainfed = final_yield_DM_t_ha_rf,
      
      # irrigated yields
      raw_yield_per_ha_irrigated   = potential_yield_DM_t_ha_ir,
      final_yield_per_ha_irrigated = final_yield_DM_t_ha_ir,
      
      # biomass multipliers
      em_bio_multiplier  = em_bio_multiplier,
      vg_bio_multiplier  = vg_bio_multiplier,
      bl_bio_multiplier  = bl_bio_multiplier,
      mt_bio_multiplier  = mt_bio_multiplier,
      
      # stress multipliers
      m_exrain_em  = em_multipliers[1],
      m_hail_em    = em_multipliers[2],
      m_fus_em     = em_multipliers[3],
      
      m_heat_vg    = vg_multipliers[1],
      m_exrain_vg  = vg_multipliers[2],
      m_hail_vg    = vg_multipliers[3],
      m_fus_vg     = vg_multipliers[4],
      m_mildew_vg  = vg_multipliers[5],
      m_thrips_vg  = vg_multipliers[6],
      
      m_heat_bl    = bl_multipliers[1],
      m_exrain_bl  = bl_multipliers[2],
      m_hail_bl    = bl_multipliers[3],
      m_fus_bl     = bl_multipliers[4],
      m_mildew_bl  = bl_multipliers[5],
      m_botrytis_bl= bl_multipliers[6],
      m_thrips_bl  = bl_multipliers[7],
      
      m_heat_mt    = mt_multipliers[1],
      m_exrain_mt  = mt_multipliers[2],
      m_hail_mt    = mt_multipliers[3],
      m_fus_mt     = mt_multipliers[4],
      m_mildew_mt     = mt_multipliers[5],
      m_botrytis_mt  = mt_multipliers[6],
      m_thrips_mt    = mt_multipliers[7],
      
      # irrigation summary
      total_irrigation_mm = sum(df$Irrig_mm, na.rm = TRUE)
    )
  }
  
  return(results)
})


# Monte Carlo simulation ----

input_variables <- read.csv("input_table_onion.csv", header = TRUE, sep = ";")

onion_mc_simulation <- mcSimulation(
  estimate          = as.estimate(input_variables),
  model_function    = onion_climate_impact,
  numberOfModelRuns = 1000,
  functionSyntax    = "plainNames"
)
