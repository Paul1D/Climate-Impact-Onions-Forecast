library(dplyr)
library(ggplot2)
library(patchwork)

# -----------------------------
# Funktion: VIP_table
# -----------------------------
VIP_table <- function (plsrResults, input_table = NULL, cut_off_line = 1, 
                       threshold = 0.8, x_axis_name = "Variable Importance in Projection", 
                       y_axis_name = NULL, legend_name = "Coefficient", 
                       legend_labels = c("Positive", "Negative"), 
                       pos_color = "cadetblue", neg_color = "firebrick", 
                       base_size = 11, ...) 
{
  assertthat::assert_that(class(plsrResults) == "mvr", 
                          msg = "plsrResults is not class 'mvr'.")
  
  VIP <- function(object) {
    if (object$method != "oscorespls") 
      stop("Only implemented for oscorespls")
    if (nrow(object$Yloadings) > 1) 
      stop("Only for single-response models")
    SS <- c(object$Yloadings)^2 * colSums(object$scores^2)
    Wnorm2 <- colSums(object$loading.weights^2)
    SSW <- sweep(object$loading.weights^2, 2, SS/Wnorm2, "*")
    sqrt(nrow(SSW) * apply(SSW, 1, cumsum)/cumsum(SS))
  }
  
  if (plsrResults$ncomp == 1) 
    vipResult <- VIP(plsrResults)
  else vipResult <- VIP(plsrResults)["Comp 1", ]
  
  coef <- plsrResults$coefficients[, , 1]
  pls_outputs <- data.frame(Variable = names(vipResult), 
                            VIP = vipResult, 
                            Coefficient = coef)
  rownames(pls_outputs) <- NULL
  
  if (!(is.null(input_table))) 
    combined_table <- dplyr::left_join(pls_outputs, input_table, 
                                       by = c(Variable = "variable"))
  else combined_table <- pls_outputs
  
  filtered_table <- dplyr::filter(combined_table, VIP > threshold)
  
  return(list(VIP_table_results = filtered_table))
}

# -----------------------------
# Labels (deutsche Namen)
# -----------------------------
labels <- c(
  "yield_reduction_seedbed_emergence"       = "Saatbett (Keimphase)",
  "yield_reduction_drought_emergence"       = "Dürre (Keimphase)",
  "yield_reduction_extreme_rain_emergence"  = "Starkregen (Keimphase)",
  "yield_reduction_hail_emergence"          = "Hagel (Keimphase)",
  "yield_reduction_fusarium_emergence"      = "Fusarium (Keimphase)",
  "yield_reduction_onion_fly_emergence"     = "Zwiebelfliege (Keimphase)",
  "yield_reduction_wireworm_emergence"      = "Drahtwurm (Keimphase)",
  
  "yield_reduction_drought_vegetative"      = "Dürre (Vegetative Phase)",
  "yield_reduction_extreme_rain_vegetative" = "Starkregen (Vegetative Phase)",
  "yield_reduction_hail_vegetative"         = "Hagel (Vegetative Phase)",
  "yield_reduction_fusarium_vegetative"     = "Fusarium (Vegetative Phase)",
  "yield_reduction_downy_mildew_vegetative" = "Falscher Mehltau (Vegetative Phase)",
  "yield_reduction_thrips_vegetative"       = "Thripse (Vegetative Phase)",
  "yield_reduction_onion_fly_vegetative"    = "Zwiebelfliege (Vegetative Phase)",
  
  "yield_reduction_drought_bulbing"         = "Dürre (Zwiebelbildung)",
  "yield_reduction_extreme_rain_bulbing"    = "Starkregen (Zwiebelbildung)",
  "yield_reduction_hail_bulbing"            = "Hagel (Zwiebelbildung)",
  "yield_reduction_botrytis_bulbing"        = "Botrytis (Zwiebelbildung)",
  "yield_reduction_fusarium_bulbing"        = "Fusarium (Zwiebelbildung)",
  "yield_reduction_downy_mildew_bulbing"    = "Falscher Mehltau (Zwiebelbildung)",
  "yield_reduction_onion_fly_bulbing"       = "Zwiebelfliege (Zwiebelbildung)",
  
  "yield_reduction_drought_maturation"      = "Dürre (Reifephase)",
  "yield_reduction_extreme_rain_maturation" = "Starkregen (Reifephase)",
  "yield_reduction_hail_maturation"         = "Hagel (Reifephase)",
  "yield_reduction_botrytis_maturation"     = "Botrytis (Reifephase)",
  "yield_reduction_fusarium_maturation"     = "Fusarium (Reifephase)",
  "yield_reduction_downy_mildew_maturation" = "Falscher Mehltau (Reifephase)",
  "yield_reduction_onion_fly_maturation"    = "Zwiebelfliege (Reifephase)"
)

# -----------------------------
# Funktion: PLSR und Plot
# -----------------------------
run_plsr_for_scenario <- function(sim_list, scenario = "historical") {
  
  # Daten vorbereiten
  sim_y <- sim_list$y %>%
    select(starts_with(paste0("weather_", scenario)))
  
  y <- sim_y[[paste0("weather_", scenario, ".final_yield_per_ha")]]
  
  x <- sim_y %>%
    select(contains("yield_reduction")) %>%
    mutate(across(everything(), ~ 1 - .))
  
  sim_obj <- list(y = as.data.frame(y), x = as.data.frame(x))
  class(sim_obj) <- c("mcSimulation", "list")
  
  # PLSR
  pls_res <- plsr.mcSimulation(
    object      = sim_obj,
    variables.x = names(sim_obj$x),
    ncomp       = 1
  )
  
  # VIP-Tabelle
  vip_tab <- VIP_table(pls_res, threshold = 0)$VIP_table_results
  
  # Präfix entfernen (weather_xxx.)
  vip_tab$CleanVar <- sub("^weather_[^.]+\\.", "", vip_tab$Variable)
  
  # Labels anwenden
  vip_tab$Variable <- ifelse(vip_tab$CleanVar %in% names(labels),
                             labels[vip_tab$CleanVar],
                             vip_tab$CleanVar)
  
  # Plot
  plot <- ggplot(vip_tab, aes(x = reorder(Variable, VIP), y = VIP)) +
    geom_col(fill = "firebrick") +
    coord_flip() +
    labs(
      title = paste("VIP Plot Ertragsmindernde Faktoren", scenario),
      x     = "Variablen",
      y     = "Variable Importance in Projection (VIP)"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(list(result = pls_res, plot = plot))
}

# -----------------------------
# Szenarien-Liste
# -----------------------------
scenarios <- c("historical", "ssp126", "ssp245", "ssp370", "ssp585")

# Alle Szenarien laufen lassen
results <- lapply(scenarios, function(scen) {
  out <- run_plsr_for_scenario(onion_mc_simulation, scen)
  out$scenario <- scen
  return(out)
})
names(results) <- scenarios

# -----------------------------
# Einzelplot (Beispiel ssp370)
# -----------------------------
results[["ssp370"]]$plot

# -----------------------------
# Vergleichsplot (alle Szenarien nebeneinander)
# -----------------------------
comparison_plot <- results[["historical"]]$plot +
  results[["ssp126"]]$plot +
  results[["ssp245"]]$plot +
  results[["ssp370"]]$plot +
  results[["ssp585"]]$plot +
  plot_layout(ncol = 2)

comparison_plot
