VIP_table <- function (plsrResults, input_table = NULL, cut_off_line = 1, 
                       threshold = 0.8, x_axis_name = "Variable Importance in Projection", 
                       y_axis_name = NULL, legend_name = "Coefficient", 
                       legend_labels = c("Positive", "Negative"), 
                       pos_color = "cadetblue", neg_color = "firebrick", 
                       base_size = 11, ...) 
{
  # allow class c("mvr", "plsr", ...)
  assertthat::assert_that(inherits(plsrResults, "mvr"),
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
  
  if (plsrResults$ncomp == 1) {
    vipResult <- VIP(plsrResults)
  } else {
    vipResult <- VIP(plsrResults)["Comp 1", ]
  }
  
  coef <- plsrResults$coefficients[, , 1]
  pls_outputs <- data.frame(
    Variable    = names(vipResult),
    VIP         = vipResult,
    Coefficient = coef
  )
  rownames(pls_outputs) <- NULL
  
  if (!is.null(input_table)) {
    combined_table <- dplyr::left_join(
      pls_outputs, input_table,
      by = c(Variable = "variable")
    )
  } else {
    combined_table <- pls_outputs
  }
  
  filtered_table <- dplyr::filter(combined_table, VIP > threshold)
  
  return(list(VIP_table_results = filtered_table))
}

labels <- c(
  # Keimphase (em)
  "m_drought_em"  = "Dürre (Keimphase)",
  "m_exrain_em"   = "Starkregen (Keimphase)",
  "m_hail_em"     = "Hagel (Keimphase)",
  "m_fus_em"      = "Fusarium (Keimphase)",
  
  # Vegetative Phase (vg)
  "m_drought_vg"      = "Dürre (Vegetative Phase)",
  "m_exrain_vg"       = "Starkregen (Vegetative Phase)",
  "m_hail_vg"         = "Hagel (Vegetative Phase)",
  "m_fus_vg"          = "Fusarium (Vegetative Phase)",
  "m_mildew_vg"       = "Falscher Mehltau (Vegetative Phase)",
  "m_thrips_vg"       = "Thripse (Vegetative Phase)",
  
  # Zwiebelbildung (bl)
  "m_drought_bl"      = "Dürre (Zwiebelbildung)",
  "m_exrain_bl"       = "Starkregen (Zwiebelbildung)",
  "m_hail_bl"         = "Hagel (Zwiebelbildung)",
  "m_botrytis_bl"     = "Botrytis (Zwiebelbildung)",
  "m_fus_bl"          = "Fusarium (Zwiebelbildung)",
  "m_mildew_bl"       = "Falscher Mehltau (Zwiebelbildung)",
  
  # Reifephase (mt)
  "m_drought_mt"      = "Dürre (Reifephase)",
  "m_exrain_mt"       = "Starkregen (Reifephase)",
  "m_hail_mt"         = "Hagel (Reifephase)",
  "m_botrytis_mt"     = "Botrytis (Reifephase)",
  "m_fus_mt"          = "Fusarium (Reifephase)",
  "m_mildew_mt"       = "Falscher Mehltau (Reifephase)"
)

run_plsr_for_scenario <- function(sim_list, scenario = "historical") {
  # ---- 1. Response: final yield ----
  y <- sim_list$y[[paste0(scenario, ".final_yield_per_ha")]]
  
  # ---- 2. Predictors: all m_* multipliers for this scenario ----
  # e.g. historical.m_drought_em, historical.m_exrain_vg, ...
  x <- sim_list$y %>%
    dplyr::select(starts_with(paste0(scenario, ".m_"))) %>%
    # convert multiplier -> yield reduction
    dplyr::mutate(across(everything(), ~ 1 - .))
  
  # ---- 3. Build minimal mcSimulation object for plsr.mcSimulation ----
  sim_obj <- list(
    y = data.frame(y = y),
    x = as.data.frame(x)
  )
  class(sim_obj) <- c("mcSimulation", "list")
  
  # ---- 4. PLSR ----
  pls_res <- plsr.mcSimulation(
    object     = sim_obj,
    resultName = "y",
    ncomp      = 1
  )
  
  # ---- 5. VIP table (no threshold, keep all) ----
  vip_tab <- VIP_table(pls_res, threshold = 0)$VIP_table_results
  
  # ---- 6. Clean variable names: remove 'scenario.' prefix ----
  vip_tab$CleanVar <- sub(paste0("^", scenario, "\\."), "", vip_tab$Variable)
  
  # ---- 7. Apply German labels where available ----
  vip_tab$Variable <- dplyr::recode(
    vip_tab$CleanVar,
    !!!labels,
    .default = vip_tab$CleanVar  # fallback: raw name
  )
  
  # ---- 8. Plot ----
  plot <- ggplot(vip_tab, aes(x = reorder(Variable, VIP), y = VIP)) +
    geom_col(fill = "firebrick") +
    coord_flip() +
    labs(
      title = paste("VIP Plot – Ertragsmindernde Faktoren", scenario),
      
      y     = "Variable Importance in Projection (VIP)"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank())
  
  return(list(result = pls_res, plot = plot, vip_table = vip_tab))
}

res_hist   <- run_plsr_for_scenario(onion_mc_simulation, "historical")
res_126    <- run_plsr_for_scenario(onion_mc_simulation, "ssp126")
res_245    <- run_plsr_for_scenario(onion_mc_simulation, "ssp245")
res_370    <- run_plsr_for_scenario(onion_mc_simulation, "ssp370")
res_585    <- run_plsr_for_scenario(onion_mc_simulation, "ssp585")

# install.packages("patchwork")
library(patchwork)

res_hist <- run_plsr_for_scenario(onion_mc_simulation, "historical")
res_126  <- run_plsr_for_scenario(onion_mc_simulation, "ssp126")
res_245  <- run_plsr_for_scenario(onion_mc_simulation, "ssp245")
res_370  <- run_plsr_for_scenario(onion_mc_simulation, "ssp370")
res_585  <- run_plsr_for_scenario(onion_mc_simulation, "ssp585")

p_hist <- res_hist$plot
p_126  <- res_126$plot
p_245  <- res_245$plot
p_370  <- res_370$plot
p_585  <- res_585$plot

# 2 rows x 3 columns (last one empty)
combined_plot <- (p_hist | p_126 | p_245) /
  (p_370 | p_585 | patchwork::plot_spacer())

combined_plot

  