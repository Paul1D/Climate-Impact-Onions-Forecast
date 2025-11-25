## Main function: build combined VIP plot for all scenarios -----------

make_onion_vip_combined_plot <- function(sim_list,
                                         scenarios = c("historical",
                                                       "ssp126",
                                                       "ssp245",
                                                       "ssp370",
                                                       "ssp585")) {
  
  # 1) Load packages --------------------------------------------------
  load_if_needed(c(
    "dplyr",
    "tidyr",
    "purrr",
    "ggplot2",
    "patchwork",
    "assertthat",
    "decisionSupport"
  )) 
  
  # 2) VIP table helper -----------------------------------------------
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
  
  # 3) German labels for stressors & phases ---------------------------
  stress_labels <- c(
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
  
  # 4) Scenario-specific PLSR helper ---------------------------------
  run_plsr_for_scenario <- function(sim_list, scenario = "historical") {
    # 1. Response: final yield
    y <- sim_list$y[[paste0(scenario, ".final_yield_per_ha")]]
    
    # 2. Predictors: all m_* multipliers for this scenario
    x <- sim_list$y %>%
      dplyr::select(dplyr::starts_with(paste0(scenario, ".m_"))) %>%
      dplyr::mutate(across(everything(), ~ 1 - .))  # multiplier -> yield reduction
    
    # 3. Build minimal mcSimulation object for plsr.mcSimulation
    sim_obj <- list(
      y = data.frame(y = y),
      x = as.data.frame(x)
    )
    class(sim_obj) <- c("mcSimulation", "list")
    
    # 4. PLSR
    pls_res <- plsr.mcSimulation(
      object     = sim_obj,
      resultName = "y",
      ncomp      = 1
    )
    
    # 5. VIP table (no threshold, keep all)
    vip_tab <- VIP_table(pls_res, threshold = 0)$VIP_table_results
    
    # 6. Clean variable names: remove 'scenario.' prefix
    vip_tab$CleanVar <- sub(paste0("^", scenario, "\\."), "", vip_tab$Variable)
    
    # 7. Apply German labels where available
    vip_tab$Variable <- dplyr::recode(
      vip_tab$CleanVar,
      !!!stress_labels,
      .default = vip_tab$CleanVar
    )
    
    # 8. Plot
    plot <- ggplot(vip_tab, aes(x = reorder(Variable, VIP), y = VIP)) +
      geom_col(fill = "firebrick") +
      coord_flip() +
      labs(
        title = paste("VIP Plot – Ertragsmindernde Faktoren", scenario),
        y     = "Variable Importance in Projection (VIP)"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title  = element_text(hjust = 0.5),
        axis.title.y = element_blank()
      )
    
    return(list(result = pls_res, plot = plot, vip_table = vip_tab))
  }
  
  # 5) Run PLSR for all scenarios & collect plots --------------------
  res_list <- purrr::map(scenarios, ~ run_plsr_for_scenario(sim_list, .x))
  plot_list <- purrr::map(res_list, "plot")
  names(plot_list) <- scenarios
  
  # Expecting 5 scenarios in this order:
  p_hist <- plot_list[["historical"]]
  p_126  <- plot_list[["ssp126"]]
  p_245  <- plot_list[["ssp245"]]
  p_370  <- plot_list[["ssp370"]]
  p_585  <- plot_list[["ssp585"]]
  
  # 6) Build combined plot (2 rows x 3 columns) ----------------------
  combined_plot <- (p_hist | p_126 | p_245) /
    (p_370 | p_585 | patchwork::plot_spacer())
  
  return(combined_plot)
}

## Display Plot ------------------------------------------------------

 combined_plot <- make_onion_vip_combined_plot(onion_mc_simulation)
 combined_plot

