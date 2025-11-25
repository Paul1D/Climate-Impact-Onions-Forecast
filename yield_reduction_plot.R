## Main function: yield reduction heatmap by phase & stressor ----

make_yield_reduction_heatmap <- function(yield_red_long) {
  
  # 1) Load packages ----------------------------------------
  load_if_needed(c(
    "dplyr",
    "tidyr",
    "ggplot2"
  ))
  
  # 2) Ordered labels (codes -> German labels) --------------
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
  
  # Desired scenario order
  scenario_order <- c("historical", "ssp126", "ssp245", "ssp370", "ssp585")
  
  # Desired phase order (top -> bottom)
  phase_order <- c(
    "Keimphase",
    "Vegetative Phase",
    "Zwiebelbildung",
    "Reifephase"
  )
  
  # 3) Build summary table with correct ordering ------------
  yield_red_summary <- yield_red_long %>%
    dplyr::mutate(
      scenario = factor(scenario, levels = scenario_order),
      
      # order stress codes according to stress_labels
      CleanVar = factor(CleanVar, levels = names(stress_labels)),
      
      # German label from named vector
      stress_lab = stress_labels[as.character(CleanVar)],
      
      # Phase = text in parentheses of the German label
      Phase = sub(".*\\((.*)\\).*", "\\1", stress_lab),
      Phase = factor(Phase, levels = phase_order)
    ) %>%
    dplyr::group_by(scenario, CleanVar, stress_lab, Phase) %>%
    dplyr::summarise(
      mean_reduction = mean(reduction, na.rm = TRUE),
      .groups = "drop"
    )
  
  # 4) Heatmap plot (Keimphase at top, grouped by phase) ----
  yield_reduction_plot <- ggplot(
    yield_red_summary,
    aes(
      x    = scenario,
      y    = stress_lab,
      fill = mean_reduction
    )
  ) +
    geom_tile(color = "white") +
    scale_fill_gradient(
      name = "mittlere\nErtragsminderung",
      low  = "white",
      high = "firebrick"
    ) +
    labs(
      title = "Ertragsminderungen nach Phase und Stressfaktor",
      x     = "Szenario",
      y     = "Stressfaktor (Phase)"
    ) +
    facet_grid(
      rows  = vars(Phase),
      scales = "free_y",
      space  = "free_y"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title        = element_text(hjust = 0.5),
      axis.text.y       = element_text(size = 9),
      axis.text.x       = element_text(angle = 45, hjust = 1),
      strip.placement   = "outside",
      strip.text.y.left = element_text(angle = 0, face = "bold")
    )
  
  return(yield_reduction_plot)
}

## Display Plot -------------------------------------------------------

yield_reduction_plot <- make_yield_reduction_heatmap(yield_red_long)
yield_reduction_plot
