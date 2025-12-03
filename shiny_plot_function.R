############################################################
## Pakete laden
############################################################
load_if_needed <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}

load_if_needed(c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "data.table",
  "decisionSupport",
  "patchwork",
  "stringr",
  "assertthat",
  "tibble"
))

############################################################
## Gemeinsame Definitionen
############################################################

# Szenariocodes und Labels (für VIP + Yield)
scenario_order_codes <- c("historical", "ssp126", "ssp245", "ssp370", "ssp585")
scenario_labels <- c(
  historical = "2020",
  ssp126     = "2075 (SSP 1-2.6)",
  ssp245     = "2075 (SSP 2-4.5)",
  ssp370     = "2075 (SSP 3-7.0)",
  ssp585     = "2075 (SSP 5-8.5)"
)

# Stressfaktoren (Codes -> deutsche Labels)
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

phase_order <- c("Keimphase", "Vegetative Phase", "Zwiebelbildung", "Reifephase")
order_labels <- unname(stress_labels)  # Reihenfolge für Y-Achse

############################################################
## 1) Daten vorbereiten: Multiplikatoren -> Schaden
##    (0 = kein Schaden, >0 = Ertragsminderung)
############################################################

# Kopie anlegen, damit das Original erhalten bleibt
onion_mc_sim_damage <- onion_mc_simulation

# alle Spalten mit ".m_" im Namen sind Stress-Multiplikatoren
stress_cols <- grep("\\.m_", names(onion_mc_sim_damage$y), value = TRUE)

# damage = 1 - multiplikator
onion_mc_sim_damage$y[stress_cols] <-
  1 - onion_mc_sim_damage$y[stress_cols]

# Ab jetzt sind die .m_-Spalten "Schaden"

############################################################
## 2) Hilfsfunktion: Variablen mit zu wenig Variation entfernen
############################################################

filter_stress_vars <- function(df, min_nonzero = 5) {
  keep <- sapply(df, function(x) sum(x != 0, na.rm = TRUE) >= min_nonzero)
  df[, keep, drop = FALSE]
}

############################################################
## 3) VIP-Hilfsfunktion (aus Spargel)
############################################################

VIP_table <- function (plsrResults, input_table = NULL, cut_off_line = 1, 
                       threshold = 0.8, x_axis_name = "Variable Importance in Projection", 
                       y_axis_name = NULL, legend_name = "Coefficient", 
                       legend_labels = c("Positive", "Negative"), 
                       pos_color = "cadetblue", neg_color = "firebrick", 
                       base_size = 11, ...) 
{
  assertthat::assert_that(inherits(plsrResults, "mvr"),
                          msg = "plsrResults ist nicht Klasse 'mvr'.")
  
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
  pls_outputs <- data.frame(
    Variable    = names(vipResult), 
    VIP         = vipResult, 
    Coefficient = coef,
    row.names   = NULL
  )
  
  if (!is.null(input_table)) {
    combined_table <- dplyr::left_join(
      pls_outputs, input_table, by = c(Variable = "variable")
    )
  } else {
    combined_table <- pls_outputs
  }
  
  dplyr::filter(combined_table, VIP > threshold) |>
    list(VIP_table_results = _)
}

############################################################
## 4) PLSR + VIP-Plot für ein Szenario (X = Schaden, mit Filter)
############################################################

run_plsr_for_scenario_onion <- function(sim_list, scen_code,
                                        min_nonzero = 5) {
  
  # Zielvariable: final_yield_per_ha für dieses Szenario
  y_col <- paste0(scen_code, ".final_yield_per_ha")
  if (!y_col %in% names(sim_list$y))
    stop("Spalte ", y_col, " nicht in sim_list$y gefunden.")
  
  y <- data.frame(yield = sim_list$y[[y_col]])
  
  # Prädiktoren: Schaden-Variablen dieses Szenarios (m_*)
  stress_base <- names(stress_labels)
  stress_cols_scen <- paste0(scen_code, ".", stress_base)
  stress_cols_scen <- stress_cols_scen[stress_cols_scen %in% names(sim_list$y)]
  
  x_raw <- sim_list$y[, stress_cols_scen, drop = FALSE]
  
  # Variablen mit zu wenig Variation rauswerfen (z.B. Hagel mit nur 1 Fall)
  x <- filter_stress_vars(x_raw, min_nonzero = min_nonzero)
  if (ncol(x) == 0) {
    stop("Für Szenario ", scen_code,
         " gibt es keine Stressvariablen mit ausreichender Variation.")
  }
  
  # Szenario-Präfix weg, damit Variablen m_drought_em heißen
  names(x) <- sub(paste0("^", scen_code, "\\."), "", names(x))
  
  sim_obj <- list(y = y, x = x)
  class(sim_obj) <- c("mcSimulation", "list")
  
  pls_res <- plsr.mcSimulation(
    object      = sim_obj,
    resultName  = "yield",
    variables.x = names(sim_obj$x),
    ncomp       = 1
  )
  
  vip_tab <- VIP_table(pls_res, threshold = 0.8)$VIP_table_results |>
    as_tibble() |>
    mutate(
      CleanVar   = Variable,
      stress_lab = stress_labels[CleanVar],
      Phase      = sub(".*\\((.*)\\).*", "\\1", stress_lab),
      Phase      = factor(Phase, levels = phase_order),
      coef_sign  = case_when(
        Coefficient >  0 ~ "positive",
        Coefficient <  0 ~ "negative",
        TRUE             ~ "zero"
      ),
      scenario   = factor(scen_code, levels = scenario_order_codes)
    )
  
  # Y-Achse: alle Stressfaktoren (auch wenn in diesem Szenario kein Punkt)
  vip_tab$stress_lab <- factor(
    vip_tab$stress_lab,
    levels  = order_labels,
    ordered = TRUE
  )
  
  # Dynamische VIP-Breaks für die Legende (3 schöne Werte im Bereich der Daten)
  vip_range  <- range(vip_tab$VIP, na.rm = TRUE)
  vip_breaks <- pretty(vip_range, n = 3)
  # Falls pretty etwas außerhalb des Datenbereichs schlägt, leicht beschneiden:
  vip_breaks <- vip_breaks[vip_breaks >= vip_range[1] & vip_breaks <= vip_range[2]]
  # Sicherheitsfallback, falls aus irgendeinem Grund nichts übrig bleibt
  if (length(vip_breaks) == 0) {
    vip_breaks <- vip_range
  }
  
  p <- ggplot(vip_tab, aes(x = scenario, y = stress_lab)) +
    geom_point(aes(size = VIP, color = coef_sign), shape = 16) +
    
    # VIP = Wichtigkeit der Variable
    # -> größere Range und dynamische Breaks für bessere Unterscheidbarkeit
    scale_size_continuous(
      range  = c(3, 10),                      # minimale & maximale Punktgröße
      breaks = vip_breaks,
      labels = round(vip_breaks, 2),
      name   = "VIP\n(Variable Importance in Projection)"
    ) +
    
    # Farben nur zur visuellen Unterscheidung – ohne eigene Legende
    scale_color_manual(
      values = c(
        negative = "firebrick",
        zero     = "grey70",
        positive = "cadetblue"
      ),
      guide = "none"
    ) +
    
    # fixe Y-Achse mit allen Stressfaktoren
    scale_y_discrete(
      limits = order_labels,
      drop   = FALSE,
      expand = expansion(mult = c(0.02, 0.06))
    ) +
    
    labs(
      title = scenario_labels[scen_code],   # Titel wie im Yield-Plot
      x     = "Szenario",
      y     = "Stressfaktor (Phase)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title         = element_text(hjust = 0.5),
      axis.text.y        = element_text(hjust = 0),
      legend.position    = "right",
      legend.title.align = 0,
      legend.box         = "vertical",
      legend.key.height  = grid::unit(0.6, "cm"),
      legend.key.width   = grid::unit(0.6, "cm"),
      legend.text        = element_text(size = 9)
    ) +
    # Legendenpunkte im VIP-Teil in Rot, mit klar sichtbarer Größe
    guides(
      size = guide_legend(
        override.aes = list(
          color = "firebrick",
          fill  = "firebrick"
        ),
        order         = 1,
        title.position = "top"
      )
    )
  
  list(result = pls_res, plot = p)
}


############################################################
## 5) VIP-Vergleichsplot über alle Szenarien
############################################################

VIP_plot_onion <- function(sim_results, min_nonzero = 5) {
  scen_codes <- scenario_order_codes
  
  results <- lapply(scen_codes, function(sc) 
    run_plsr_for_scenario_onion(sim_results, sc, min_nonzero = min_nonzero)
  )
  names(results) <- scen_codes
  
  no_x <- theme(axis.title.x = element_blank(),
                axis.text.x  = element_blank(),
                axis.ticks.x = element_blank())
  no_y <- theme(axis.title.y = element_blank(),
                axis.text.y  = element_blank(),
                axis.ticks.y = element_blank())
  
  comparison_plot <-
    # erster Plot: Y-Beschriftungen behalten, aber LEGENDE aus
    (results[[1]]$plot + no_x + theme(legend.position = "none")) +
    (results[[2]]$plot + no_x + no_y + theme(legend.position = "none")) +
    (results[[3]]$plot + no_x + no_y + theme(legend.position = "none")) +
    (results[[4]]$plot + no_x + no_y + theme(legend.position = "none")) +
    # letzter Plot: ohne Y-Beschriftung, aber mit Legende
    (results[[5]]$plot + no_x + no_y) +
    plot_layout(ncol = 5)
  
  comparison_plot
}

############################################################
## 6) Ertrags-Boxplots (Yield-Plot)
############################################################

onion_plots <- function(onion_mc_sim_damage) {
  
  results_yield_onion <- rbind(
    data.frame(
      Ertrag               = onion_mc_sim_damage$y$historical.raw_yield_per_ha,
      vermarktbarer_Ertrag = onion_mc_sim_damage$y$historical.final_yield_per_ha,
      id                   = seq_len(nrow(onion_mc_sim_damage$y)),
      scenario_code        = "historical"
    ),
    data.frame(
      Ertrag               = onion_mc_sim_damage$y$ssp126.raw_yield_per_ha,
      vermarktbarer_Ertrag = onion_mc_sim_damage$y$ssp126.final_yield_per_ha,
      id                   = seq_len(nrow(onion_mc_sim_damage$y)),
      scenario_code        = "ssp126"
    ),
    data.frame(
      Ertrag               = onion_mc_sim_damage$y$ssp245.raw_yield_per_ha,
      vermarktbarer_Ertrag = onion_mc_sim_damage$y$ssp245.final_yield_per_ha,
      id                   = seq_len(nrow(onion_mc_sim_damage$y)),
      scenario_code        = "ssp245"
    ),
    data.frame(
      Ertrag               = onion_mc_sim_damage$y$ssp370.raw_yield_per_ha,
      vermarktbarer_Ertrag = onion_mc_sim_damage$y$ssp370.final_yield_per_ha,
      id                   = seq_len(nrow(onion_mc_sim_damage$y)),
      scenario_code        = "ssp370"
    ),
    data.frame(
      Ertrag               = onion_mc_sim_damage$y$ssp585.raw_yield_per_ha,
      vermarktbarer_Ertrag = onion_mc_sim_damage$y$ssp585.final_yield_per_ha,
      id                   = seq_len(nrow(onion_mc_sim_damage$y)),
      scenario_code        = "ssp585"
    )
  )
  
  # Szenario-Faktor mit gleichen Labels wie im VIP-Titel
  results_yield_onion$scenario <- factor(
    results_yield_onion$scenario_code,
    levels = scenario_order_codes,
    labels = scenario_labels[scenario_order_codes]
  )
  
  names(results_yield_onion)[1:3] <- c("Ertrag", "vermarktbarer_Ertrag", "id")
  
  results_yield_onion_longer <- results_yield_onion |>
    tidyr::pivot_longer(cols = c(Ertrag, vermarktbarer_Ertrag),
                        names_to = "name", values_to = "value") |>
    mutate(
      name   = factor(name,
                      levels = c("Ertrag", "vermarktbarer_Ertrag"),
                      labels = c("Potentieller Ertrag ohne Schäden",
                                 "Tatsächlicher Ertrag")),
      period = ifelse(scenario == "2020", "Historical", "Future")
    )
  
  loss_df <- results_yield_onion_longer |>
    group_by(scenario, name) |>
    summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = name, values_from = mean_value) |>
    mutate(
      loss_percent =
        (`Potentieller Ertrag ohne Schäden` - `Tatsächlicher Ertrag`) /
        `Potentieller Ertrag ohne Schäden` * 100
    )
  
  max_y <- max(results_yield_onion_longer$value, na.rm = TRUE)
  
  # Baseline & Text für die Legende
  baseline_hist  <- 45.7
  baseline_label <- paste0(
    "Durchschnittsertrag Regierungsbezirk\nKöln (", baseline_hist, " t/ha)"
  )
  
  onion_yield_plot <- ggplot(results_yield_onion_longer,
                             aes(x = scenario, y = value, fill = name)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +
    
    # Prozentzahlen über den Boxplots
    geom_text(
      data = loss_df,
      aes(x = scenario, y = max_y * 1.03, label = paste0(round(loss_percent), "%")),
      inherit.aes = FALSE,
      size = 4,
      fontface = "bold"
    ) +
    
    # Baseline als eigene Legenden-Zeile
    geom_hline(
      data = data.frame(yint = baseline_hist),
      aes(yintercept = yint, linetype = "baseline"),
      linewidth = 0.5,
      colour = "black",
      inherit.aes = FALSE
    ) +
    scale_linetype_manual(
      name   = NULL,
      values = c(baseline = "dashed"),
      labels = c(baseline = baseline_label)
    ) +
    
    theme(
      legend.title = element_blank(),
      legend.position = "right",
      strip.background = element_rect(fill = "lightgrey"),
      strip.text = element_text(size = 12, face = "bold")
    ) +
    scale_x_discrete(name = "Klimaszenario") +
    scale_y_continuous(name = "Zwiebel-Ertrag [t/ha]") +
    
    # Reihenfolge der Legenden: erst Linie, dann Füllfarben
    guides(
      linetype = guide_legend(order = 1),
      fill     = guide_legend(order = 2)
    )
  
  onion_vip_plot <- VIP_plot_onion(onion_mc_sim_damage)
  
  list(
    vip_plot   = onion_vip_plot,
    yield_plot = onion_yield_plot
  )
}


############################################################
## 7) Heatmap: mittlere Ertragsminderung (optional)
############################################################

build_stress_damage_long <- function(onion_mc_sim_damage) {
  onion_mc_sim_damage$y |>
    dplyr::select(dplyr::matches("^(historical|ssp126|ssp245|ssp370|ssp585)\\.m_")) |>
    tidyr::pivot_longer(cols = dplyr::everything(),
                        names_to = "full_name",
                        values_to = "damage") |>
    tidyr::separate(full_name, into = c("scenario", "CleanVar"),
                    sep = "\\.", extra = "merge") |>
    dplyr::mutate(
      scenario = factor(scenario, levels = scenario_order_codes,
                        labels = scenario_labels[scenario_order_codes])
    )
}

make_yield_reduction_heatmap <- function(yield_damage_long) {
  
  # 1) Zusammenfassung: mittlere Ertragsminderung in PROZENT ----------------
  # damage wird als Anteil (0–1) angenommen -> *100 für Prozent
  
  yield_red_summary <- yield_damage_long |>
    dplyr::mutate(
      CleanVar   = factor(CleanVar, levels = names(stress_labels)),
      stress_lab = stress_labels[as.character(CleanVar)],
      # Phase aus dem Klammerinhalt des deutschen Labels ziehen
      Phase      = sub(".*\\((.*)\\).*", "\\1", stress_lab),
      Phase      = factor(Phase, levels = phase_order)
    ) |>
    dplyr::group_by(scenario, CleanVar, stress_lab, Phase) |>
    dplyr::summarise(
      mean_reduction = mean(damage, na.rm = TRUE) * 100,  # jetzt in Prozent!
      .groups        = "drop"
    )
  
  # 2) Heatmap-Plot ---------------------------------------------------------
  
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
      name   = "mittlere\nErtragsminderung [%]",
      low    = "white",
      high   = "firebrick"
      # Optional: eigene Breaks setzen, z.B.:
      # breaks = c(0, 5, 10, 20),
      # limits = c(0, 20)
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


############################################################
## VERWENDUNG
############################################################
## VIP + Ertrags-Boxplots:
  res <- onion_plots(onion_mc_sim_damage)
   res$vip_plot
   res$yield_plot
##
## Heatmap (optional):
   dmg_long  <- build_stress_damage_long(onion_mc_sim_damage)
   heat_plot <- make_yield_reduction_heatmap(dmg_long)
   heat_plot
############################################################
