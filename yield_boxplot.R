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
  "decisionSupport",
  "tidyr",
  "data.table"
))


# -----------------------------
# Prepare onion yield data (same structure as asparagus)
# -----------------------------

results_yield_onion <- rbind(
  data.frame(
    Ertrag = onion_mc_simulation$y$historical.raw_yield_per_ha,
    vermarktbarer_Ertrag = onion_mc_simulation$y$historical.final_yield_per_ha,
    id = seq_along(onion_mc_simulation$y$historical.raw_yield_per_ha),
    scenario = "Historical Weather"
  ),
  data.frame(
    Ertrag = onion_mc_simulation$y$ssp126.raw_yield_per_ha,
    vermarktbarer_Ertrag = onion_mc_simulation$y$ssp126.final_yield_per_ha,
    id = seq_along(onion_mc_simulation$y$ssp126.raw_yield_per_ha),
    scenario = "SSP 1-2.6"
  ),
  data.frame(
    Ertrag = onion_mc_simulation$y$ssp245.raw_yield_per_ha,
    vermarktbarer_Ertrag = onion_mc_simulation$y$ssp245.final_yield_per_ha,
    id = seq_along(onion_mc_simulation$y$ssp245.raw_yield_per_ha),
    scenario = "SSP 2-4.5"
  ),
  data.frame(
    Ertrag = onion_mc_simulation$y$ssp370.raw_yield_per_ha,
    vermarktbarer_Ertrag = onion_mc_simulation$y$ssp370.final_yield_per_ha,
    id = seq_along(onion_mc_simulation$y$ssp370.raw_yield_per_ha),
    scenario = "SSP 3-7.0"
  ),
  data.frame(
    Ertrag = onion_mc_simulation$y$ssp585.raw_yield_per_ha,
    vermarktbarer_Ertrag = onion_mc_simulation$y$ssp585.final_yield_per_ha,
    id = seq_along(onion_mc_simulation$y$ssp585.raw_yield_per_ha),
    scenario = "SSP 5-8.5"
  )
)

# Rename like asparagus code
names(results_yield_onion) <- c("Ertrag", "vermarktbarer_Ertrag", "id", "scenario")

# Long table: two yield types
results_yield_onion_longer <- results_yield_onion %>%
  pivot_longer(cols = c(Ertrag, vermarktbarer_Ertrag)) 

results_yield_onion_longer$name <- factor(
  results_yield_onion_longer$name,
  levels = c("Ertrag", "vermarktbarer_Ertrag")
)

# Historical vs future grouping (optional)
results_yield_onion_longer <- results_yield_onion_longer %>%
  mutate(period = ifelse(scenario == "Historical Weather", "Historical", "Future"))

results_yield_onion_longer$period   <- as.factor(results_yield_onion_longer$period)
results_yield_onion_longer$scenario <- as.factor(results_yield_onion_longer$scenario)

# -----------------------------
# Summary stats for % labels and reference line
# -----------------------------

summary_df_onion <- results_yield_onion_longer %>%
  group_by(scenario, name) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    q25 = quantile(value, 0.25, na.rm = TRUE),
    q75 = quantile(value, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from  = name,
    values_from = c(mean_value, q25, q75)
  ) %>%
  mutate(
    percent = (mean_value_vermarktbarer_Ertrag / mean_value_Ertrag) * 100,
    IQR_marketable_yield = q75_vermarktbarer_Ertrag - q25_vermarktbarer_Ertrag
  )

# Baseline reference = historical final yield mean
baseline_hist <- summary_df_onion %>%
  filter(scenario == "Historical Weather") %>%
  pull(mean_value_vermarktbarer_Ertrag)

# -----------------------------
# Plot (styled like asparagus)
# -----------------------------

onion_yield_plot <- ggplot(
  results_yield_onion_longer,
  aes(x = scenario, y = value, fill = name)
) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_text(
    data = summary_df_onion,
    aes(
      x = scenario,
      y = max(results_yield_onion_longer$value, na.rm = TRUE) * 1.02,
      label = paste0(round(percent), "%")
    ),
    inherit.aes = FALSE,
    size = 4
  ) +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    strip.background = element_rect(fill = "lightgrey"),
    strip.text = element_text(size = 12, face = "bold") 
  ) +
  scale_x_discrete(name = "Klimaszenario") +
  scale_y_continuous(name = "Zwiebel-Ertrag [t/ha]") +
  geom_hline(yintercept = 47, linetype = "dashed") +
  annotate(
    geom = "text",
    x = -Inf,
    y = baseline_hist * 0.97,
    label = paste0(round(47, 1), " t/ha\nhistorischer Ertrag"),
    color = "black",
    fontface = "plain",
    hjust = 0,
    vjust = 1,
    angle = 0,
    size = 3
  )

onion_yield_plot
