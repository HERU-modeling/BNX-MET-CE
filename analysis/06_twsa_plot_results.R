rm(list = ls()) # to clean the workspace

library(dplyr) # to manipulate data
library(reshape2) # to transform data
library(ggplot2) # for nice looking plots
library(tidyverse)
library(data.table)
library(tidyr)
library(RColorBrewer)
library(Rmisc)
library(grid)
library(gridExtra)
library(lattice)

# Call model setup functions
source("R/input_parameter_functions.R")
source("R/outcome_functions.R")

# Load parameters
source("Analysis/00_load_parameters.R")

load(file = "outputs/dsa/twsa/incremental_twsa_itt_ps_scaled.RData")
load(file = "outputs/dsa/twsa/incremental_twsa_pp_scaled.RData")

# Load DSA parameters
### BNX threshold SA ###
# ITT-PS
df_dsa_twsa_itt_ps_labels <- read.csv(file = "data/dsa/two-way/twsa_itt_ps.csv", header = TRUE, colClasses = c(NA, NA, "NULL", "NULL", "NULL", "NULL"))
# PP
df_dsa_twsa_pp_labels <- read.csv(file = "data/dsa/two-way/twsa_pp.csv", header = TRUE, colClasses = c(NA, NA, "NULL", "NULL", "NULL", "NULL"))

# Combine outcomes into data frame
df_incremental_temp_itt_ps <- cbind(df_dsa_twsa_itt_ps_labels, df_incremental_twsa_itt_ps_scaled)
df_incremental_temp_pp <- cbind(df_dsa_twsa_pp_labels, df_incremental_twsa_pp_scaled)

df_twsa_itt_ps <- df_incremental_temp_itt_ps %>%
  select("perc_improvement_tx", "perc_improvement_death", "n_inc_qalys_adj_2020_scaled", "n_inc_odf_adj_2020_scaled", "n_inc_odn_adj_2020_scaled")
df_twsa_pp <- df_incremental_temp_pp %>%
  select("perc_improvement_tx", "perc_improvement_death", "n_inc_qalys_adj_2020_scaled", "n_inc_odf_adj_2020_scaled", "n_inc_odn_adj_2020_scaled")

# Cohort scaling factor
# n_pop_cohort <- l_params_bnx_itt$n_pop_oat

# Scale outputs
df_twsa_itt_ps_scaled <- df_twsa_itt_ps %>%
  as_tibble() %>%
  mutate(
    n_inc_qalys_adj_total_max_scaled = n_inc_qalys_adj_2020_scaled,
    n_inc_odf_adj_max_scaled = n_inc_odf_adj_2020_scaled,
    n_inc_odn_adj_max_scaled = n_inc_odn_adj_2020_scaled
  )

df_twsa_pp_scaled <- df_twsa_pp %>%
  as_tibble() %>%
  mutate(
    n_inc_qalys_adj_total_max_scaled = n_inc_qalys_adj_2020_scaled,
    n_inc_odf_adj_max_scaled = n_inc_odf_adj_2020_scaled,
    n_inc_odn_adj_max_scaled = n_inc_odn_adj_2020_scaled
  )

plot_twsa_itt_ps_ly <- ggplot(df_twsa_itt_ps_scaled, aes(x = perc_improvement_tx, y = perc_improvement_death, fill = n_inc_qalys_adj_total_max_scaled)) +
  theme_bw() +
  geom_tile() +
  scale_fill_gradientn(
    colours = c("#FF0000", "#FFFFCC", "#075AFF"),
    limits = c(-4000, 4000),
    breaks = c(-4000, 0, 4000),
    labels = c("Favours methadone", "Scenarios equal", "Favours BNX"),
    guide = guide_colourbar(
      direction = "horizontal",
      title = NULL
    )
  ) +
  geom_text(aes(label = round(n_inc_qalys_adj_total_max_scaled, 0)), color = "black", size = 1.8, fontface = "bold") +
  scale_x_continuous(breaks = c(-.40, -.30, -.20, -.10, 0, .10, .20, .30, .40), labels = c("+40%", "+30%", "+20%", "+10%", "Base", "-10%", "-20%", "-30%", "-40%")) +
  scale_y_continuous(breaks = c(-.40, -.30, -.20, -.10, 0, .10, .20, .30, .40), labels = c("+40%", "+30%", "+20%", "+10%", "Base", "-10%", "-20%", "-30%", "-40%")) +
  labs(
    x = "Change in hazard ratio on risk of treatment discontinuation (BNX vs. methadone)",
    y = "Change in hazard ratio on mortality risk in treatment (BNX vs. methadone)"
  ) +
  theme(
    axis.title.x = element_text(size = 8, hjust = 0.5),
    axis.title.y = element_text(size = 8, hjust = 0.5)
  ) +
  theme(
    axis.text = element_text(size = 7, color = "black")
  ) +
  theme(plot.margin = unit(
    c(0, 0.35, 0, 0.35),
    "cm"
  )) +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(1.5, "cm"),
    legend.text = element_text(size = 8)
  ) +
  coord_fixed()

ggsave(plot_twsa_itt_ps_ly,
  filename = "plots/dsa/twsa/twsa_itt_ps_ly.png",
  width = 6, height = 6, dpi = 600
)
