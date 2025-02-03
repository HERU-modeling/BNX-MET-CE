rm(list = ls()) # to clean the workspace

library(dplyr)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(parallel)
library(foreach)
library(doParallel)
library(tidyr)
library(gridExtra)

# Call model setup functions
source("R/input_parameter_functions.R")
source("R/generate_psa_parameters.R")
source("R/model_setup_functions.R")
source("R/mortality_function.R")
source("R/overdose_prob_function.R")
source("R/calibration_functions.R")
source("R/cali_scenario_comp.R")

# Load parameters
source("Analysis/00_load_parameters.R")

# Set number of simulations
n_sim <- 10000

# Load PSA inputs
load(file = "outputs/psa/inputs/df_psa_params_itt.RData")
load(file = "outputs/psa/inputs/df_psa_params_pp.RData")

m_calib_post_itt <- l_imis_output_itt$m_calib_post
m_calib_post_pp <- l_imis_output_pp$m_calib_post

#### Load calibration targets ####
df_cali_targets <- read.csv(file = "data/calibration/cali_targets.csv", header = TRUE)

l_cali_targets <- list(
  df_odf = subset(df_cali_targets, target == "fatal overdoses"),
  df_dno = subset(df_cali_targets, target == "non-overdose deaths"),
  df_odn = subset(df_cali_targets, target == "non-fatal overdoses"),
  df_acm = subset(df_cali_targets, target == "all-cause mortality")
)

# Set number of cores
n_cores <- detectCores()
cl <- makeCluster(n_cores - 1)
clusterExport(cl, c("calibration_out", "update_param_list", "markov_model", "overdose", "mort")) # load all functions used on clusters
registerDoParallel(cl)

# Run PSA blockwise
n_runs <- n_sim # n_sim to run entire PSA
n_block_size <- 1000 # size of block for each loop
n_blocks <- if_else(n_runs / n_block_size < 1, 1, n_runs / n_block_size) # to run entire set
n_start <- 0 # set to 0 if running full PSA

###############
### Run PSA ###
###############
Sys.time()
l_met_output_comp_itt <- cali_output_comp(
  scenario = "met_itt",
  l_cali_targets = l_cali_targets
)
Sys.time()
l_bnx_output_comp_itt <- cali_output_comp(
  scenario = "bnx_itt",
  l_cali_targets = l_cali_targets
)
Sys.time()

stopImplicitCluster()

# Save model outputs
save(
  l_met_output_comp_itt,
  l_bnx_output_comp_itt,
  file = "outputs/psa/model_output_target.RData"
)

# Load model outputs
load(file = "outputs/psa/model_output_target.RData")

# Add all-cause mortality to outputs (sum of ODF and DNO)
# MET
m_outcomes_acm_met_itt <- l_met_output_comp_itt$m_outcomes_odf + l_met_output_comp_itt$m_outcomes_dno
# BNX
m_outcomes_acm_bnx_itt <- l_bnx_output_comp_itt$m_outcomes_odf + l_bnx_output_comp_itt$m_outcomes_dno

# Model outputs
# MET (ITT)
# Fatal overdoses
m_outcomes_odf_comb_stats_met_itt <- cbind(
  matrixStats::colQuantiles(l_met_output_comp_itt$m_outcomes_odf_comb,
    probs = c(0.025, 0.5, 0.975)
  ),
  matrixStats::colMeans2(l_met_output_comp_itt$m_outcomes_odf_comb)
)
# Non-overdose deaths
m_outcomes_dno_comb_stats_met_itt <- cbind(
  matrixStats::colQuantiles(l_met_output_comp_itt$m_outcomes_dno_comb,
    probs = c(0.025, 0.5, 0.975)
  ),
  matrixStats::colMeans2(l_met_output_comp_itt$m_outcomes_dno_comb)
)
# All-cause mortality
m_outcomes_acm_stats_met_itt <- cbind(
  matrixStats::colQuantiles(m_outcomes_acm_met_itt,
    probs = c(0.025, 0.5, 0.975)
  ),
  matrixStats::colMeans2(m_outcomes_acm_met_itt)
)
# BNX (ITT)
# Fatal overdoses
m_outcomes_odf_comb_stats_bnx_itt <- cbind(
  matrixStats::colQuantiles(l_bnx_output_comp_itt$m_outcomes_odf_comb,
    probs = c(0.025, 0.5, 0.975)
  ),
  matrixStats::colMeans2(l_bnx_output_comp_itt$m_outcomes_odf_comb)
)
# Non-overdose deaths
m_outcomes_dno_comb_stats_bnx_itt <- cbind(
  matrixStats::colQuantiles(l_bnx_output_comp_itt$m_outcomes_dno_comb,
    probs = c(0.025, 0.5, 0.975)
  ),
  matrixStats::colMeans2(l_bnx_output_comp_itt$m_outcomes_dno_comb)
)
# All-cause mortality
m_outcomes_acm_stats_bnx_itt <- cbind(
  matrixStats::colQuantiles(m_outcomes_acm_bnx_itt,
    probs = c(0.025, 0.5, 0.975)
  ),
  matrixStats::colMeans2(m_outcomes_acm_bnx_itt)
)

m_time <- matrix(l_cali_targets$df_odf$time)
m_pop <- matrix(l_cali_targets$df_odf$pop)

# Combine model outputs with time and population
# MET (ITT)
m_outcomes_odf_comb_fit_met_itt <- cbind(m_outcomes_odf_comb_stats_met_itt, m_time, m_pop)
m_outcomes_acm_comb_fit_met_itt <- cbind(m_outcomes_acm_stats_met_itt, m_time, m_pop)
# BNX (ITT)
m_outcomes_odf_comb_fit_bnx_itt <- cbind(m_outcomes_odf_comb_stats_bnx_itt, m_time, m_pop)
m_outcomes_acm_comb_fit_bnx_itt <- cbind(m_outcomes_acm_stats_bnx_itt, m_time, m_pop)

# Create dataframes
# Fatal overdoses
# MET (ITT)
df_model_targets_odf_fit_met_itt <- m_outcomes_odf_comb_fit_met_itt %>%
  as_tibble() %>%
  setNames(c("ci_low", "median", "ci_high", "pe", "time", "pop")) %>%
  mutate(
    group = "Methadone-only (counterfactual)",
    num = pe * m_pop,
    low = ci_low * m_pop,
    high = ci_high * m_pop
  ) %>%
  select(time, group, num, low, high)
# BNX (ITT)
df_model_targets_odf_fit_bnx_itt <- m_outcomes_odf_comb_fit_bnx_itt %>%
  as_tibble() %>%
  setNames(c("ci_low", "median", "ci_high", "pe", "time", "pop")) %>%
  mutate(
    group = "BNX-only (counterfactual)",
    num = pe * m_pop,
    low = ci_low * m_pop,
    high = ci_high * m_pop
  ) %>%
  select(time, group, num, low, high)

# All-cause mortality
# MET (ITT)
df_model_targets_acm_fit_met_itt <- m_outcomes_acm_comb_fit_met_itt %>%
  as_tibble() %>%
  setNames(c("ci_low", "median", "ci_high", "pe", "time", "pop")) %>%
  mutate(
    group = "Methadone-only (counterfactual)",
    num = pe * m_pop,
    low = ci_low * m_pop,
    high = ci_high * m_pop
  ) %>%
  select(time, group, num, low, high)
# BNX (ITT)
df_model_targets_acm_fit_bnx_itt <- m_outcomes_acm_comb_fit_bnx_itt %>%
  as_tibble() %>%
  setNames(c("ci_low", "median", "ci_high", "pe", "time", "pop")) %>%
  mutate(
    group = "BNX-only (counterfactual)",
    num = pe * m_pop,
    low = ci_low * m_pop,
    high = ci_high * m_pop
  ) %>%
  select(time, group, num, low, high)

# Load targets
# Fatal overdoses
df_targets_odf <- l_cali_targets$df_odf %>%
  as_tibble() %>%
  mutate(
    group = "Observed data",
    num = pe * m_pop,
    low = pe * m_pop,
    high = pe * m_pop
  ) %>%
  select(time, group, num, low, high)
# All-cause mortality
df_targets_acm <- l_cali_targets$df_acm %>%
  as_tibble() %>%
  mutate(
    group = "Observed data",
    num = pe * m_pop,
    low = pe * m_pop,
    high = pe * m_pop
  ) %>%
  select(time, group, num, low, high)

# Combine
df_fit_odf_itt <- bind_rows(df_targets_odf, df_model_targets_odf_fit_met_itt, df_model_targets_odf_fit_bnx_itt)
df_fit_acm_itt <- bind_rows(df_targets_acm, df_model_targets_acm_fit_met_itt, df_model_targets_acm_fit_bnx_itt)

# Plot fit vs. targets
## Fatal overdose
# ITT
p_temp_odf_itt <- ggplot(df_fit_odf_itt) +
  geom_pointrange(
    aes(x = time, y = num, ymin = low, ymax = high, color = group),
    position = position_dodge(width = 25),
    linewidth = .75
  )
plot_fit_odf_itt <- p_temp_odf_itt + labs(title = NULL, x = "Year", y = "Fatal overdoses") +
  scale_color_manual(values = c("#4B92DB", "#0C2340", "#C8102E")) +
  scale_x_continuous(
    breaks = l_cali_targets$df_odf$time,
    labels = c("2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020*")
  ) +
  scale_y_continuous(
    breaks = c(100, 200, 300, 400, 500),
    limits = c(0, 600)
  ) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.02, vjust = -7),
    legend.position = "none",
    legend.title = element_blank(),
    text = element_text(size = 15)
  )

## All-cause mortality
# ITT
p_temp_acm_itt <- ggplot(df_fit_acm_itt) +
  geom_pointrange(
    aes(x = time, y = num, ymin = low, ymax = high, color = group),
    position = position_dodge(width = 25),
    linewidth = .75
  )
plot_fit_acm_itt <- p_temp_acm_itt + labs(title = NULL, x = "Year", y = "All-cause deaths") +
  scale_color_manual(values = c("#4B92DB", "#0C2340", "#C8102E")) +
  scale_x_continuous(
    breaks = l_cali_targets$df_acm$time,
    labels = c("2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020*")
  ) +
  scale_y_continuous(
    breaks = c(100, 200, 300, 400, 500),
    limits = c(0, 600)
  ) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.02, vjust = -7),
    legend.position = "none",
    legend.title = element_blank(),
    text = element_text(size = 15)
  )
# ggsave(plot_fit_acm_itt,
#   filename = "plots/psa/model_output_target_acm_itt.png",
#   width = 6, height = 6, dpi = 350
# )

# Plot for extracting legend only
plot_fit_odf_leg <- p_temp_odf_itt + labs(title = NULL, x = "Year", y = "Fatal overdoses") +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom"
  ) +
  scale_color_manual(values = c("#4B92DB", "#0C2340", "#C8102E")) +
  scale_x_continuous(
    breaks = l_cali_targets$df_odf$time,
    labels = c("2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020*")
  )
# Code to extract legend from plots
g_legend <- function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
mylegend <- g_legend(plot_fit_odf_leg)
# mylegend2 <- g_legend(plot_cali_fit_acm_itt)

# Combine plots (ITT)
plot_fit_comb_itt <- grid.arrange(arrangeGrob(plot_fit_odf_itt, plot_fit_acm_itt, nrow = 2),
  mylegend,
  nrow = 2, heights = c(6, .5), widths = 6
)
ggsave(plot_fit_comb_itt,
  filename = "plots/psa/model_output_target_itt.png",
  width = 6, height = 6
)
