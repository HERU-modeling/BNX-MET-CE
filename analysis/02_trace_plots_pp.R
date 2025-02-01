rm(list = ls()) # to clean the workspace

library(dplyr) # to manipulate data
library(reshape2) # to transform data
library(ggplot2) # for nice looking plots
library(tidyverse)
library(Rmisc)
library(grid)
library(gridExtra)
library(lattice)

# Call model setup functions
source("R/input_parameter_functions.R")
source("R/model_setup_functions.R")
source("R/plot_functions.R")
source("R/mortality_function.R")
source("R/overdose_prob_function.R")

# Load parameters
source("Analysis/00_load_parameters.R")
# load(file = "outputs/calibration/pp/imis_output_pp.RData")

l_params_bnx_pp <- update_param_list(l_params_all = l_params_bnx_pp, params_updated = l_imis_output_pp$v_calib_post_mean)
l_params_met_pp <- update_param_list(l_params_all = l_params_met_pp, params_updated = l_imis_output_pp$v_calib_post_mean)
l_params_cali_pp <- update_param_list(l_params_all = l_params_cali_pp, params_updated = l_imis_output_pp$v_calib_post_mean)

# Run model
l_out_markov_bnx_pp <- markov_model(l_params_all = l_params_bnx_pp, err_stop = FALSE, verbose = TRUE, checks = FALSE, time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only")
l_out_markov_met_pp <- markov_model(l_params_all = l_params_met_pp, err_stop = FALSE, verbose = TRUE, checks = TRUE, time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only")
l_out_markov_cali_pp <- markov_model(l_params_all = l_params_cali_pp, err_stop = FALSE, verbose = TRUE, checks = FALSE, time_horizon = "cali", ce_est = "cali", analytic_cohort = "cali")

# Output Markov trace
write.csv(l_out_markov_bnx_pp$m_m_agg_trace, "outputs/trace/pp/agg_trace_bnx_pp.csv", row.names = TRUE)
write.csv(l_out_markov_met_pp$m_m_agg_trace, "outputs/trace/pp/agg_trace_met_pp.csv", row.names = TRUE)
write.csv(l_out_markov_cali_pp$m_m_agg_trace, "outputs/trace/pp/agg_trace_cali_pp.csv", row.names = TRUE)

#### Create plots ####
l_trace_bnx_pp <- trace_plots(outcomes = l_out_markov_bnx_pp, cali = FALSE)
l_trace_met_pp <- trace_plots(outcomes = l_out_markov_met_pp, cali = FALSE)
l_trace_cali_pp <- trace_plots(outcomes = l_out_markov_cali_pp, cali = TRUE)

### Outputs ###
# BNX
plot_full_trace_bnx_pp <- grid.arrange(l_trace_bnx_pp[[1]], l_trace_bnx_pp[[2]], nrow = 2, ncol = 1, heights = 2:1)
ggsave(plot_full_trace_bnx_pp,
  filename = "Plots/Markov Trace/full_trace_bnx_pp.png",
  width = 8, height = 9, dpi = 600
)

# MET
plot_full_trace_met_pp <- grid.arrange(l_trace_met_pp[[1]], l_trace_met_pp[[2]], nrow = 2, ncol = 1, heights = 2:1)
ggsave(plot_full_trace_met_pp,
  filename = "Plots/Markov Trace/full_trace_met_pp.png",
  width = 8, height = 9, dpi = 600
)

# CALI
plot_full_trace_cali_pp <- grid.arrange(l_trace_cali_pp[[1]], l_trace_cali_pp[[2]], nrow = 2, ncol = 1, heights = 2:1)
ggsave(plot_full_trace_cali_pp,
  filename = "Plots/Markov Trace/full_trace_cali_pp.png",
  width = 8, height = 9, dpi = 600
)
