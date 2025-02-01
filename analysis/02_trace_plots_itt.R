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

l_params_bnx_itt <- update_param_list(l_params_all = l_params_bnx_itt, params_updated = l_imis_output_itt$v_calib_post_mean)
l_params_met_itt <- update_param_list(l_params_all = l_params_met_itt, params_updated = l_imis_output_itt$v_calib_post_mean)
l_params_cali_itt <- update_param_list(l_params_all = l_params_cali_itt, params_updated = l_imis_output_itt$v_calib_post_mean)

# Run model
l_out_markov_bnx_itt_ps <- markov_model(l_params_all = l_params_bnx_itt, err_stop = FALSE, verbose = TRUE, checks = FALSE, time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only")
l_out_markov_met_itt_ps <- markov_model(l_params_all = l_params_met_itt, err_stop = FALSE, verbose = TRUE, checks = FALSE, time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only")
l_out_markov_cali_itt <- markov_model(l_params_all = l_params_cali_itt, err_stop = FALSE, verbose = TRUE, checks = FALSE, time_horizon = "cali", ce_est = "cali", analytic_cohort = "cali")

# Output Markov trace
write.csv(l_out_markov_cali_itt$m_m_agg_trace, "outputs/trace/itt/agg_trace_cali_itt.csv", row.names = TRUE)

#### Create plots ####
l_trace_bnx_itt_ps <- trace_plots(outcomes = l_out_markov_bnx_itt_ps, analysis = "full")
l_trace_met_itt_ps <- trace_plots(outcomes = l_out_markov_met_itt_ps, analysis = "full")
l_trace_cali_itt <- trace_plots(outcomes = l_out_markov_cali_itt, analysis = "cali")
# l_trace_cohort_balance_itt <- trace_plots(outcomes = l_out_markov_met_itt_ps, analysis = "cohort_balance")

### Outputs ###
# BNX
plot_full_trace_bnx_itt_ps <- grid.arrange(l_trace_bnx_itt_ps[[1]], l_trace_bnx_itt_ps[[2]], nrow = 2, ncol = 1, heights = 2:1)
ggsave(plot_full_trace_bnx_itt_ps,
  filename = "Plots/Markov Trace/full_trace_bnx_itt_ps.png",
  width = 8, height = 9, dpi = 600
)

# MET
plot_full_trace_met_itt_ps <- grid.arrange(l_trace_met_itt_ps[[1]], l_trace_met_itt_ps[[2]], nrow = 2, ncol = 1, heights = 2:1)
ggsave(plot_full_trace_met_itt_ps,
  filename = "Plots/Markov Trace/full_trace_met_itt_ps.png",
  width = 8, height = 9, dpi = 600
)

# CALI
plot_full_trace_cali_itt <- grid.arrange(l_trace_cali_itt[[1]], l_trace_cali_itt[[2]], nrow = 2, ncol = 1, heights = 2:1)
ggsave(plot_full_trace_cali_itt,
  filename = "Plots/Markov Trace/full_trace_cali_itt.png",
  width = 8, height = 9, dpi = 600
)

# Cohort balance
# plot_cohort_balance_itt <- l_trace_cohort_balance_itt[[1]] # , l_trace_cohort_balance_itt[[2]], nrow = 2, ncol = 1, heights = 2:1)
# ggsave(plot_cohort_balance_itt,
#   filename = "Plots/Markov Trace/cohort_balance_itt.png",
#   width = 8, height = 9, dpi = 600
# )
