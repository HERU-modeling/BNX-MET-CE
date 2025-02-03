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

# Run model
l_out_markov_bnx_itt_ps <- markov_model(l_params_all = l_params_bnx_itt, err_stop = FALSE, verbose = TRUE, checks = FALSE, time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only")
l_out_markov_met_itt_ps <- markov_model(l_params_all = l_params_met_itt, err_stop = FALSE, verbose = TRUE, checks = FALSE, time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only")

#### Create plots ####
l_trace_bnx_itt_ps <- trace_plots(outcomes = l_out_markov_bnx_itt_ps, analysis = "full")
l_trace_met_itt_ps <- trace_plots(outcomes = l_out_markov_met_itt_ps, analysis = "full")

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
