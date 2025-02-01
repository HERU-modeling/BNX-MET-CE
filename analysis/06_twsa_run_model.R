rm(list = ls()) # to clean the workspace

library(dplyr)
library(ggplot2)
library(tidyr)
library(parallel)
library(foreach)
library(doParallel)

# Source functions
source("R/input_parameter_functions.R")
source("R/model_setup_functions.R")
source("R/outcome_functions.R")
source("R/mortality_function.R")
source("R/overdose_prob_function.R")

# Load parameters
source("Analysis/00_load_parameters.R")

# Load DSA parameters
# ITT-PS
df_dsa_twsa_itt_ps <- read.csv(file = "data/dsa/two-way/twsa_itt_ps.csv", header = TRUE, colClasses = c("NULL", "NULL", NA, NA, NA, NA))
df_dsa_twsa_itt_ps_labels <- read.csv(file = "data/dsa/two-way/twsa_itt_ps.csv", header = TRUE, colClasses = c(NA, NA, "NULL", "NULL", "NULL", "NULL"))
# pp
df_dsa_twsa_pp <- read.csv(file = "data/dsa/two-way/twsa_pp.csv", header = TRUE, colClasses = c("NULL", "NULL", NA, NA, NA, NA))
df_dsa_twsa_pp_labels <- read.csv(file = "data/dsa/two-way/twsa_pp.csv", header = TRUE, colClasses = c(NA, NA, "NULL", "NULL", "NULL", "NULL"))

# Initialize matrices
# ITT-PS
v_dsa_twsa_itt_ps_names <- colnames(df_dsa_twsa_itt_ps)
v_dsa_twsa_itt_ps_rownames <- rownames(df_dsa_twsa_itt_ps)
# PP
v_dsa_twsa_pp_names <- colnames(df_dsa_twsa_pp)
v_dsa_twsa_pp_rownames <- rownames(df_dsa_twsa_pp)

# ITT-PS
m_dsa_twsa_itt_ps <- array(0,
  dim = c(nrow(df_dsa_twsa_itt_ps), length(df_dsa_twsa_itt_ps)),
  dimnames = list(v_dsa_twsa_itt_ps_rownames, v_dsa_twsa_itt_ps_names)
)
# PP
m_dsa_twsa_pp <- array(0,
  dim = c(nrow(df_dsa_twsa_pp), length(df_dsa_twsa_pp)),
  dimnames = list(v_dsa_twsa_pp_rownames, v_dsa_twsa_pp_names)
)

## Threshold SA ##
# ITT-PS
for (i in seq_len(nrow(df_dsa_twsa_itt_ps))) {
  m_dsa_twsa_itt_ps[i, ] <- unlist(df_dsa_twsa_itt_ps[i, ])
}
# PP
for (i in seq_len(nrow(df_dsa_twsa_pp))) {
  m_dsa_twsa_pp[i, ] <- unlist(df_dsa_twsa_pp[i, ])
}

# Initialize lists
l_outcomes_bnx <- list()
l_outcomes_met <- list()
l_inc_outcomes <- list()

# Set number of cores
n_cores <- detectCores() # number of cores is dictated more by memory than cpu
cl <- makeCluster(n_cores - 1)
clusterExport(cl, c("markov_model", "update_param_list", "outcomes", "inc_outcomes")) # load all functions used on clusters
registerDoParallel(cl)

# Run model
# Intention to treat - PS
Sys.time()
l_twsa_itt_ps <- foreach(i = seq_len(nrow(m_dsa_twsa_itt_ps)), .combine = combine_custom_twsa, .packages = "tidyr") %dopar% {
  l_outcomes_bnx <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = m_dsa_twsa_itt_ps[i, ], time_horizon = "full", ce_est = "dsa", analytic_cohort = "bnx_only", checks = FALSE)
  l_outcomes_met <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = m_dsa_twsa_itt_ps[i, ], time_horizon = "full", ce_est = "dsa", analytic_cohort = "met_only", checks = FALSE)
  l_inc_outcomes <- inc_outcomes(outcomes_comp = l_outcomes_met, outcomes_int = l_outcomes_bnx)

  df_outcomes_bnx_twsa <- l_outcomes_bnx$df_outcomes
  df_outcomes_met_twsa <- l_outcomes_met$df_outcomes
  df_incremental_twsa <- l_inc_outcomes$df_incremental
  df_incremental_twsa_scaled <- l_inc_outcomes$df_incremental_scaled

  return(list(
    df_outcomes_bnx_twsa = df_outcomes_bnx_twsa,
    df_outcomes_met_twsa = df_outcomes_met_twsa,
    df_incremental_twsa = df_incremental_twsa,
    df_incremental_twsa_scaled = df_incremental_twsa_scaled
  ))
}
Sys.time() # Takes ~4 min with 19 cores

# Per-protocol
l_twsa_pp <- foreach(i = seq_len(nrow(m_dsa_twsa_pp)), .combine = combine_custom_twsa, .packages = "tidyr") %dopar% {
  l_outcomes_bnx <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = m_dsa_twsa_pp[i, ], time_horizon = "full", ce_est = "dsa", analytic_cohort = "bnx_only", checks = FALSE)
  l_outcomes_met <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = m_dsa_twsa_pp[i, ], time_horizon = "full", ce_est = "dsa", analytic_cohort = "met_only", checks = FALSE)
  l_inc_outcomes <- inc_outcomes(outcomes_comp = l_outcomes_met, outcomes_int = l_outcomes_bnx)

  df_outcomes_bnx_twsa <- l_outcomes_bnx$df_outcomes
  df_outcomes_met_twsa <- l_outcomes_met$df_outcomes
  df_incremental_twsa <- l_inc_outcomes$df_incremental
  df_incremental_twsa_scaled <- l_inc_outcomes$df_incremental_scaled

  return(list(
    df_outcomes_bnx_twsa = df_outcomes_bnx_twsa,
    df_outcomes_met_twsa = df_outcomes_met_twsa,
    df_incremental_twsa = df_incremental_twsa,
    df_incremental_twsa_scaled = df_incremental_twsa_scaled
  ))
}
stopImplicitCluster()

# ITT-PS
df_outcomes_bnx_twsa_itt_ps <- l_twsa_itt_ps$df_outcomes_bnx_twsa
df_outcomes_met_twsa_itt_ps <- l_twsa_itt_ps$df_outcomes_met_twsa
df_incremental_twsa_itt_ps <- l_twsa_itt_ps$df_incremental_twsa
df_incremental_twsa_itt_ps_scaled <- l_twsa_itt_ps$df_incremental_twsa_scaled
# PP
df_outcomes_bnx_twsa_pp <- l_twsa_pp$df_outcomes_bnx_twsa
df_outcomes_met_twsa_pp <- l_twsa_pp$df_outcomes_met_twsa
df_incremental_twsa_pp <- l_twsa_pp$df_incremental_twsa
df_incremental_twsa_pp_scaled <- l_twsa_pp$df_incremental_twsa_scaled

# Save outputs
## As .RData ##
save(df_outcomes_bnx_twsa_itt_ps,
  file = "outputs/dsa/twsa/outcomes_bnx_twsa_itt_ps.RData"
)
save(df_outcomes_met_twsa_itt_ps,
  file = "outputs/dsa/twsa/outcomes_met_twsa_itt_ps.RData"
)
save(df_incremental_twsa_itt_ps,
  file = "outputs/dsa/twsa/incremental_twsa_itt_ps.RData"
)
save(df_incremental_twsa_itt_ps_scaled,
  file = "outputs/dsa/twsa/incremental_twsa_itt_ps_scaled.RData"
)
save(df_outcomes_bnx_twsa_pp,
  file = "outputs/dsa/twsa/outcomes_bnx_twsa_pp.RData"
)
save(df_outcomes_met_twsa_pp,
  file = "outputs/dsa/twsa/outcomes_met_twsa_pp.RData"
)
save(df_incremental_twsa_pp,
  file = "outputs/dsa/twsa/incremental_twsa_pp.RData"
)
save(df_incremental_twsa_pp_scaled,
  file = "outputs/dsa/twsa/incremental_twsa_pp_scaled.RData"
)
