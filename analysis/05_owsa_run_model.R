rm(list = ls()) # to clean the workspace

library(dplyr) # to manipulate data
library(ggplot2) # for nice looking plots
library(tidyr)
library(tibble)

# Source functions
source("R/input_parameter_functions.R")
source("R/model_setup_functions.R")
source("R/outcome_functions.R")
source("R/mortality_function.R")
source("R/overdose_prob_function.R")
source("R/owsa_read_data_function.R")

# Load parameters
source("Analysis/00_load_parameters.R")

# Load DSA parameters
# DSA data
df_dsa_owsa <- read.csv(file = "data/dsa/owsa/owsa_ranges.csv", row.names = 1, header = TRUE)
# 25% Inc; 75% Prev
v_init_dist_met_25_inc <- read.csv(file = "data/dsa/owsa/init_dist_met_25_inc.csv", row.names = 1, header = TRUE)
v_init_dist_bnx_25_inc <- read.csv(file = "data/dsa/owsa/init_dist_bnx_25_inc.csv", row.names = 1, header = TRUE)
# 50% Inc; 50% Prev
v_init_dist_met_50_inc <- read.csv(file = "data/dsa/owsa/init_dist_met_50_inc.csv", row.names = 1, header = TRUE)
v_init_dist_bnx_50_inc <- read.csv(file = "data/dsa/owsa/init_dist_bnx_50_inc.csv", row.names = 1, header = TRUE)
# 75% Inc; 25% Prev
v_init_dist_met_75_inc <- read.csv(file = "data/dsa/owsa/init_dist_met_75_inc.csv", row.names = 1, header = TRUE)
v_init_dist_bnx_75_inc <- read.csv(file = "data/dsa/owsa/init_dist_bnx_75_inc.csv", row.names = 1, header = TRUE)

# Read data into vector
l_owsa_met_itt_ps <- owsa_params(
  df_dsa_owsa = df_dsa_owsa,
  v_init_dist_met_25_inc = v_init_dist_met_25_inc,
  v_init_dist_met_50_inc = v_init_dist_met_50_inc,
  v_init_dist_met_75_inc = v_init_dist_met_75_inc,
  v_init_dist_bnx_25_inc = v_init_dist_bnx_25_inc,
  v_init_dist_bnx_50_inc = v_init_dist_bnx_50_inc,
  v_init_dist_bnx_75_inc = v_init_dist_bnx_75_inc,
  tx = "met",
  ce_est = "itt_ps"
)
l_owsa_bnx_itt_ps <- owsa_params(
  df_dsa_owsa = df_dsa_owsa,
  v_init_dist_met_25_inc = v_init_dist_met_25_inc,
  v_init_dist_met_50_inc = v_init_dist_met_50_inc,
  v_init_dist_met_75_inc = v_init_dist_met_75_inc,
  v_init_dist_bnx_25_inc = v_init_dist_bnx_25_inc,
  v_init_dist_bnx_50_inc = v_init_dist_bnx_50_inc,
  v_init_dist_bnx_75_inc = v_init_dist_bnx_75_inc,
  tx = "bnx",
  ce_est = "itt_ps"
)
l_owsa_met_pp <- owsa_params(
  df_dsa_owsa = df_dsa_owsa,
  v_init_dist_met_25_inc = v_init_dist_met_25_inc,
  v_init_dist_met_50_inc = v_init_dist_met_50_inc,
  v_init_dist_met_75_inc = v_init_dist_met_75_inc,
  v_init_dist_bnx_25_inc = v_init_dist_bnx_25_inc,
  v_init_dist_bnx_50_inc = v_init_dist_bnx_50_inc,
  v_init_dist_bnx_75_inc = v_init_dist_bnx_75_inc,
  tx = "met",
  ce_est = "pp"
)
l_owsa_bnx_pp <- owsa_params(
  df_dsa_owsa = df_dsa_owsa,
  v_init_dist_met_25_inc = v_init_dist_met_25_inc,
  v_init_dist_met_50_inc = v_init_dist_met_50_inc,
  v_init_dist_met_75_inc = v_init_dist_met_75_inc,
  v_init_dist_bnx_25_inc = v_init_dist_bnx_25_inc,
  v_init_dist_bnx_50_inc = v_init_dist_bnx_50_inc,
  v_init_dist_bnx_75_inc = v_init_dist_bnx_75_inc,
  tx = "bnx",
  ce_est = "pp"
)

# MET
# ITT
v_owsa_met_itt_ps_low <- l_owsa_met_itt_ps$v_owsa_low
v_owsa_met_itt_ps_high <- l_owsa_met_itt_ps$v_owsa_high
# PP
v_owsa_met_pp_low <- l_owsa_met_pp$v_owsa_low
v_owsa_met_pp_high <- l_owsa_met_pp$v_owsa_high
# Cohort balance (same for ITT and PP)
v_owsa_met_cohort_balance_25_inc <- l_owsa_met_itt_ps$v_cohort_balance_25_inc
v_owsa_met_cohort_balance_50_inc <- l_owsa_met_itt_ps$v_cohort_balance_50_inc
v_owsa_met_cohort_balance_75_inc <- l_owsa_met_itt_ps$v_cohort_balance_75_inc

# BNX
# ITT
v_owsa_bnx_itt_ps_low <- l_owsa_bnx_itt_ps$v_owsa_low
v_owsa_bnx_itt_ps_high <- l_owsa_bnx_itt_ps$v_owsa_high
# PP
v_owsa_bnx_pp_low <- l_owsa_bnx_pp$v_owsa_low
v_owsa_bnx_pp_high <- l_owsa_bnx_pp$v_owsa_high
# Cohort balance (same for ITT and PP)
v_owsa_bnx_cohort_balance_25_inc <- l_owsa_bnx_itt_ps$v_cohort_balance_25_inc
v_owsa_bnx_cohort_balance_50_inc <- l_owsa_bnx_itt_ps$v_cohort_balance_50_inc
v_owsa_bnx_cohort_balance_75_inc <- l_owsa_bnx_itt_ps$v_cohort_balance_75_inc

#### Run model ####

########################
#### Initial params ####
########################
# Base case deterministic estimates
## ITT-PS
l_out_bnx_base_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = NULL, time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_base_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = NULL, time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_base_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_base_itt_ps, outcomes_int = l_out_bnx_base_itt_ps)
# Combine
n_base_itt_ps <- l_inc_base_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled
## PP
l_out_bnx_base_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = NULL, time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_base_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = NULL, time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_base_pp <- inc_outcomes(outcomes_comp = l_out_met_base_pp, outcomes_int = l_out_bnx_base_pp)
# Combine
n_base_pp <- c(l_inc_base_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)

########################
#### Cohort balance ####
########################
## ITT-PS
# 25% Incident
l_out_bnx_cohort_balance_25_inc_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_cohort_balance_25_inc, time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_cohort_balance_25_inc_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_cohort_balance_25_inc, time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_cohort_balance_25_inc_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_cohort_balance_25_inc_itt_ps, outcomes_int = l_out_bnx_cohort_balance_25_inc_itt_ps)
# 50% Incident
# l_out_bnx_cohort_balance_50_inc_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_cohort_balance_50_inc, time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
# l_out_met_cohort_balance_50_inc_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_cohort_balance_50_inc, time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
# l_inc_cohort_balance_50_inc_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_cohort_balance_50_inc_itt_ps, outcomes_int = l_out_bnx_cohort_balance_50_inc_itt_ps)
# 75% Incident
l_out_bnx_cohort_balance_75_inc_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_cohort_balance_75_inc, time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_cohort_balance_75_inc_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_cohort_balance_75_inc, time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_cohort_balance_75_inc_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_cohort_balance_75_inc_itt_ps, outcomes_int = l_out_bnx_cohort_balance_75_inc_itt_ps)
# Combine
v_cohort_balance_itt_ps <- c(l_inc_cohort_balance_75_inc_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_cohort_balance_25_inc_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)

## PP
# 25% Incident
l_out_bnx_cohort_balance_25_inc_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_cohort_balance_25_inc, time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_cohort_balance_25_inc_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_cohort_balance_25_inc, time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_cohort_balance_25_inc_pp <- inc_outcomes(outcomes_comp = l_out_met_cohort_balance_25_inc_pp, outcomes_int = l_out_bnx_cohort_balance_25_inc_pp)
# 50% Incident
# l_out_bnx_cohort_balance_50_inc_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_cohort_balance_50_inc, time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
# l_out_met_cohort_balance_50_inc_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_cohort_balance_50_inc, time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
# l_inc_cohort_balance_50_inc_pp <- inc_outcomes(outcomes_comp = l_out_met_cohort_balance_50_inc_pp, outcomes_int = l_out_bnx_cohort_balance_50_inc_pp)
# 75% Incident
l_out_bnx_cohort_balance_75_inc_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_cohort_balance_75_inc, time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_cohort_balance_75_inc_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_cohort_balance_75_inc, time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_cohort_balance_75_inc_pp <- inc_outcomes(outcomes_comp = l_out_met_cohort_balance_75_inc_pp, outcomes_int = l_out_bnx_cohort_balance_75_inc_pp)
# Combine
v_cohort_balance_pp <- c(l_inc_cohort_balance_75_inc_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_cohort_balance_25_inc_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)

############################
#### Non-overdose death ####
############################
# Non-overdose mortality (out-of-treatment vs. in-treatment)
## ITT-PS
# Low
l_out_bnx_dno_ou_oat_low_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_low["hr_dno_ou_oat"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_dno_ou_oat_low_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_low["hr_dno_ou_oat"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_dno_ou_oat_low_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_dno_ou_oat_low_itt_ps, outcomes_int = l_out_bnx_dno_ou_oat_low_itt_ps)
# High
l_out_bnx_dno_ou_oat_high_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_high["hr_dno_ou_oat"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_dno_ou_oat_high_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_high["hr_dno_ou_oat"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_dno_ou_oat_high_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_dno_ou_oat_high_itt_ps, outcomes_int = l_out_bnx_dno_ou_oat_high_itt_ps)
# Combine
v_dno_ou_oat_itt_ps <- c(l_inc_dno_ou_oat_low_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_dno_ou_oat_high_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)
## PP
# Low
l_out_bnx_dno_ou_oat_low_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_low["hr_dno_ou_oat"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_dno_ou_oat_low_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_low["hr_dno_ou_oat"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_dno_ou_oat_low_pp <- inc_outcomes(outcomes_comp = l_out_met_dno_ou_oat_low_pp, outcomes_int = l_out_bnx_dno_ou_oat_low_pp)
# High
l_out_bnx_dno_ou_oat_high_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_high["hr_dno_ou_oat"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_dno_ou_oat_high_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_high["hr_dno_ou_oat"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_dno_ou_oat_high_pp <- inc_outcomes(outcomes_comp = l_out_met_dno_ou_oat_high_pp, outcomes_int = l_out_bnx_dno_ou_oat_high_pp)
# Combine
v_dno_ou_oat_pp <- c(l_inc_dno_ou_oat_low_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_dno_ou_oat_high_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)

###################################
#### Overdose (non-calibrated) ####
###################################
# Overdose rate in abstinence
## ITT-PS
# Low
l_out_bnx_abs_od_low_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_low["n_abs_od"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_abs_od_low_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_low["n_abs_od"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_abs_od_low_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_abs_od_low_itt_ps, outcomes_int = l_out_bnx_abs_od_low_itt_ps)
# High
l_out_bnx_abs_od_high_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_high["n_abs_od"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_abs_od_high_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_high["n_abs_od"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_abs_od_high_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_abs_od_high_itt_ps, outcomes_int = l_out_bnx_abs_od_high_itt_ps)
# Combine
v_abs_od_itt_ps <- c(l_inc_abs_od_low_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_abs_od_high_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)
## PP
# Low
l_out_bnx_abs_od_low_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_low["n_abs_od"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_abs_od_low_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_low["n_abs_od"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_abs_od_low_pp <- inc_outcomes(outcomes_comp = l_out_met_abs_od_low_pp, outcomes_int = l_out_bnx_abs_od_low_pp)
# High
l_out_bnx_abs_od_high_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_high["n_abs_od"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_abs_od_high_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_high["n_abs_od"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_abs_od_high_pp <- inc_outcomes(outcomes_comp = l_out_met_abs_od_high_pp, outcomes_int = l_out_bnx_abs_od_high_pp)
# Combine
v_abs_od_pp <- c(l_inc_abs_od_low_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_abs_od_high_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)

# Multiplier for first week out of treatment
## ITT-PS
# Low
l_out_bnx_ou_od_mult_low_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_low["n_ou_od_mult"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_ou_od_mult_low_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_low["n_ou_od_mult"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_ou_od_mult_low_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_ou_od_mult_low_itt_ps, outcomes_int = l_out_bnx_ou_od_mult_low_itt_ps)
# High
l_out_bnx_ou_od_mult_high_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_high["n_ou_od_mult"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_ou_od_mult_high_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_high["n_ou_od_mult"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_ou_od_mult_high_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_ou_od_mult_high_itt_ps, outcomes_int = l_out_bnx_ou_od_mult_high_itt_ps)
# Combine
v_ou_od_mult_itt_ps <- c(l_inc_ou_od_mult_low_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_ou_od_mult_high_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)
## PP
# Low
l_out_bnx_ou_od_mult_low_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_low["n_ou_od_mult"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_ou_od_mult_low_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_low["n_ou_od_mult"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_ou_od_mult_low_pp <- inc_outcomes(outcomes_comp = l_out_met_ou_od_mult_low_pp, outcomes_int = l_out_bnx_ou_od_mult_low_pp)
# High
l_out_bnx_ou_od_mult_high_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_high["n_ou_od_mult"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_ou_od_mult_high_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_high["n_ou_od_mult"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_ou_od_mult_high_pp <- inc_outcomes(outcomes_comp = l_out_met_ou_od_mult_high_pp, outcomes_int = l_out_bnx_ou_od_mult_high_pp)
# Combine
v_ou_od_mult_pp <- c(l_inc_ou_od_mult_low_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_ou_od_mult_high_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)

# Non-fatal overdose rate multiplier for out-of-treatment vs. OAT
## ITT-PS
# Low
l_out_bnx_ou_oat_odn_mult_low_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_low["n_ou_oat_odn_mult"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_ou_oat_odn_mult_low_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_low["n_ou_oat_odn_mult"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_ou_oat_odn_mult_low_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_ou_oat_odn_mult_low_itt_ps, outcomes_int = l_out_bnx_ou_oat_odn_mult_low_itt_ps)
# High
l_out_bnx_ou_oat_odn_mult_high_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_high["n_ou_oat_odn_mult"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_ou_oat_odn_mult_high_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_high["n_ou_oat_odn_mult"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_ou_oat_odn_mult_high_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_ou_oat_odn_mult_high_itt_ps, outcomes_int = l_out_bnx_ou_oat_odn_mult_high_itt_ps)
# Combine
v_ou_oat_odn_mult_itt_ps <- c(l_inc_ou_oat_odn_mult_low_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_ou_oat_odn_mult_high_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)
## PP
# Low
l_out_bnx_ou_oat_odn_mult_low_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_low["n_ou_oat_odn_mult"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_ou_oat_odn_mult_low_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_low["n_ou_oat_odn_mult"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_ou_oat_odn_mult_low_pp <- inc_outcomes(outcomes_comp = l_out_met_ou_oat_odn_mult_low_pp, outcomes_int = l_out_bnx_ou_oat_odn_mult_low_pp)
# High
l_out_bnx_ou_oat_odn_mult_high_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_high["n_ou_oat_odn_mult"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_ou_oat_odn_mult_high_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_high["n_ou_oat_odn_mult"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_ou_oat_odn_mult_high_pp <- inc_outcomes(outcomes_comp = l_out_met_ou_oat_odn_mult_high_pp, outcomes_int = l_out_bnx_ou_oat_odn_mult_high_pp)
# Combine
v_ou_oat_odn_mult_pp <- c(l_inc_ou_oat_odn_mult_low_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_ou_oat_odn_mult_high_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)

# Fatal overdose rate multiplier for out-of-treatment vs. OAT
## ITT-PS
# Low
l_out_bnx_ou_oat_odf_mult_low_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_low["n_ou_oat_odf_mult"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_ou_oat_odf_mult_low_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_low["n_ou_oat_odf_mult"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_ou_oat_odf_mult_low_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_ou_oat_odf_mult_low_itt_ps, outcomes_int = l_out_bnx_ou_oat_odf_mult_low_itt_ps)
# High
l_out_bnx_ou_oat_odf_mult_high_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_high["n_ou_oat_odf_mult"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_ou_oat_odf_mult_high_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_high["n_ou_oat_odf_mult"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_ou_oat_odf_mult_high_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_ou_oat_odf_mult_high_itt_ps, outcomes_int = l_out_bnx_ou_oat_odf_mult_high_itt_ps)
# Combine
v_ou_oat_odf_mult_itt_ps <- c(l_inc_ou_oat_odf_mult_low_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_ou_oat_odf_mult_high_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)
## PP
# Low
l_out_bnx_ou_oat_odf_mult_low_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_low["n_ou_oat_odf_mult"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_ou_oat_odf_mult_low_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_low["n_ou_oat_odf_mult"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_ou_oat_odf_mult_low_pp <- inc_outcomes(outcomes_comp = l_out_met_ou_oat_odf_mult_low_pp, outcomes_int = l_out_bnx_ou_oat_odf_mult_low_pp)
# High
l_out_bnx_ou_oat_odf_mult_high_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_high["n_ou_oat_odf_mult"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_ou_oat_odf_mult_high_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_high["n_ou_oat_odf_mult"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_ou_oat_odf_mult_high_pp <- inc_outcomes(outcomes_comp = l_out_met_ou_oat_odf_mult_high_pp, outcomes_int = l_out_bnx_ou_oat_odf_mult_high_pp)
# Combine
v_ou_oat_odf_mult_pp <- c(l_inc_ou_oat_odf_mult_low_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_ou_oat_odf_mult_high_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)

#############################################
#### Comparative effectiveness estimates ####
#############################################
# Treatment retention - Incident users
## ITT-PS
# Low
l_out_bnx_hr_tx_inc_low_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_low["hr_tx_itt_ps_inc"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_hr_tx_inc_low_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_low["hr_tx_itt_ps_inc"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_hr_tx_inc_low_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_hr_tx_inc_low_itt_ps, outcomes_int = l_out_bnx_hr_tx_inc_low_itt_ps)
# High
l_out_bnx_hr_tx_inc_high_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_high["hr_tx_itt_ps_inc"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_hr_tx_inc_high_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_high["hr_tx_itt_ps_inc"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_hr_tx_inc_high_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_hr_tx_inc_high_itt_ps, outcomes_int = l_out_bnx_hr_tx_inc_high_itt_ps)
# Combine
v_hr_tx_inc_itt_ps <- c(l_inc_hr_tx_inc_low_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_hr_tx_inc_high_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)
## PP
# Low
l_out_bnx_hr_tx_inc_low_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_low["hr_tx_pp_inc"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_hr_tx_inc_low_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_low["hr_tx_pp_inc"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_hr_tx_inc_low_pp <- inc_outcomes(outcomes_comp = l_out_met_hr_tx_inc_low_pp, outcomes_int = l_out_bnx_hr_tx_inc_low_pp)
# High
l_out_bnx_hr_tx_inc_high_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_high["hr_tx_pp_inc"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_hr_tx_inc_high_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_high["hr_tx_pp_inc"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_hr_tx_inc_high_pp <- inc_outcomes(outcomes_comp = l_out_met_hr_tx_inc_high_pp, outcomes_int = l_out_bnx_hr_tx_inc_high_pp)
# Combine
v_hr_tx_inc_pp <- c(l_inc_hr_tx_inc_low_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_hr_tx_inc_high_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)
# Treatment retention - Prevalent users
## ITT-PS
# Low
l_out_bnx_hr_tx_prev_low_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_low["hr_tx_itt_ps_prev"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_hr_tx_prev_low_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_low["hr_tx_itt_ps_prev"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_hr_tx_prev_low_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_hr_tx_prev_low_itt_ps, outcomes_int = l_out_bnx_hr_tx_prev_low_itt_ps)
# High
l_out_bnx_hr_tx_prev_high_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_high["hr_tx_itt_ps_prev"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_hr_tx_prev_high_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_high["hr_tx_itt_ps_prev"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_hr_tx_prev_high_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_hr_tx_prev_high_itt_ps, outcomes_int = l_out_bnx_hr_tx_prev_high_itt_ps)
# Combine
v_hr_tx_prev_itt_ps <- c(l_inc_hr_tx_prev_low_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_hr_tx_prev_high_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)
## PP
# Low
l_out_bnx_hr_tx_prev_low_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_low["hr_tx_pp_prev"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_hr_tx_prev_low_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_low["hr_tx_pp_prev"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_hr_tx_prev_low_pp <- inc_outcomes(outcomes_comp = l_out_met_hr_tx_prev_low_pp, outcomes_int = l_out_bnx_hr_tx_prev_low_pp)
# High
l_out_bnx_hr_tx_prev_high_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_high["hr_tx_pp_prev"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_hr_tx_prev_high_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_high["hr_tx_pp_prev"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_hr_tx_prev_high_pp <- inc_outcomes(outcomes_comp = l_out_met_hr_tx_prev_high_pp, outcomes_int = l_out_bnx_hr_tx_prev_high_pp)
# Combine
v_hr_tx_prev_pp <- c(l_inc_hr_tx_prev_low_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_hr_tx_prev_high_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)

# Mortality risk - Incident users
## ITT-PS
# Low
l_out_bnx_hr_death_inc_low_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_low["hr_death_itt_ps_inc"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_hr_death_inc_low_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_low["hr_death_itt_ps_inc"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_hr_death_inc_low_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_hr_death_inc_low_itt_ps, outcomes_int = l_out_bnx_hr_death_inc_low_itt_ps)
# High
l_out_bnx_hr_death_inc_high_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_high["hr_death_itt_ps_inc"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_hr_death_inc_high_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_high["hr_death_itt_ps_inc"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_hr_death_inc_high_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_hr_death_inc_high_itt_ps, outcomes_int = l_out_bnx_hr_death_inc_high_itt_ps)
# Combine
v_hr_death_inc_itt_ps <- c(l_inc_hr_death_inc_low_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_hr_death_inc_high_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)
## PP
# Low
l_out_bnx_hr_death_inc_low_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_low["hr_death_pp_inc"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_hr_death_inc_low_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_low["hr_death_pp_inc"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_hr_death_inc_low_pp <- inc_outcomes(outcomes_comp = l_out_met_hr_death_inc_low_pp, outcomes_int = l_out_bnx_hr_death_inc_low_pp)
# High
l_out_bnx_hr_death_inc_high_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_high["hr_death_pp_inc"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_hr_death_inc_high_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_high["hr_death_pp_inc"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_hr_death_inc_high_pp <- inc_outcomes(outcomes_comp = l_out_met_hr_death_inc_high_pp, outcomes_int = l_out_bnx_hr_death_inc_high_pp)
# Combine
v_hr_death_inc_pp <- c(l_inc_hr_death_inc_low_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_hr_death_inc_high_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)
# Mortality risk - Prevalent users
## ITT-PS
# Low
l_out_bnx_hr_death_prev_low_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_low["hr_death_itt_ps_prev"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_hr_death_prev_low_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_low["hr_death_itt_ps_prev"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_hr_death_prev_low_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_hr_death_prev_low_itt_ps, outcomes_int = l_out_bnx_hr_death_prev_low_itt_ps)
# High
l_out_bnx_hr_death_prev_high_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_high["hr_death_itt_ps_prev"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_hr_death_prev_high_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_high["hr_death_itt_ps_prev"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_hr_death_prev_high_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_hr_death_prev_high_itt_ps, outcomes_int = l_out_bnx_hr_death_prev_high_itt_ps)
# Combine
v_hr_death_prev_itt_ps <- c(l_inc_hr_death_prev_low_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_hr_death_prev_high_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)
## PP
# Low
l_out_bnx_hr_death_prev_low_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_low["hr_death_pp_prev"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_hr_death_prev_low_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_low["hr_death_pp_prev"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_hr_death_prev_low_pp <- inc_outcomes(outcomes_comp = l_out_met_hr_death_prev_low_pp, outcomes_int = l_out_bnx_hr_death_prev_low_pp)
# High
l_out_bnx_hr_death_prev_high_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_high["hr_death_pp_prev"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_hr_death_prev_high_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_high["hr_death_pp_prev"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_hr_death_prev_high_pp <- inc_outcomes(outcomes_comp = l_out_met_hr_death_prev_high_pp, outcomes_int = l_out_bnx_hr_death_prev_high_pp)
# Combine
v_hr_death_prev_pp <- c(l_inc_hr_death_prev_low_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_hr_death_prev_high_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)

###############################
#### Calibrated parameters ####
###############################
# Overdose rate in OAT
## ITT-PS
# Low
l_out_bnx_oat_od_low_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_low["n_oat_od"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_oat_od_low_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_low["n_oat_od"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_oat_od_low_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_oat_od_low_itt_ps, outcomes_int = l_out_bnx_oat_od_low_itt_ps)
# High
l_out_bnx_oat_od_high_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_high["n_oat_od"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_oat_od_high_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_high["n_oat_od"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_oat_od_high_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_oat_od_high_itt_ps, outcomes_int = l_out_bnx_oat_od_high_itt_ps)
# Combine
v_oat_od_itt_ps <- c(l_inc_oat_od_low_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_oat_od_high_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)
## PP
# Low
l_out_bnx_oat_od_low_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_low["n_oat_od"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_oat_od_low_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_low["n_oat_od"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_oat_od_low_pp <- inc_outcomes(outcomes_comp = l_out_met_oat_od_low_pp, outcomes_int = l_out_bnx_oat_od_low_pp)
# High
l_out_bnx_oat_od_high_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_high["n_oat_od"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_oat_od_high_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_high["n_oat_od"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_oat_od_high_pp <- inc_outcomes(outcomes_comp = l_out_met_oat_od_high_pp, outcomes_int = l_out_bnx_oat_od_high_pp)
# Combine
v_oat_od_pp <- c(l_inc_oat_od_low_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_oat_od_high_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)

# Overdose rate multiplier on fentanyl prevalence
## ITT-PS
# Low
l_out_bnx_fent_prev_od_mult_low_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_low["n_fent_prev_od_mult"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_fent_prev_od_mult_low_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_low["n_fent_prev_od_mult"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_fent_prev_od_mult_low_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_fent_prev_od_mult_low_itt_ps, outcomes_int = l_out_bnx_fent_prev_od_mult_low_itt_ps)
# High
l_out_bnx_fent_prev_od_mult_high_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_high["n_fent_prev_od_mult"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_fent_prev_od_mult_high_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_high["n_fent_prev_od_mult"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_fent_prev_od_mult_high_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_fent_prev_od_mult_high_itt_ps, outcomes_int = l_out_bnx_fent_prev_od_mult_high_itt_ps)
# Combine
v_fent_prev_od_mult_itt_ps <- c(l_inc_fent_prev_od_mult_low_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_fent_prev_od_mult_high_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)
## PP
# Low
l_out_bnx_fent_prev_od_mult_low_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_low["n_fent_prev_od_mult"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_fent_prev_od_mult_low_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_low["n_fent_prev_od_mult"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_fent_prev_od_mult_low_pp <- inc_outcomes(outcomes_comp = l_out_met_fent_prev_od_mult_low_pp, outcomes_int = l_out_bnx_fent_prev_od_mult_low_pp)
# High
l_out_bnx_fent_prev_od_mult_high_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_high["n_fent_prev_od_mult"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_fent_prev_od_mult_high_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_high["n_fent_prev_od_mult"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_fent_prev_od_mult_high_pp <- inc_outcomes(outcomes_comp = l_out_met_fent_prev_od_mult_high_pp, outcomes_int = l_out_bnx_fent_prev_od_mult_high_pp)
# Combine
v_fent_prev_od_mult_pp <- c(l_inc_fent_prev_od_mult_low_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_fent_prev_od_mult_high_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)

# Overdose rate multiplier on change in fentanyl prevalence
## ITT-PS
# Low
l_out_bnx_fent_delta_od_mult_low_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_low["n_fent_delta_od_mult"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_fent_delta_od_mult_low_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_low["n_fent_delta_od_mult"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_fent_delta_od_mult_low_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_fent_delta_od_mult_low_itt_ps, outcomes_int = l_out_bnx_fent_delta_od_mult_low_itt_ps)
# High
l_out_bnx_fent_delta_od_mult_high_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_high["n_fent_delta_od_mult"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_fent_delta_od_mult_high_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_high["n_fent_delta_od_mult"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_fent_delta_od_mult_high_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_fent_delta_od_mult_high_itt_ps, outcomes_int = l_out_bnx_fent_delta_od_mult_high_itt_ps)
# Combine
v_fent_delta_od_mult_itt_ps <- c(l_inc_fent_delta_od_mult_low_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_fent_delta_od_mult_high_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)
## PP
# Low
l_out_bnx_fent_delta_od_mult_low_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_low["n_fent_delta_od_mult"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_fent_delta_od_mult_low_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_low["n_fent_delta_od_mult"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_fent_delta_od_mult_low_pp <- inc_outcomes(outcomes_comp = l_out_met_fent_delta_od_mult_low_pp, outcomes_int = l_out_bnx_fent_delta_od_mult_low_pp)
# High
l_out_bnx_fent_delta_od_mult_high_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_high["n_fent_delta_od_mult"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_fent_delta_od_mult_high_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_high["n_fent_delta_od_mult"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_fent_delta_od_mult_high_pp <- inc_outcomes(outcomes_comp = l_out_met_fent_delta_od_mult_high_pp, outcomes_int = l_out_bnx_fent_delta_od_mult_high_pp)
# Combine
v_fent_delta_od_mult_pp <- c(l_inc_fent_delta_od_mult_low_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_fent_delta_od_mult_high_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)

# Fatal overdose rate in OAT
## ITT-PS
# Low
l_out_bnx_fatal_od_oat_low_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_low["n_fatal_od_oat"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_fatal_od_oat_low_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_low["n_fatal_od_oat"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_fatal_od_oat_low_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_fatal_od_oat_low_itt_ps, outcomes_int = l_out_bnx_fatal_od_oat_low_itt_ps)
# High
l_out_bnx_fatal_od_oat_high_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_high["n_fatal_od_oat"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_fatal_od_oat_high_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_high["n_fatal_od_oat"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_fatal_od_oat_high_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_fatal_od_oat_high_itt_ps, outcomes_int = l_out_bnx_fatal_od_oat_high_itt_ps)
# Combine
v_fatal_od_oat_itt_ps <- c(l_inc_fatal_od_oat_low_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_fatal_od_oat_high_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)
## PP
# Low
l_out_bnx_fatal_od_oat_low_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_low["n_fatal_od_oat"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_fatal_od_oat_low_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_low["n_fatal_od_oat"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_fatal_od_oat_low_pp <- inc_outcomes(outcomes_comp = l_out_met_fatal_od_oat_low_pp, outcomes_int = l_out_bnx_fatal_od_oat_low_pp)
# High
l_out_bnx_fatal_od_oat_high_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_high["n_fatal_od_oat"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_fatal_od_oat_high_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_high["n_fatal_od_oat"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_fatal_od_oat_high_pp <- inc_outcomes(outcomes_comp = l_out_met_fatal_od_oat_high_pp, outcomes_int = l_out_bnx_fatal_od_oat_high_pp)
# Combine
v_fatal_od_oat_pp <- c(l_inc_fatal_od_oat_low_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_fatal_od_oat_high_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)

# Non-overdose mortality rate in OAT
## ITT-PS
# Low
l_out_bnx_hr_oat_low_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_low["hr_oat"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_hr_oat_low_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_low["hr_oat"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_hr_oat_low_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_hr_oat_low_itt_ps, outcomes_int = l_out_bnx_hr_oat_low_itt_ps)
# High
l_out_bnx_hr_oat_high_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_high["hr_oat"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_hr_oat_high_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_high["hr_oat"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_hr_oat_high_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_hr_oat_high_itt_ps, outcomes_int = l_out_bnx_hr_oat_high_itt_ps)
# Combine
v_hr_oat_itt_ps <- c(l_inc_hr_oat_low_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_hr_oat_high_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)
## PP
# Low
l_out_bnx_hr_oat_low_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_low["hr_oat"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_hr_oat_low_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_low["hr_oat"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_hr_oat_low_pp <- inc_outcomes(outcomes_comp = l_out_met_hr_oat_low_pp, outcomes_int = l_out_bnx_hr_oat_low_pp)
# High
l_out_bnx_hr_oat_high_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_high["hr_oat"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_hr_oat_high_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_high["hr_oat"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_hr_oat_high_pp <- inc_outcomes(outcomes_comp = l_out_met_hr_oat_high_pp, outcomes_int = l_out_bnx_hr_oat_high_pp)
# Combine
v_hr_oat_pp <- c(l_inc_hr_oat_low_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_hr_oat_high_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)

# Weibull scale (out-of-treatment)
## ITT-PS
# Low
l_out_bnx_weibull_scale_ou_low_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_low["p_weibull_scale_ou"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_weibull_scale_ou_low_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_low["p_weibull_scale_ou"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_weibull_scale_ou_low_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_weibull_scale_ou_low_itt_ps, outcomes_int = l_out_bnx_weibull_scale_ou_low_itt_ps)
# High
l_out_bnx_weibull_scale_ou_high_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_high["p_weibull_scale_ou"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_weibull_scale_ou_high_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_high["p_weibull_scale_ou"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_weibull_scale_ou_high_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_weibull_scale_ou_high_itt_ps, outcomes_int = l_out_bnx_weibull_scale_ou_high_itt_ps)
# Combine
v_weibull_scale_ou_itt_ps <- c(l_inc_weibull_scale_ou_low_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_weibull_scale_ou_high_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)
## PP
# Low
l_out_bnx_weibull_scale_ou_low_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_low["p_weibull_scale_ou"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_weibull_scale_ou_low_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_low["p_weibull_scale_ou"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_weibull_scale_ou_low_pp <- inc_outcomes(outcomes_comp = l_out_met_weibull_scale_ou_low_pp, outcomes_int = l_out_bnx_weibull_scale_ou_low_pp)
# High
l_out_bnx_weibull_scale_ou_high_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_high["p_weibull_scale_ou"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_weibull_scale_ou_high_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_high["p_weibull_scale_ou"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_weibull_scale_ou_high_pp <- inc_outcomes(outcomes_comp = l_out_met_weibull_scale_ou_high_pp, outcomes_int = l_out_bnx_weibull_scale_ou_high_pp)
# Combine
v_weibull_scale_ou_pp <- c(l_inc_weibull_scale_ou_low_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_weibull_scale_ou_high_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)

# Weibull shape (out-of-treatment)
## ITT-PS
# Low
l_out_bnx_weibull_shape_ou_low_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_low["p_weibull_shape_ou"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_weibull_shape_ou_low_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_low["p_weibull_shape_ou"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_weibull_shape_ou_low_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_weibull_shape_ou_low_itt_ps, outcomes_int = l_out_bnx_weibull_shape_ou_low_itt_ps)
# High
l_out_bnx_weibull_shape_ou_high_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_high["p_weibull_shape_ou"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_weibull_shape_ou_high_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_high["p_weibull_shape_ou"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_weibull_shape_ou_high_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_weibull_shape_ou_high_itt_ps, outcomes_int = l_out_bnx_weibull_shape_ou_high_itt_ps)
# Combine
v_weibull_shape_ou_itt_ps <- c(l_inc_weibull_shape_ou_low_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_weibull_shape_ou_high_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)
## PP
# Low
l_out_bnx_weibull_shape_ou_low_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_low["p_weibull_shape_ou"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_weibull_shape_ou_low_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_low["p_weibull_shape_ou"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_weibull_shape_ou_low_pp <- inc_outcomes(outcomes_comp = l_out_met_weibull_shape_ou_low_pp, outcomes_int = l_out_bnx_weibull_shape_ou_low_pp)
# High
l_out_bnx_weibull_shape_ou_high_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_high["p_weibull_shape_ou"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_weibull_shape_ou_high_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_high["p_weibull_shape_ou"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_weibull_shape_ou_high_pp <- inc_outcomes(outcomes_comp = l_out_met_weibull_shape_ou_high_pp, outcomes_int = l_out_bnx_weibull_shape_ou_high_pp)
# Combine
v_weibull_shape_ou_pp <- c(l_inc_weibull_shape_ou_low_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_weibull_shape_ou_high_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)

# Weibull scale (abstinence)
## ITT-PS
# Low
l_out_bnx_weibull_scale_abs_low_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_low["p_weibull_scale_abs"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_weibull_scale_abs_low_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_low["p_weibull_scale_abs"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_weibull_scale_abs_low_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_weibull_scale_abs_low_itt_ps, outcomes_int = l_out_bnx_weibull_scale_abs_low_itt_ps)
# High
l_out_bnx_weibull_scale_abs_high_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_high["p_weibull_scale_abs"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_weibull_scale_abs_high_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_high["p_weibull_scale_abs"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_weibull_scale_abs_high_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_weibull_scale_abs_high_itt_ps, outcomes_int = l_out_bnx_weibull_scale_abs_high_itt_ps)
# Combine
v_weibull_scale_abs_itt_ps <- c(l_inc_weibull_scale_abs_low_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_weibull_scale_abs_high_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)
## PP
# Low
l_out_bnx_weibull_scale_abs_low_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_low["p_weibull_scale_abs"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_weibull_scale_abs_low_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_low["p_weibull_scale_abs"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_weibull_scale_abs_low_pp <- inc_outcomes(outcomes_comp = l_out_met_weibull_scale_abs_low_pp, outcomes_int = l_out_bnx_weibull_scale_abs_low_pp)
# High
l_out_bnx_weibull_scale_abs_high_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_high["p_weibull_scale_abs"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_weibull_scale_abs_high_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_high["p_weibull_scale_abs"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_weibull_scale_abs_high_pp <- inc_outcomes(outcomes_comp = l_out_met_weibull_scale_abs_high_pp, outcomes_int = l_out_bnx_weibull_scale_abs_high_pp)
# Combine
v_weibull_scale_abs_pp <- c(l_inc_weibull_scale_abs_low_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_weibull_scale_abs_high_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)

# Weibull shape (abstinence)
## ITT-PS
# Low
l_out_bnx_weibull_shape_abs_low_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_low["p_weibull_shape_abs"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_weibull_shape_abs_low_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_low["p_weibull_shape_abs"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_weibull_shape_abs_low_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_weibull_shape_abs_low_itt_ps, outcomes_int = l_out_bnx_weibull_shape_abs_low_itt_ps)
# High
l_out_bnx_weibull_shape_abs_high_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_high["p_weibull_shape_abs"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_weibull_shape_abs_high_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_high["p_weibull_shape_abs"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_weibull_shape_abs_high_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_weibull_shape_abs_high_itt_ps, outcomes_int = l_out_bnx_weibull_shape_abs_high_itt_ps)
# Combine
v_weibull_shape_abs_itt_ps <- c(l_inc_weibull_shape_abs_low_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_weibull_shape_abs_high_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)
## PP
# Low
l_out_bnx_weibull_shape_abs_low_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_low["p_weibull_shape_abs"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_weibull_shape_abs_low_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_low["p_weibull_shape_abs"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_weibull_shape_abs_low_pp <- inc_outcomes(outcomes_comp = l_out_met_weibull_shape_abs_low_pp, outcomes_int = l_out_bnx_weibull_shape_abs_low_pp)
# High
l_out_bnx_weibull_shape_abs_high_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_high["p_weibull_shape_abs"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_weibull_shape_abs_high_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_high["p_weibull_shape_abs"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_weibull_shape_abs_high_pp <- inc_outcomes(outcomes_comp = l_out_met_weibull_shape_abs_high_pp, outcomes_int = l_out_bnx_weibull_shape_abs_high_pp)
# Combine
v_weibull_shape_abs_pp <- c(l_inc_weibull_shape_abs_low_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_weibull_shape_abs_high_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)

##################################
#### Weibull methadone remain ####
##################################
# Weibull scale (methadone - incident)
## ITT-PS
# Low
l_out_bnx_weibull_scale_met_inc_low_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_low["p_weibull_scale_met_inc"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_weibull_scale_met_inc_low_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_low["p_weibull_scale_met_inc"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_weibull_scale_met_inc_low_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_weibull_scale_met_inc_low_itt_ps, outcomes_int = l_out_bnx_weibull_scale_met_inc_low_itt_ps)
# High
l_out_bnx_weibull_scale_met_inc_high_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_high["p_weibull_scale_met_inc"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_weibull_scale_met_inc_high_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_high["p_weibull_scale_met_inc"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_weibull_scale_met_inc_high_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_weibull_scale_met_inc_high_itt_ps, outcomes_int = l_out_bnx_weibull_scale_met_inc_high_itt_ps)
# Combine
v_weibull_scale_met_inc_itt_ps <- c(l_inc_weibull_scale_met_inc_low_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_weibull_scale_met_inc_high_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)
## PP
# Low
l_out_bnx_weibull_scale_met_inc_low_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_low["p_weibull_scale_met_inc"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_weibull_scale_met_inc_low_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_low["p_weibull_scale_met_inc"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_weibull_scale_met_inc_low_pp <- inc_outcomes(outcomes_comp = l_out_met_weibull_scale_met_inc_low_pp, outcomes_int = l_out_bnx_weibull_scale_met_inc_low_pp)
# High
l_out_bnx_weibull_scale_met_inc_high_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_high["p_weibull_scale_met_inc"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_weibull_scale_met_inc_high_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_high["p_weibull_scale_met_inc"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_weibull_scale_met_inc_high_pp <- inc_outcomes(outcomes_comp = l_out_met_weibull_scale_met_inc_high_pp, outcomes_int = l_out_bnx_weibull_scale_met_inc_high_pp)
# Combine
v_weibull_scale_met_inc_pp <- c(l_inc_weibull_scale_met_inc_low_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_weibull_scale_met_inc_high_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)

# Weibull shape (methadone - incident)
## ITT-PS
# Low
l_out_bnx_weibull_shape_met_inc_low_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_low["p_weibull_shape_met_inc"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_weibull_shape_met_inc_low_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_low["p_weibull_shape_met_inc"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_weibull_shape_met_inc_low_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_weibull_shape_met_inc_low_itt_ps, outcomes_int = l_out_bnx_weibull_shape_met_inc_low_itt_ps)
# High
l_out_bnx_weibull_shape_met_inc_high_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_high["p_weibull_shape_met_inc"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_weibull_shape_met_inc_high_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_high["p_weibull_shape_met_inc"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_weibull_shape_met_inc_high_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_weibull_shape_met_inc_high_itt_ps, outcomes_int = l_out_bnx_weibull_shape_met_inc_high_itt_ps)
# Combine
v_weibull_shape_met_inc_itt_ps <- c(l_inc_weibull_shape_met_inc_low_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_weibull_shape_met_inc_high_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)
## PP
# Low
l_out_bnx_weibull_shape_met_inc_low_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_low["p_weibull_shape_met_inc"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_weibull_shape_met_inc_low_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_low["p_weibull_shape_met_inc"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_weibull_shape_met_inc_low_pp <- inc_outcomes(outcomes_comp = l_out_met_weibull_shape_met_inc_low_pp, outcomes_int = l_out_bnx_weibull_shape_met_inc_low_pp)
# High
l_out_bnx_weibull_shape_met_inc_high_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_high["p_weibull_shape_met_inc"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_weibull_shape_met_inc_high_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_high["p_weibull_shape_met_inc"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_weibull_shape_met_inc_high_pp <- inc_outcomes(outcomes_comp = l_out_met_weibull_shape_met_inc_high_pp, outcomes_int = l_out_bnx_weibull_shape_met_inc_high_pp)
# Combine
v_weibull_shape_met_inc_pp <- c(l_inc_weibull_shape_met_inc_low_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_weibull_shape_met_inc_high_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)

# Weibull scale (methadone - incident)
## ITT-PS
# Low
l_out_bnx_weibull_scale_met_inc_low_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_low["p_weibull_scale_met_inc"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_weibull_scale_met_inc_low_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_low["p_weibull_scale_met_inc"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_weibull_scale_met_inc_low_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_weibull_scale_met_inc_low_itt_ps, outcomes_int = l_out_bnx_weibull_scale_met_inc_low_itt_ps)
# High
l_out_bnx_weibull_scale_met_inc_high_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_high["p_weibull_scale_met_inc"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_weibull_scale_met_inc_high_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_high["p_weibull_scale_met_inc"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_weibull_scale_met_inc_high_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_weibull_scale_met_inc_high_itt_ps, outcomes_int = l_out_bnx_weibull_scale_met_inc_high_itt_ps)
# Combine
v_weibull_scale_met_inc_itt_ps <- c(l_inc_weibull_scale_met_inc_low_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_weibull_scale_met_inc_high_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)
## PP
# Low
l_out_bnx_weibull_scale_met_inc_low_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_low["p_weibull_scale_met_inc"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_weibull_scale_met_inc_low_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_low["p_weibull_scale_met_inc"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_weibull_scale_met_inc_low_pp <- inc_outcomes(outcomes_comp = l_out_met_weibull_scale_met_inc_low_pp, outcomes_int = l_out_bnx_weibull_scale_met_inc_low_pp)
# High
l_out_bnx_weibull_scale_met_inc_high_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_high["p_weibull_scale_met_inc"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_weibull_scale_met_inc_high_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_high["p_weibull_scale_met_inc"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_weibull_scale_met_inc_high_pp <- inc_outcomes(outcomes_comp = l_out_met_weibull_scale_met_inc_high_pp, outcomes_int = l_out_bnx_weibull_scale_met_inc_high_pp)
# Combine
v_weibull_scale_met_inc_pp <- c(l_inc_weibull_scale_met_inc_low_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_weibull_scale_met_inc_high_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)

# Weibull shape (methadone - prevalent)
## ITT-PS
# Low
l_out_bnx_weibull_shape_met_prev_low_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_low["p_weibull_shape_met_prev"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_weibull_shape_met_prev_low_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_low["p_weibull_shape_met_prev"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_weibull_shape_met_prev_low_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_weibull_shape_met_prev_low_itt_ps, outcomes_int = l_out_bnx_weibull_shape_met_prev_low_itt_ps)
# High
l_out_bnx_weibull_shape_met_prev_high_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_high["p_weibull_shape_met_prev"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_weibull_shape_met_prev_high_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_high["p_weibull_shape_met_prev"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_weibull_shape_met_prev_high_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_weibull_shape_met_prev_high_itt_ps, outcomes_int = l_out_bnx_weibull_shape_met_prev_high_itt_ps)
# Combine
v_weibull_shape_met_prev_itt_ps <- c(l_inc_weibull_shape_met_prev_low_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_weibull_shape_met_prev_high_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)
## PP
# Low
l_out_bnx_weibull_shape_met_prev_low_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_low["p_weibull_shape_met_prev"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_weibull_shape_met_prev_low_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_low["p_weibull_shape_met_prev"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_weibull_shape_met_prev_low_pp <- inc_outcomes(outcomes_comp = l_out_met_weibull_shape_met_prev_low_pp, outcomes_int = l_out_bnx_weibull_shape_met_prev_low_pp)
# High
l_out_bnx_weibull_shape_met_prev_high_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_high["p_weibull_shape_met_prev"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_weibull_shape_met_prev_high_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_high["p_weibull_shape_met_prev"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_weibull_shape_met_prev_high_pp <- inc_outcomes(outcomes_comp = l_out_met_weibull_shape_met_prev_high_pp, outcomes_int = l_out_bnx_weibull_shape_met_prev_high_pp)
# Combine
v_weibull_shape_met_prev_pp <- c(l_inc_weibull_shape_met_prev_low_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_weibull_shape_met_prev_high_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)

# Weibull scale (methadone - prevalent)
## ITT-PS
# Low
l_out_bnx_weibull_scale_met_prev_low_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_low["p_weibull_scale_met_prev"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_weibull_scale_met_prev_low_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_low["p_weibull_scale_met_prev"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_weibull_scale_met_prev_low_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_weibull_scale_met_prev_low_itt_ps, outcomes_int = l_out_bnx_weibull_scale_met_prev_low_itt_ps)
# High
l_out_bnx_weibull_scale_met_prev_high_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_bnx_itt_ps_high["p_weibull_scale_met_prev"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_weibull_scale_met_prev_high_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_mean, v_params_dsa = v_owsa_met_itt_ps_high["p_weibull_scale_met_prev"], time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only", checks = FALSE)
l_inc_weibull_scale_met_prev_high_itt_ps <- inc_outcomes(outcomes_comp = l_out_met_weibull_scale_met_prev_high_itt_ps, outcomes_int = l_out_bnx_weibull_scale_met_prev_high_itt_ps)
# Combine
v_weibull_scale_met_prev_itt_ps <- c(l_inc_weibull_scale_met_prev_low_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_weibull_scale_met_prev_high_itt_ps$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)
## PP
# Low
l_out_bnx_weibull_scale_met_prev_low_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_low["p_weibull_scale_met_prev"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_weibull_scale_met_prev_low_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_low["p_weibull_scale_met_prev"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_weibull_scale_met_prev_low_pp <- inc_outcomes(outcomes_comp = l_out_met_weibull_scale_met_prev_low_pp, outcomes_int = l_out_bnx_weibull_scale_met_prev_low_pp)
# High
l_out_bnx_weibull_scale_met_prev_high_pp <- outcomes(l_params_all = l_params_bnx_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_bnx_pp_high["p_weibull_scale_met_prev"], time_horizon = "full", ce_est = "pp", analytic_cohort = "bnx_only", checks = FALSE)
l_out_met_weibull_scale_met_prev_high_pp <- outcomes(l_params_all = l_params_met_pp, v_params_calib = l_imis_output_pp$v_calib_post_mean, v_params_dsa = v_owsa_met_pp_high["p_weibull_scale_met_prev"], time_horizon = "full", ce_est = "pp", analytic_cohort = "met_only", checks = FALSE)
l_inc_weibull_scale_met_prev_high_pp <- inc_outcomes(outcomes_comp = l_out_met_weibull_scale_met_prev_high_pp, outcomes_int = l_out_bnx_weibull_scale_met_prev_high_pp)
# Combine
v_weibull_scale_met_prev_pp <- c(l_inc_weibull_scale_met_prev_low_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled, l_inc_weibull_scale_met_prev_high_pp$df_incremental_scaled$n_inc_qalys_adj_2020_scaled)

# Combine into matrix
# ITT-PS
m_owsa_itt_ps <- rbind(
  v_cohort_balance_itt_ps,
  v_dno_ou_oat_itt_ps, v_abs_od_itt_ps, v_ou_od_mult_itt_ps, v_ou_oat_odn_mult_itt_ps, v_ou_oat_odf_mult_itt_ps,
  v_hr_tx_inc_itt_ps, v_hr_death_inc_itt_ps, v_hr_tx_prev_itt_ps, v_hr_death_prev_itt_ps,
  v_oat_od_itt_ps, v_fent_prev_od_mult_itt_ps, v_fent_delta_od_mult_itt_ps, v_fatal_od_oat_itt_ps, v_hr_oat_itt_ps, v_weibull_scale_ou_itt_ps, v_weibull_shape_ou_itt_ps, v_weibull_scale_abs_itt_ps, v_weibull_shape_abs_itt_ps,
  v_weibull_scale_met_inc_itt_ps, v_weibull_shape_met_inc_itt_ps, v_weibull_scale_met_prev_itt_ps, v_weibull_shape_met_prev_itt_ps
)
# PP
m_owsa_pp <- rbind(
  v_cohort_balance_pp,
  v_dno_ou_oat_pp, v_abs_od_pp, v_ou_od_mult_pp, v_ou_oat_odn_mult_pp, v_ou_oat_odf_mult_pp,
  v_hr_tx_inc_pp, v_hr_death_inc_pp, v_hr_tx_prev_pp, v_hr_death_prev_pp,
  v_oat_od_pp, v_fent_prev_od_mult_pp, v_fent_delta_od_mult_pp, v_fatal_od_oat_pp, v_hr_oat_pp, v_weibull_scale_ou_pp, v_weibull_shape_ou_pp, v_weibull_scale_abs_pp, v_weibull_shape_abs_pp,
  v_weibull_scale_met_inc_pp, v_weibull_shape_met_inc_pp, v_weibull_scale_met_prev_pp, v_weibull_shape_met_prev_pp
)

df_owsa_itt_ps <- as.data.frame(m_owsa_itt_ps)
df_owsa_pp <- as.data.frame(m_owsa_pp)
colnames(df_owsa_itt_ps) <- colnames(df_owsa_pp) <- c("Lower", "Upper")
df_owsa_itt_ps <- as_tibble(df_owsa_itt_ps) %>% add_column(var_name = c(
  "Cohort balance (percentage incident users)",
  "RR non-overdose death (out-of-treatment vs. OAT)",
  "Overdose rate (abstinence)",
  "RR overdose (1st week out-of-treatment vs. week 2+)",
  "RR overdose (out-of-treatment vs. OAT)",
  "RR fatal overdose (out-of-treatment vs. OAT)",
  "CE-HR treatment discontinuation (incident users)",
  "CE-HR mortality (incident users)",
  "CE-HR treatment discontinuation (experienced users)",
  "CE-HR mortality (experienced users)",
  "Overdose rate (OAT)",
  "Overdose rate multiplier (fentanyl prevalence)",
  "Overdose rate multiplier (fentanyl delta)",
  "Fatal overdose rate (OAT)",
  "RR non-overdose death (OAT vs. abstinence)",
  "Weibull scale (out-of-treatment)",
  "Weibull shape (out-of-treatment)",
  "Weibull scale (abstinence)",
  "Weibull shape (abstinence)",
  "Weibull scale (methadone - incident users)",
  "Weibull shape (methadone - incident users)",
  "Weibull scale (methadone - experienced users)",
  "Weibull shape (methadone - experienced users)"
))
df_owsa_pp <- as_tibble(df_owsa_pp) %>% add_column(var_name = c(
  "Cohort balance (percentage incident users)",
  "RR non-overdose death (out-of-treatment vs. OAT)",
  "Overdose rate (abstinence)",
  "RR overdose (1st week out-of-treatment vs. week 2+)",
  "RR overdose (out-of-treatment vs. OAT)",
  "RR fatal overdose (out-of-treatment vs. OAT)",
  "CE-HR treatment discontinuation (incident users)",
  "CE-HR mortality (incident users)",
  "CE-HR treatment discontinuation (experienced users)",
  "CE-HR mortality (experienced users)",
  "Overdose rate (OAT)",
  "Overdose rate multiplier (fentanyl prevalence)",
  "Overdose rate multiplier (fentanyl delta)",
  "Fatal overdose rate (OAT)",
  "RR non-overdose death (OAT vs. abstinence)",
  "Weibull scale (out-of-treatment)",
  "Weibull shape (out-of-treatment)",
  "Weibull scale (abstinence)",
  "Weibull shape (abstinence)",
  "Weibull scale (methadone - incident users)",
  "Weibull shape (methadone - incident users)",
  "Weibull scale (methadone - experienced users)",
  "Weibull shape (methadone - experienced users)"
))

l_owsa_itt_ps <- list(
  df_owsa_itt_ps = df_owsa_itt_ps,
  n_base_itt_ps = n_base_itt_ps
)
l_owsa_pp <- list(
  df_owsa_pp = df_owsa_pp,
  n_base_pp = n_base_pp
)
# Save outputs
## As .RData ##
save(l_owsa_itt_ps,
  file = "outputs/dsa/owsa/outcomes_owsa_itt_ps.RData"
)
save(l_owsa_pp,
  file = "outputs/dsa/owsa/outcomes_owsa_pp.RData"
)
