rm(list = ls()) # to clean the workspace

library(dplyr)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(parallel)
library(foreach)
library(doParallel)
library(tidyr)

# Call model setup functions
source("R/input_parameter_functions.R")
source("R/generate_psa_parameters.R")
source("R/model_setup_functions.R")
source("R/outcome_functions.R")
source("R/psa_analysis.R")
source("R/mortality_function.R")
source("R/overdose_prob_function.R")

# Load parameters
source("Analysis/00_load_parameters.R")

# Set population size for dirichlet draws
n_pop_cohort <- l_params_bnx_itt$n_pop_oat # update with cohort size
# n_pop_est <- l_params_bnx_itt$n_pop_est
n_sim <- 10000 # just to test function (will be set as n_sim)

### PSA model outputs
### Run Markov model for PSA draws and return outputs ###
# Generate PSA parameter draws
# Intent-to-treat
df_psa_params_itt <- generate_psa_params(
  n_sim = n_sim, seed = 3730687, n_pop = n_pop_cohort,
  file.death_hr = "data/death_hr.csv",
  file.weibull = "data/weibull_itt.csv",
  file.ce_tx = "data/ce_tx.csv",
  file.ce_death = "data/ce_death.csv",
  file.overdose = "data/overdose.csv",
  file.fentanyl = "data/fentanyl.csv",
  file.naloxone = "data/naloxone.csv",
  file.imis_output = "outputs/Calibration/itt/imis_output_itt.RData"
)

df_psa_params_pp <- generate_psa_params(
  n_sim = n_sim, seed = 3730687, n_pop = n_pop_cohort,
  file.death_hr = "data/death_hr.csv",
  file.weibull = "data/weibull_pp.csv",
  file.ce_tx = "data/ce_tx.csv",
  file.ce_death = "data/ce_death.csv",
  file.overdose = "data/overdose.csv",
  file.fentanyl = "data/fentanyl.csv",
  file.naloxone = "data/naloxone.csv",
  file.imis_output = "outputs/Calibration/pp/imis_output_pp.RData"
)

# Per-protocol

# Output data
## As .RData
save(df_psa_params_itt,
  file = "outputs/psa/inputs/df_psa_params_itt.RData"
)
save(df_psa_params_pp,
  file = "outputs/psa/inputs/df_psa_params_pp.RData"
)

## As .csv
write.csv(df_psa_params_itt, "outputs/psa/inputs/input_psa_itt.csv",
  row.names = TRUE
)
write.csv(df_psa_params_pp, "outputs/psa/inputs/input_psa_pp.csv",
  row.names = TRUE
)

# Load PSA inputs
# load(file = "outputs/psa/inputs/df_psa_params_itt.RData")
load(file = "outputs/psa/inputs/df_psa_params_pp.RData")

# Set number of cores
n_cores <- detectCores()
cl <- makeCluster(n_cores - 1)
clusterExport(cl, c("markov_model", "update_param_list", "outcomes", "inc_outcomes", "overdose", "mort")) # load all functions used on clusters
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
l_psa_itt_ps <- psa_analysis(ce_est = "itt_ps")
Sys.time()
l_psa_pp <- psa_analysis(ce_est = "pp")
Sys.time()
l_psa_pp_hd <- psa_analysis(ce_est = "pp_hd")
Sys.time()

stopImplicitCluster()

### Output results
## ITT - Propensity score
df_outcomes_met_psa_itt_ps <- l_psa_itt_ps$df_outcomes_met_psa_comb
df_outcomes_bnx_psa_itt_ps <- l_psa_itt_ps$df_outcomes_bnx_psa_comb
df_incremental_psa_itt_ps <- l_psa_itt_ps$df_incremental_psa_comb
df_incremental_psa_itt_ps_scaled <- l_psa_itt_ps$df_incremental_psa_scaled_comb

## PP - Per protocol
df_outcomes_met_psa_pp <- l_psa_pp$df_outcomes_met_psa_comb
df_outcomes_bnx_psa_pp <- l_psa_pp$df_outcomes_bnx_psa_comb
df_incremental_psa_pp <- l_psa_pp$df_incremental_psa_comb
df_incremental_psa_pp_scaled <- l_psa_pp$df_incremental_psa_scaled_comb

## PP - Per protocol (high dose)
df_outcomes_met_psa_pp_hd <- l_psa_pp_hd$df_outcomes_met_psa_comb
df_outcomes_bnx_psa_pp_hd <- l_psa_pp_hd$df_outcomes_bnx_psa_comb
df_incremental_psa_pp_hd <- l_psa_pp_hd$df_incremental_psa_comb
df_incremental_psa_pp_hd_scaled <- l_psa_pp_hd$df_incremental_psa_scaled_comb


## As .RData
save(df_outcomes_met_psa_itt_ps,
  df_outcomes_bnx_psa_itt_ps,
  df_incremental_psa_itt_ps,
  df_incremental_psa_itt_ps_scaled,
  df_outcomes_met_psa_pp,
  df_outcomes_bnx_psa_pp,
  df_incremental_psa_pp,
  df_incremental_psa_pp_scaled,
  df_outcomes_met_psa_pp_hd,
  df_outcomes_bnx_psa_pp_hd,
  df_incremental_psa_pp_hd,
  df_incremental_psa_pp_hd_scaled,
  file = "outputs/psa/outcomes_psa.RData"
)
