library(dplyr) # to manipulate data
library(EDCimport)

# Call model setup functions
# To-do: Move into package eventually
source("R/input_parameter_functions.R")

# Load parameters
# Calibrated parameter values
l_imis_output_itt <- EDCimport::load_as_list(file = "outputs/calibration/itt/imis_output_itt.RData")
l_imis_output_pp <- EDCimport::load_as_list(file = "outputs/calibration/pp/imis_output_pp.RData")

# BNX scenario
l_params_bnx_itt <- load_all_params(
  file.init = "data/init_params.csv",
  file.pop_scaling = "data/pop_scaling.csv",
  file.init_dist = "data/init_dist_bnx.csv",
  file.cohort_balance = "data/cohort_balance.csv",
  file.mort = "data/all_cause_mortality.csv",
  file.death_hr = "data/death_hr.csv",
  file.weibull = "data/weibull_itt.csv",
  file.ce_tx = "data/ce_tx.csv",
  file.ce_death = "data/ce_death.csv",
  file.exit_dest = "data/unconditional_bnx_itt.csv",
  file.overdose = "data/overdose.csv",
  file.cali_params = "data/calibration/cali_priors.csv",
  file.cali_targets = "data/calibration/cali_targets.csv",
  file.fentanyl = "data/fentanyl.csv",
  file.naloxone = "data/naloxone.csv",
  file.qalys = "data/qalys.csv"
)

l_params_bnx_pp <- load_all_params(
  file.init = "data/init_params.csv",
  file.pop_scaling = "data/pop_scaling.csv",
  file.init_dist = "data/init_dist_bnx.csv",
  file.cohort_balance = "data/cohort_balance.csv",
  file.mort = "data/all_cause_mortality.csv",
  file.death_hr = "data/death_hr.csv",
  file.weibull = "data/weibull_pp.csv",
  file.ce_tx = "data/ce_tx.csv",
  file.ce_death = "data/ce_death.csv",
  file.exit_dest = "data/unconditional_bnx_pp.csv",
  file.overdose = "data/overdose.csv",
  file.cali_params = "data/calibration/cali_priors.csv",
  file.cali_targets = "data/calibration/cali_targets.csv",
  file.fentanyl = "data/fentanyl.csv",
  file.naloxone = "data/naloxone.csv",
  file.qalys = "data/qalys.csv"
)

# Methadone scenario
l_params_met_itt <- load_all_params(
  file.init = "data/init_params.csv",
  file.pop_scaling = "data/pop_scaling.csv",
  file.init_dist = "data/init_dist_met.csv",
  file.cohort_balance = "data/cohort_balance.csv",
  file.mort = "data/all_cause_mortality.csv",
  file.death_hr = "data/death_hr.csv",
  file.weibull = "data/weibull_itt.csv",
  file.ce_tx = "data/ce_tx.csv",
  file.ce_death = "data/ce_death.csv",
  file.exit_dest = "data/unconditional_met_itt.csv",
  file.overdose = "data/overdose.csv",
  file.cali_params = "data/calibration/cali_priors.csv",
  file.cali_targets = "data/calibration/cali_targets.csv",
  file.fentanyl = "data/fentanyl.csv",
  file.naloxone = "data/naloxone.csv",
  file.qalys = "data/qalys.csv"
)

l_params_met_pp <- load_all_params(
  file.init = "data/init_params.csv",
  file.pop_scaling = "data/pop_scaling.csv",
  file.init_dist = "data/init_dist_met.csv",
  file.cohort_balance = "data/cohort_balance.csv",
  file.mort = "data/all_cause_mortality.csv",
  file.death_hr = "data/death_hr.csv",
  file.weibull = "data/weibull_pp.csv",
  file.ce_tx = "data/ce_tx.csv",
  file.ce_death = "data/ce_death.csv",
  file.exit_dest = "data/unconditional_met_pp.csv",
  file.overdose = "data/overdose.csv",
  file.cali_params = "data/calibration/cali_priors.csv",
  file.cali_targets = "data/calibration/cali_targets.csv",
  file.fentanyl = "data/fentanyl.csv",
  file.naloxone = "data/naloxone.csv",
  file.qalys = "data/qalys.csv"
)

l_params_cali_pp <- load_all_params(
  file.init = "data/init_params.csv",
  file.pop_scaling = "data/pop_scaling.csv",
  file.init_dist = "data/calibration/init_dist_cali.csv",
  file.cohort_balance = "data/cohort_balance.csv",
  file.mort = "data/all_cause_mortality.csv",
  file.death_hr = "data/death_hr.csv",
  file.weibull = "data/weibull_pp.csv",
  file.ce_tx = "data/ce_tx.csv",
  file.ce_death = "data/ce_death.csv",
  file.exit_dest = "data/calibration/unconditional_cali_pp.csv",
  file.overdose = "data/overdose.csv",
  file.cali_params = "data/calibration/cali_priors.csv",
  file.cali_targets = "data/calibration/cali_targets.csv",
  file.fentanyl = "data/fentanyl.csv",
  file.naloxone = "data/naloxone.csv",
  file.qalys = "data/qalys.csv"
)

l_params_cali_itt <- load_all_params(
  file.init = "data/init_params.csv",
  file.pop_scaling = "data/pop_scaling.csv",
  file.init_dist = "data/calibration/init_dist_cali.csv",
  file.cohort_balance = "data/cohort_balance.csv",
  file.mort = "data/all_cause_mortality.csv",
  file.death_hr = "data/death_hr.csv",
  file.weibull = "data/weibull_itt.csv",
  file.ce_tx = "data/ce_tx.csv",
  file.ce_death = "data/ce_death.csv",
  file.exit_dest = "data/calibration/unconditional_cali_itt.csv",
  file.overdose = "data/overdose.csv",
  file.cali_params = "data/calibration/cali_priors.csv",
  file.cali_targets = "data/calibration/cali_targets.csv",
  file.fentanyl = "data/fentanyl.csv",
  file.naloxone = "data/naloxone.csv",
  file.qalys = "data/qalys.csv"
)
