rm(list = ls()) # to clean the workspace

library(dplyr) # to manipulate data
library(reshape2) # to transform data
library(ggplot2) # for nice looking plots
library(tidyverse)
# library(formattable)

# Call model setup functions
source("R/input_parameter_functions.R")
source("R/model_setup_functions.R")
source("R/mortality_function.R")
source("R/overdose_prob_function.R")
source("R/outcome_functions.R")

# Load parameters
source("Analysis/00_load_parameters.R") # load all model parameters for each scenario + calibrated parameters

#### Produce model outputs ####
#### Intent to treat ####
#### Propensity score ####
l_outcomes_bnx_itt_ps <- outcomes(l_params_all = l_params_bnx_itt, v_params_calib = l_imis_output_itt$v_calib_post_map, time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only")
l_outcomes_met_itt_ps <- outcomes(l_params_all = l_params_met_itt, v_params_calib = l_imis_output_itt$v_calib_post_map, time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only")

df_outcomes_itt_ps <- rbind(l_outcomes_bnx_itt_ps$df_outcomes, l_outcomes_met_itt_ps$df_outcomes)
rownames(df_outcomes_itt_ps) <- c("BNX", "Methadone")

# Incremental outcomes
l_inc_outcomes_itt_ps <- inc_outcomes(outcomes_comp = l_outcomes_met_itt_ps, outcomes_int = l_outcomes_bnx_itt_ps)

# Full model trace
write.csv(l_outcomes_met_itt_ps$m_m_trace, "outputs/trace/itt/full_trace_met_itt_ps.csv", row.names = TRUE)
write.csv(l_outcomes_bnx_itt_ps$m_m_trace, "outputs/trace/itt/full_trace_bnx_itt_ps.csv", row.names = TRUE)

# Aggregate trace
write.csv(l_outcomes_met_itt_ps$m_m_agg_trace, "outputs/trace/itt/agg_trace_met.csv", row.names = TRUE)
write.csv(l_outcomes_bnx_itt_ps$m_m_agg_trace, "outputs/trace/itt/agg_trace_bnx.csv", row.names = TRUE)

# Outcomes
# Totals by scenario
df_outcomes_itt_ps <- rbind(l_outcomes_bnx_itt_ps$df_outcomes, l_outcomes_met_itt_ps$df_outcomes)
rownames(df_outcomes_itt_ps) <- c("BNX", "Methadone")

# Incremental life-years
# Unscaled
df_incremental_itt_ps <- l_inc_outcomes_itt_ps$df_incremental
rownames(df_incremental_itt_ps) <- c("BNX vs. Methadone")
# Scaled
df_incremental_itt_ps_scaled <- l_inc_outcomes_itt_ps$df_incremental_scaled
rownames(df_incremental_itt_ps_scaled) <- c("BNX vs. Methadone")

# Output
save(df_incremental_itt_ps,
  file = "outputs/outcomes/itt/incremental_det_itt_ps.RData"
)
save(df_incremental_itt_ps_scaled,
  file = "outputs/outcomes/itt/incremental_det_itt_ps_scaled.RData"
)

write.csv(df_outcomes_itt_ps, "outputs/outcomes/itt/main_output_det_itt_ps.csv", row.names = TRUE)
write.csv(df_incremental_itt_ps, "outputs/outcomes/itt/incremental_det_itt_ps.csv", row.names = TRUE)
write.csv(df_incremental_itt_ps_scaled, "outputs/outcomes/itt/incremental_det_itt_ps_scaled.csv", row.names = TRUE)
