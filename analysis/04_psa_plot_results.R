library(dplyr)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(tidyr)

# Call model setup functions
source("R/input_parameter_functions.R")

# Load parameters
source("Analysis/00_load_parameters.R")

# Set population size for dirichlet draws
n_pop_cohort <- l_params_bnx_itt$n_pop_oat # update with cohort size
# n_pop_est <- l_params_bnx_itt$n_pop_est
# n_sim <- 10000 # just to test function (will be set as n_sim)

### Process PSA results
load(file = "outputs/psa/outcomes_psa_itt_ps.RData")
# load(file = "outputs/psa/outcomes_psa_itt_rr_sa.RData")
load(file = "outputs/psa/outcomes_psa_itt_rr_100_inc_sa.RData")
load(file = "outputs/psa/outcomes_psa_itt_rr_0_inc_sa.RData")
load(file = "outputs/psa/outcomes_psa_pp.RData")
load(file = "outputs/psa/outcomes_psa_pp_hd.RData")

### Summary stats ###
# Incremental
##############
### ITT-PS ###
##############
# Life years
tbl_df_summ_inc_ly_psa_itt_ps <- df_incremental_psa_itt_ps_scaled %>%
  as_tibble() %>%
  select(
    n_inc_qalys_adj_2010_scaled,
    n_inc_qalys_adj_2012_scaled,
    n_inc_qalys_adj_2014_scaled,
    n_inc_qalys_adj_2016_scaled,
    n_inc_qalys_adj_2018_scaled,
    n_inc_qalys_adj_2020_scaled
  ) %>%
  gather(
    key = "variable",
    value = "value"
  ) %>%
  group_by(variable) %>%
  dplyr::summarize(
    mean = mean(value),
    sd = sd(value),
    q50 = quantile(value, probs = .5),
    q025 = quantile(value, probs = .025),
    q975 = quantile(value, probs = .975),
    min = min(value),
    max = max(value)
  ) %>%
  mutate(
    time = ifelse(variable == "n_inc_qalys_adj_2010_scaled", 2010,
      ifelse(variable == "n_inc_qalys_adj_2012_scaled", 2012,
        ifelse(variable == "n_inc_qalys_adj_2014_scaled", 2014,
          ifelse(variable == "n_inc_qalys_adj_2016_scaled", 2016,
            ifelse(variable == "n_inc_qalys_adj_2018_scaled", 2018,
              ifelse(variable == "n_inc_qalys_adj_2020_scaled", 2020, NA)
            )
          )
        )
      )
    ),
    scenario = "Initiator analysis (base-case)"
  ) %>%
  arrange(time) %>%
  select(scenario, time, mean, q025, q975)

# Life years (BNX-ITT)
tbl_df_summ_bnx_ly_psa_itt_ps <- df_outcomes_bnx_itt_ps_scaled %>%
  as_tibble() %>%
  select(
    n_qalys_2010_ann_scaled,
    n_qalys_2012_ann_scaled,
    n_qalys_2014_ann_scaled,
    n_qalys_2016_ann_scaled,
    n_qalys_2018_ann_scaled,
    n_qalys_2020_ann_scaled
  ) %>%
  gather("variable", "value") %>%
  group_by(variable) %>%
  dplyr::summarize(
    mean = mean(value),
    sd = sd(value),
    q50 = quantile(value, probs = .5),
    q025 = quantile(value, probs = .025),
    q975 = quantile(value, probs = .975),
    min = min(value),
    max = max(value)
  ) %>%
  mutate(
    time = ifelse(variable == "n_qalys_2010_ann_scaled", 2010,
      ifelse(variable == "n_qalys_2012_ann_scaled", 2012,
        ifelse(variable == "n_qalys_2014_ann_scaled", 2014,
          ifelse(variable == "n_qalys_2016_ann_scaled", 2016,
            ifelse(variable == "n_qalys_2018_ann_scaled", 2018,
              ifelse(variable == "n_qalys_2020_ann_scaled", 2020, NA)
            )
          )
        )
      )
    ),
    scenario = "Initiator analysis (base-case)"
  ) %>%
  arrange(time) %>%
  select(scenario, time, mean, q025, q975)

# Life years (MET-ITT)
tbl_df_summ_met_ly_psa_itt_ps <- df_outcomes_met_itt_ps_scaled %>%
  as_tibble() %>%
  select(
    n_qalys_2010_ann_scaled,
    n_qalys_2012_ann_scaled,
    n_qalys_2014_ann_scaled,
    n_qalys_2016_ann_scaled,
    n_qalys_2018_ann_scaled,
    n_qalys_2020_ann_scaled
  ) %>%
  gather("variable", "value") %>%
  group_by(variable) %>%
  dplyr::summarize(
    mean = mean(value),
    sd = sd(value),
    q50 = quantile(value, probs = .5),
    q025 = quantile(value, probs = .025),
    q975 = quantile(value, probs = .975),
    min = min(value),
    max = max(value)
  ) %>%
  mutate(
    time = ifelse(variable == "n_qalys_2010_ann_scaled", 2010,
      ifelse(variable == "n_qalys_2012_ann_scaled", 2012,
        ifelse(variable == "n_qalys_2014_ann_scaled", 2014,
          ifelse(variable == "n_qalys_2016_ann_scaled", 2016,
            ifelse(variable == "n_qalys_2018_ann_scaled", 2018,
              ifelse(variable == "n_qalys_2020_ann_scaled", 2020, NA)
            )
          )
        )
      )
    ),
    scenario = "Initiator analysis (base-case)"
  ) %>%
  arrange(time) %>%
  select(scenario, time, mean, q025, q975)

# Fatal overdoses
tbl_df_summ_inc_odf_psa_itt_ps <- df_incremental_psa_itt_ps_scaled %>%
  as_tibble() %>%
  select(
    n_inc_odf_adj_2010_scaled,
    n_inc_odf_adj_2012_scaled,
    n_inc_odf_adj_2014_scaled,
    n_inc_odf_adj_2016_scaled,
    n_inc_odf_adj_2018_scaled,
    n_inc_odf_adj_2020_scaled
  ) %>%
  gather("variable", "value") %>%
  group_by(variable) %>%
  dplyr::summarize(
    mean = mean(value),
    sd = sd(value),
    q50 = quantile(value, probs = .5),
    q025 = quantile(value, probs = .025),
    q975 = quantile(value, probs = .975),
    min = min(value),
    max = max(value)
  ) %>%
  mutate(
    time = ifelse(variable == "n_inc_odf_adj_2010_scaled", 2010,
      ifelse(variable == "n_inc_odf_adj_2012_scaled", 2012,
        ifelse(variable == "n_inc_odf_adj_2014_scaled", 2014,
          ifelse(variable == "n_inc_odf_adj_2016_scaled", 2016,
            ifelse(variable == "n_inc_odf_adj_2018_scaled", 2018,
              ifelse(variable == "n_inc_odf_adj_2020_scaled", 2020, NA)
            )
          )
        )
      )
    ),
    scenario = "Initiator analysis (base-case)"
  ) %>%
  arrange(time) %>%
  select(scenario, time, mean, q025, q975)

# Non-fatal overdoses
tbl_df_summ_inc_odn_psa_itt_ps <- df_incremental_psa_itt_ps_scaled %>%
  as_tibble() %>%
  select(
    n_inc_odn_adj_2010_scaled,
    n_inc_odn_adj_2012_scaled,
    n_inc_odn_adj_2014_scaled,
    n_inc_odn_adj_2016_scaled,
    n_inc_odn_adj_2018_scaled,
    n_inc_odn_adj_2020_scaled
  ) %>%
  gather("variable", "value") %>%
  group_by(variable) %>%
  dplyr::summarize(
    mean = mean(value),
    sd = sd(value),
    q50 = quantile(value, probs = .5),
    q025 = quantile(value, probs = .025),
    q975 = quantile(value, probs = .975),
    min = min(value),
    max = max(value)
  ) %>%
  mutate(
    time = ifelse(variable == "n_inc_odn_adj_2010_scaled", 2010,
      ifelse(variable == "n_inc_odn_adj_2012_scaled", 2012,
        ifelse(variable == "n_inc_odn_adj_2014_scaled", 2014,
          ifelse(variable == "n_inc_odn_adj_2016_scaled", 2016,
            ifelse(variable == "n_inc_odn_adj_2018_scaled", 2018,
              ifelse(variable == "n_inc_odn_adj_2020_scaled", 2020, NA)
            )
          )
        )
      )
    ),
    scenario = "Initiator analysis (base-case)"
  ) %>%
  arrange(time) %>%
  select(scenario, time, mean, q025, q975)

# All-cause deaths
tbl_df_summ_inc_acm_psa_itt_ps <- df_incremental_psa_itt_ps_scaled %>%
  as_tibble() %>%
  select(
    n_inc_acm_adj_2010_scaled,
    n_inc_acm_adj_2012_scaled,
    n_inc_acm_adj_2014_scaled,
    n_inc_acm_adj_2016_scaled,
    n_inc_acm_adj_2018_scaled,
    n_inc_acm_adj_2020_scaled
  ) %>%
  gather("variable", "value") %>%
  group_by(variable) %>%
  dplyr::summarize(
    mean = mean(value),
    sd = sd(value),
    q50 = quantile(value, probs = .5),
    q025 = quantile(value, probs = .025),
    q975 = quantile(value, probs = .975),
    min = min(value),
    max = max(value)
  ) %>%
  mutate(
    time = ifelse(variable == "n_inc_acm_adj_2010_scaled", 2010,
      ifelse(variable == "n_inc_acm_adj_2012_scaled", 2012,
        ifelse(variable == "n_inc_acm_adj_2014_scaled", 2014,
          ifelse(variable == "n_inc_acm_adj_2016_scaled", 2016,
            ifelse(variable == "n_inc_acm_adj_2018_scaled", 2018,
              ifelse(variable == "n_inc_acm_adj_2020_scaled", 2020, NA)
            )
          )
        )
      )
    ),
    scenario = "Initiator analysis (base-case)"
  ) %>%
  arrange(time) %>%
  select(scenario, time, mean, q025, q975)

##################
### ITT-R&R SA ###
##################
# 20K runs
# Life years
# tbl_df_summ_inc_ly_psa_itt_rr_sa <- df_incremental_psa_itt_rr_sa_scaled %>%
#   as_tibble() %>%
#   select(
#     n_inc_qalys_adj_2010_scaled,
#     n_inc_qalys_adj_2012_scaled,
#     n_inc_qalys_adj_2014_scaled,
#     n_inc_qalys_adj_2016_scaled,
#     n_inc_qalys_adj_2018_scaled,
#     n_inc_qalys_adj_2020_scaled
#   ) %>%
#   gather("variable", "value") %>%
#   group_by(variable) %>%
#   dplyr::summarize(
#     mean = mean(value),
#     sd = sd(value),
#     q50 = quantile(value, probs = .5),
#     q025 = quantile(value, probs = .025),
#     q975 = quantile(value, probs = .975),
#     min = min(value),
#     max = max(value)
#   ) %>%
#   mutate(
#     time = ifelse(variable == "n_inc_qalys_adj_2010_scaled", 2010,
#       ifelse(variable == "n_inc_qalys_adj_2012_scaled", 2012,
#         ifelse(variable == "n_inc_qalys_adj_2014_scaled", 2014,
#           ifelse(variable == "n_inc_qalys_adj_2016_scaled", 2016,
#             ifelse(variable == "n_inc_qalys_adj_2018_scaled", 2018,
#               ifelse(variable == "n_inc_qalys_adj_2020_scaled", 2020, NA)
#             )
#           )
#         )
#       )
#     ),
#     scenario = "Initiator analysis (20,000 run SA)"
#   ) %>%
#   arrange(time) %>%
#   select(scenario, time, mean, q025, q975)

# 100% incident OAT
# Life years
tbl_df_summ_inc_ly_psa_itt_rr_100_inc_sa <- df_incremental_psa_itt_rr_100_inc_sa_scaled %>%
  as_tibble() %>%
  select(
    n_inc_qalys_adj_2010_scaled,
    n_inc_qalys_adj_2012_scaled,
    n_inc_qalys_adj_2014_scaled,
    n_inc_qalys_adj_2016_scaled,
    n_inc_qalys_adj_2018_scaled,
    n_inc_qalys_adj_2020_scaled
  ) %>%
  gather("variable", "value") %>%
  group_by(variable) %>%
  dplyr::summarize(
    mean = mean(value),
    sd = sd(value),
    q50 = quantile(value, probs = .5),
    q025 = quantile(value, probs = .025),
    q975 = quantile(value, probs = .975),
    min = min(value),
    max = max(value)
  ) %>%
  mutate(
    time = ifelse(variable == "n_inc_qalys_adj_2010_scaled", 2010,
      ifelse(variable == "n_inc_qalys_adj_2012_scaled", 2012,
        ifelse(variable == "n_inc_qalys_adj_2014_scaled", 2014,
          ifelse(variable == "n_inc_qalys_adj_2016_scaled", 2016,
            ifelse(variable == "n_inc_qalys_adj_2018_scaled", 2018,
              ifelse(variable == "n_inc_qalys_adj_2020_scaled", 2020, NA)
            )
          )
        )
      )
    ),
    scenario = "Initiator analysis (100% incident OAT SA)"
  ) %>%
  arrange(time) %>%
  select(scenario, time, mean, q025, q975)

# 0% incident OAT
# Life years
tbl_df_summ_inc_ly_psa_itt_rr_0_inc_sa <- df_incremental_psa_itt_rr_0_inc_sa_scaled %>%
  as_tibble() %>%
  select(
    n_inc_qalys_adj_2010_scaled,
    n_inc_qalys_adj_2012_scaled,
    n_inc_qalys_adj_2014_scaled,
    n_inc_qalys_adj_2016_scaled,
    n_inc_qalys_adj_2018_scaled,
    n_inc_qalys_adj_2020_scaled
  ) %>%
  gather("variable", "value") %>%
  group_by(variable) %>%
  dplyr::summarize(
    mean = mean(value),
    sd = sd(value),
    q50 = quantile(value, probs = .5),
    q025 = quantile(value, probs = .025),
    q975 = quantile(value, probs = .975),
    min = min(value),
    max = max(value)
  ) %>%
  mutate(
    time = ifelse(variable == "n_inc_qalys_adj_2010_scaled", 2010,
      ifelse(variable == "n_inc_qalys_adj_2012_scaled", 2012,
        ifelse(variable == "n_inc_qalys_adj_2014_scaled", 2014,
          ifelse(variable == "n_inc_qalys_adj_2016_scaled", 2016,
            ifelse(variable == "n_inc_qalys_adj_2018_scaled", 2018,
              ifelse(variable == "n_inc_qalys_adj_2020_scaled", 2020, NA)
            )
          )
        )
      )
    ),
    scenario = "Initiator analysis (0% incident OAT SA)"
  ) %>%
  arrange(time) %>%
  select(scenario, time, mean, q025, q975)

##########
### PP ###
##########
# Life years
tbl_df_summ_inc_ly_psa_pp <- df_incremental_psa_pp_scaled %>%
  as_tibble() %>%
  select(
    n_inc_qalys_adj_2010_scaled,
    n_inc_qalys_adj_2012_scaled,
    n_inc_qalys_adj_2014_scaled,
    n_inc_qalys_adj_2016_scaled,
    n_inc_qalys_adj_2018_scaled,
    n_inc_qalys_adj_2020_scaled
  ) %>%
  gather("variable", "value") %>%
  group_by(variable) %>%
  dplyr::summarize(
    mean = mean(value),
    sd = sd(value),
    q50 = quantile(value, probs = .5),
    q025 = quantile(value, probs = .025),
    q975 = quantile(value, probs = .975),
    min = min(value),
    max = max(value)
  ) %>%
  mutate(
    time = ifelse(variable == "n_inc_qalys_adj_2010_scaled", 2010,
      ifelse(variable == "n_inc_qalys_adj_2012_scaled", 2012,
        ifelse(variable == "n_inc_qalys_adj_2014_scaled", 2014,
          ifelse(variable == "n_inc_qalys_adj_2016_scaled", 2016,
            ifelse(variable == "n_inc_qalys_adj_2018_scaled", 2018,
              ifelse(variable == "n_inc_qalys_adj_2020_scaled", 2020, NA)
            )
          )
        )
      )
    ),
    scenario = "Per-protocol analysis (alternate-case)"
  ) %>%
  arrange(time) %>%
  select(scenario, time, mean, q025, q975)

# Fatal overdoses
tbl_df_summ_inc_odf_psa_pp <- df_incremental_psa_pp_scaled %>%
  as_tibble() %>%
  select(
    n_inc_odf_adj_2010_scaled,
    n_inc_odf_adj_2012_scaled,
    n_inc_odf_adj_2014_scaled,
    n_inc_odf_adj_2016_scaled,
    n_inc_odf_adj_2018_scaled,
    n_inc_odf_adj_2020_scaled
  ) %>%
  gather("variable", "value") %>%
  group_by(variable) %>%
  dplyr::summarize(
    mean = mean(value),
    sd = sd(value),
    q50 = quantile(value, probs = .5),
    q025 = quantile(value, probs = .025),
    q975 = quantile(value, probs = .975),
    min = min(value),
    max = max(value)
  ) %>%
  mutate(
    time = ifelse(variable == "n_inc_odf_adj_2010_scaled", 2010,
      ifelse(variable == "n_inc_odf_adj_2012_scaled", 2012,
        ifelse(variable == "n_inc_odf_adj_2014_scaled", 2014,
          ifelse(variable == "n_inc_odf_adj_2016_scaled", 2016,
            ifelse(variable == "n_inc_odf_adj_2018_scaled", 2018,
              ifelse(variable == "n_inc_odf_adj_2020_scaled", 2020, NA)
            )
          )
        )
      )
    ),
    scenario = "Per-protocol analysis (alternate-case)"
  ) %>%
  arrange(time) %>%
  select(scenario, time, mean, q025, q975)

# Non-fatal overdoses
tbl_df_summ_inc_odn_psa_pp <- df_incremental_psa_pp_scaled %>%
  as_tibble() %>%
  select(
    n_inc_odn_adj_2010_scaled,
    n_inc_odn_adj_2012_scaled,
    n_inc_odn_adj_2014_scaled,
    n_inc_odn_adj_2016_scaled,
    n_inc_odn_adj_2018_scaled,
    n_inc_odn_adj_2020_scaled
  ) %>%
  gather("variable", "value") %>%
  group_by(variable) %>%
  dplyr::summarize(
    mean = mean(value),
    sd = sd(value),
    q50 = quantile(value, probs = .5),
    q025 = quantile(value, probs = .025),
    q975 = quantile(value, probs = .975),
    min = min(value),
    max = max(value)
  ) %>%
  mutate(
    time = ifelse(variable == "n_inc_odn_adj_2010_scaled", 2010,
      ifelse(variable == "n_inc_odn_adj_2012_scaled", 2012,
        ifelse(variable == "n_inc_odn_adj_2014_scaled", 2014,
          ifelse(variable == "n_inc_odn_adj_2016_scaled", 2016,
            ifelse(variable == "n_inc_odn_adj_2018_scaled", 2018,
              ifelse(variable == "n_inc_odn_adj_2020_scaled", 2020, NA)
            )
          )
        )
      )
    ),
    scenario = "Per-protocol analysis (alternate-case)"
  ) %>%
  arrange(time) %>%
  select(scenario, time, mean, q025, q975)

# All-cause deaths
tbl_df_summ_inc_acm_psa_pp <- df_incremental_psa_pp_scaled %>%
  as_tibble() %>%
  select(
    n_inc_acm_adj_2010_scaled,
    n_inc_acm_adj_2012_scaled,
    n_inc_acm_adj_2014_scaled,
    n_inc_acm_adj_2016_scaled,
    n_inc_acm_adj_2018_scaled,
    n_inc_acm_adj_2020_scaled
  ) %>%
  gather("variable", "value") %>%
  group_by(variable) %>%
  dplyr::summarize(
    mean = mean(value),
    sd = sd(value),
    q50 = quantile(value, probs = .5),
    q025 = quantile(value, probs = .025),
    q975 = quantile(value, probs = .975),
    min = min(value),
    max = max(value)
  ) %>%
  mutate(
    time = ifelse(variable == "n_inc_acm_adj_2010_scaled", 2010,
      ifelse(variable == "n_inc_acm_adj_2012_scaled", 2012,
        ifelse(variable == "n_inc_acm_adj_2014_scaled", 2014,
          ifelse(variable == "n_inc_acm_adj_2016_scaled", 2016,
            ifelse(variable == "n_inc_acm_adj_2018_scaled", 2018,
              ifelse(variable == "n_inc_acm_adj_2020_scaled", 2020, NA)
            )
          )
        )
      )
    ),
    scenario = "Per-protocol analysis (alternate-case)"
  ) %>%
  arrange(time) %>%
  select(scenario, time, mean, q025, q975)

######################
### PP (high dose) ###
######################
# Life years
tbl_df_summ_inc_ly_psa_pp_hd <- df_incremental_psa_pp_hd_scaled %>%
  as_tibble() %>%
  select(
    n_inc_qalys_adj_2010_scaled,
    n_inc_qalys_adj_2012_scaled,
    n_inc_qalys_adj_2014_scaled,
    n_inc_qalys_adj_2016_scaled,
    n_inc_qalys_adj_2018_scaled,
    n_inc_qalys_adj_2020_scaled
  ) %>%
  gather("variable", "value") %>%
  group_by(variable) %>%
  dplyr::summarize(
    mean = mean(value),
    sd = sd(value),
    q50 = quantile(value, probs = .5),
    q025 = quantile(value, probs = .025),
    q975 = quantile(value, probs = .975),
    min = min(value),
    max = max(value)
  ) %>%
  mutate(
    time = ifelse(variable == "n_inc_qalys_adj_2010_scaled", 2010,
      ifelse(variable == "n_inc_qalys_adj_2012_scaled", 2012,
        ifelse(variable == "n_inc_qalys_adj_2014_scaled", 2014,
          ifelse(variable == "n_inc_qalys_adj_2016_scaled", 2016,
            ifelse(variable == "n_inc_qalys_adj_2018_scaled", 2018,
              ifelse(variable == "n_inc_qalys_adj_2020_scaled", 2020, NA)
            )
          )
        )
      )
    ),
    scenario = "High-dose (sensitivity analysis)"
  ) %>%
  arrange(time) %>%
  select(scenario, time, mean, q025, q975)

# Fatal overdoses
tbl_df_summ_inc_odf_psa_pp_hd <- df_incremental_psa_pp_hd_scaled %>%
  as_tibble() %>%
  select(
    n_inc_odf_adj_2010_scaled,
    n_inc_odf_adj_2012_scaled,
    n_inc_odf_adj_2014_scaled,
    n_inc_odf_adj_2016_scaled,
    n_inc_odf_adj_2018_scaled,
    n_inc_odf_adj_2020_scaled
  ) %>%
  gather("variable", "value") %>%
  group_by(variable) %>%
  dplyr::summarize(
    mean = mean(value),
    sd = sd(value),
    q50 = quantile(value, probs = .5),
    q025 = quantile(value, probs = .025),
    q975 = quantile(value, probs = .975),
    min = min(value),
    max = max(value)
  ) %>%
  mutate(
    time = ifelse(variable == "n_inc_odf_adj_2010_scaled", 2010,
      ifelse(variable == "n_inc_odf_adj_2012_scaled", 2012,
        ifelse(variable == "n_inc_odf_adj_2014_scaled", 2014,
          ifelse(variable == "n_inc_odf_adj_2016_scaled", 2016,
            ifelse(variable == "n_inc_odf_adj_2018_scaled", 2018,
              ifelse(variable == "n_inc_odf_adj_2020_scaled", 2020, NA)
            )
          )
        )
      )
    ),
    scenario = "High-dose (sensitivity analysis)"
  ) %>%
  arrange(time) %>%
  select(scenario, time, mean, q025, q975)

# Non-fatal overdoses
tbl_df_summ_inc_odn_psa_pp_hd <- df_incremental_psa_pp_hd_scaled %>%
  as_tibble() %>%
  select(
    n_inc_odn_adj_2010_scaled,
    n_inc_odn_adj_2012_scaled,
    n_inc_odn_adj_2014_scaled,
    n_inc_odn_adj_2016_scaled,
    n_inc_odn_adj_2018_scaled,
    n_inc_odn_adj_2020_scaled
  ) %>%
  gather("variable", "value") %>%
  group_by(variable) %>%
  dplyr::summarize(
    mean = mean(value),
    sd = sd(value),
    q50 = quantile(value, probs = .5),
    q025 = quantile(value, probs = .025),
    q975 = quantile(value, probs = .975),
    min = min(value),
    max = max(value)
  ) %>%
  mutate(
    time = ifelse(variable == "n_inc_odn_adj_2010_scaled", 2010,
      ifelse(variable == "n_inc_odn_adj_2012_scaled", 2012,
        ifelse(variable == "n_inc_odn_adj_2014_scaled", 2014,
          ifelse(variable == "n_inc_odn_adj_2016_scaled", 2016,
            ifelse(variable == "n_inc_odn_adj_2018_scaled", 2018,
              ifelse(variable == "n_inc_odn_adj_2020_scaled", 2020, NA)
            )
          )
        )
      )
    ),
    scenario = "High-dose (sensitivity analysis)"
  ) %>%
  arrange(time) %>%
  select(scenario, time, mean, q025, q975)

# All-cause deaths
tbl_df_summ_inc_acm_psa_pp_hd <- df_incremental_psa_pp_hd_scaled %>%
  as_tibble() %>%
  select(
    n_inc_acm_adj_2010_scaled,
    n_inc_acm_adj_2012_scaled,
    n_inc_acm_adj_2014_scaled,
    n_inc_acm_adj_2016_scaled,
    n_inc_acm_adj_2018_scaled,
    n_inc_acm_adj_2020_scaled
  ) %>%
  gather("variable", "value") %>%
  group_by(variable) %>%
  dplyr::summarize(
    mean = mean(value),
    sd = sd(value),
    q50 = quantile(value, probs = .5),
    q025 = quantile(value, probs = .025),
    q975 = quantile(value, probs = .975),
    min = min(value),
    max = max(value)
  ) %>%
  mutate(
    time = ifelse(variable == "n_inc_acm_adj_2010_scaled", 2010,
      ifelse(variable == "n_inc_acm_adj_2012_scaled", 2012,
        ifelse(variable == "n_inc_acm_adj_2014_scaled", 2014,
          ifelse(variable == "n_inc_acm_adj_2016_scaled", 2016,
            ifelse(variable == "n_inc_acm_adj_2018_scaled", 2018,
              ifelse(variable == "n_inc_acm_adj_2020_scaled", 2020, NA)
            )
          )
        )
      )
    ),
    scenario = "High-dose (sensitivity analysis)"
  ) %>%
  arrange(time) %>%
  select(scenario, time, mean, q025, q975)

# Combine
tbl_df_summ_inc_ly_psa_comb <- rbind(tbl_df_summ_inc_ly_psa_itt_ps, tbl_df_summ_inc_ly_psa_pp, tbl_df_summ_inc_ly_psa_pp_hd)
tbl_df_summ_inc_odf_psa_comb <- rbind(tbl_df_summ_inc_odf_psa_itt_ps, tbl_df_summ_inc_odf_psa_pp, tbl_df_summ_inc_odf_psa_pp_hd)
tbl_df_summ_inc_odn_psa_comb <- rbind(tbl_df_summ_inc_odn_psa_itt_ps, tbl_df_summ_inc_odn_psa_pp, tbl_df_summ_inc_odn_psa_pp_hd)
tbl_df_summ_inc_acm_psa_comb <- rbind(tbl_df_summ_inc_acm_psa_itt_ps, tbl_df_summ_inc_acm_psa_pp, tbl_df_summ_inc_acm_psa_pp_hd)

tbl_df_summ_inc_ly_psa_comb_sa <- rbind(tbl_df_summ_inc_ly_psa_itt_ps, tbl_df_summ_inc_ly_psa_itt_rr_sa, tbl_df_summ_inc_ly_psa_pp, tbl_df_summ_inc_ly_psa_pp_hd)
tbl_df_summ_inc_odf_psa_comb_sa <- rbind(tbl_df_summ_inc_odf_psa_itt_ps, tbl_df_summ_inc_odf_psa_itt_rr_sa, tbl_df_summ_inc_odf_psa_pp, tbl_df_summ_inc_odf_psa_pp_hd)
tbl_df_summ_inc_odn_psa_comb_sa <- rbind(tbl_df_summ_inc_odn_psa_itt_ps, tbl_df_summ_inc_odn_psa_itt_rr_sa, tbl_df_summ_inc_odn_psa_pp, tbl_df_summ_inc_odn_psa_pp_hd)
tbl_df_summ_inc_acm_psa_comb_sa <- rbind(tbl_df_summ_inc_acm_psa_itt_ps, tbl_df_summ_inc_acm_psa_itt_rr_sa, tbl_df_summ_inc_acm_psa_pp, tbl_df_summ_inc_acm_psa_pp_hd)

tbl_df_summ_ly_psa_comb_sa <- rbind(tbl_df_summ_bnx_ly_psa_itt_ps, tbl_df_summ_met_ly_psa_itt_ps)

tbl_df_summ_ly_psa_comb_inc_prev_sa <- rbind(tbl_df_summ_inc_ly_psa_itt_ps, tbl_df_summ_inc_ly_psa_itt_rr_100_inc_sa, tbl_df_summ_inc_ly_psa_itt_rr_0_inc_sa)

## As .csv ####
write.csv(tbl_df_summ_inc_ly_psa_comb,
  file = "outputs/psa/summary_incremental_ly_psa_combined.csv",
  row.names = FALSE
)
write.csv(tbl_df_summ_inc_odf_psa_comb,
  file = "outputs/psa/summary_incremental_odf_psa_combined.csv",
  row.names = FALSE
)
write.csv(tbl_df_summ_inc_odn_psa_comb,
  file = "outputs/psa/summary_incremental_odn_psa_combined.csv",
  row.names = FALSE
)
write.csv(tbl_df_summ_inc_acm_psa_comb,
  file = "outputs/psa/summary_incremental_acm_psa_combined.csv",
  row.names = FALSE
)

############################
### Life-years-lost plot ###
############################
# Output scaled up to cohort
tbl_df_summ_inc_ly_psa_comb$scenario <- factor(
  tbl_df_summ_inc_ly_psa_comb$scenario,
  levels = c(
    "Initiator analysis (base-case)",
    "Per-protocol analysis (alternate-case)",
    "High-dose (sensitivity analysis)"
  )
)
# Main plot
plot_psa_ly_scaled <- ggplot(tbl_df_summ_inc_ly_psa_comb) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 0.75) + # NEW
  geom_pointrange(
    aes(x = time, y = mean, ymin = q025, ymax = q975, color = scenario, linetype = scenario),
    position = position_dodge(width = .75),
    linewidth = .75
  ) +
  scale_color_manual(values = c("#FC4C02", "#005778", "#008E97")) +
  scale_linetype_manual(values = c(
    "Initiator (primary analysis)" = "solid",
    "Per-protocol (sensitivity analysis)" = "solid",
    "High-dose (sensitivity analysis)" = "solid"
  )) +
  labs(y = "Incremental life years (BNX vs. methadone)", x = "Year") +
  scale_x_continuous(breaks = c(2010, 2012, 2014, 2016, 2018, 2020), labels = c("2010", "2012", "2014", "2016", "2018", "2020*"), limits = c(2009, 2021)) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),
    legend.key = element_rect(fill = "transparent", colour = "transparent"),
    plot.title = element_text(hjust = 0.02, vjust = -7),
    legend.position = "bottom",
    legend.title = element_blank(),
    text = element_text(size = 15)
  ) +
  guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))

ggsave(plot_psa_ly_scaled,
  filename = "plots/psa/psa-life-years-lost_scaled.png",
  width = 8, height = 6, dpi = 350
)

# Sensitivity analysis plot
# tbl_df_summ_inc_ly_psa_comb_sa$scenario <- factor(
#   tbl_df_summ_inc_ly_psa_comb_sa$scenario,
#   levels = c(
#     "Initiator analysis (base-case)",
#     "Initiator analysis (20,000 run SA)",
#     "Per-protocol analysis (alternate-case)",
#     "High-dose (sensitivity analysis)"
#   )
# )
# plot_psa_ly_scaled_sa <- ggplot(tbl_df_summ_inc_ly_psa_comb_sa) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 0.75) + # NEW
#   geom_pointrange(
#     aes(x = time, y = mean, ymin = q025, ymax = q975, color = scenario, linetype = scenario),
#     position = position_dodge(width = .75),
#     linewidth = .75
#   ) +
#   scale_color_manual(values = c("#FC4C02", "#C8102E", "#005778", "#008E97")) +
#   scale_linetype_manual(values = c(
#     "Initiator analysis (base-case)" = "solid",
#     "Initiator analysis (20,000 run SA)" = "solid",
#     "Per-protocol analysis (alternate-case)" = "solid",
#     "High-dose (sensitivity analysis)" = "solid"
#   )) +
#   labs(y = "Incremental life years (BNX vs. methadone)", x = "Year") +
#   scale_x_continuous(breaks = c(2010, 2012, 2014, 2016, 2018, 2020), labels = c("2010", "2012", "2014", "2016", "2018", "2020*"), limits = c(2009, 2021)) +
#   theme(
#     panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),
#     legend.key = element_rect(fill = "transparent", colour = "transparent"),
#     plot.title = element_text(hjust = 0.02, vjust = -7),
#     legend.position = "bottom",
#     legend.title = element_blank(),
#     text = element_text(size = 15)
#   ) +
#   guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))

# ggsave(plot_psa_ly_scaled_sa,
#   filename = "plots/psa/psa-life-years-lost_scaled-SA.png",
#   width = 8, height = 6, dpi = 350
# )

# Incident prevalence sensitivity analysis plot
tbl_df_summ_inc_ly_psa_comb_inc_prev_sa$scenario <- factor(
  tbl_df_summ_inc_ly_psa_comb_inc_prev_sa$scenario,
  levels = c(
    "Initiator (primary analysis)",
    "Initiator (100% Incident OAT SA)",
    "Initiator (100% Experienced OAT SA)"
  )
)
plot_psa_ly_scaled_inc_prev_sa <- ggplot(tbl_df_summ_inc_ly_psa_comb_inc_prev_sa) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 0.75) + # NEW
  geom_pointrange(
    aes(x = time, y = mean, ymin = q025, ymax = q975, color = scenario, linetype = scenario),
    position = position_dodge(width = .75),
    linewidth = .75
  ) +
  scale_color_manual(values = c("#FC4C02", "#005778", "#008E97")) +
  scale_linetype_manual(values = c(
    "Initiator (primary analysis)" = "solid",
    "Initiator (100% Incident OAT SA)" = "solid",
    "Initiator (100% Experienced OAT SA)" = "solid"
  )) +
  labs(y = "Incremental life years (BNX vs. methadone)", x = "Year") +
  scale_x_continuous(breaks = c(2010, 2012, 2014, 2016, 2018, 2020), labels = c("2010", "2012", "2014", "2016", "2018", "2020*"), limits = c(2009, 2021)) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),
    legend.key = element_rect(fill = "transparent", colour = "transparent"),
    plot.title = element_text(hjust = 0.02, vjust = -7),
    legend.position = "bottom",
    legend.title = element_blank(),
    text = element_text(size = 15)
  ) +
  guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))

ggsave(plot_psa_ly_scaled_inc_prev_sa,
  filename = "plots/psa/psa-life-years-lost_scaled-inc-prev-SA.png",
  width = 8, height = 6, dpi = 350
)

# Total life years plot
tbl_df_summ_ly_psa_comb_sa$scenario <- factor(
  tbl_df_summ_ly_psa_comb_sa$scenario,
  levels = c(
    "BNX",
    "Methadone"
  )
)
plot_psa_ly_total_scaled_sa <- ggplot(tbl_df_summ_ly_psa_comb_sa) +
  geom_pointrange(
    aes(x = time, y = mean, ymin = q025, ymax = q975, color = scenario, linetype = scenario),
    position = position_dodge(width = .75),
    linewidth = .75
  ) +
  scale_color_manual(values = c("#008E97", "#005778")) +
  scale_linetype_manual(values = c(
    "BNX" = "solid",
    "Methadone" = "solid"
  )) +
  labs(y = "Total life years", x = "Year") +
  scale_x_continuous(breaks = c(2010, 2012, 2014, 2016, 2018, 2020), labels = c("2010", "2012", "2014", "2016", "2018", "2020*"), limits = c(2009, 2021)) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),
    legend.key = element_rect(fill = "transparent", colour = "transparent"),
    plot.title = element_text(hjust = 0.02, vjust = -7),
    legend.position = "bottom",
    legend.title = element_blank(),
    text = element_text(size = 15)
  ) +
  guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))

###########################
### Fatal overdose plot ###
###########################
# Output scaled up to cohort
plot_psa_odf_scaled <- ggplot(tbl_df_summ_inc_odf_psa_comb) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = .75) + # NEW
  geom_pointrange(
    aes(x = time, y = mean, ymin = q025, ymax = q975, color = scenario),
    position = position_dodge(width = .75),
    linewidth = .75
  ) +
  scale_color_manual(values = c("#FC4C02", "#005778", "#008E97")) + # NEW
  labs(y = "Fatal overdoses (BNX vs. methadone)", x = "Year") +
  scale_x_continuous(breaks = c(2010, 2012, 2014, 2016, 2018, 2020), labels = c("2010", "2012", "2014", "2016", "2018", "2020*"), limits = c(2009, 2021)) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),
    legend.key = element_rect(fill = "transparent", colour = "transparent"),
    plot.title = element_text(hjust = 0.02, vjust = -7),
    legend.position = "right",
    legend.title = element_blank(),
    text = element_text(size = 15)
  )

ggsave(plot_psa_odf_scaled,
  filename = "plots/psa/psa-fatal-overdose_scaled.png",
  width = 8, height = 6, dpi = 350
)

###############################
### Non-fatal overdose plot ###
###############################
# Output scaled up to cohort
plot_psa_odn_scaled <- ggplot(tbl_df_summ_inc_odn_psa_comb) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = .75) + # NEW
  geom_pointrange(
    aes(x = time, y = mean, ymin = q025, ymax = q975, color = scenario),
    position = position_dodge(width = .75),
    linewidth = .75
  ) +
  scale_color_manual(values = c("#FC4C02", "#005778", "#008E97")) + # NEW
  labs(y = "Non-fatal overdoses (BNX vs. methadone)", x = "Year") +
  scale_x_continuous(breaks = c(2010, 2012, 2014, 2016, 2018, 2020), labels = c("2010", "2012", "2014", "2016", "2018", "2020*"), limits = c(2009, 2021)) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),
    legend.key = element_rect(fill = "transparent", colour = "transparent"),
    plot.title = element_text(hjust = 0.02, vjust = -7),
    legend.position = "right",
    legend.title = element_blank(),
    text = element_text(size = 15)
  )

ggsave(plot_psa_odn_scaled,
  filename = "plots/psa/psa-non-fatal-overdose_scaled.png",
  width = 8, height = 6, dpi = 350
)

#############################
### All-cause deaths plot ###
#############################
# Output scaled up to cohort
plot_psa_acm_scaled <- ggplot(tbl_df_summ_inc_acm_psa_comb) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = .75) + # NEW
  geom_pointrange(
    aes(x = time, y = mean, ymin = q025, ymax = q975, color = scenario),
    position = position_dodge(width = .75),
    linewidth = .75
  ) +
  scale_color_manual(values = c("#FC4C02", "#005778", "#008E97")) + # NEW
  labs(y = "All-cause deaths (BNX vs. methadone)", x = "Year") +
  scale_x_continuous(breaks = c(2010, 2012, 2014, 2016, 2018, 2020), labels = c("2010", "2012", "2014", "2016", "2018", "2020*"), limits = c(2009, 2021)) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),
    legend.key = element_rect(fill = "transparent", colour = "transparent"),
    plot.title = element_text(hjust = 0.02, vjust = -7),
    legend.position = "right",
    legend.title = element_blank(),
    text = element_text(size = 15)
  )

ggsave(plot_psa_acm_scaled,
  filename = "plots/psa/psa-all-cause-deaths_scaled.png",
  width = 8, height = 6, dpi = 350
)
