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
# n_pop_cohort <- l_params_bnx_itt$n_pop_oat # update with cohort size
# n_pop_est <- l_params_bnx_itt$n_pop_est
n_sim <- 10000 # just to test function (will be set as n_sim)

### Load PSA results
load(file = "outputs/psa/outcomes_psa.RData")

##############
### ITT-PS ###
##############
# Life years
tbl_df_inc_ly_psa_itt_ps <- df_incremental_psa_itt_ps_scaled %>%
  as_tibble() %>%
  mutate(
    n_ly_scaled_2010 = n_inc_qalys_adj_2010_scaled,
    n_ly_scaled_2012 = n_inc_qalys_adj_2012_scaled,
    n_ly_scaled_2014 = n_inc_qalys_adj_2014_scaled,
    n_ly_scaled_2016 = n_inc_qalys_adj_2016_scaled,
    n_ly_scaled_2018 = n_inc_qalys_adj_2018_scaled,
    n_ly_scaled_2020 = n_inc_qalys_adj_2020_scaled
  ) %>%
  mutate(
    sim_ly_0_2020 = ifelse((n_ly_scaled_2020 <= 0), 1, 0),
    sim_ly_500_2020 = ifelse((n_ly_scaled_2020 <= -500), 1, 0),
    sim_ly_1000_2020 = ifelse((n_ly_scaled_2020 <= -1000), 1, 0),
    sim_ly_1500_2020 = ifelse((n_ly_scaled_2020 <= -1500), 1, 0),
    sim_ly_2000_2020 = ifelse((n_ly_scaled_2020 <= -2000), 1, 0),
    sim_ly_2500_2020 = ifelse((n_ly_scaled_2020 <= -2500), 1, 0),
    sim_ly_3000_2020 = ifelse((n_ly_scaled_2020 <= -3000), 1, 0)
  ) %>%
  select(sim_ly_0_2020, sim_ly_500_2020, sim_ly_1000_2020, sim_ly_1500_2020, sim_ly_2000_2020, sim_ly_2500_2020, sim_ly_3000_2020) %>%
  summarise_all(mean)

##########
### PP ###
##########
# Life years
tbl_df_inc_ly_psa_pp <- df_incremental_psa_pp_scaled %>%
  as_tibble() %>%
  mutate(
    n_ly_scaled_2010 = n_inc_qalys_adj_2010_scaled,
    n_ly_scaled_2012 = n_inc_qalys_adj_2012_scaled,
    n_ly_scaled_2014 = n_inc_qalys_adj_2014_scaled,
    n_ly_scaled_2016 = n_inc_qalys_adj_2016_scaled,
    n_ly_scaled_2018 = n_inc_qalys_adj_2018_scaled,
    n_ly_scaled_2020 = n_inc_qalys_adj_2020_scaled
  ) %>%
  mutate(
    sim_ly_0_2020 = ifelse((n_ly_scaled_2020 <= 0), 1, 0),
    sim_ly_500_2020 = ifelse((n_ly_scaled_2020 <= -500), 1, 0),
    sim_ly_1000_2020 = ifelse((n_ly_scaled_2020 <= -1000), 1, 0),
    sim_ly_1500_2020 = ifelse((n_ly_scaled_2020 <= -1500), 1, 0),
    sim_ly_2000_2020 = ifelse((n_ly_scaled_2020 <= -2000), 1, 0),
    sim_ly_2500_2020 = ifelse((n_ly_scaled_2020 <= -2500), 1, 0),
    sim_ly_3000_2020 = ifelse((n_ly_scaled_2020 <= -3000), 1, 0)
  ) %>%
  select(sim_ly_0_2020, sim_ly_500_2020, sim_ly_1000_2020, sim_ly_1500_2020, sim_ly_2000_2020, sim_ly_2500_2020, sim_ly_3000_2020) %>%
  summarise_all(mean)

######################
### PP (high dose) ###
######################
# Life years
tbl_df_inc_ly_psa_pp_hd <- df_incremental_psa_pp_hd_scaled %>%
  as_tibble() %>%
  mutate(
    n_ly_scaled_2010 = n_inc_qalys_adj_2010_scaled,
    n_ly_scaled_2012 = n_inc_qalys_adj_2012_scaled,
    n_ly_scaled_2014 = n_inc_qalys_adj_2014_scaled,
    n_ly_scaled_2016 = n_inc_qalys_adj_2016_scaled,
    n_ly_scaled_2018 = n_inc_qalys_adj_2018_scaled,
    n_ly_scaled_2020 = n_inc_qalys_adj_2020_scaled
  ) %>%
  mutate(
    sim_ly_0_2020 = ifelse((n_ly_scaled_2020 <= 0), 1, 0),
    sim_ly_500_2020 = ifelse((n_ly_scaled_2020 <= -500), 1, 0),
    sim_ly_1000_2020 = ifelse((n_ly_scaled_2020 <= -1000), 1, 0),
    sim_ly_1500_2020 = ifelse((n_ly_scaled_2020 <= -1500), 1, 0),
    sim_ly_2000_2020 = ifelse((n_ly_scaled_2020 <= -2000), 1, 0),
    sim_ly_2500_2020 = ifelse((n_ly_scaled_2020 <= -2500), 1, 0),
    sim_ly_3000_2020 = ifelse((n_ly_scaled_2020 <= -3000), 1, 0)
  ) %>%
  select(sim_ly_0_2020, sim_ly_500_2020, sim_ly_1000_2020, sim_ly_1500_2020, sim_ly_2000_2020, sim_ly_2500_2020, sim_ly_3000_2020) %>%
  summarise_all(mean)

tbl_df_inc_ly_psa <- rbind(tbl_df_inc_ly_psa_itt_ps, tbl_df_inc_ly_psa_pp, tbl_df_inc_ly_psa_pp_hd)
row_names <- c("ITT-PS", "PP", "PP-HD")
rownames(tbl_df_inc_ly_psa) <- row_names
write.csv(tbl_df_inc_ly_psa,
  file = "outputs/psa/perc_sim_ly.csv",
  row.names = TRUE
)
