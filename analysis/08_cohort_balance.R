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

#### Load targets ####
df_cohort_balance <- read.csv(file = "data/calibration/cohort_balance.csv", header = TRUE)

# Rename columns in df_cohort_balance
df_cohort_balance <- df_cohort_balance %>%
  dplyr::rename("Incident user (obs)" = p_inc, "Experienced user (obs)" = p_prev) %>%
  pivot_longer(cols = c("Incident user (obs)", "Experienced user (obs)"), names_to = "state", values_to = "proportion")

# df_cohort_balance <- df_cohort_balance# %>%
# rename("Incident user (obs)" = p_inc, "Experienced user (obs)" = p_prev) # %>%
# pivot_longer(cols = c("Incident user (obs)", "Experienced user (obs)"), names_to = "state", values_to = "proportion")

# l_params_bnx_itt <- update_param_list(l_params_all = l_params_bnx_itt, params_updated = l_imis_output_itt$v_calib_post_mean)
l_params_met_itt <- update_param_list(l_params_all = l_params_met_itt, params_updated = l_imis_output_itt$v_calib_post_mean)
# l_params_cali_itt <- update_param_list(l_params_all = l_params_cali_itt, params_updated = l_imis_output_itt$v_calib_post_mean)

# Run model
# l_out_markov_bnx_itt_ps <- markov_model(l_params_all = l_params_bnx_itt, err_stop = FALSE, verbose = TRUE, checks = FALSE, time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "bnx_only")
l_out_markov_met_itt_ps <- markov_model(l_params_all = l_params_met_itt, err_stop = FALSE, verbose = TRUE, checks = FALSE, time_horizon = "full", ce_est = "itt_ps", analytic_cohort = "met_only")
# l_out_markov_cali_itt <- markov_model(l_params_all = l_params_cali_itt, err_stop = FALSE, verbose = TRUE, checks = FALSE, time_horizon = "cali", ce_est = "cali", analytic_cohort = "cali")

## Cohort balance
df_m_agg_cohort_balance_trace <- as.data.frame(l_out_markov_met_itt_ps$m_m_cohort_balance_trace)
df_m_agg_cohort_balance_trace$time <- as.numeric(rownames(df_m_agg_cohort_balance_trace))
df_m_agg_cohort_balance_trace$year <- df_m_agg_cohort_balance_trace$time / 52
df_m_agg_cohort_balance_trace$total <- df_m_agg_cohort_balance_trace$`Incident user` + df_m_agg_cohort_balance_trace$`Experienced user`
df_m_agg_cohort_balance_trace <- df_m_agg_cohort_balance_trace %>%
  mutate("Incident user (model)" = `Incident user` / total, "Experienced user (model)" = `Experienced user` / total)

df_m_agg_cohort_balance_trace_plot <- df_m_agg_cohort_balance_trace %>% gather(state, proportion, "Incident user (model)", "Experienced user (model)") # health states to plot
df_m_agg_cohort_balance_trace_plot <- df_m_agg_cohort_balance_trace_plot %>%
  filter(year == as.integer(year)) %>%
  # mutate(state = ifelse(state == "Incident user", "Incident user (model)", "Experienced user (model)")) %>%
  select(year, state, proportion)

# Merge the data frames on 'year'
df_merged <- rbind(df_m_agg_cohort_balance_trace_plot, df_cohort_balance)

# Plot
p_cohort_balance <- ggplot(data = df_merged, aes(x = year, y = proportion, color = state)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(
    values = c(
      "Incident user (model)" = "#4575b4",
      "Incident user (obs)" = "#91bfdb",
      "Experienced user (model)" = "#d73027",
      "Experienced user (obs)" = "#fc8d59"
    ),
    name = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  ) +
  guides(color = guide_legend(nrow = 2)) +
  labs(x = "Year", y = "Proportion") +
  scale_x_continuous(
    breaks = seq(0, 10, 1),
    labels = c("2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020*")
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.1),
    labels = scales::percent
  )

ggsave(p_cohort_balance,
  filename = "plots/markov trace/cohort_balance.png",
  width = 8, height = 6, units = "in", dpi = 350
)
