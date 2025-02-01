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
n_sim <- 10000 # just to test function (will be set as n_sim)

### Process PSA results
load(file = "outputs/psa/outcomes_psa.RData")

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
    scenario = "Initiator analysis"
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
    scenario = "Initiator analysis"
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
    scenario = "Initiator analysis"
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
    scenario = "Initiator analysis"
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
    scenario = "Per-protocol sensitivity analysis"
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
    scenario = "Per-protocol sensitivity analysis"
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
    scenario = "Per-protocol sensitivity analysis"
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
    scenario = "Per-protocol sensitivity analysis"
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
    scenario = "Per-protocol sensitivity analysis (high-dose)"
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
    scenario = "Per-protocol sensitivity analysis (high-dose)"
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
    scenario = "Per-protocol sensitivity analysis (high-dose)"
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
    scenario = "Per-protocol sensitivity analysis (high-dose)"
  ) %>%
  arrange(time) %>%
  select(scenario, time, mean, q025, q975)

# Combine
tbl_df_summ_inc_ly_psa_comb <- rbind(tbl_df_summ_inc_ly_psa_itt_ps, tbl_df_summ_inc_ly_psa_pp, tbl_df_summ_inc_ly_psa_pp_hd)
tbl_df_summ_inc_odf_psa_comb <- rbind(tbl_df_summ_inc_odf_psa_itt_ps, tbl_df_summ_inc_odf_psa_pp, tbl_df_summ_inc_odf_psa_pp_hd)
tbl_df_summ_inc_odn_psa_comb <- rbind(tbl_df_summ_inc_odn_psa_itt_ps, tbl_df_summ_inc_odn_psa_pp, tbl_df_summ_inc_odn_psa_pp_hd)
tbl_df_summ_inc_acm_psa_comb <- rbind(tbl_df_summ_inc_acm_psa_itt_ps, tbl_df_summ_inc_acm_psa_pp, tbl_df_summ_inc_acm_psa_pp_hd)

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
plot_psa_ly_scaled <- ggplot(tbl_df_summ_inc_ly_psa_comb) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 0.75) + # NEW
  geom_pointrange(
    aes(x = time, y = mean, ymin = q025, ymax = q975, color = scenario, linetype = scenario),
    position = position_dodge(width = .75),
    linewidth = .75
  ) +
  scale_color_manual(values = c("#FC4C02", "#005778", "#008E97")) +
  # scale_color_manual(values = c("Initiator analysis (full population)" = "#081d58", "Per-protocol SA (full population)" = "#41b6c4", "Per-protocol SA (high-dose)" = "#c7e9b4")) +
  scale_linetype_manual(values = c("Initiator analysis" = "solid", "Per-protocol sensitivity analysis" = "solid", "Per-protocol sensitivity analysis (high-dose)" = "solid")) +
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
