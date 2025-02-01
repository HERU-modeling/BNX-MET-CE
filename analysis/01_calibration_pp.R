rm(list = ls()) # to clean the workspace

#### Load packages, data and functions ####
#### Load packages and functions ####
library(dplyr) # to manipulate data
library(reshape2) # to transform data
library(ggplot2) # for nice looking plots
library(ggridges) # specialized ridge plots
library(tidyverse)
library(lhs)
library(IMIS)
library(grid)
library(gridExtra)
library(lattice)
library(parallel)
library(foreach)
library(doParallel)

source("R/input_parameter_functions.R")
source("R/model_setup_functions.R")
source("R/calibration_functions.R")
source("R/mortality_function.R")
source("R/overdose_prob_function.R")

# Load model inputs #
l_params_all <- load_all_params(
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

# Load calibration inputs #
v_cali_param_names <- c(
  "'Non-fentanyl overdose rate (OAT)'",
  "'Fentanyl overdose mult (prevalence)'",
  "'Fentanyl overdose mult (delta)'",
  "'Fatal overdose rate (OAT)'",
  "'Non-overdose mortality mult (OAT)'",
  "'Remain out-of-treatment (scale)'",
  "'Remain out-of-treatment (shape)'",
  "'Remain long-term abstinence (scale)'",
  "'Remain long-term abstinence (shape)'"
)

v_par1 <- c(
  n_oat_od_shape = l_params_all$n_oat_od_shape,
  n_fent_prev_od_mult_low = l_params_all$n_fent_prev_od_mult_low,
  n_fent_delta_od_mult_low = l_params_all$n_fent_delta_od_mult_low,
  n_fatal_od_oat_shape = l_params_all$n_fatal_od_oat_shape,
  hr_oat_shape = l_params_all$hr_oat_shape,
  p_weibull_scale_ou_shape = l_params_all$p_weibull_scale_ou_shape,
  p_weibull_shape_ou_shape = l_params_all$p_weibull_shape_ou_shape,
  p_weibull_scale_abs_shape = l_params_all$p_weibull_scale_abs_shape,
  p_weibull_shape_abs_shape = l_params_all$p_weibull_shape_abs_shape
)
v_par2 <- c(
  n_oat_od_scale = l_params_all$n_oat_od_scale,
  n_fent_prev_od_mult_high = l_params_all$n_fent_prev_od_mult_high,
  n_fent_delta_od_mult_high = l_params_all$n_fent_delta_od_mult_high,
  n_fatal_od_oat_scale = l_params_all$n_fatal_od_oat_scale,
  hr_oat_scale = l_params_all$hr_oat_scale,
  p_weibull_scale_ou_scale = l_params_all$p_weibull_scale_ou_scale,
  p_weibull_shape_ou_scale = l_params_all$p_weibull_shape_ou_scale,
  p_weibull_scale_abs_scale = l_params_all$p_weibull_scale_abs_scale,
  p_weibull_shape_abs_scale = l_params_all$p_weibull_shape_abs_scale
)

#### Load calibration targets ####
df_cali_targets <- read.csv(file = "data/calibration/cali_targets.csv", header = TRUE)

l_cali_targets <- list(
  df_odf = subset(df_cali_targets, target == "fatal overdoses"),
  df_dno = subset(df_cali_targets, target == "non-overdose deaths"),
  df_odn = subset(df_cali_targets, target == "non-fatal overdoses"),
  df_oot = subset(df_cali_targets, target == "time out of treatment (pp)"),
  df_abs = subset(df_cali_targets, target == "time in long-term abstinence")
)

### Number of calibration targets
v_target_names <- c("Fatal overdoses", "Non-overdose deaths", "Non-fatal overdoses", "Time out of treatment", "Time in long-term abstinence") # "Non-Overdose Deaths"
n_target <- length(v_target_names)

#### Specify calibration parameters ####
### Set seed
set.seed(3730687)

### Number of random samples to obtain from the posterior distribution
n_resamp <- 10000 # to match number of PSA draws

# Functions required by IMIS: prior(x), likelihood(x), sample.prior(n)
t_cali_start <- Sys.time()
#### Run IMIS algorithm ####
l_fit_imis <- IMIS::IMIS(
  B = 100, # n_samp = B*10 (was 100 incremental sample size at each iteration of IMIS)
  B.re = n_resamp, # "n_resamp" desired posterior sample size
  number_k = 2000, # maximum number of iterations in IMIS (originally 10)
  D = 0 # D originally set to zero, but try with other numbers (Raftery suggested 10)
) # originally 0
t_imis_end <- Sys.time()

# Unique parameter sets
n_unique <- length(unique(l_fit_imis$resample[, 1])) # 6392
# Effective sample size
n_ess <- round(sum(table(l_fit_imis$resample[, 1]))^2 / sum(table(l_fit_imis$resample[, 1])^2), 0) # 4496
# Max weight
n_max_wt <- max(table(l_fit_imis$resample[, 1])) / sum(table(l_fit_imis$resample[, 1])) # 0.001
# Calibration stats
df_cali_stats <- data.frame(n_unique, n_ess, n_max_wt)

### Obtain posterior
m_calib_post <- l_fit_imis$resample

### Compute posterior mean
v_calib_post_mean <- colMeans(m_calib_post)
### Compute posterior median and 95% credible interval
m_calib_post_95cr <- matrixStats::colQuantiles(m_calib_post,
  probs = c(0.025, 0.5, 0.975)
)
### Compute posterior values for draw
v_calib_post <- exp(log_post(m_calib_post))
### Compute maximum-a-posteriori (MAP) as the mode of the sampled values
v_calib_post_map <- m_calib_post[which.max(v_calib_post), ]
# Summary statistics
df_posterior_summ <- data.frame(
  Parameter = v_cali_param_names,
  Mean = v_calib_post_mean,
  m_calib_post_95cr,
  MAP = v_calib_post_map,
  check.names = FALSE
)

### Save calibration stats
## As .RData
# PP
save(df_cali_stats,
  m_calib_post,
  v_calib_post,
  v_calib_post_map,
  v_calib_post_mean,
  df_posterior_summ,
  file = "outputs/calibration/pp/imis_output_pp.RData"
)

## As .csv
write.csv(df_cali_stats,
  file = "outputs/calibration/pp/cali_stats_pp.csv",
  row.names = FALSE
)
write.csv(df_posterior_summ,
  file = "outputs/calibration/pp/summary_posterior_pp.csv",
  row.names = FALSE
)

#### Visualization of posterior distribution ####
### Rescale posterior to plot density of plots
v_calib_alpha <- scales::rescale(v_calib_post)

### Plot 10000 draws from the posterior with marginal histograms
png("plots/calibration/pp/posterior_distribution_marginal.png",
  width = 8, height = 6, units = "in", res = 300
)
psych::pairs.panels(m_calib_post)
dev.off()

#### START PLOTTING CODE HERE ####
#### Plot prior vs. posterior distribution for calibration parameters ####
# Load posterior
load(file = "outputs/calibration/pp/imis_output_pp.RData")

# Draw sample prior
m_calib_prior <- sample.prior(n_samp = n_resamp)

# Prepare data
df_samp_prior <- melt(
  cbind(
    Distribution = "Prior",
    as.data.frame(m_calib_prior[1:n_resamp, ])
  ),
  variable.name = "Parameter"
)
df_samp_post_imis <- melt(
  cbind(
    Distribution = "Posterior",
    as.data.frame(m_calib_post[1:n_resamp, ])
  ),
  variable.name = "Parameter"
)
df_calib_prior_post <- rbind(df_samp_prior, df_samp_post_imis)
df_calib_prior_post$Distribution <- ordered(df_calib_prior_post$Distribution,
  levels = c(
    "Prior",
    "Posterior"
  )
)
df_calib_prior_post$Parameter <- factor(df_calib_prior_post$Parameter,
  levels = levels(df_calib_prior_post$Parameter),
  ordered = TRUE,
  labels = v_cali_param_names
)

### Plot priors and IMIS posteriors
# TO-DO: Add vertical lines for prior mean and MAP
prior_v_posterior <- ggplot(
  df_calib_prior_post,
  aes(x = value, y = after_stat(density), fill = Distribution)
) +
  facet_wrap(~Parameter,
    scales = "free",
    ncol = 3,
    labeller = label_parsed
  ) +
  geom_density(alpha = 0.5) +
  theme_bw(base_size = 16) +
  # scale_color_manual(values = c("#0080C6", "#FFC20E")) +
  guides(
    fill = guide_legend(title = "", order = 1),
    linetype = guide_legend(title = "", order = 2),
    color = guide_legend(title = "", order = 2)
  ) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background = element_rect(
      fill = "white",
      color = "white"
    ),
    strip.text = element_text(hjust = 0)
  )
ggsave(prior_v_posterior,
  filename = "plots/calibration/pp/prior-v-posterior.png",
  width = 10, height = 7
)
t_cali_end <- Sys.time()
#### Plot model fit against calibration targets ####
# Set number of cores
n_cores <- parallel::detectCores() - 1
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)

# Run blockwise
n_runs <- n_resamp
n_block_size <- 1000 # size of block for each loop
n_blocks <- n_runs / n_block_size # to run entire set
n_start <- 0

l_cali_outcomes_odf <- list()
l_cali_outcomes_dno <- list()
l_cali_outcomes_odn <- list()
l_cali_outcomes_oot <- list()
l_cali_outcomes_abs <- list()

for (j in (0:(n_blocks - 1))) {
  l_cali_target_fit <- foreach(i = (n_start + 1 + j * n_block_size):(n_start + (j + 1) * n_block_size), .combine = combine_custom_cali, .packages = "tidyr") %dopar% {
    l_model_target_fit <- calibration_out(
      v_params_calib = m_calib_post[i, ],
      l_params_all = l_params_all
    )
    m_model_targets_odf <- l_model_target_fit$fatal_overdose
    m_model_targets_dno <- l_model_target_fit$death_non_od
    m_model_targets_odn <- l_model_target_fit$overdose
    m_model_targets_oot <- l_model_target_fit$time_oot
    m_model_targets_abs <- l_model_target_fit$time_abs

    return(list(
      m_model_targets_odf = m_model_targets_odf,
      m_model_targets_dno = m_model_targets_dno,
      m_model_targets_odn = m_model_targets_odn,
      m_model_targets_oot = m_model_targets_oot,
      m_model_targets_abs = m_model_targets_abs
    ))
  }

  m_cali_outcomes_odf <- l_cali_target_fit$m_model_targets_odf
  m_cali_outcomes_dno <- l_cali_target_fit$m_model_targets_dno
  m_cali_outcomes_odn <- l_cali_target_fit$m_model_targets_odn
  m_cali_outcomes_oot <- l_cali_target_fit$m_model_targets_oot
  m_cali_outcomes_abs <- l_cali_target_fit$m_model_targets_abs

  l_cali_outcomes_odf[[j + 1]] <- m_cali_outcomes_odf
  l_cali_outcomes_dno[[j + 1]] <- m_cali_outcomes_dno
  l_cali_outcomes_odn[[j + 1]] <- m_cali_outcomes_odn
  l_cali_outcomes_oot[[j + 1]] <- m_cali_outcomes_oot
  l_cali_outcomes_abs[[j + 1]] <- m_cali_outcomes_abs

  out <- paste0("Block ", (j + 1), " of ", n_blocks, " complete.") # Status
  print(out)
}

# combine output data sets
# initialize empty matrices
m_outcomes_odf_comb <- m_outcomes_dno_comb <- m_outcomes_odn_comb <- m_outcomes_oot_comb <- m_outcomes_abs_comb <- matrix(0, nrow = 0, ncol = 9)

for (i in 1:n_blocks) {
  m_temp <- l_cali_outcomes_odf[[i]]
  m_outcomes_odf_comb <- rbind(m_outcomes_odf_comb, m_temp)
}
for (i in 1:n_blocks) {
  m_temp <- l_cali_outcomes_dno[[i]]
  m_outcomes_dno_comb <- rbind(m_outcomes_dno_comb, m_temp)
}
for (i in 1:n_blocks) {
  m_temp <- l_cali_outcomes_odn[[i]]
  m_outcomes_odn_comb <- rbind(m_outcomes_odn_comb, m_temp)
}
for (i in 1:n_blocks) {
  m_temp <- l_cali_outcomes_oot[[i]]
  m_outcomes_oot_comb <- rbind(m_outcomes_oot_comb, m_temp)
}
for (i in 1:n_blocks) {
  m_temp <- l_cali_outcomes_abs[[i]]
  m_outcomes_abs_comb <- rbind(m_outcomes_abs_comb, m_temp)
}

stopImplicitCluster()

Sys.time()
## As .RData
save(m_outcomes_odf_comb,
  file = "outputs/calibration/pp/model_targets_odf_pp.RData"
)
save(m_outcomes_dno_comb,
  file = "outputs/calibration/pp/model_targets_dno_pp.RData"
)
save(m_outcomes_odn_comb,
  file = "outputs/calibration/pp/model_targets_odn_pp.RData"
)
save(m_outcomes_oot_comb,
  file = "outputs/calibration/pp/model_targets_oot_pp.RData"
)
save(m_outcomes_abs_comb,
  file = "outputs/calibration/pp/model_targets_abs_pp.RData"
)

# Load data
load(file = "outputs/calibration/pp/model_targets_odf_pp.RData")
load(file = "outputs/calibration/pp/model_targets_dno_pp.RData")
load(file = "outputs/calibration/pp/model_targets_odn_pp.RData")
load(file = "outputs/calibration/pp/model_targets_oot_pp.RData")
load(file = "outputs/calibration/pp/model_targets_abs_pp.RData")

# Model outputs
m_outcomes_odf_comb_stats <- cbind(
  matrixStats::colQuantiles(m_outcomes_odf_comb,
    probs = c(0.025, 0.5, 0.975)
  ),
  matrixStats::colMeans2(m_outcomes_odf_comb)
)
m_outcomes_dno_comb_stats <- cbind(
  matrixStats::colQuantiles(m_outcomes_dno_comb,
    probs = c(0.025, 0.5, 0.975)
  ),
  matrixStats::colMeans2(m_outcomes_dno_comb)
)
m_outcomes_odn_comb_stats <- cbind(
  matrixStats::colQuantiles(m_outcomes_odn_comb,
    probs = c(0.025, 0.5, 0.975)
  ),
  matrixStats::colMeans2(m_outcomes_odn_comb)
)
m_outcomes_oot_comb_stats <- cbind(
  matrixStats::colQuantiles(m_outcomes_oot_comb,
    probs = c(0.025, 0.5, 0.975)
  ),
  matrixStats::colMeans2(m_outcomes_oot_comb)
)
m_outcomes_abs_comb_stats <- cbind(
  matrixStats::colQuantiles(m_outcomes_abs_comb,
    probs = c(0.025, 0.5, 0.975)
  ),
  matrixStats::colMeans2(m_outcomes_abs_comb)
)

m_time <- matrix(l_cali_targets$df_odf$time)
# m_year <- matrix(l_cali_targets$df_odf$year)
m_pop <- matrix(l_cali_targets$df_odf$pop)
m_outcomes_odf_comb_fit <- cbind(m_outcomes_odf_comb_stats, m_time, m_pop)
m_outcomes_dno_comb_fit <- cbind(m_outcomes_dno_comb_stats, m_time, m_pop)
m_outcomes_odn_comb_fit <- cbind(m_outcomes_odn_comb_stats, m_time, m_pop)
m_outcomes_oot_comb_fit <- cbind(m_outcomes_oot_comb_stats, m_time, m_pop)
m_outcomes_abs_comb_fit <- cbind(m_outcomes_abs_comb_stats, m_time, m_pop)

df_model_targets_odf_fit <- m_outcomes_odf_comb_fit %>%
  as_tibble() %>%
  setNames(c("ci_low", "median", "ci_high", "pe", "time", "pop")) %>%
  mutate(
    group = "Model output (95% CI)",
    num = pe * m_pop,
    low = ci_low * m_pop,
    high = ci_high * m_pop
  ) %>%
  select(time, group, num, low, high)

df_model_targets_dno_fit <- m_outcomes_dno_comb_fit %>%
  as_tibble() %>%
  setNames(c("ci_low", "median", "ci_high", "pe", "time", "pop")) %>%
  mutate(
    group = "Model output (95% CI)",
    num = pe * m_pop,
    low = ci_low * m_pop,
    high = ci_high * m_pop
  ) %>%
  select(time, group, num, low, high)

df_model_targets_odn_fit <- m_outcomes_odn_comb_fit %>%
  as_tibble() %>%
  setNames(c("ci_low", "median", "ci_high", "pe", "time", "pop")) %>%
  mutate(
    group = "Model output (95% CI)",
    num = pe * m_pop,
    low = ci_low * m_pop,
    high = ci_high * m_pop
  ) %>%
  select(time, group, num, low, high)

df_model_targets_oot_fit <- m_outcomes_oot_comb_fit %>%
  as_tibble() %>%
  setNames(c("ci_low", "median", "ci_high", "pe", "time", "pop")) %>%
  mutate(
    group = "Model output (95% CI)",
    num = pe,
    low = ci_low,
    high = ci_high
  ) %>%
  select(time, group, num, low, high)

df_model_targets_abs_fit <- m_outcomes_abs_comb_fit %>%
  as_tibble() %>%
  setNames(c("ci_low", "median", "ci_high", "pe", "time", "pop")) %>%
  mutate(
    group = "Model output (95% CI)",
    num = pe,
    low = ci_low,
    high = ci_high
  ) %>%
  select(time, group, num, low, high)

# Targets
df_targets_odf <- l_cali_targets$df_odf %>%
  as_tibble() %>%
  mutate(
    group = "Cali target (95% CI)",
    num = pe * m_pop,
    low = low * m_pop,
    high = high * m_pop
  ) %>%
  select(time, group, num, low, high)

df_targets_dno <- l_cali_targets$df_dno %>%
  as_tibble() %>%
  mutate(
    group = "Cali target (95% CI)",
    num = pe * m_pop,
    low = low * m_pop,
    high = high * m_pop
  ) %>%
  select(time, group, num, low, high)

df_targets_odn <- l_cali_targets$df_odn %>%
  as_tibble() %>%
  mutate(
    group = "Cali target (95% CI)",
    num = pe * m_pop,
    low = low * m_pop,
    high = high * m_pop
  ) %>%
  select(time, group, num, low, high)

df_targets_oot <- l_cali_targets$df_oot %>%
  as_tibble() %>%
  mutate(
    group = "Cali target (95% CI)",
    num = pe,
    low = low,
    high = high
  ) %>%
  select(time, group, num, low, high)

df_targets_abs <- l_cali_targets$df_abs %>%
  as_tibble() %>%
  mutate(
    group = "Cali target (95% CI)",
    num = pe,
    low = low,
    high = high
  ) %>%
  select(time, group, num, low, high)

# Combine
df_fit_odf <- bind_rows(df_targets_odf, df_model_targets_odf_fit)
df_fit_dno <- bind_rows(df_targets_dno, df_model_targets_dno_fit)
df_fit_odn <- bind_rows(df_targets_odn, df_model_targets_odn_fit)
df_fit_oot <- bind_rows(df_targets_oot, df_model_targets_oot_fit)
df_fit_abs <- bind_rows(df_targets_abs, df_model_targets_abs_fit)

# Plot fit vs. targets
# Fatal overdose
p_temp_odf <- ggplot(df_fit_odf, aes(x = time, y = num, group = group, color = group)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = low, ymax = high),
    width = .5,
    position = position_dodge(0.05)
  )

plot_fit_odf <- p_temp_odf + labs(title = NULL, x = "Year", y = "Fatal overdoses") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(legend.title = element_blank()) +
  theme(
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 7, color = "black")
  ) +
  scale_color_manual(values = c("#696158", "#A6192E")) +
  scale_x_continuous(
    breaks = l_cali_targets$df_odf$time,
    labels = c("2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020*")
  ) +
  ylim(0, 210)

# Non-overdose deaths
p_temp_dno <- ggplot(df_fit_dno, aes(x = time, y = num, group = group, color = group)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = low, ymax = high),
    width = .5,
    position = position_dodge(0.05)
  )

plot_fit_dno <- p_temp_dno + labs(title = NULL, x = "Year", y = "Non-Overdose Deaths") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(legend.title = element_blank()) +
  theme(
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 7, color = "black")
  ) +
  scale_color_manual(values = c("#696158", "#A6192E")) +
  scale_x_continuous(
    breaks = l_cali_targets$df_dno$time,
    labels = c("2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020*")
  ) +
  ylim(0, 125)

# Non-fatal overdose
p_temp_odn <- ggplot(df_fit_odn, aes(x = time, y = num, group = group, color = group)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = low, ymax = high),
    width = .5,
    position = position_dodge(0.05)
  )

plot_fit_odn <- p_temp_odn + labs(title = NULL, x = "Year", y = "Non-fatal overdoses") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(legend.title = element_blank()) +
  theme(
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 7, color = "black")
  ) +
  scale_color_manual(values = c("#696158", "#A6192E")) +
  scale_x_continuous(
    breaks = l_cali_targets$df_odf$time,
    labels = c("2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020*")
  ) +
  ylim(0, 7500)

# Time out of treatment
p_temp_oot <- ggplot(df_fit_oot, aes(x = time, y = num, group = group, color = group)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = low, ymax = high),
    width = .5,
    position = position_dodge(0.05)
  )

plot_fit_oot <- p_temp_oot + labs(title = NULL, x = "Year", y = "Time spent out of treatment") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(legend.title = element_blank()) +
  theme(
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 7, color = "black")
  ) +
  scale_color_manual(values = c("#696158", "#A6192E")) +
  scale_x_continuous(
    breaks = l_cali_targets$df_odf$time,
    labels = c("2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020*")
  ) +
  ylim(0.25, 0.75)

# Time in abstinence
p_temp_abs <- ggplot(df_fit_abs, aes(x = time, y = num, group = group, color = group)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = low, ymax = high),
    width = .5,
    position = position_dodge(0.05)
  )

plot_fit_abs <- p_temp_abs + labs(title = NULL, x = "Year", y = "Time long-term abstinence") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(legend.title = element_blank()) +
  theme(
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 7, color = "black")
  ) +
  scale_color_manual(values = c("#696158", "#A6192E")) +
  scale_x_continuous(
    breaks = l_cali_targets$df_odf$time,
    labels = c("2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020*")
  ) +
  ylim(0, 0.1)

# plot_fit_odn
# Plot for extracting legend only
plot_fit_odf_leg <- p_temp_odf + labs(title = NULL, x = "Year", y = "Fatal overdoses") +
  theme_classic() +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = c("#696158", "#A6192E")) +
  scale_x_continuous(
    breaks = l_cali_targets$df_odf$time,
    labels = c("2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020*")
  )
# Code to extract legend from plots
g_legend <- function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

mylegend <- g_legend(plot_fit_odf_leg)

# plot_fit_comb <- grid.arrange(arrangeGrob(plot_fit_odf, plot_fit_dno, plot_fit_oot, nrow = 3),
#   mylegend,
#   nrow = 2, heights = c(6, .5)
# )

plot_fit_comb_all <- grid.arrange(arrangeGrob(plot_fit_odf, plot_fit_dno, plot_fit_odn, plot_fit_oot, plot_fit_abs, mylegend, nrow = 3, ncol = 2),
  # mylegend,
  # nrow = 2,
  heights = 6, widths = 7
)

# plot_fit_comb
# Outputs
# ggsave(plot_fit_odf,
#   filename = "plots/calibration/pp/target-fit-odf.png",
#   width = 4, height = 4
# )
# ggsave(plot_fit_dno,
#   filename = "plots/calibration/pp/target-fit-dno.png",
#   width = 4, height = 4
# )
# ggsave(plot_fit_odn,
#   filename = "plots/calibration/pp/target-fit-odn.png",
#   width = 4, height = 4
# )
# ggsave(plot_fit_oot,
#   filename = "plots/calibration/pp/target-fit-oot.png",
#   width = 4, height = 4
# )
# ggsave(plot_fit_abs,
#   filename = "plots/calibration/pp/target-fit-abs.png",
#   width = 4, height = 4
# )
# ggsave(plot_fit_comb,
#   filename = "plots/calibration/pp/target-fit-comb.png",
#   width = 6, height = 8
# )
ggsave(plot_fit_comb_all,
  filename = "plots/calibration/pp/target-fit-comb-all.png",
  width = 6, height = 6
)
t_cali_fit_end <- Sys.time()

t_imis_run <- difftime(t_imis_end, t_cali_start, units = "hours")
t_cali_run <- difftime(t_cali_end, t_cali_start, units = "hours")
t_cali_fit_run <- difftime(t_cali_fit_end, t_cali_end, units = "hours")
t_full_run <- difftime(t_cali_fit_end, t_cali_start, units = "hours")

df_run_time <- rbind(t_imis_run, t_cali_run, t_cali_fit_run, t_full_run)
colnames(df_run_time) <- c("time (hr)")

write.csv(df_run_time,
  file = "outputs/calibration/pp/cali_run_time_pp.csv",
  row.names = TRUE
)
