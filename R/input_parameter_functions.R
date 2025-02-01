#' Load mortality data
#'
#' \code{load_mort_data} is used to load age-specific mortality from .csv file
#' into vector.
#'
#' @param file String with the location and name of the file with mortality
#' data.
#' @return
#' A vector with mortality by age.
#' @export
load_mort_params <- function(file.mort = NULL, n_male) {
  df_lt_bc <- read.csv(file = file.mort)
  v_r_mort_by_age_male <- df_lt_bc %>%
    select(male) %>%
    as.matrix()
  v_r_mort_by_age_female <- df_lt_bc %>%
    select(female) %>%
    as.matrix()
  v_r_mort_by_age <- (v_r_mort_by_age_male * n_male) + (v_r_mort_by_age_female * (1 - n_male)) # weighted mortality
  return(v_r_mort_by_age)
}

#' Load all parameters
#'
#' \code{load_all_params} loads all parameters for the decision model from multiple sources and creates a list.
#'
#' @param file.init String with the location and name of the file with initial set of parameters
#' @param file.init_dist String with the location and name of the file with initial distributions
#' @param file.mort String with the location and name of the file with mortality data
#' @param file.death_hr String with the location and name of death hazard ratios
#' @param file.weibull_scale String with the location and name of the file with weibull scale
#' @param file.weibull_shape String with the location and name of the file with weibull shape
#' @param file.exit_dest String with the location and name of the file with empirical destination states
#' @param file.overdose String with the location and name of the file with overdose/fentanyl-related parameters
#' @param file.cali_params String with the location and name of the file with calibrated parameters
#' @param file.fentanyl String with the location and name of the file with fentanyl exposure parameters
#' @param file.naloxone String with the location and name of the file with naloxone parameters
#'
#' @return
#' A list of all parameters used for the decision model.
#' @export
load_all_params <- function(
    file.init = NULL,
    file.pop_scaling = NULL,
    file.init_dist = NULL,
    file.cohort_balance = NULL,
    file.mort = NULL,
    file.death_hr = NULL,
    file.weibull = NULL,
    file.ce_tx = NULL,
    file.ce_death = NULL,
    file.exit_dest = NULL,
    file.overdose = NULL,
    file.cali_params = NULL,
    file.cali_targets = NULL,
    file.fentanyl = NULL,
    file.naloxone = NULL,
    file.qalys = NULL) {
  # Load files of all baseline model parameters
  df_init_params <- read.csv(file = file.init, row.names = 1, header = TRUE) # Initial parameter values
  df_pop_scaling <- read.csv(file = file.pop_scaling, row.names = 1, header = TRUE) # Annual population scaling factor
  df_init_dist <- read.csv(file = file.init_dist, row.names = 1, header = TRUE) # Initial parameter values
  df_cohort_balance <- read.csv(file = file.cohort_balance, row.names = 1, header = TRUE)
  df_death_hr <- read.csv(file = file.death_hr, row.names = 1, header = TRUE) # Mortality hazard ratios
  df_weibull <- read.csv(file = file.weibull, row.names = 1, header = TRUE) # Weibull params
  df_ce_tx <- read.csv(file = file.ce_tx, row.names = 1, header = TRUE) # Comparative effectiveness parameters (treatment retention)
  df_ce_death <- read.csv(file = file.ce_death, row.names = 1, header = TRUE) # Comparative effectiveness parameters (mortality)
  df_exit_dest <- read.csv(file = file.exit_dest, row.names = 1, header = TRUE) # Unconditional transition probs
  df_overdose <- read.csv(file = file.overdose, row.names = 1, header = TRUE) # Overdose-fentanyl parameters
  df_cali_params <- read.csv(file = file.cali_params, row.names = 1, header = TRUE) # Calibrated parameters
  df_cali_targets <- read.csv(file = file.cali_targets, header = TRUE) # Calibration targets
  df_fentanyl <- read.csv(file = file.fentanyl, row.names = 1, header = TRUE) # Fentanyl exposure parameters
  df_naloxone <- read.csv(file = file.naloxone, row.names = 1, header = TRUE) # Time-varying naloxone parameters for calibration
  df_qalys <- read.csv(file = file.qalys, row.names = 1, header = TRUE) # QALY weights

  l_params_all <- list(

    #### Initial parameters ####
    # n_pop_oat = df_init_params["pe", "pop_oat"],
    # n_pop_est = df_init_params["pe", "pop_est_2020"],
    n_age_init = df_init_params["pe", "age_init"], # age at baseline
    n_age_max_full = df_init_params["pe", "age_max_full"], # maximum age of follow up
    # n_age_max_fent = df_init_params["pe", "age_max_fent"], # maximum age of follow up
    n_per = df_init_params["pe", "period_yr"], # periods per year (e.g. 12-months)
    # n_discount = df_init_params["pe", "discount"], # discount rate
    n_cali_per = df_init_params["pe", "cali_per"], # number of calibration periods
    p_male = df_init_params["pe", "male_prop"], # % male
    # p_inj = df_init_params["pe", "inj_prop"], # % injection
    n_per_cali = df_init_params["pe", "n_per_cali"],
    n_per_full = df_init_params["pe", "n_per_full"],

    #### Annual OAT population scaling ####
    # Total individuals in each year
    n_pop_scaling_2010 = df_pop_scaling["2010", "n"],
    n_pop_scaling_2011 = df_pop_scaling["2011", "n"],
    n_pop_scaling_2012 = df_pop_scaling["2012", "n"],
    n_pop_scaling_2013 = df_pop_scaling["2013", "n"],
    n_pop_scaling_2014 = df_pop_scaling["2014", "n"],
    n_pop_scaling_2015 = df_pop_scaling["2015", "n"],
    n_pop_scaling_2016 = df_pop_scaling["2016", "n"],
    n_pop_scaling_2017 = df_pop_scaling["2017", "n"],
    n_pop_scaling_2018 = df_pop_scaling["2018", "n"],
    n_pop_scaling_2019 = df_pop_scaling["2019", "n"],
    n_pop_scaling_2020 = df_pop_scaling["2020", "n"],

    #### Initial state distribution ####
    #v_init_dist = as.vector(df_init_dist["pe", ]),
    bnx_inc = df_init_dist["pe", "bnx_inc"],
    met_inc = df_init_dist["pe", "met_inc"],
    abs_inc = df_init_dist["pe", "abs_inc"],
    oum_inc = df_init_dist["pe", "oum_inc"],
    oub_inc = df_init_dist["pe", "oub_inc"],
    ouo_inc = df_init_dist["pe", "ouo_inc"],
    odn_inc = df_init_dist["pe", "odn_inc"],
    odf_inc = df_init_dist["pe", "odf_inc"],
    abs_inc = df_init_dist["pe", "abs_inc"],
    bnx_prev = df_init_dist["pe", "bnx_prev"],
    met_prev = df_init_dist["pe", "met_prev"],
    abs_prev = df_init_dist["pe", "abs_prev"],
    oum_prev = df_init_dist["pe", "oum_prev"],
    oub_prev = df_init_dist["pe", "oub_prev"],
    ouo_prev = df_init_dist["pe", "ouo_prev"],
    odn_prev = df_init_dist["pe", "odn_prev"],
    odf_prev = df_init_dist["pe", "odf_prev"],
    abs_prev = df_init_dist["pe", "abs_prev"],

    #### Cohort balance ####
    # Proportion rebalance from incident to prevalent (net)
    p_inc_prev_2010 = df_cohort_balance["2010", "inc_prev"],
    p_inc_prev_2011 = df_cohort_balance["2011", "inc_prev"],
    p_inc_prev_2012 = df_cohort_balance["2012", "inc_prev"],
    p_inc_prev_2013 = df_cohort_balance["2013", "inc_prev"],
    p_inc_prev_2014 = df_cohort_balance["2014", "inc_prev"],
    p_inc_prev_2015 = df_cohort_balance["2015", "inc_prev"],
    p_inc_prev_2016 = df_cohort_balance["2016", "inc_prev"],
    p_inc_prev_2017 = df_cohort_balance["2017", "inc_prev"],
    p_inc_prev_2018 = df_cohort_balance["2018", "inc_prev"],
    p_inc_prev_2019 = df_cohort_balance["2019", "inc_prev"],
    p_inc_prev_2020 = df_cohort_balance["2020", "inc_prev"],
    # Proportion rebalance from prevalent to incident (net)
    p_prev_inc_2010 = df_cohort_balance["2010", "prev_inc"],
    p_prev_inc_2011 = df_cohort_balance["2011", "prev_inc"],
    p_prev_inc_2012 = df_cohort_balance["2012", "prev_inc"],
    p_prev_inc_2013 = df_cohort_balance["2013", "prev_inc"],
    p_prev_inc_2014 = df_cohort_balance["2014", "prev_inc"],
    p_prev_inc_2015 = df_cohort_balance["2015", "prev_inc"],
    p_prev_inc_2016 = df_cohort_balance["2016", "prev_inc"],
    p_prev_inc_2017 = df_cohort_balance["2017", "prev_inc"],
    p_prev_inc_2018 = df_cohort_balance["2018", "prev_inc"],
    p_prev_inc_2019 = df_cohort_balance["2019", "prev_inc"],
    p_prev_inc_2020 = df_cohort_balance["2020", "prev_inc"],

    # Adjustment for new entries into treatment
    p_inc_entry_adj_2010 = df_cohort_balance["2010", "inc_entry"],
    p_inc_entry_adj_2011 = df_cohort_balance["2011", "inc_entry"],
    p_inc_entry_adj_2012 = df_cohort_balance["2012", "inc_entry"],
    p_inc_entry_adj_2013 = df_cohort_balance["2013", "inc_entry"],
    p_inc_entry_adj_2014 = df_cohort_balance["2014", "inc_entry"],
    p_inc_entry_adj_2015 = df_cohort_balance["2015", "inc_entry"],
    p_inc_entry_adj_2016 = df_cohort_balance["2016", "inc_entry"],
    p_inc_entry_adj_2017 = df_cohort_balance["2017", "inc_entry"],
    p_inc_entry_adj_2018 = df_cohort_balance["2018", "inc_entry"],
    p_inc_entry_adj_2019 = df_cohort_balance["2019", "inc_entry"],
    p_inc_entry_adj_2020 = df_cohort_balance["2020", "inc_entry"],

    # Time-points for calibration
    n_per_2012 = subset(df_cali_targets, year == 2012 & target == "fatal overdoses")$time,
    n_per_2013 = subset(df_cali_targets, year == 2013 & target == "fatal overdoses")$time,
    n_per_2014 = subset(df_cali_targets, year == 2014 & target == "fatal overdoses")$time,
    n_per_2015 = subset(df_cali_targets, year == 2015 & target == "fatal overdoses")$time,
    n_per_2016 = subset(df_cali_targets, year == 2016 & target == "fatal overdoses")$time,
    n_per_2017 = subset(df_cali_targets, year == 2017 & target == "fatal overdoses")$time,
    n_per_2018 = subset(df_cali_targets, year == 2018 & target == "fatal overdoses")$time,
    n_per_2019 = subset(df_cali_targets, year == 2019 & target == "fatal overdoses")$time,
    n_per_2020 = subset(df_cali_targets, year == 2020 & target == "fatal overdoses")$time,

    #### Mortality ####

    v_r_mort_by_age = load_mort_params(file = file.mort, n_male = df_init_params["pe", "male_prop"]), # vector of age-specific mortality

    #### Hazard ratios for non-overdose mortality ####
    hr_abs = df_death_hr["pe", "abs"],
    hr_dno_ou_oat = df_death_hr["pe", "dno_ou_oat"],

    #### BNX vs. methadone comparative effectiveness ####
    # Discontinuation (inc)
    hr_tx_itt_ps_inc = df_ce_tx["pe", "itt_ps_inc"],
    hr_tx_pp_inc = df_ce_tx["pe", "pp_inc"],
    hr_tx_pp_hd_inc = df_ce_tx["pe", "pp_hd_inc"],
    # Discontinuation (prev)
    hr_tx_itt_ps_prev = df_ce_tx["pe", "itt_ps_prev"],
    hr_tx_pp_prev = df_ce_tx["pe", "pp_prev"],
    hr_tx_pp_hd_prev = df_ce_tx["pe", "pp_hd_prev"],
    # Mortality (inc)
    hr_death_itt_ps_inc = df_ce_death["pe", "itt_ps_inc"],
    hr_death_pp_inc = df_ce_death["pe", "pp_inc"],
    hr_death_pp_hd_inc = df_ce_death["pe", "pp_hd_inc"],
    # Mortality (prev)
    hr_death_itt_ps_prev = df_ce_death["pe", "itt_ps_prev"],
    hr_death_pp_prev = df_ce_death["pe", "pp_prev"],
    hr_death_pp_hd_prev = df_ce_death["pe", "pp_hd_prev"],

    #### Load weibull parameters ####
    # Incident new users (for episode 1)
    # Shape
    p_weibull_shape_bnx_inc = df_weibull["pe", "bnx_shape_inc"],
    p_weibull_shape_met_inc = df_weibull["pe", "met_shape_inc"],
    p_weibull_shape_oat_inc = df_weibull["pe", "oat_shape_inc"],

    # scale
    p_weibull_scale_bnx_inc = df_weibull["pe", "bnx_scale_inc"],
    p_weibull_scale_met_inc = df_weibull["pe", "met_scale_inc"],
    p_weibull_scale_oat_inc = df_weibull["pe", "oat_scale_inc"],

    # Prevalent users (for episode 2+)
    # Shape
    p_weibull_shape_bnx_prev = df_weibull["pe", "bnx_shape_prev"],
    p_weibull_shape_met_prev = df_weibull["pe", "met_shape_prev"],
    p_weibull_shape_oat_prev = df_weibull["pe", "oat_shape_prev"],

    # scale
    p_weibull_scale_bnx_prev = df_weibull["pe", "bnx_scale_prev"],
    p_weibull_scale_met_prev = df_weibull["pe", "met_scale_prev"],
    p_weibull_scale_oat_prev = df_weibull["pe", "oat_scale_prev"],

    #### Transition probabilities on state exit ####
    # From BNX
    p_bnx_met_inc = df_exit_dest["bnx_inc", "met_inc"],
    p_bnx_abs_inc = df_exit_dest["bnx_inc", "abs_inc"],
    p_bnx_oum_inc = df_exit_dest["bnx_inc", "oum_inc"],
    p_bnx_oub_inc = df_exit_dest["bnx_inc", "oub_inc"],
    p_bnx_ouo_inc = df_exit_dest["bnx_inc", "ouo_inc"],
    p_bnx_met_prev = df_exit_dest["bnx_prev", "met_prev"],
    p_bnx_abs_prev = df_exit_dest["bnx_prev", "abs_prev"],
    p_bnx_oum_prev = df_exit_dest["bnx_prev", "oum_prev"],
    p_bnx_oub_prev = df_exit_dest["bnx_prev", "oub_prev"],
    p_bnx_ouo_prev = df_exit_dest["bnx_prev", "ouo_prev"],
    # From MET
    p_met_bnx_inc = df_exit_dest["met_inc", "bnx_inc"],
    p_met_abs_inc = df_exit_dest["met_inc", "abs_inc"],
    p_met_oum_inc = df_exit_dest["met_inc", "oum_inc"],
    p_met_oub_inc = df_exit_dest["met_inc", "oub_inc"],
    p_met_ouo_inc = df_exit_dest["met_inc", "ouo_inc"],
    p_met_bnx_prev = df_exit_dest["met_prev", "bnx_prev"],
    p_met_abs_prev = df_exit_dest["met_prev", "abs_prev"],
    p_met_oum_prev = df_exit_dest["met_prev", "oum_prev"],
    p_met_oub_prev = df_exit_dest["met_prev", "oub_prev"],
    p_met_ouo_prev = df_exit_dest["met_prev", "ouo_prev"],
    # From ABS
    p_abs_oum_inc = df_exit_dest["abs_inc", "oum_inc"],
    p_abs_oub_inc = df_exit_dest["abs_inc", "oub_inc"],
    p_abs_ouo_inc = df_exit_dest["abs_inc", "ouo_inc"],
    p_abs_met_inc = df_exit_dest["abs_inc", "met_inc"],
    p_abs_bnx_inc = df_exit_dest["abs_inc", "bnx_inc"],
    p_abs_oum_prev = df_exit_dest["abs_prev", "oum_prev"],
    p_abs_oub_prev = df_exit_dest["abs_prev", "oub_prev"],
    p_abs_ouo_prev = df_exit_dest["abs_prev", "ouo_prev"],
    p_abs_met_prev = df_exit_dest["abs_prev", "met_prev"],
    p_abs_bnx_prev = df_exit_dest["abs_prev", "bnx_prev"],
    # From OUM
    p_oum_met_inc = df_exit_dest["oum_inc", "met_inc"],
    p_oum_bnx_inc = df_exit_dest["oum_inc", "bnx_inc"],
    p_oum_abs_inc = df_exit_dest["oum_inc", "abs_inc"],
    p_oum_oub_inc = df_exit_dest["oum_inc", "oub_inc"],
    p_oum_ouo_inc = df_exit_dest["oum_inc", "ouo_inc"],
    p_oum_met_prev = df_exit_dest["oum_prev", "met_prev"],
    p_oum_bnx_prev = df_exit_dest["oum_prev", "bnx_prev"],
    p_oum_abs_prev = df_exit_dest["oum_prev", "abs_prev"],
    p_oum_oub_prev = df_exit_dest["oum_prev", "oub_prev"],
    p_oum_ouo_prev = df_exit_dest["oum_prev", "ouo_prev"],
    # From OUB
    p_oub_met_inc = df_exit_dest["oub_inc", "met_inc"],
    p_oub_bnx_inc = df_exit_dest["oub_inc", "bnx_inc"],
    p_oub_abs_inc = df_exit_dest["oub_inc", "abs_inc"],
    p_oub_oum_inc = df_exit_dest["oub_inc", "oum_inc"],
    p_oub_ouo_inc = df_exit_dest["oub_inc", "ouo_inc"],
    p_oub_met_prev = df_exit_dest["oub_prev", "met_prev"],
    p_oub_bnx_prev = df_exit_dest["oub_prev", "bnx_prev"],
    p_oub_abs_prev = df_exit_dest["oub_prev", "abs_prev"],
    p_oub_oum_prev = df_exit_dest["oub_prev", "oum_prev"],
    p_oub_ouo_prev = df_exit_dest["oub_prev", "ouo_prev"],
    # From OUO
    p_ouo_met_inc = df_exit_dest["ouo_inc", "met_inc"],
    p_ouo_bnx_inc = df_exit_dest["ouo_inc", "bnx_inc"],
    p_ouo_abs_inc = df_exit_dest["ouo_inc", "abs_inc"],
    p_ouo_oum_inc = df_exit_dest["ouo_inc", "oum_inc"],
    p_ouo_oub_inc = df_exit_dest["ouo_inc", "oub_inc"],
    p_ouo_met_prev = df_exit_dest["ouo_prev", "met_prev"],
    p_ouo_bnx_prev = df_exit_dest["ouo_prev", "bnx_prev"],
    p_ouo_abs_prev = df_exit_dest["ouo_prev", "abs_prev"],
    p_ouo_oum_prev = df_exit_dest["ouo_prev", "oum_prev"],
    p_ouo_oub_prev = df_exit_dest["ouo_prev", "oub_prev"],
    # From OD
    p_odn_met_inc = df_exit_dest["odn_inc", "met_inc"],
    p_odn_bnx_inc = df_exit_dest["odn_inc", "bnx_inc"],
    p_odn_abs_inc = df_exit_dest["odn_inc", "abs_inc"],
    p_odn_oum_inc = df_exit_dest["odn_inc", "oum_inc"],
    p_odn_oub_inc = df_exit_dest["odn_inc", "oub_inc"],
    p_odn_ouo_inc = df_exit_dest["odn_inc", "ouo_inc"],
    p_odn_met_prev = df_exit_dest["odn_prev", "met_prev"],
    p_odn_bnx_prev = df_exit_dest["odn_prev", "bnx_prev"],
    p_odn_abs_prev = df_exit_dest["odn_prev", "abs_prev"],
    p_odn_oum_prev = df_exit_dest["odn_prev", "oum_prev"],
    p_odn_oub_prev = df_exit_dest["odn_prev", "oub_prev"],
    p_odn_ouo_prev = df_exit_dest["odn_prev", "ouo_prev"],

    #### Overdose ####
    # Overdose parameters loaded as rates and converted to probabilities in model
    # Includes additional calibration related parameters
    # Base overdose params
    # Calibrated params
    # Mean estimates
    n_oat_od = df_cali_params["pe", "n_oat_od"],
    n_fent_prev_od_mult = df_cali_params["pe", "n_fent_prev_od_mult"],
    n_fent_delta_od_mult = df_cali_params["pe", "n_fent_delta_od_mult"],
    n_fatal_od_oat = df_cali_params["pe", "n_fatal_od_oat"],
    hr_oat = df_cali_params["pe", "hr_oat"],
    p_weibull_scale_ou = df_cali_params["pe", "p_weibull_scale_ou"],
    p_weibull_shape_ou = df_cali_params["pe", "p_weibull_shape_ou"],
    p_weibull_scale_abs = df_cali_params["pe", "p_weibull_scale_abs"],
    p_weibull_shape_abs = df_cali_params["pe", "p_weibull_shape_abs"],

    # Gamma shape parameter (prior)
    n_oat_od_shape = df_cali_params["shape", "n_oat_od"],
    n_fent_prev_od_mult_low = df_cali_params["low", "n_fent_prev_od_mult"],
    n_fent_delta_od_mult_low = df_cali_params["low", "n_fent_delta_od_mult"],
    n_fatal_od_oat_shape = df_cali_params["shape", "n_fatal_od_oat"],
    hr_oat_shape = df_cali_params["shape", "hr_oat"],
    p_weibull_scale_ou_shape = df_cali_params["shape", "p_weibull_scale_ou"],
    p_weibull_shape_ou_shape = df_cali_params["shape", "p_weibull_shape_ou"],
    p_weibull_scale_abs_shape = df_cali_params["shape", "p_weibull_scale_abs"],
    p_weibull_shape_abs_shape = df_cali_params["shape", "p_weibull_shape_abs"],

    # Gamma scale parameter (prior)
    n_oat_od_scale = df_cali_params["scale", "n_oat_od"],
    n_fent_prev_od_mult_high = df_cali_params["high", "n_fent_prev_od_mult"],
    n_fent_delta_od_mult_high = df_cali_params["high", "n_fent_delta_od_mult"],
    n_fatal_od_oat_scale = df_cali_params["scale", "n_fatal_od_oat"],
    hr_oat_scale = df_cali_params["scale", "hr_oat"],
    p_weibull_scale_ou_scale = df_cali_params["scale", "p_weibull_scale_ou"],
    p_weibull_shape_ou_scale = df_cali_params["scale", "p_weibull_shape_ou"],
    p_weibull_scale_abs_scale = df_cali_params["scale", "p_weibull_scale_abs"],
    p_weibull_shape_abs_scale = df_cali_params["scale", "p_weibull_shape_abs"],

    # Non-calibrated overdose params
    # ABS
    n_abs_od = df_overdose["pe", "abs_od"],
    # BNX
    n_bnx_od_mult = df_overdose["pe", "bnx_od_mult"],
    # MET
    n_met_od_mult = df_overdose["pe", "met_od_mult"],

    # Opioid use
    n_ou_od_mult = df_overdose["pe", "ou_od_mult"],
    # n_ou_od_mult_shape = df_overdose["shape", "ou_od_mult"],
    # n_ou_od_mult_scale = df_overdose["scale", "ou_od_mult"],
    # n_oum_od_mult = df_overdose["pe", "ou_od_mult"],
    # n_oub_od_mult = df_overdose["pe", "ou_od_mult"],
    # n_ouo_od_mult = df_overdose["pe", "ou_od_mult"],

    # Abstinence
    n_abs_od_mult = df_overdose["pe", "abs_od_mult"],
    n_abs_od_mult_shape = df_overdose["shape", "abs_od_mult"],
    n_abs_od_mult_scale = df_overdose["scale", "abs_od_mult"],

    # Relative multipliers for out-of-treatment (opioid use) vs. in treatment (OAT)
    # Non-fatal overdose
    n_ou_oat_odn_mult = df_overdose["pe", "ou_oat_odn_mult"],
    n_ou_oat_odn_mult_shape = df_overdose["shape", "ou_oat_odn_mult"],
    n_ou_oat_odn_mult_scale = df_overdose["scale", "ou_oat_odn_mult"],

    # Fatal overdose
    n_ou_oat_odf_mult = df_overdose["pe", "ou_oat_odf_mult"],
    n_ou_oat_odf_mult_shape = df_overdose["shape", "ou_oat_odf_mult"],
    n_ou_oat_odf_mult_scale = df_overdose["scale", "ou_oat_odf_mult"],

    # Fentanyl and naloxone
    v_fent_prev = as.vector(df_fentanyl[, "prev"]),
    v_fent_delta = as.vector(df_fentanyl[, "delta"]),
    v_nx_rev = as.vector(df_naloxone[, "pe"]),

    # QALY weights (all set to 1 for alive)
    u_bnx = df_qalys["pe", "bnx"],
    u_met = df_qalys["pe", "met"],
    u_oum = df_qalys["pe", "oum"],
    u_oub = df_qalys["pe", "oub"],
    u_ouo = df_qalys["pe", "ouo"],
    u_odn = df_qalys["pe", "odn"],
    u_odf = df_qalys["pe", "odf"],
    u_abs = df_qalys["pe", "abs"]
  ) # Close list
  return(l_params_all) # Return full parameter list
}

#' Update parameters
#'
#' \code{update_param_list} is used to update list of all parameters with new
#' values for specific parameters.
#'
#' @param l_params_all List with all parameters of decision model
#' @param params_updated Parameters for which values need to be updated
#' @return
#' A modifed list with all parameters updated.
#' @export

update_param_list <- function(l_params_all, params_updated) {
  if (typeof(params_updated) != "list") {
    params_updated <- split(unname(params_updated), names(params_updated)) # convert the named vector to a list
  }
  l_params_all <- modifyList(l_params_all, params_updated) # update the values
  return(l_params_all)
}
