#' One-Way Sensitivity Analysis Parameters
#'
#' This function extracts and names the parameters for a one-way sensitivity analysis (OWSA) from a given data frame.
#'
#' @param df_dsa_owsa A data frame containing the OWSA parameters. It should have rows for "low" and "high" values of each parameter.
#' @param ce_est select which comparative effectiveness specification to use.
#'
#' @return A list containing the named parameters for the OWSA.
#' @export
#'
#' @examples
#' df_dsa_owsa <- data.frame(
#'   low = c(male_prop = 0.4, hr_dno_ou_oat = 0.8, n_abs_od = 0.1, n_ou_od_mult = 1.5),
#'   high = c(male_prop = 0.6, hr_dno_ou_oat = 1.2, n_abs_od = 0.3, n_ou_od_mult = 2.0)
#' )
#' ce_est <- NULL
#' params <- owsa_params(df_dsa_owsa, ce_est)
#' print(params)
owsa_params <- function(
    df_dsa_owsa = df_dsa_owsa,
    v_init_dist_met_25_inc = v_init_dist_met_25_inc,
    v_init_dist_met_50_inc = v_init_dist_met_50_inc,
    v_init_dist_met_75_inc = v_init_dist_met_75_inc,
    v_init_dist_bnx_25_inc = v_init_dist_bnx_25_inc,
    v_init_dist_bnx_50_inc = v_init_dist_bnx_50_inc,
    v_init_dist_bnx_75_inc = v_init_dist_bnx_75_inc,
    tx = tx,
    ce_est = ce_est) {
  ## General params (not specific to intention-to-treat or per-protocol)
  ########################
  #### Initial params ####
  ########################
  # male proportion of population
  p_male_low <- unlist(df_dsa_owsa["low", "male_prop"])
  p_male_high <- unlist(df_dsa_owsa["high", "male_prop"])
  names(p_male_low) <- c("male_prop")
  names(p_male_high) <- c("male_prop")

  ########################
  #### Cohort balance ####
  ########################
  p_inc_prev_2010 <- p_inc_prev_2011 <- p_inc_prev_2012 <- p_inc_prev_2013 <- p_inc_prev_2014 <- 0
  p_inc_prev_2015 <- p_inc_prev_2016 <- p_inc_prev_2017 <- p_inc_prev_2018 <- p_inc_prev_2019 <- p_inc_prev_2020 <- 0
  v_inc_prev <- c(
    p_inc_prev_2010, p_inc_prev_2011, p_inc_prev_2012, p_inc_prev_2013, p_inc_prev_2014, p_inc_prev_2015,
    p_inc_prev_2016, p_inc_prev_2017, p_inc_prev_2018, p_inc_prev_2019, p_inc_prev_2020
  )
  names(v_inc_prev) <- c(
    "p_inc_prev_2010", "p_inc_prev_2011", "p_inc_prev_2012", "p_inc_prev_2013", "p_inc_prev_2014", "p_inc_prev_2015",
    "p_inc_prev_2016", "p_inc_prev_2017", "p_inc_prev_2018", "p_inc_prev_2019", "p_inc_prev_2020"
  )
  # Create combined vector of parameters for cohort balance (set transitions from inc to prev to zero, set baseline to fixed %)
  if (tx == "met") {
    v_cohort_balance_25_inc <- c(v_init_dist_met_25_inc, v_inc_prev)
    v_cohort_balance_50_inc <- c(v_init_dist_met_50_inc, v_inc_prev)
    v_cohort_balance_75_inc <- c(v_init_dist_met_75_inc, v_inc_prev)
  } else if (tx == "bnx") {
    v_cohort_balance_25_inc <- c(v_init_dist_bnx_25_inc, v_inc_prev)
    v_cohort_balance_50_inc <- c(v_init_dist_bnx_50_inc, v_inc_prev)
    v_cohort_balance_75_inc <- c(v_init_dist_bnx_75_inc, v_inc_prev)
  } else {
    stop("tx must be either 'met' or 'bnx'")
  }

  ############################
  #### Non-overdose death ####
  ############################
  # Opioid use vs. OAT
  hr_dno_ou_oat_dsa_low <- unlist(df_dsa_owsa["low", "hr_dno_ou_oat"])
  hr_dno_ou_oat_dsa_high <- unlist(df_dsa_owsa["high", "hr_dno_ou_oat"])
  names(hr_dno_ou_oat_dsa_low) <- c("hr_dno_ou_oat")
  names(hr_dno_ou_oat_dsa_high) <- c("hr_dno_ou_oat")

  ###################################
  #### Overdose (non-calibrated) ####
  ###################################
  # Overdose rate in abstinence
  n_abs_od_low <- unlist(df_dsa_owsa["low", "n_abs_od"])
  n_abs_od_high <- unlist(df_dsa_owsa["high", "n_abs_od"])
  names(n_abs_od_low) <- c("n_abs_od")
  names(n_abs_od_high) <- c("n_abs_od")

  # Multiplier for first week out of treatment
  n_ou_od_mult_low <- unlist(df_dsa_owsa["low", "n_ou_od_mult"])
  n_ou_od_mult_high <- unlist(df_dsa_owsa["high", "n_ou_od_mult"])
  names(n_ou_od_mult_low) <- c("n_ou_od_mult")
  names(n_ou_od_mult_high) <- c("n_ou_od_mult")

  # Non-fatal overdose rate multiplier for out-of-treatment vs. OAT
  n_ou_oat_odn_mult_low <- unlist(df_dsa_owsa["low", "n_ou_oat_odn_mult"])
  n_ou_oat_odn_mult_high <- unlist(df_dsa_owsa["high", "n_ou_oat_odn_mult"])
  names(n_ou_oat_odn_mult_low) <- c("n_ou_oat_odn_mult")
  names(n_ou_oat_odn_mult_high) <- c("n_ou_oat_odn_mult")

  # Fatal overdose rate multiplier for out-of-treatment vs. OAT
  n_ou_oat_odf_mult_low <- unlist(df_dsa_owsa["low", "n_ou_oat_odf_mult"])
  n_ou_oat_odf_mult_high <- unlist(df_dsa_owsa["high", "n_ou_oat_odf_mult"])
  names(n_ou_oat_odf_mult_low) <- c("n_ou_oat_odf_mult")
  names(n_ou_oat_odf_mult_high) <- c("n_ou_oat_odf_mult")

  #############################################
  #### Comparative effectiveness estimates ####
  #############################################
  if (ce_est == "itt_ps") {
    # ITT-PS
    # TX - Inc
    hr_tx_itt_ps_inc_dsa_low <- unlist(df_dsa_owsa["low", "hr_tx_itt_ps_inc"])
    hr_tx_itt_ps_inc_dsa_high <- unlist(df_dsa_owsa["high", "hr_tx_itt_ps_inc"])
    names(hr_tx_itt_ps_inc_dsa_low) <- c("hr_tx_itt_ps_inc")
    names(hr_tx_itt_ps_inc_dsa_high) <- c("hr_tx_itt_ps_inc")
    # TX - Prev
    hr_tx_itt_ps_prev_dsa_low <- unlist(df_dsa_owsa["low", "hr_tx_itt_ps_prev"])
    hr_tx_itt_ps_prev_dsa_high <- unlist(df_dsa_owsa["high", "hr_tx_itt_ps_prev"])
    names(hr_tx_itt_ps_prev_dsa_low) <- c("hr_tx_itt_ps_prev")
    names(hr_tx_itt_ps_prev_dsa_high) <- c("hr_tx_itt_ps_prev")
    # Death - Inc
    hr_death_itt_ps_inc_dsa_low <- unlist(df_dsa_owsa["low", "hr_death_itt_ps_inc"])
    hr_death_itt_ps_inc_dsa_high <- unlist(df_dsa_owsa["high", "hr_death_itt_ps_inc"])
    names(hr_death_itt_ps_inc_dsa_low) <- c("hr_death_itt_ps_inc")
    names(hr_death_itt_ps_inc_dsa_high) <- c("hr_death_itt_ps_inc")
    # Death - Prev
    hr_death_itt_ps_prev_dsa_low <- unlist(df_dsa_owsa["low", "hr_death_itt_ps_prev"])
    hr_death_itt_ps_prev_dsa_high <- unlist(df_dsa_owsa["high", "hr_death_itt_ps_prev"])
    names(hr_death_itt_ps_prev_dsa_low) <- c("hr_death_itt_ps_prev")
    names(hr_death_itt_ps_prev_dsa_high) <- c("hr_death_itt_ps_prev")

    ###############################
    #### Calibrated parameters ####
    ###############################
    # Overdose rate in OAT
    n_oat_od_itt_low <- unlist(df_dsa_owsa["low", "n_oat_od_itt"])
    n_oat_od_itt_high <- unlist(df_dsa_owsa["high", "n_oat_od_itt"])
    names(n_oat_od_itt_low) <- c("n_oat_od")
    names(n_oat_od_itt_high) <- c("n_oat_od")

    # Overdose rate multiplier on fentanyl prevalence
    n_fent_prev_od_mult_itt_low <- unlist(df_dsa_owsa["low", "n_fent_prev_od_mult_itt"])
    n_fent_prev_od_mult_itt_high <- unlist(df_dsa_owsa["high", "n_fent_prev_od_mult_itt"])
    names(n_fent_prev_od_mult_itt_low) <- c("n_fent_prev_od_mult")
    names(n_fent_prev_od_mult_itt_high) <- c("n_fent_prev_od_mult")

    # Overdose rate multiplier on change in fentanyl prevalence
    n_fent_delta_od_mult_itt_low <- unlist(df_dsa_owsa["low", "n_fent_delta_od_mult_itt"])
    n_fent_delta_od_mult_itt_high <- unlist(df_dsa_owsa["high", "n_fent_delta_od_mult_itt"])
    names(n_fent_delta_od_mult_itt_low) <- c("n_fent_delta_od_mult")
    names(n_fent_delta_od_mult_itt_high) <- c("n_fent_delta_od_mult")

    # Fatal overdose rate in OAT
    n_fatal_od_oat_itt_low <- unlist(df_dsa_owsa["low", "n_fatal_od_oat_itt"])
    n_fatal_od_oat_itt_high <- unlist(df_dsa_owsa["high", "n_fatal_od_oat_itt"])
    names(n_fatal_od_oat_itt_low) <- c("n_fatal_od_oat")
    names(n_fatal_od_oat_itt_high) <- c("n_fatal_od_oat")

    # Non-overdose mortality rate in OAT
    hr_oat_itt_low <- unlist(df_dsa_owsa["low", "hr_oat_itt"])
    hr_oat_itt_high <- unlist(df_dsa_owsa["high", "hr_oat_itt"])
    names(hr_oat_itt_low) <- c("hr_oat")
    names(hr_oat_itt_high) <- c("hr_oat")

    # Weibull scale (out-of-treatment)
    p_weibull_scale_ou_itt_low <- unlist(df_dsa_owsa["low", "p_weibull_scale_ou_itt"])
    p_weibull_scale_ou_itt_high <- unlist(df_dsa_owsa["high", "p_weibull_scale_ou_itt"])
    names(p_weibull_scale_ou_itt_low) <- c("p_weibull_scale_ou")
    names(p_weibull_scale_ou_itt_high) <- c("p_weibull_scale_ou")

    # Weibull shape (out-of-treatment)
    p_weibull_shape_ou_itt_low <- unlist(df_dsa_owsa["low", "p_weibull_shape_ou_itt"])
    p_weibull_shape_ou_itt_high <- unlist(df_dsa_owsa["high", "p_weibull_shape_ou_itt"])
    names(p_weibull_shape_ou_itt_low) <- c("p_weibull_shape_ou")
    names(p_weibull_shape_ou_itt_high) <- c("p_weibull_shape_ou")

    # Weibull scale (abstinence)
    p_weibull_scale_abs_itt_low <- unlist(df_dsa_owsa["low", "p_weibull_scale_abs_itt"])
    p_weibull_scale_abs_itt_high <- unlist(df_dsa_owsa["high", "p_weibull_scale_abs_itt"])
    names(p_weibull_scale_abs_itt_low) <- c("p_weibull_scale_abs")
    names(p_weibull_scale_abs_itt_high) <- c("p_weibull_scale_abs")

    # Weibull shape (abstinence)
    p_weibull_shape_abs_itt_low <- unlist(df_dsa_owsa["low", "p_weibull_shape_abs_itt"])
    p_weibull_shape_abs_itt_high <- unlist(df_dsa_owsa["high", "p_weibull_shape_abs_itt"])
    names(p_weibull_shape_abs_itt_low) <- c("p_weibull_shape_abs")
    names(p_weibull_shape_abs_itt_high) <- c("p_weibull_shape_abs")

    ##################################
    #### Weibull methadone remain ####
    ##################################
    # Weibull scale methadone (inc)
    p_weibull_scale_met_inc_itt_low <- unlist(df_dsa_owsa["low", "p_weibull_scale_met_inc_itt"])
    p_weibull_scale_met_inc_itt_high <- unlist(df_dsa_owsa["high", "p_weibull_scale_met_inc_itt"])
    names(p_weibull_scale_met_inc_itt_low) <- c("p_weibull_scale_met_inc")
    names(p_weibull_scale_met_inc_itt_high) <- c("p_weibull_scale_met_inc")

    # Weibull shape methadone (inc)
    p_weibull_shape_met_inc_itt_low <- unlist(df_dsa_owsa["low", "p_weibull_shape_met_inc_itt"])
    p_weibull_shape_met_inc_itt_high <- unlist(df_dsa_owsa["high", "p_weibull_shape_met_inc_itt"])
    names(p_weibull_shape_met_inc_itt_low) <- c("p_weibull_shape_met_inc")
    names(p_weibull_shape_met_inc_itt_high) <- c("p_weibull_shape_met_inc")

    # Weibull scale methadone (prev)
    p_weibull_scale_met_prev_itt_low <- unlist(df_dsa_owsa["low", "p_weibull_scale_met_prev_itt"])
    p_weibull_scale_met_prev_itt_high <- unlist(df_dsa_owsa["high", "p_weibull_scale_met_prev_itt"])
    names(p_weibull_scale_met_prev_itt_low) <- c("p_weibull_scale_met_prev")
    names(p_weibull_scale_met_prev_itt_high) <- c("p_weibull_scale_met_prev")

    # Weibull shape (abstinence)
    p_weibull_shape_met_prev_itt_low <- unlist(df_dsa_owsa["low", "p_weibull_shape_met_prev_itt"])
    p_weibull_shape_met_prev_itt_high <- unlist(df_dsa_owsa["high", "p_weibull_shape_met_prev_itt"])
    names(p_weibull_shape_met_prev_itt_low) <- c("p_weibull_shape_met_prev")
    names(p_weibull_shape_met_prev_itt_high) <- c("p_weibull_shape_met_prev")

    # Combine into vector
    # Low values
    v_owsa_low <- c(
      # Non-specific
      p_male_low, hr_dno_ou_oat_dsa_low, n_abs_od_low, n_ou_od_mult_low, n_ou_oat_odn_mult_low, n_ou_oat_odf_mult_low,
      # ITT-PS
      hr_tx_itt_ps_inc_dsa_low, hr_tx_itt_ps_prev_dsa_low, hr_death_itt_ps_inc_dsa_low, hr_death_itt_ps_prev_dsa_low,
      # ITT
      n_oat_od_itt_low, n_fent_prev_od_mult_itt_low, n_fent_delta_od_mult_itt_low, n_fatal_od_oat_itt_low, hr_oat_itt_low,
      p_weibull_scale_ou_itt_low, p_weibull_shape_ou_itt_low, p_weibull_scale_abs_itt_low, p_weibull_shape_abs_itt_low,
      p_weibull_scale_met_inc_itt_low, p_weibull_shape_met_inc_itt_low, p_weibull_scale_met_prev_itt_low, p_weibull_shape_met_prev_itt_low
    )
    # High values
    v_owsa_high <- c(
      # Non-specific
      p_male_high, hr_dno_ou_oat_dsa_high, n_abs_od_high, n_ou_od_mult_high, n_ou_oat_odn_mult_high, n_ou_oat_odf_mult_high,
      # ITT-PS
      hr_tx_itt_ps_inc_dsa_high, hr_tx_itt_ps_prev_dsa_high, hr_death_itt_ps_inc_dsa_high, hr_death_itt_ps_prev_dsa_high,
      # ITT
      n_oat_od_itt_high, n_fent_prev_od_mult_itt_high, n_fent_delta_od_mult_itt_high, n_fatal_od_oat_itt_high, hr_oat_itt_high,
      p_weibull_scale_ou_itt_high, p_weibull_shape_ou_itt_high, p_weibull_scale_abs_itt_high, p_weibull_shape_abs_itt_high,
      p_weibull_scale_met_inc_itt_high, p_weibull_shape_met_inc_itt_high, p_weibull_scale_met_prev_itt_high, p_weibull_shape_met_prev_itt_high
    )
    ## Per-protocol OWSA parameters
  } else if (ce_est == "pp") {
    # PP
    # TX - Inc
    hr_tx_pp_inc_dsa_low <- unlist(df_dsa_owsa["low", "hr_tx_pp_inc"])
    hr_tx_pp_inc_dsa_high <- unlist(df_dsa_owsa["high", "hr_tx_pp_inc"])
    names(hr_tx_pp_inc_dsa_low) <- c("hr_tx_pp_inc")
    names(hr_tx_pp_inc_dsa_high) <- c("hr_tx_pp_inc")
    # TX - Prev
    hr_tx_pp_prev_dsa_low <- unlist(df_dsa_owsa["low", "hr_tx_pp_prev"])
    hr_tx_pp_prev_dsa_high <- unlist(df_dsa_owsa["high", "hr_tx_pp_prev"])
    names(hr_tx_pp_prev_dsa_low) <- c("hr_tx_pp_prev")
    names(hr_tx_pp_prev_dsa_high) <- c("hr_tx_pp_prev")
    # Death - Inc
    hr_death_pp_inc_dsa_low <- unlist(df_dsa_owsa["low", "hr_death_pp_inc"])
    hr_death_pp_inc_dsa_high <- unlist(df_dsa_owsa["high", "hr_death_pp_inc"])
    names(hr_death_pp_inc_dsa_low) <- c("hr_death_pp_inc")
    names(hr_death_pp_inc_dsa_high) <- c("hr_death_pp_inc")
    # Death - Prev
    hr_death_pp_prev_dsa_low <- unlist(df_dsa_owsa["low", "hr_death_pp_prev"])
    hr_death_pp_prev_dsa_high <- unlist(df_dsa_owsa["high", "hr_death_pp_prev"])
    names(hr_death_pp_prev_dsa_low) <- c("hr_death_pp_prev")
    names(hr_death_pp_prev_dsa_high) <- c("hr_death_pp_prev")

    ###############################
    #### Calibrated parameters ####
    ###############################
    # Overdose rate in OAT
    n_oat_od_pp_low <- unlist(df_dsa_owsa["low", "n_oat_od_pp"])
    n_oat_od_pp_high <- unlist(df_dsa_owsa["high", "n_oat_od_pp"])
    names(n_oat_od_pp_low) <- c("n_oat_od")
    names(n_oat_od_pp_high) <- c("n_oat_od")

    # Overdose rate multiplier on fentanyl prevalence
    n_fent_prev_od_mult_pp_low <- unlist(df_dsa_owsa["low", "n_fent_prev_od_mult_pp"])
    n_fent_prev_od_mult_pp_high <- unlist(df_dsa_owsa["high", "n_fent_prev_od_mult_pp"])
    names(n_fent_prev_od_mult_pp_low) <- c("n_fent_prev_od_mult")
    names(n_fent_prev_od_mult_pp_high) <- c("n_fent_prev_od_mult")

    # Overdose rate multiplier on change in fentanyl prevalence
    n_fent_delta_od_mult_pp_low <- unlist(df_dsa_owsa["low", "n_fent_delta_od_mult_pp"])
    n_fent_delta_od_mult_pp_high <- unlist(df_dsa_owsa["high", "n_fent_delta_od_mult_pp"])
    names(n_fent_delta_od_mult_pp_low) <- c("n_fent_delta_od_mult")
    names(n_fent_delta_od_mult_pp_high) <- c("n_fent_delta_od_mult")

    # Fatal overdose rate in OAT
    n_fatal_od_oat_pp_low <- unlist(df_dsa_owsa["low", "n_fatal_od_oat_pp"])
    n_fatal_od_oat_pp_high <- unlist(df_dsa_owsa["high", "n_fatal_od_oat_pp"])
    names(n_fatal_od_oat_pp_low) <- c("n_fatal_od_oat")
    names(n_fatal_od_oat_pp_high) <- c("n_fatal_od_oat")

    # Non-overdose mortality rate in OAT
    hr_oat_pp_low <- unlist(df_dsa_owsa["low", "hr_oat_pp"])
    hr_oat_pp_high <- unlist(df_dsa_owsa["high", "hr_oat_pp"])
    names(hr_oat_pp_low) <- c("hr_oat")
    names(hr_oat_pp_high) <- c("hr_oat")

    # Weibull scale (out-of-treatment)
    p_weibull_scale_ou_pp_low <- unlist(df_dsa_owsa["low", "p_weibull_scale_ou_pp"])
    p_weibull_scale_ou_pp_high <- unlist(df_dsa_owsa["high", "p_weibull_scale_ou_pp"])
    names(p_weibull_scale_ou_pp_low) <- c("p_weibull_scale_ou")
    names(p_weibull_scale_ou_pp_high) <- c("p_weibull_scale_ou")

    # Weibull shape (out-of-treatment)
    p_weibull_shape_ou_pp_low <- unlist(df_dsa_owsa["low", "p_weibull_shape_ou_pp"])
    p_weibull_shape_ou_pp_high <- unlist(df_dsa_owsa["high", "p_weibull_shape_ou_pp"])
    names(p_weibull_shape_ou_pp_low) <- c("p_weibull_shape_ou")
    names(p_weibull_shape_ou_pp_high) <- c("p_weibull_shape_ou")

    # Weibull scale (abstinence)
    p_weibull_scale_abs_pp_low <- unlist(df_dsa_owsa["low", "p_weibull_scale_abs_pp"])
    p_weibull_scale_abs_pp_high <- unlist(df_dsa_owsa["high", "p_weibull_scale_abs_pp"])
    names(p_weibull_scale_abs_pp_low) <- c("p_weibull_scale_abs")
    names(p_weibull_scale_abs_pp_high) <- c("p_weibull_scale_abs")

    # Weibull shape (abstinence)
    p_weibull_shape_abs_pp_low <- unlist(df_dsa_owsa["low", "p_weibull_shape_abs_pp"])
    p_weibull_shape_abs_pp_high <- unlist(df_dsa_owsa["high", "p_weibull_shape_abs_pp"])
    names(p_weibull_shape_abs_pp_low) <- c("p_weibull_shape_abs")
    names(p_weibull_shape_abs_pp_high) <- c("p_weibull_shape_abs")

    ##################################
    #### Weibull methadone remain ####
    ##################################
    # Weibull scale methadone (inc)
    p_weibull_scale_met_inc_pp_low <- unlist(df_dsa_owsa["low", "p_weibull_scale_met_inc_pp"])
    p_weibull_scale_met_inc_pp_high <- unlist(df_dsa_owsa["high", "p_weibull_scale_met_inc_pp"])
    names(p_weibull_scale_met_inc_pp_low) <- c("p_weibull_scale_met_inc")
    names(p_weibull_scale_met_inc_pp_high) <- c("p_weibull_scale_met_inc")

    # Weibull shape methadone (inc)
    p_weibull_shape_met_inc_pp_low <- unlist(df_dsa_owsa["low", "p_weibull_shape_met_inc_pp"])
    p_weibull_shape_met_inc_pp_high <- unlist(df_dsa_owsa["high", "p_weibull_shape_met_inc_pp"])
    names(p_weibull_shape_met_inc_pp_low) <- c("p_weibull_shape_met_inc")
    names(p_weibull_shape_met_inc_pp_high) <- c("p_weibull_shape_met_inc")

    # Weibull scale methadone (prev)
    p_weibull_scale_met_prev_pp_low <- unlist(df_dsa_owsa["low", "p_weibull_scale_met_prev_pp"])
    p_weibull_scale_met_prev_pp_high <- unlist(df_dsa_owsa["high", "p_weibull_scale_met_prev_pp"])
    names(p_weibull_scale_met_prev_pp_low) <- c("p_weibull_scale_met_prev")
    names(p_weibull_scale_met_prev_pp_high) <- c("p_weibull_scale_met_prev")

    # Weibull shape (abstinence)
    p_weibull_shape_met_prev_pp_low <- unlist(df_dsa_owsa["low", "p_weibull_shape_met_prev_pp"])
    p_weibull_shape_met_prev_pp_high <- unlist(df_dsa_owsa["high", "p_weibull_shape_met_prev_pp"])
    names(p_weibull_shape_met_prev_pp_low) <- c("p_weibull_shape_met_prev")
    names(p_weibull_shape_met_prev_pp_high) <- c("p_weibull_shape_met_prev")

    # Combine into vector
    # Vector of low values
    v_owsa_low <- c(
      # Non-specific
      p_male_low, hr_dno_ou_oat_dsa_low, n_abs_od_low, n_ou_od_mult_low, n_ou_oat_odn_mult_low, n_ou_oat_odf_mult_low,
      # PP
      hr_tx_pp_inc_dsa_low, hr_tx_pp_prev_dsa_low, hr_death_pp_inc_dsa_low, hr_death_pp_prev_dsa_low,
      # PP
      n_oat_od_pp_low, n_fent_prev_od_mult_pp_low, n_fent_delta_od_mult_pp_low, n_fatal_od_oat_pp_low, hr_oat_pp_low,
      p_weibull_scale_ou_pp_low, p_weibull_shape_ou_pp_low, p_weibull_scale_abs_pp_low, p_weibull_shape_abs_pp_low,
      p_weibull_scale_met_inc_pp_low, p_weibull_shape_met_inc_pp_low, p_weibull_scale_met_prev_pp_low, p_weibull_shape_met_prev_pp_low
    )
    # Vector of high values
    v_owsa_high <- c(
      # Non-specific
      p_male_high, hr_dno_ou_oat_dsa_high, n_abs_od_high, n_ou_od_mult_high, n_ou_oat_odn_mult_high, n_ou_oat_odf_mult_high,
      # PP
      hr_tx_pp_inc_dsa_high, hr_tx_pp_prev_dsa_high, hr_death_pp_inc_dsa_high, hr_death_pp_prev_dsa_high,
      # PP
      n_oat_od_pp_high, n_fent_prev_od_mult_pp_high, n_fent_delta_od_mult_pp_high, n_fatal_od_oat_pp_high, hr_oat_pp_high,
      p_weibull_scale_ou_pp_high, p_weibull_shape_ou_pp_high, p_weibull_scale_abs_pp_high, p_weibull_shape_abs_pp_high,
      p_weibull_scale_met_inc_pp_high, p_weibull_shape_met_inc_pp_high, p_weibull_scale_met_prev_pp_high, p_weibull_shape_met_prev_pp_high
    )
  } else {
    stop("ce_est must be either 'itt_ps' or 'pp'")
  }

  return(list(
    v_owsa_low = v_owsa_low,
    v_owsa_high = v_owsa_high,
    v_cohort_balance_25_inc = v_cohort_balance_25_inc,
    v_cohort_balance_50_inc = v_cohort_balance_50_inc,
    v_cohort_balance_75_inc = v_cohort_balance_75_inc
  ))
}
