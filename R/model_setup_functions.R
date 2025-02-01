#' Markov model
#'
#' \code{markov_model} implements the main model functions to calculate Markov trace.
#'
#' @param l_params_all List with all parameters
#' @param err_stop Logical variable to stop model run if transition array is invalid, if TRUE. Default = FALSE.
#' @param verbose Logical variable to indicate print out of messages. Default = FALSE
#' @param checks Logical variable to indicate output of visual checks (e.g. slices of transition array)
#' @param cali Logical variable to adjust model cutoff to only calibration period
#' @param ce_est Variable to select comparative effectiveness estimate to use (options: ITT_PS; ITT_IV; PP)
#' @return
#' a_TDP: Transition probability array
#' m_M_trace: Fully stratified markov cohort trace
#' m_M_agg_trace: Aggregated markov trace over base health states
#' m_M_agg_trace_death: State-specific mortality from each health state
#' @export
markov_model <- function(l_params_all,
                         err_stop = FALSE,
                         verbose = FALSE,
                         checks = FALSE,
                         # cali = FALSE,
                         time_horizon = time_horizon,
                         ce_est = ce_est,
                         analytic_cohort = analytic_cohort) {
  ### Definition:
  ##   Markov model implementation function
  ### Prefixes:
  ##   l_* denotes list
  ##   n_* denotes number
  ##   a_* denotes 3-D array
  ##   m_* denotes 2-D matrix
  ##   v_* denotes vector
  ##   df_* denotes data frame
  ##   p_* denotes transition parameters
  ##   c_* denotes costs
  ##   u_* denotes utilities
  ### Arguments:
  ##   l_params_all: List with all parameters
  ##   verbose: Logical variable to indicate print out of messages
  ### Returns:
  ##   a_TDP: Transition probability array.
  ##   m_M_trace: Fully disaggregated matrix cohort trace.
  ##   m_M_agg_trace: Aggregated trace over selected base health states.
  ##   m_M_agg_trace_death: State-specific mortality from each health state.
  ##
  with(as.list(l_params_all), {
    #### Set up model states ####
    l_dim_s <- list() # list of health states

    # Base health states
    base <- l_dim_s[[1]] <- c("met", "bnx", "abs", "oum", "oub", "ouo", "odn", "odf")
    # Episodes (1, 2+; for incident and prevalent new users)
    ep <- l_dim_s[[2]] <- c("1", "2")
    # Set model periods
    if (time_horizon == "cali") {
      # Calibration periods
      # n_t <- round(n_cali_per * n_per)
      n_t <- n_per_cali
    } else if (time_horizon == "full") {
      # Maximum model periods(regular)
      # n_t <- round((n_age_max_full - n_age_init) * n_per) # convert years
      n_t <- n_per_full
    } else if (time_horizon == "fent_era") {
      # Maximum model periods(regular)
      n_t <- round((n_age_max_fent - n_age_init) * n_per) # convert years
    } else {
      print("No time horizon selected")
    }

    df_flat <- expand.grid(l_dim_s) # combine all elements together into vector of health states
    df_flat <- dplyr::rename(df_flat,
      base = Var1,
      ep = Var2
    )

    # Create index of states to populate transition matrices
    # All treatment
    oat <- df_flat$base == "bnx" | df_flat$base == "met"
    # All out-of-treatment (incl ABS)
    oot <- df_flat$base == "oum" | df_flat$base == "oub" | df_flat$base == "ouo" | df_flat$base == "odn" #| df_flat$base == "odf"
    # Buprenorphine-naloxone
    bnx <- df_flat$base == "bnx"
    # Methadone
    met <- df_flat$base == "met"
    # Opioid use (following MET discontinuation)
    oum <- df_flat$base == "oum"
    # Opioid use (following BNX discontinuation)
    oub <- df_flat$base == "oub"
    # Opioid use (overall)
    ouo <- df_flat$base == "ouo"
    # All opioid use
    ou <- df_flat$base == "oum" | df_flat$base == "oub" | df_flat$base == "ouo"
    # Overdose
    all_od <- df_flat$base == "odn" | df_flat$base == "odf"
    non_od <- df_flat$base != "odn" & df_flat$base != "odf"
    odn <- df_flat$base == "odn" # non-fatal overdose
    odf <- df_flat$base == "odf" # fatal overdose
    # Abstinence
    abs <- df_flat$base == "abs"
    # Episodes
    ep1 <- df_flat$ep == "1"
    ep2 <- df_flat$ep == "2"

    df_n <- tidyr::unite(df_flat, newCol) # combine columns into one data frame of all health states
    v_n_states <- df_n[, 1] # convert df into vector
    n_states <- length(v_n_states) # total number of health states
    l_index_s <- list(
      oat = oat, oot = oot,
      bnx = bnx,
      met = met,
      oum = oum,
      oub = oub,
      ouo = ouo,
      ou = ou,
      all_od = all_od, odn = odn, odf = odf,
      abs = abs,
      ep1 = ep1, ep2 = ep2
    )

    #### Choose estimate for BNX vs methadone hazard ratio ####
    if (ce_est == "itt_ps") {
      # Treatment discontinuation
      hr_ce_tx_inc <- hr_tx_itt_ps_inc
      hr_ce_tx_prev <- hr_tx_itt_ps_prev
      # Mortality
      hr_ce_death_inc <- hr_death_itt_ps_inc
      hr_ce_death_prev <- hr_death_itt_ps_prev
    } else if (ce_est == "pp_hd") {
      # Treatment discontinuation
      hr_ce_tx_inc <- hr_tx_pp_hd_inc
      hr_ce_tx_prev <- hr_tx_pp_hd_prev
      # Mortality
      hr_ce_death_inc <- hr_death_pp_hd_inc
      hr_ce_death_prev <- hr_death_pp_hd_prev
    } else if (ce_est == "pp") {
      # Treatment discontinuation
      hr_ce_tx_inc <- hr_tx_pp_inc
      hr_ce_tx_prev <- hr_tx_pp_prev
      # Mortality
      hr_ce_death_inc <- hr_death_pp_inc
      hr_ce_death_prev <- hr_death_pp_prev
    } else if (ce_est == "dsa") {
      # Treatment discontinuation
      hr_ce_tx_inc <- hr_tx_dsa_inc
      hr_ce_tx_prev <- hr_tx_dsa_prev
      # Mortality
      hr_ce_death_inc <- hr_death_dsa_inc
      hr_ce_death_prev <- hr_death_dsa_prev
    } else if (ce_est == "cali") {
      # print("Calibration selected")
    } else {
      stop("No comparative effectiveness estimate selected")
    }

    # Module to calculate probability of overdose from states
    # Four separate matrices to account for state-time (first month vs. second+), and model-time (changing fentanyl prevalence, etc.)
    #### Time-dependent overdose probabilities ####
    # Time periods
    time_periods <- length(v_fent_prev)

    # Empty 2-D matrix
    m_odn <- m_odn_first <- m_odf <- m_odf_first <- array(0,
      dim = c(n_states, time_periods),
      dimnames = list(v_n_states, 1:time_periods)
    )

    for (i in 1:time_periods) {
      # Probability of overdose
      # Call user-defined overdose probability function ("overdose_prob_function.R")
      # Non-fatal (first week)
      m_odn_first[bnx, i] <- overdose(l_params_all = l_params_all, rate = n_oat_od, rate_fatal = n_fatal_od_oat, multiplier = n_bnx_od_mult, fent_mult = n_fent_prev_od_mult, fent_delta_mult = n_fent_delta_od_mult, ce_fatal_od_mult = NULL, time = i, first_week = TRUE, fatal = FALSE)
      m_odn_first[met, i] <- overdose(l_params_all = l_params_all, rate = n_oat_od, rate_fatal = n_fatal_od_oat, multiplier = n_met_od_mult, fent_mult = n_fent_prev_od_mult, fent_delta_mult = n_fent_delta_od_mult, ce_fatal_od_mult = NULL, time = i, first_week = TRUE, fatal = FALSE)
      m_odn_first[oum, i] <- overdose(l_params_all = l_params_all, rate = (n_oat_od * n_ou_oat_odn_mult), rate_fatal = (n_fatal_od_oat * n_ou_oat_odf_mult), multiplier = n_ou_od_mult, fent_mult = n_fent_prev_od_mult, fent_delta_mult = n_fent_delta_od_mult, ce_fatal_od_mult = NULL, time = i, first_week = TRUE, fatal = FALSE)
      m_odn_first[oub, i] <- overdose(l_params_all = l_params_all, rate = (n_oat_od * n_ou_oat_odn_mult), rate_fatal = (n_fatal_od_oat * n_ou_oat_odf_mult), multiplier = n_ou_od_mult, fent_mult = n_fent_prev_od_mult, fent_delta_mult = n_fent_delta_od_mult, ce_fatal_od_mult = NULL, time = i, first_week = TRUE, fatal = FALSE)
      m_odn_first[ouo, i] <- overdose(l_params_all = l_params_all, rate = (n_oat_od * n_ou_oat_odn_mult), rate_fatal = (n_fatal_od_oat * n_ou_oat_odf_mult), multiplier = n_ou_od_mult, fent_mult = n_fent_prev_od_mult, fent_delta_mult = n_fent_delta_od_mult, ce_fatal_od_mult = NULL, time = i, first_week = TRUE, fatal = FALSE)
      m_odn_first[odn, i] <- overdose(l_params_all = l_params_all, rate = (n_oat_od * n_ou_oat_odn_mult), rate_fatal = (n_fatal_od_oat * n_ou_oat_odf_mult), multiplier = n_ou_od_mult, fent_mult = n_fent_prev_od_mult, fent_delta_mult = n_fent_delta_od_mult, ce_fatal_od_mult = NULL, time = i, first_week = TRUE, fatal = FALSE)
      m_odn_first[abs, i] <- overdose(l_params_all = l_params_all, rate = n_abs_od, rate_fatal = (n_fatal_od_oat * n_ou_oat_odf_mult), multiplier = n_abs_od_mult, fent_mult = n_fent_prev_od_mult, fent_delta_mult = n_fent_delta_od_mult, ce_fatal_od_mult = NULL, time = i, first_week = TRUE, fatal = FALSE)

      # Fatal (first week)
      if (time_horizon == "cali") {
        m_odf_first[bnx, i] <- overdose(l_params_all = l_params_all, rate = n_oat_od, rate_fatal = n_fatal_od_oat, multiplier = n_bnx_od_mult, fent_mult = n_fent_prev_od_mult, fent_delta_mult = n_fent_delta_od_mult, ce_fatal_od_mult = NULL, time = i, first_week = TRUE, fatal = TRUE)
      } else {
        m_odf_first[bnx & ep1, i] <- overdose(l_params_all = l_params_all, rate = n_oat_od, rate_fatal = n_fatal_od_oat, multiplier = n_bnx_od_mult, fent_mult = n_fent_prev_od_mult, fent_delta_mult = n_fent_delta_od_mult, ce_fatal_od_mult = hr_ce_death_inc, time = i, first_week = TRUE, fatal = TRUE, bnx = TRUE)
        m_odf_first[bnx & ep2, i] <- overdose(l_params_all = l_params_all, rate = n_oat_od, rate_fatal = n_fatal_od_oat, multiplier = n_bnx_od_mult, fent_mult = n_fent_prev_od_mult, fent_delta_mult = n_fent_delta_od_mult, ce_fatal_od_mult = hr_ce_death_prev, time = i, first_week = TRUE, fatal = TRUE, bnx = TRUE)
      }
      # First week of opioid use following dropout from BNX (not applied to week 2+)
      if (time_horizon == "cali") {
        m_odf_first[oub, i] <- overdose(l_params_all = l_params_all, rate = (n_oat_od * n_ou_oat_odn_mult), rate_fatal = (n_fatal_od_oat * n_ou_oat_odf_mult), multiplier = n_ou_od_mult, fent_mult = n_fent_prev_od_mult, fent_delta_mult = n_fent_delta_od_mult, ce_fatal_od_mult = NULL, time = i, first_week = TRUE, fatal = TRUE)
      } else {
        m_odf_first[oub & ep1, i] <- overdose(l_params_all = l_params_all, rate = (n_oat_od * n_ou_oat_odn_mult), rate_fatal = (n_fatal_od_oat * n_ou_oat_odf_mult), multiplier = n_ou_od_mult, fent_mult = n_fent_prev_od_mult, fent_delta_mult = n_fent_delta_od_mult, ce_fatal_od_mult = hr_ce_death_inc, time = i, first_week = TRUE, fatal = TRUE, bnx = TRUE)
        m_odf_first[oub & ep2, i] <- overdose(l_params_all = l_params_all, rate = (n_oat_od * n_ou_oat_odn_mult), rate_fatal = (n_fatal_od_oat * n_ou_oat_odf_mult), multiplier = n_ou_od_mult, fent_mult = n_fent_prev_od_mult, fent_delta_mult = n_fent_delta_od_mult, ce_fatal_od_mult = hr_ce_death_prev, time = i, first_week = TRUE, fatal = TRUE, bnx = TRUE)
      }

      m_odf_first[met, i] <- overdose(l_params_all = l_params_all, rate = n_oat_od, rate_fatal = n_fatal_od_oat, multiplier = n_met_od_mult, fent_mult = n_fent_prev_od_mult, fent_delta_mult = n_fent_delta_od_mult, ce_fatal_od_mult = NULL, time = i, first_week = TRUE, fatal = TRUE)
      m_odf_first[oum, i] <- overdose(l_params_all = l_params_all, rate = (n_oat_od * n_ou_oat_odn_mult), rate_fatal = (n_fatal_od_oat * n_ou_oat_odf_mult), multiplier = n_ou_od_mult, fent_mult = n_fent_prev_od_mult, fent_delta_mult = n_fent_delta_od_mult, ce_fatal_od_mult = NULL, time = i, first_week = TRUE, fatal = TRUE)
      m_odf_first[ouo, i] <- overdose(l_params_all = l_params_all, rate = (n_oat_od * n_ou_oat_odn_mult), rate_fatal = (n_fatal_od_oat * n_ou_oat_odf_mult), multiplier = n_ou_od_mult, fent_mult = n_fent_prev_od_mult, fent_delta_mult = n_fent_delta_od_mult, ce_fatal_od_mult = NULL, time = i, first_week = TRUE, fatal = TRUE)
      m_odf_first[odn, i] <- overdose(l_params_all = l_params_all, rate = (n_oat_od * n_ou_oat_odn_mult), rate_fatal = (n_fatal_od_oat * n_ou_oat_odf_mult), multiplier = n_ou_od_mult, fent_mult = n_fent_prev_od_mult, fent_delta_mult = n_fent_delta_od_mult, ce_fatal_od_mult = NULL, time = i, first_week = TRUE, fatal = TRUE)
      m_odf_first[abs, i] <- overdose(l_params_all = l_params_all, rate = n_abs_od, rate_fatal = (n_fatal_od_oat * n_ou_oat_odf_mult), multiplier = n_abs_od_mult, fent_mult = n_fent_prev_od_mult, fent_delta_mult = n_fent_delta_od_mult, ce_fatal_od_mult = NULL, time = i, first_week = TRUE, fatal = TRUE)

      # Non-fatal (week 2+)
      m_odn[bnx, i] <- overdose(l_params_all = l_params_all, rate = n_oat_od, rate_fatal = n_fatal_od_oat, multiplier = n_bnx_od_mult, fent_mult = n_fent_prev_od_mult, fent_delta_mult = n_fent_delta_od_mult, ce_fatal_od_mult = NULL, time = i, first_week = FALSE, fatal = FALSE)
      m_odn[met, i] <- overdose(l_params_all = l_params_all, rate = n_oat_od, rate_fatal = n_fatal_od_oat, multiplier = n_met_od_mult, fent_mult = n_fent_prev_od_mult, fent_delta_mult = n_fent_delta_od_mult, ce_fatal_od_mult = NULL, time = i, first_week = FALSE, fatal = FALSE)
      m_odn[oum, i] <- overdose(l_params_all = l_params_all, rate = (n_oat_od * n_ou_oat_odn_mult), rate_fatal = (n_fatal_od_oat * n_ou_oat_odf_mult), multiplier = n_ou_od_mult, fent_mult = n_fent_prev_od_mult, fent_delta_mult = n_fent_delta_od_mult, ce_fatal_od_mult = NULL, time = i, first_week = FALSE, fatal = FALSE)
      m_odn[oub, i] <- overdose(l_params_all = l_params_all, rate = (n_oat_od * n_ou_oat_odn_mult), rate_fatal = (n_fatal_od_oat * n_ou_oat_odf_mult), multiplier = n_ou_od_mult, fent_mult = n_fent_prev_od_mult, fent_delta_mult = n_fent_delta_od_mult, ce_fatal_od_mult = NULL, time = i, first_week = FALSE, fatal = FALSE)
      m_odn[ouo, i] <- overdose(l_params_all = l_params_all, rate = (n_oat_od * n_ou_oat_odn_mult), rate_fatal = (n_fatal_od_oat * n_ou_oat_odf_mult), multiplier = n_ou_od_mult, fent_mult = n_fent_prev_od_mult, fent_delta_mult = n_fent_delta_od_mult, ce_fatal_od_mult = NULL, time = i, first_week = FALSE, fatal = FALSE)
      m_odn[odn, i] <- overdose(l_params_all = l_params_all, rate = (n_oat_od * n_ou_oat_odn_mult), rate_fatal = (n_fatal_od_oat * n_ou_oat_odf_mult), multiplier = n_ou_od_mult, fent_mult = n_fent_prev_od_mult, fent_delta_mult = n_fent_delta_od_mult, ce_fatal_od_mult = NULL, time = i, first_week = FALSE, fatal = FALSE)
      m_odn[abs, i] <- overdose(l_params_all = l_params_all, rate = n_abs_od, rate_fatal = (n_fatal_od_oat * n_ou_oat_odf_mult), multiplier = n_abs_od_mult, fent_mult = n_fent_prev_od_mult, fent_delta_mult = n_fent_delta_od_mult, ce_fatal_od_mult = NULL, time = i, first_week = FALSE, fatal = FALSE)

      # Fatal (week 2+)
      if (time_horizon == "cali") {
        m_odf[bnx, i] <- overdose(l_params_all = l_params_all, rate = n_oat_od, rate_fatal = n_fatal_od_oat, multiplier = n_bnx_od_mult, fent_mult = n_fent_prev_od_mult, fent_delta_mult = n_fent_delta_od_mult, ce_fatal_od_mult = NULL, time = i, first_week = FALSE, fatal = TRUE)
      } else {
        m_odf[bnx & ep1, i] <- overdose(l_params_all = l_params_all, rate = n_oat_od, rate_fatal = n_fatal_od_oat, multiplier = n_bnx_od_mult, fent_mult = n_fent_prev_od_mult, fent_delta_mult = n_fent_delta_od_mult, ce_fatal_od_mult = hr_ce_death_inc, time = i, first_week = FALSE, fatal = TRUE, bnx = TRUE)
        m_odf[bnx & ep2, i] <- overdose(l_params_all = l_params_all, rate = n_oat_od, rate_fatal = n_fatal_od_oat, multiplier = n_bnx_od_mult, fent_mult = n_fent_prev_od_mult, fent_delta_mult = n_fent_delta_od_mult, ce_fatal_od_mult = hr_ce_death_prev, time = i, first_week = FALSE, fatal = TRUE, bnx = TRUE)
      }

      m_odf[met, i] <- overdose(l_params_all = l_params_all, rate = n_oat_od, rate_fatal = n_fatal_od_oat, multiplier = n_met_od_mult, fent_mult = n_fent_prev_od_mult, fent_delta_mult = n_fent_delta_od_mult, ce_fatal_od_mult = NULL, time = i, first_week = FALSE, fatal = TRUE)
      m_odf[oum, i] <- overdose(l_params_all = l_params_all, rate = (n_oat_od * n_ou_oat_odn_mult), rate_fatal = (n_fatal_od_oat * n_ou_oat_odf_mult), multiplier = n_ou_od_mult, fent_mult = n_fent_prev_od_mult, fent_delta_mult = n_fent_delta_od_mult, ce_fatal_od_mult = NULL, time = i, first_week = FALSE, fatal = TRUE)
      m_odf[oub, i] <- overdose(l_params_all = l_params_all, rate = (n_oat_od * n_ou_oat_odn_mult), rate_fatal = (n_fatal_od_oat * n_ou_oat_odf_mult), multiplier = n_ou_od_mult, fent_mult = n_fent_prev_od_mult, fent_delta_mult = n_fent_delta_od_mult, ce_fatal_od_mult = NULL, time = i, first_week = FALSE, fatal = TRUE)
      m_odf[ouo, i] <- overdose(l_params_all = l_params_all, rate = (n_oat_od * n_ou_oat_odn_mult), rate_fatal = (n_fatal_od_oat * n_ou_oat_odf_mult), multiplier = n_ou_od_mult, fent_mult = n_fent_prev_od_mult, fent_delta_mult = n_fent_delta_od_mult, ce_fatal_od_mult = NULL, time = i, first_week = FALSE, fatal = TRUE)
      m_odf[odn, i] <- overdose(l_params_all = l_params_all, rate = (n_oat_od * n_ou_oat_odn_mult), rate_fatal = (n_fatal_od_oat * n_ou_oat_odf_mult), multiplier = n_ou_od_mult, fent_mult = n_fent_prev_od_mult, fent_delta_mult = n_fent_delta_od_mult, ce_fatal_od_mult = NULL, time = i, first_week = FALSE, fatal = TRUE)
      m_odf[abs, i] <- overdose(l_params_all = l_params_all, rate = n_abs_od, rate_fatal = (n_fatal_od_oat * n_ou_oat_odf_mult), multiplier = n_abs_od_mult, fent_mult = n_fent_prev_od_mult, fent_delta_mult = n_fent_delta_od_mult, ce_fatal_od_mult = NULL, time = i, first_week = FALSE, fatal = TRUE)
    }

    if (checks) {
      write.csv(m_odn_first, "checks/overdose/m_odn_first.csv", row.names = TRUE)
      write.csv(m_odf_first, "checks/overdose/m_odf_first.csv", row.names = TRUE)
      write.csv(m_odn, "checks/overdose/m_odn.csv", row.names = TRUE)
      write.csv(m_odf, "checks/overdose/m_odf.csv", row.names = TRUE)
    } else {

    }

    # Probability of non-overdose
    m_non_od_first <- 1 - (m_odn_first + m_odf_first)
    m_non_od <- 1 - (m_odn + m_odf)

    #### Time-dependent remain probabilities ####
    a_remain <- array(0,
      dim = c(n_states, time_periods, n_t),
      dimnames = list(v_n_states, 1:time_periods, 1:n_t)
    )

    # Probability of remaining in health state
    # All remain in fatal overdose, remain probability = 1
    # Remove frailty terms for episode (separate weibull estimates for incident (EP1) and prevalent users (EP2))
    for (j in 1:time_periods) {
      for (i in 1:n_t) {
        # Incident users
        # Integration of comparative effectiveness estimate for BNX v MET treatment discontinuation
        if (time_horizon == "cali") {
          a_remain[ep1 & bnx, j, i] <- 0
          a_remain[ep1 & met, j, i] <- as.vector(exp(p_weibull_scale_oat_inc * (((i - 1)^p_weibull_shape_oat_inc) - (i^p_weibull_shape_oat_inc)))) # "MET" becomes "all-OAT" state for calibration
        } else {
          a_remain[ep1 & bnx, j, i] <- as.vector(exp(hr_ce_tx_inc * p_weibull_scale_met_inc * (((i - 1)^p_weibull_shape_met_inc) - (i^p_weibull_shape_met_inc))))
          a_remain[ep1 & met, j, i] <- as.vector(exp(p_weibull_scale_met_inc * (((i - 1)^p_weibull_shape_met_inc) - (i^p_weibull_shape_met_inc))))
        }
        a_remain[ep1 & abs, j, i] <- as.vector(exp(p_weibull_scale_abs * (((i - 1)^p_weibull_shape_abs) - (i^p_weibull_shape_abs))))
        a_remain[ep1 & oum, j, i] <- as.vector(exp(p_weibull_scale_ou * (((i - 1)^p_weibull_shape_ou) - (i^p_weibull_shape_ou))))
        a_remain[ep1 & oub, j, i] <- as.vector(exp(p_weibull_scale_ou * (((i - 1)^p_weibull_shape_ou) - (i^p_weibull_shape_ou))))
        a_remain[ep1 & ouo, j, i] <- as.vector(exp(p_weibull_scale_ou * (((i - 1)^p_weibull_shape_ou) - (i^p_weibull_shape_ou))))

        a_remain[ep1 & odn, j, i] <- 0
        a_remain[ep1 & odf, j, i] <- 1

        # Prevalent users
        # Integration of comparative effectiveness estimate for BNX v MET treatment discontinuation
        if (time_horizon == "cali") {
          a_remain[ep2 & bnx, j, i] <- 0
          a_remain[ep2 & met, j, i] <- as.vector(exp(p_weibull_scale_oat_prev * (((i - 1)^p_weibull_shape_oat_prev) - (i^p_weibull_shape_oat_prev))))
        } else {
          a_remain[ep2 & bnx, j, i] <- as.vector(exp(hr_ce_tx_prev * p_weibull_scale_met_prev * (((i - 1)^p_weibull_shape_met_prev) - (i^p_weibull_shape_met_prev))))
          a_remain[ep2 & met, j, i] <- as.vector(exp(p_weibull_scale_met_prev * (((i - 1)^p_weibull_shape_met_prev) - (i^p_weibull_shape_met_prev))))
        }
        a_remain[ep2 & abs, j, i] <- as.vector(exp(p_weibull_scale_abs * (((i - 1)^p_weibull_shape_abs) - (i^p_weibull_shape_abs))))
        a_remain[ep2 & oum, j, i] <- as.vector(exp(p_weibull_scale_ou * (((i - 1)^p_weibull_shape_ou) - (i^p_weibull_shape_ou))))
        a_remain[ep2 & oub, j, i] <- as.vector(exp(p_weibull_scale_ou * (((i - 1)^p_weibull_shape_ou) - (i^p_weibull_shape_ou))))
        a_remain[ep2 & ouo, j, i] <- as.vector(exp(p_weibull_scale_ou * (((i - 1)^p_weibull_shape_ou) - (i^p_weibull_shape_ou))))

        a_remain[ep2 & odn, j, i] <- 0
        a_remain[ep2 & odf, j, i] <- 1
      }
    }
    # Modify TDP for non-overdose
    # Probability of state-exit
    m_remain_1 <- a_remain[, 1, ] # 2010
    m_remain_2 <- a_remain[, 2, ] # 2011
    m_remain_3 <- a_remain[, 3, ] # 2012
    m_remain_4 <- a_remain[, 4, ] # 2013
    m_remain_5 <- a_remain[, 5, ] # 2014
    m_remain_6 <- a_remain[, 6, ] # 2015
    m_remain_7 <- a_remain[, 7, ] # 2016
    m_remain_8 <- a_remain[, 8, ] # Calibration 2017
    m_remain_9 <- a_remain[, 9, ] # Calibration 2018
    m_remain_10 <- a_remain[, 10, ] # Calibration 2019
    m_remain_11 <- a_remain[, 11, ] # Calibration 2020

    # Adjust remain probability in opioid use to account for new treatment initiations
    m_remain_1[ou, ] <- m_remain_1[ou, ] * (1 - p_inc_entry_adj_2010)
    m_remain_2[ou, ] <- m_remain_2[ou, ] * (1 - p_inc_entry_adj_2011)
    m_remain_3[ou, ] <- m_remain_3[ou, ] * (1 - p_inc_entry_adj_2012)
    m_remain_4[ou, ] <- m_remain_4[ou, ] * (1 - p_inc_entry_adj_2013)
    m_remain_5[ou, ] <- m_remain_5[ou, ] * (1 - p_inc_entry_adj_2014)
    m_remain_6[ou, ] <- m_remain_6[ou, ] * (1 - p_inc_entry_adj_2015)
    m_remain_7[ou, ] <- m_remain_7[ou, ] * (1 - p_inc_entry_adj_2016)
    m_remain_8[ou, ] <- m_remain_8[ou, ] * (1 - p_inc_entry_adj_2017)
    m_remain_9[ou, ] <- m_remain_9[ou, ] * (1 - p_inc_entry_adj_2018)
    m_remain_10[ou, ] <- m_remain_10[ou, ] * (1 - p_inc_entry_adj_2019)
    m_remain_11[ou, ] <- m_remain_11[ou, ] * (1 - p_inc_entry_adj_2020)

    m_leave_1 <- 1 - m_remain_1 # 2010
    m_leave_2 <- 1 - m_remain_2 # 2011
    m_leave_3 <- 1 - m_remain_3 # 2012
    m_leave_4 <- 1 - m_remain_4 # 2013
    m_leave_5 <- 1 - m_remain_5 # 2014
    m_leave_6 <- 1 - m_remain_6 # 2015
    m_leave_7 <- 1 - m_remain_7 # 2016
    m_leave_8 <- 1 - m_remain_8 # 2017
    m_leave_9 <- 1 - m_remain_9 # 2018
    m_leave_10 <- 1 - m_remain_10 # 2019
    m_leave_11 <- 1 - m_remain_11 # 2020

    if (checks) {
      # Time dependent state-exit probabilities (from weibull estimates)
      write.csv(m_remain_1, "checks/transitions/m_tdp_1.csv", row.names = TRUE)
      write.csv(m_remain_11, "checks/transitions/m_tdp_11.csv", row.names = TRUE)

      write.csv(m_leave_1, "checks/transitions/m_leave_1.csv", row.names = TRUE)
      write.csv(m_leave_11, "checks/transitions/m_leave_11.csv", row.names = TRUE)
    } else {

    }

    # Mortality CE estimate
    # Call user-defined mortality calculation function ("mortality_function.R")
    # BNX and methadone equivalent for calibration
    if (time_horizon == "cali") {
      v_mort_bnx_inc <- mort(l_params_all = l_params_all, hr = hr_oat, per = n_per)
      v_mort_bnx_prev <- mort(l_params_all = l_params_all, hr = hr_oat, per = n_per)
      # Mortality adjustment for BNX in analytic runs
    } else {
      v_mort_bnx_inc <- mort(l_params_all = l_params_all, hr = (hr_ce_death_inc * hr_oat), per = n_per)
      v_mort_bnx_prev <- mort(l_params_all = l_params_all, hr = (hr_ce_death_prev * hr_oat), per = n_per)
    }
    # Background mortality
    v_mort_met <- mort(l_params_all = l_params_all, hr = hr_oat, per = n_per)
    v_mort_oum <- mort(l_params_all = l_params_all, hr = (hr_oat * hr_dno_ou_oat), per = n_per)
    v_mort_oub <- mort(l_params_all = l_params_all, hr = (hr_oat * hr_dno_ou_oat), per = n_per)
    v_mort_ouo <- mort(l_params_all = l_params_all, hr = (hr_oat * hr_dno_ou_oat), per = n_per)
    v_mort_odn <- mort(l_params_all = l_params_all, hr = (hr_oat * hr_dno_ou_oat), per = n_per) # Background mortality equal for REL/ODN (fatal overdoses counted separately)
    v_mort_odf <- rep(0, n_t) # stay in ODF, not tracked in "death"
    v_mort_abs <- mort(l_params_all = l_params_all, hr = hr_abs, per = n_per)

    # Create empty mortality matrix
    m_mort <- array(0,
      dim = c(n_states, n_t),
      dimnames = list(v_n_states, 1:n_t)
    )
    # Populate mortality matrix (death probability from each state)
    for (i in 1:n_t) {
      m_mort[ep1 & bnx, i] <- v_mort_bnx_inc[i]
      m_mort[ep2 & bnx, i] <- v_mort_bnx_prev[i]
      m_mort[met, i] <- v_mort_met[i]
      m_mort[oum, i] <- v_mort_ouo[i]
      m_mort[oub, i] <- v_mort_ouo[i]
      m_mort[ouo, i] <- v_mort_ouo[i]
      m_mort[odn, i] <- v_mort_odn[i] # using background excess mortality for relapse in non-fatal overdose
      m_mort[odf, i] <- v_mort_odf[i] # transition to death = 0
      m_mort[abs, i] <- v_mort_abs[i]
    }

    if (checks) {
      # Mortality matrix
      write.csv(m_mort, "checks/mortality/m_mort.csv", row.names = TRUE)
    } else {

    }

    # Alive probability in each period
    m_alive <- 1 - m_mort

    #### Unconditional transition probabilities ####
    # Create as array (different probabilities for model-time varying overdose)
    a_up <- a_up_first <- array(0,
      dim = c(n_states, n_states, time_periods),
      dimnames = list(v_n_states, v_n_states, 1:time_periods)
    )

    # Overdose probability populated first, accounting for higher probability of overdose transition in first month
    for (i in 1:time_periods) {
      # From BNX
      # First week
      a_up_first[bnx & ep1, met & ep1, i] <- p_bnx_met_inc
      a_up_first[bnx & ep1, abs & ep1, i] <- p_bnx_abs_inc
      a_up_first[bnx & ep1, oub & ep1, i] <- p_bnx_oub_inc
      # Week 2+
      a_up[bnx & ep1, met & ep1, i] <- p_bnx_met_inc
      a_up[bnx & ep1, abs & ep1, i] <- p_bnx_abs_inc
      a_up[bnx & ep1, oub & ep1, i] <- p_bnx_oub_inc
      # First week
      a_up_first[bnx & ep2, met & ep2, i] <- p_bnx_met_prev
      a_up_first[bnx & ep2, abs & ep2, i] <- p_bnx_abs_prev
      a_up_first[bnx & ep2, oub & ep2, i] <- p_bnx_oub_prev
      # Week 2+
      a_up[bnx & ep2, met & ep2, i] <- p_bnx_met_prev
      a_up[bnx & ep2, abs & ep2, i] <- p_bnx_abs_prev
      a_up[bnx & ep2, oub & ep2, i] <- p_bnx_oub_prev

      # From MET
      # First week
      a_up_first[met & ep1, bnx & ep1, i] <- p_met_bnx_inc
      a_up_first[met & ep1, abs & ep1, i] <- p_met_abs_inc
      a_up_first[met & ep1, oum & ep1, i] <- p_met_oum_inc
      # Week 2+
      a_up[met & ep1, bnx & ep1, i] <- p_met_bnx_inc
      a_up[met & ep1, abs & ep1, i] <- p_met_abs_inc
      a_up[met & ep1, oum & ep1, i] <- p_met_oum_inc
      # First week
      a_up_first[met & ep2, bnx & ep2, i] <- p_met_bnx_prev
      a_up_first[met & ep2, abs & ep2, i] <- p_met_abs_prev
      a_up_first[met & ep2, oum & ep2, i] <- p_met_oum_prev
      # Week 2+
      a_up[met & ep2, bnx & ep2, i] <- p_met_bnx_prev
      a_up[met & ep2, abs & ep2, i] <- p_met_abs_prev
      a_up[met & ep2, oum & ep2, i] <- p_met_oum_prev

      # From ABS
      # First week
      a_up_first[abs & ep1, ouo & ep1, i] <- p_abs_ouo_inc
      a_up_first[abs & ep1, met & ep1, i] <- p_abs_met_inc
      a_up_first[abs & ep1, bnx & ep1, i] <- p_abs_bnx_inc
      # Week 2+
      a_up[abs & ep1, ouo & ep1, i] <- p_abs_ouo_inc
      a_up[abs & ep1, met & ep1, i] <- p_abs_met_inc
      a_up[abs & ep1, bnx & ep1, i] <- p_abs_bnx_inc
      # First week
      a_up_first[abs & ep2, ouo & ep2, i] <- p_abs_ouo_prev
      a_up_first[abs & ep2, met & ep2, i] <- p_abs_met_prev
      a_up_first[abs & ep2, bnx & ep2, i] <- p_abs_bnx_prev
      # Week 2+
      a_up[abs & ep2, ouo & ep2, i] <- p_abs_ouo_prev
      a_up[abs & ep2, met & ep2, i] <- p_abs_met_prev
      a_up[abs & ep2, bnx & ep2, i] <- p_abs_bnx_prev

      # From OUM
      # First week
      a_up_first[oum & ep1, met & ep1, i] <- p_oum_met_inc
      a_up_first[oum & ep1, bnx & ep1, i] <- p_oum_bnx_inc
      a_up_first[oum & ep1, abs & ep1, i] <- p_oum_abs_inc
      # Week 2+
      a_up[oum & ep1, met & ep1, i] <- p_oum_met_inc
      a_up[oum & ep1, bnx & ep1, i] <- p_oum_bnx_inc
      a_up[oum & ep1, abs & ep1, i] <- p_oum_abs_inc
      # First week
      a_up_first[oum & ep2, met & ep2, i] <- p_oum_met_prev
      a_up_first[oum & ep2, bnx & ep2, i] <- p_oum_bnx_prev
      a_up_first[oum & ep2, abs & ep2, i] <- p_oum_abs_prev
      # Week 2+
      a_up[oum & ep2, met & ep2, i] <- p_oum_met_prev
      a_up[oum & ep2, bnx & ep2, i] <- p_oum_bnx_prev
      a_up[oum & ep2, abs & ep2, i] <- p_oum_abs_prev

      # From OUB
      # First week
      a_up_first[oub & ep1, met & ep1, i] <- p_oub_met_inc
      a_up_first[oub & ep1, bnx & ep1, i] <- p_oub_bnx_inc
      a_up_first[oub & ep1, abs & ep1, i] <- p_oub_abs_inc
      # Week 2+
      a_up[oub & ep1, met & ep1, i] <- p_oub_met_inc
      a_up[oub & ep1, bnx & ep1, i] <- p_oub_bnx_inc
      a_up[oub & ep1, abs & ep1, i] <- p_oub_abs_inc
      # First week
      a_up_first[oub & ep2, met & ep2, i] <- p_oub_met_prev
      a_up_first[oub & ep2, bnx & ep2, i] <- p_oub_bnx_prev
      a_up_first[oub & ep2, abs & ep2, i] <- p_oub_abs_prev
      # Week 2+
      a_up[oub & ep2, met & ep2, i] <- p_oub_met_prev
      a_up[oub & ep2, bnx & ep2, i] <- p_oub_bnx_prev
      a_up[oub & ep2, abs & ep2, i] <- p_oub_abs_prev

      # From OUO
      # First week
      a_up_first[ouo & ep1, met & ep1, i] <- p_ouo_met_inc
      a_up_first[ouo & ep1, bnx & ep1, i] <- p_ouo_bnx_inc
      a_up_first[ouo & ep1, abs & ep1, i] <- p_ouo_abs_inc
      # Week 2+
      a_up[ouo & ep1, met & ep1, i] <- p_ouo_met_inc
      a_up[ouo & ep1, bnx & ep1, i] <- p_ouo_bnx_inc
      a_up[ouo & ep1, abs & ep1, i] <- p_ouo_abs_inc
      # First week
      a_up_first[ouo & ep2, met & ep2, i] <- p_ouo_met_prev
      a_up_first[ouo & ep2, bnx & ep2, i] <- p_ouo_bnx_prev
      a_up_first[ouo & ep2, abs & ep2, i] <- p_ouo_abs_prev
      # Week 2+
      a_up[ouo & ep2, met & ep2, i] <- p_ouo_met_prev
      a_up[ouo & ep2, bnx & ep2, i] <- p_ouo_bnx_prev
      a_up[ouo & ep2, abs & ep2, i] <- p_ouo_abs_prev

      # From OD (all time periods same)
      a_up[odn & ep1, met & ep1, i] <- a_up_first[odn & ep1, met & ep1, i] <- p_odn_met_inc
      a_up[odn & ep1, bnx & ep1, i] <- a_up_first[odn & ep1, bnx & ep1, i] <- p_odn_bnx_inc
      a_up[odn & ep1, abs & ep1, i] <- a_up_first[odn & ep1, abs & ep1, i] <- p_odn_abs_inc
      a_up[odn & ep1, ouo & ep1, i] <- a_up_first[odn & ep1, ouo & ep1, i] <- p_odn_ouo_inc
      a_up[odn & ep2, met & ep2, i] <- a_up_first[odn & ep2, met & ep2, i] <- p_odn_met_prev
      a_up[odn & ep2, bnx & ep2, i] <- a_up_first[odn & ep2, bnx & ep2, i] <- p_odn_bnx_prev
      a_up[odn & ep2, abs & ep2, i] <- a_up_first[odn & ep2, abs & ep2, i] <- p_odn_abs_prev
      a_up[odn & ep2, ouo & ep2, i] <- a_up_first[odn & ep2, ouo & ep2, i] <- p_odn_ouo_prev
    }

    if (checks) {
      # Time dependent state-exit probabilities (from weibull estimates)
      write.csv(a_up_first[, , 1], "checks/transitions/a_up_first_1.csv", row.names = TRUE)
      write.csv(a_up_first[, , 11], "checks/transitions/a_up_first_11.csv", row.names = TRUE)
      write.csv(a_up[, , 1], "checks/transitions/a_up_1.csv", row.names = TRUE)
      write.csv(a_up[, , 11], "checks/transitions/a_up_11.csv", row.names = TRUE)
    } else {

    }

    #### Create full time-dependent transition array ####
    # Empty 3-D array
    # Model set up for 10 years of calendar time variation
    a_tdp_1 <- a_tdp_2 <- a_tdp_3 <- a_tdp_4 <- a_tdp_5 <- a_tdp_6 <- a_tdp_7 <- a_tdp_8 <- a_tdp_9 <- a_tdp_10 <- a_tdp_11 <- array(0,
      dim = c(n_states, n_states, n_t),
      dimnames = list(v_n_states, v_n_states, 1:n_t)
    )

    # Add transitions conditional on state-exit (m_leave = 1 - remain)
    # Multiple arrays to account for model-time-varying overdose parameters
    # Modified transitions for first month (state-time)
    for (i in 1) {
      a_tdp_1[, , i] <- a_up_first[, , 1] * m_leave_1[, i]
      a_tdp_2[, , i] <- a_up_first[, , 2] * m_leave_2[, i]
      a_tdp_3[, , i] <- a_up_first[, , 3] * m_leave_3[, i]
      a_tdp_4[, , i] <- a_up_first[, , 4] * m_leave_4[, i]
      a_tdp_5[, , i] <- a_up_first[, , 5] * m_leave_5[, i]
      a_tdp_6[, , i] <- a_up_first[, , 6] * m_leave_6[, i]
      a_tdp_7[, , i] <- a_up_first[, , 7] * m_leave_7[, i]
      a_tdp_8[, , i] <- a_up_first[, , 8] * m_leave_8[, i]
      a_tdp_9[, , i] <- a_up_first[, , 9] * m_leave_9[, i]
      a_tdp_10[, , i] <- a_up_first[, , 10] * m_leave_10[, i]
      a_tdp_11[, , i] <- a_up_first[, , 11] * m_leave_11[, i]
    }
    ## All transitions 2+ months
    for (i in 2:n_t) {
      a_tdp_1[, , i] <- a_up[, , 1] * m_leave_1[, i]
      a_tdp_2[, , i] <- a_up[, , 2] * m_leave_2[, i]
      a_tdp_3[, , i] <- a_up[, , 3] * m_leave_3[, i]
      a_tdp_4[, , i] <- a_up[, , 4] * m_leave_4[, i]
      a_tdp_5[, , i] <- a_up[, , 5] * m_leave_5[, i]
      a_tdp_6[, , i] <- a_up[, , 6] * m_leave_6[, i]
      a_tdp_7[, , i] <- a_up[, , 7] * m_leave_7[, i]
      a_tdp_8[, , i] <- a_up[, , 8] * m_leave_8[, i]
      a_tdp_9[, , i] <- a_up[, , 9] * m_leave_9[, i]
      a_tdp_10[, , i] <- a_up[, , 10] * m_leave_10[, i]
      a_tdp_11[, , i] <- a_up[, , 11] * m_leave_11[, i]
    }

    # Add time-dependent remain probabilities
    # Note change here for cohort rebalancing
    for (i in 1:n_t) {
      # Incidence user
      # BNX
      a_tdp_1[ep1 & bnx, ep1 & bnx, i] <- m_remain_1[ep1 & bnx, i]
      a_tdp_2[ep1 & bnx, ep1 & bnx, i] <- m_remain_2[ep1 & bnx, i]
      a_tdp_3[ep1 & bnx, ep1 & bnx, i] <- m_remain_3[ep1 & bnx, i]
      a_tdp_4[ep1 & bnx, ep1 & bnx, i] <- m_remain_4[ep1 & bnx, i]
      a_tdp_5[ep1 & bnx, ep1 & bnx, i] <- m_remain_5[ep1 & bnx, i]
      a_tdp_6[ep1 & bnx, ep1 & bnx, i] <- m_remain_6[ep1 & bnx, i]
      a_tdp_7[ep1 & bnx, ep1 & bnx, i] <- m_remain_7[ep1 & bnx, i]
      a_tdp_8[ep1 & bnx, ep1 & bnx, i] <- m_remain_8[ep1 & bnx, i]
      a_tdp_9[ep1 & bnx, ep1 & bnx, i] <- m_remain_9[ep1 & bnx, i]
      a_tdp_10[ep1 & bnx, ep1 & bnx, i] <- m_remain_10[ep1 & bnx, i]
      a_tdp_11[ep1 & bnx, ep1 & bnx, i] <- m_remain_11[ep1 & bnx, i]
      # MET
      a_tdp_1[ep1 & met, ep1 & met, i] <- m_remain_1[ep1 & met, i]
      a_tdp_2[ep1 & met, ep1 & met, i] <- m_remain_2[ep1 & met, i]
      a_tdp_3[ep1 & met, ep1 & met, i] <- m_remain_3[ep1 & met, i]
      a_tdp_4[ep1 & met, ep1 & met, i] <- m_remain_4[ep1 & met, i]
      a_tdp_5[ep1 & met, ep1 & met, i] <- m_remain_5[ep1 & met, i]
      a_tdp_6[ep1 & met, ep1 & met, i] <- m_remain_6[ep1 & met, i]
      a_tdp_7[ep1 & met, ep1 & met, i] <- m_remain_7[ep1 & met, i]
      a_tdp_8[ep1 & met, ep1 & met, i] <- m_remain_8[ep1 & met, i]
      a_tdp_9[ep1 & met, ep1 & met, i] <- m_remain_9[ep1 & met, i]
      a_tdp_10[ep1 & met, ep1 & met, i] <- m_remain_10[ep1 & met, i]
      a_tdp_11[ep1 & met, ep1 & met, i] <- m_remain_11[ep1 & met, i]
      # ABS
      a_tdp_1[ep1 & abs, ep1 & abs, i] <- m_remain_1[ep1 & abs, i]
      a_tdp_2[ep1 & abs, ep1 & abs, i] <- m_remain_2[ep1 & abs, i]
      a_tdp_3[ep1 & abs, ep1 & abs, i] <- m_remain_3[ep1 & abs, i]
      a_tdp_4[ep1 & abs, ep1 & abs, i] <- m_remain_4[ep1 & abs, i]
      a_tdp_5[ep1 & abs, ep1 & abs, i] <- m_remain_5[ep1 & abs, i]
      a_tdp_6[ep1 & abs, ep1 & abs, i] <- m_remain_6[ep1 & abs, i]
      a_tdp_7[ep1 & abs, ep1 & abs, i] <- m_remain_7[ep1 & abs, i]
      a_tdp_8[ep1 & abs, ep1 & abs, i] <- m_remain_8[ep1 & abs, i]
      a_tdp_9[ep1 & abs, ep1 & abs, i] <- m_remain_9[ep1 & abs, i]
      a_tdp_10[ep1 & abs, ep1 & abs, i] <- m_remain_10[ep1 & abs, i]
      a_tdp_11[ep1 & abs, ep1 & abs, i] <- m_remain_11[ep1 & abs, i]
      # OUM
      a_tdp_1[ep1 & oum, ep1 & oum, i] <- m_remain_1[ep1 & oum, i]
      a_tdp_2[ep1 & oum, ep1 & oum, i] <- m_remain_2[ep1 & oum, i]
      a_tdp_3[ep1 & oum, ep1 & oum, i] <- m_remain_3[ep1 & oum, i]
      a_tdp_4[ep1 & oum, ep1 & oum, i] <- m_remain_4[ep1 & oum, i]
      a_tdp_5[ep1 & oum, ep1 & oum, i] <- m_remain_5[ep1 & oum, i]
      a_tdp_6[ep1 & oum, ep1 & oum, i] <- m_remain_6[ep1 & oum, i]
      a_tdp_7[ep1 & oum, ep1 & oum, i] <- m_remain_7[ep1 & oum, i]
      a_tdp_8[ep1 & oum, ep1 & oum, i] <- m_remain_8[ep1 & oum, i]
      a_tdp_9[ep1 & oum, ep1 & oum, i] <- m_remain_9[ep1 & oum, i]
      a_tdp_10[ep1 & oum, ep1 & oum, i] <- m_remain_10[ep1 & oum, i]
      a_tdp_11[ep1 & oum, ep1 & oum, i] <- m_remain_11[ep1 & oum, i]
      # OUB
      a_tdp_1[ep1 & oub, ep1 & oub, i] <- m_remain_1[ep1 & oub, i]
      a_tdp_2[ep1 & oub, ep1 & oub, i] <- m_remain_2[ep1 & oub, i]
      a_tdp_3[ep1 & oub, ep1 & oub, i] <- m_remain_3[ep1 & oub, i]
      a_tdp_4[ep1 & oub, ep1 & oub, i] <- m_remain_4[ep1 & oub, i]
      a_tdp_5[ep1 & oub, ep1 & oub, i] <- m_remain_5[ep1 & oub, i]
      a_tdp_6[ep1 & oub, ep1 & oub, i] <- m_remain_6[ep1 & oub, i]
      a_tdp_7[ep1 & oub, ep1 & oub, i] <- m_remain_7[ep1 & oub, i]
      a_tdp_8[ep1 & oub, ep1 & oub, i] <- m_remain_8[ep1 & oub, i]
      a_tdp_9[ep1 & oub, ep1 & oub, i] <- m_remain_9[ep1 & oub, i]
      a_tdp_10[ep1 & oub, ep1 & oub, i] <- m_remain_10[ep1 & oub, i]
      a_tdp_11[ep1 & oub, ep1 & oub, i] <- m_remain_11[ep1 & oub, i]
      # OUO
      a_tdp_1[ep1 & ouo, ep1 & ouo, i] <- m_remain_1[ep1 & ouo, i]
      a_tdp_2[ep1 & ouo, ep1 & ouo, i] <- m_remain_2[ep1 & ouo, i]
      a_tdp_3[ep1 & ouo, ep1 & ouo, i] <- m_remain_3[ep1 & ouo, i]
      a_tdp_4[ep1 & ouo, ep1 & ouo, i] <- m_remain_4[ep1 & ouo, i]
      a_tdp_5[ep1 & ouo, ep1 & ouo, i] <- m_remain_5[ep1 & ouo, i]
      a_tdp_6[ep1 & ouo, ep1 & ouo, i] <- m_remain_6[ep1 & ouo, i]
      a_tdp_7[ep1 & ouo, ep1 & ouo, i] <- m_remain_7[ep1 & ouo, i]
      a_tdp_8[ep1 & ouo, ep1 & ouo, i] <- m_remain_8[ep1 & ouo, i]
      a_tdp_9[ep1 & ouo, ep1 & ouo, i] <- m_remain_9[ep1 & ouo, i]
      a_tdp_10[ep1 & ouo, ep1 & ouo, i] <- m_remain_10[ep1 & ouo, i]
      a_tdp_11[ep1 & ouo, ep1 & ouo, i] <- m_remain_11[ep1 & ouo, i]
      # ODF
      a_tdp_1[ep1 & odf, ep1 & odf, i] <- m_remain_1[ep1 & odf, i]
      a_tdp_2[ep1 & odf, ep1 & odf, i] <- m_remain_2[ep1 & odf, i]
      a_tdp_3[ep1 & odf, ep1 & odf, i] <- m_remain_3[ep1 & odf, i]
      a_tdp_4[ep1 & odf, ep1 & odf, i] <- m_remain_4[ep1 & odf, i]
      a_tdp_5[ep1 & odf, ep1 & odf, i] <- m_remain_5[ep1 & odf, i]
      a_tdp_6[ep1 & odf, ep1 & odf, i] <- m_remain_6[ep1 & odf, i]
      a_tdp_7[ep1 & odf, ep1 & odf, i] <- m_remain_7[ep1 & odf, i]
      a_tdp_8[ep1 & odf, ep1 & odf, i] <- m_remain_8[ep1 & odf, i]
      a_tdp_9[ep1 & odf, ep1 & odf, i] <- m_remain_9[ep1 & odf, i]
      a_tdp_10[ep1 & odf, ep1 & odf, i] <- m_remain_10[ep1 & odf, i]
      a_tdp_11[ep1 & odf, ep1 & odf, i] <- m_remain_11[ep1 & odf, i]

      # Prevalent user
      # BNX
      a_tdp_1[ep2 & bnx, ep2 & bnx, i] <- m_remain_1[ep2 & bnx, i]
      a_tdp_2[ep2 & bnx, ep2 & bnx, i] <- m_remain_2[ep2 & bnx, i]
      a_tdp_3[ep2 & bnx, ep2 & bnx, i] <- m_remain_3[ep2 & bnx, i]
      a_tdp_4[ep2 & bnx, ep2 & bnx, i] <- m_remain_4[ep2 & bnx, i]
      a_tdp_5[ep2 & bnx, ep2 & bnx, i] <- m_remain_5[ep2 & bnx, i]
      a_tdp_6[ep2 & bnx, ep2 & bnx, i] <- m_remain_6[ep2 & bnx, i]
      a_tdp_7[ep2 & bnx, ep2 & bnx, i] <- m_remain_7[ep2 & bnx, i]
      a_tdp_8[ep2 & bnx, ep2 & bnx, i] <- m_remain_8[ep2 & bnx, i]
      a_tdp_9[ep2 & bnx, ep2 & bnx, i] <- m_remain_9[ep2 & bnx, i]
      a_tdp_10[ep2 & bnx, ep2 & bnx, i] <- m_remain_10[ep2 & bnx, i]
      a_tdp_11[ep2 & bnx, ep2 & bnx, i] <- m_remain_11[ep2 & bnx, i]
      # MET
      a_tdp_1[ep2 & met, ep2 & met, i] <- m_remain_1[ep2 & met, i]
      a_tdp_2[ep2 & met, ep2 & met, i] <- m_remain_2[ep2 & met, i]
      a_tdp_3[ep2 & met, ep2 & met, i] <- m_remain_3[ep2 & met, i]
      a_tdp_4[ep2 & met, ep2 & met, i] <- m_remain_4[ep2 & met, i]
      a_tdp_5[ep2 & met, ep2 & met, i] <- m_remain_5[ep2 & met, i]
      a_tdp_6[ep2 & met, ep2 & met, i] <- m_remain_6[ep2 & met, i]
      a_tdp_7[ep2 & met, ep2 & met, i] <- m_remain_7[ep2 & met, i]
      a_tdp_8[ep2 & met, ep2 & met, i] <- m_remain_8[ep2 & met, i]
      a_tdp_9[ep2 & met, ep2 & met, i] <- m_remain_9[ep2 & met, i]
      a_tdp_10[ep2 & met, ep2 & met, i] <- m_remain_10[ep2 & met, i]
      a_tdp_11[ep2 & met, ep2 & met, i] <- m_remain_11[ep2 & met, i]
      # ABS
      a_tdp_1[ep2 & abs, ep2 & abs, i] <- m_remain_1[ep2 & abs, i]
      a_tdp_2[ep2 & abs, ep2 & abs, i] <- m_remain_2[ep2 & abs, i]
      a_tdp_3[ep2 & abs, ep2 & abs, i] <- m_remain_3[ep2 & abs, i]
      a_tdp_4[ep2 & abs, ep2 & abs, i] <- m_remain_4[ep2 & abs, i]
      a_tdp_5[ep2 & abs, ep2 & abs, i] <- m_remain_5[ep2 & abs, i]
      a_tdp_6[ep2 & abs, ep2 & abs, i] <- m_remain_6[ep2 & abs, i]
      a_tdp_7[ep2 & abs, ep2 & abs, i] <- m_remain_7[ep2 & abs, i]
      a_tdp_8[ep2 & abs, ep2 & abs, i] <- m_remain_8[ep2 & abs, i]
      a_tdp_9[ep2 & abs, ep2 & abs, i] <- m_remain_9[ep2 & abs, i]
      a_tdp_10[ep2 & abs, ep2 & abs, i] <- m_remain_10[ep2 & abs, i]
      a_tdp_11[ep2 & abs, ep2 & abs, i] <- m_remain_11[ep2 & abs, i]
      # OUM
      a_tdp_1[ep2 & oum, ep2 & oum, i] <- m_remain_1[ep2 & oum, i]
      a_tdp_2[ep2 & oum, ep2 & oum, i] <- m_remain_2[ep2 & oum, i]
      a_tdp_3[ep2 & oum, ep2 & oum, i] <- m_remain_3[ep2 & oum, i]
      a_tdp_4[ep2 & oum, ep2 & oum, i] <- m_remain_4[ep2 & oum, i]
      a_tdp_5[ep2 & oum, ep2 & oum, i] <- m_remain_5[ep2 & oum, i]
      a_tdp_6[ep2 & oum, ep2 & oum, i] <- m_remain_6[ep2 & oum, i]
      a_tdp_7[ep2 & oum, ep2 & oum, i] <- m_remain_7[ep2 & oum, i]
      a_tdp_8[ep2 & oum, ep2 & oum, i] <- m_remain_8[ep2 & oum, i]
      a_tdp_9[ep2 & oum, ep2 & oum, i] <- m_remain_9[ep2 & oum, i]
      a_tdp_10[ep2 & oum, ep2 & oum, i] <- m_remain_10[ep2 & oum, i]
      a_tdp_11[ep2 & oum, ep2 & oum, i] <- m_remain_11[ep2 & oum, i]
      # OUB
      a_tdp_1[ep2 & oub, ep2 & oub, i] <- m_remain_1[ep2 & oub, i]
      a_tdp_2[ep2 & oub, ep2 & oub, i] <- m_remain_2[ep2 & oub, i]
      a_tdp_3[ep2 & oub, ep2 & oub, i] <- m_remain_3[ep2 & oub, i]
      a_tdp_4[ep2 & oub, ep2 & oub, i] <- m_remain_4[ep2 & oub, i]
      a_tdp_5[ep2 & oub, ep2 & oub, i] <- m_remain_5[ep2 & oub, i]
      a_tdp_6[ep2 & oub, ep2 & oub, i] <- m_remain_6[ep2 & oub, i]
      a_tdp_7[ep2 & oub, ep2 & oub, i] <- m_remain_7[ep2 & oub, i]
      a_tdp_8[ep2 & oub, ep2 & oub, i] <- m_remain_8[ep2 & oub, i]
      a_tdp_9[ep2 & oub, ep2 & oub, i] <- m_remain_9[ep2 & oub, i]
      a_tdp_10[ep2 & oub, ep2 & oub, i] <- m_remain_10[ep2 & oub, i]
      a_tdp_11[ep2 & oub, ep2 & oub, i] <- m_remain_11[ep2 & oub, i]
      # OUO
      a_tdp_1[ep2 & ouo, ep2 & ouo, i] <- m_remain_1[ep2 & ouo, i]
      a_tdp_2[ep2 & ouo, ep2 & ouo, i] <- m_remain_2[ep2 & ouo, i]
      a_tdp_3[ep2 & ouo, ep2 & ouo, i] <- m_remain_3[ep2 & ouo, i]
      a_tdp_4[ep2 & ouo, ep2 & ouo, i] <- m_remain_4[ep2 & ouo, i]
      a_tdp_5[ep2 & ouo, ep2 & ouo, i] <- m_remain_5[ep2 & ouo, i]
      a_tdp_6[ep2 & ouo, ep2 & ouo, i] <- m_remain_6[ep2 & ouo, i]
      a_tdp_7[ep2 & ouo, ep2 & ouo, i] <- m_remain_7[ep2 & ouo, i]
      a_tdp_8[ep2 & ouo, ep2 & ouo, i] <- m_remain_8[ep2 & ouo, i]
      a_tdp_9[ep2 & ouo, ep2 & ouo, i] <- m_remain_9[ep2 & ouo, i]
      a_tdp_10[ep2 & ouo, ep2 & ouo, i] <- m_remain_10[ep2 & ouo, i]
      a_tdp_11[ep2 & ouo, ep2 & ouo, i] <- m_remain_11[ep2 & ouo, i]
      # ODF
      a_tdp_1[ep2 & odf, ep2 & odf, i] <- m_remain_1[ep2 & odf, i]
      a_tdp_2[ep2 & odf, ep2 & odf, i] <- m_remain_2[ep2 & odf, i]
      a_tdp_3[ep2 & odf, ep2 & odf, i] <- m_remain_3[ep2 & odf, i]
      a_tdp_4[ep2 & odf, ep2 & odf, i] <- m_remain_4[ep2 & odf, i]
      a_tdp_5[ep2 & odf, ep2 & odf, i] <- m_remain_5[ep2 & odf, i]
      a_tdp_6[ep2 & odf, ep2 & odf, i] <- m_remain_6[ep2 & odf, i]
      a_tdp_7[ep2 & odf, ep2 & odf, i] <- m_remain_7[ep2 & odf, i]
      a_tdp_8[ep2 & odf, ep2 & odf, i] <- m_remain_8[ep2 & odf, i]
      a_tdp_9[ep2 & odf, ep2 & odf, i] <- m_remain_9[ep2 & odf, i]
      a_tdp_10[ep2 & odf, ep2 & odf, i] <- m_remain_10[ep2 & odf, i]
      a_tdp_11[ep2 & odf, ep2 & odf, i] <- m_remain_11[ep2 & odf, i]
    }

    if (checks) {
      # State transitions pre-overdose
      # First month (state-time)
      write.csv(a_tdp_1[, , 1], "checks/transitions/a_tdp_pre_od_2010.csv")
      write.csv(a_tdp_11[, , 1], "checks/transitions/a_tdp_pre_od_2020.csv")
    }

    # Modify array for non-overdose prob
    for (i in 1) {
      a_tdp_1[, , i] <- a_tdp_1[, , i] * m_non_od_first[, 1]
      a_tdp_2[, , i] <- a_tdp_2[, , i] * m_non_od_first[, 2]
      a_tdp_3[, , i] <- a_tdp_3[, , i] * m_non_od_first[, 3]
      a_tdp_4[, , i] <- a_tdp_4[, , i] * m_non_od_first[, 4]
      a_tdp_5[, , i] <- a_tdp_5[, , i] * m_non_od_first[, 5]
      a_tdp_6[, , i] <- a_tdp_6[, , i] * m_non_od_first[, 6]
      a_tdp_7[, , i] <- a_tdp_7[, , i] * m_non_od_first[, 7]
      a_tdp_8[, , i] <- a_tdp_8[, , i] * m_non_od_first[, 8]
      a_tdp_9[, , i] <- a_tdp_9[, , i] * m_non_od_first[, 9]
      a_tdp_10[, , i] <- a_tdp_10[, , i] * m_non_od_first[, 10]
      a_tdp_11[, , i] <- a_tdp_11[, , i] * m_non_od_first[, 11]
    }
    ## All transitions 2+ months
    for (i in 2:n_t) {
      a_tdp_1[, , i] <- a_tdp_1[, , i] * m_non_od[, 1]
      a_tdp_2[, , i] <- a_tdp_2[, , i] * m_non_od[, 2]
      a_tdp_3[, , i] <- a_tdp_3[, , i] * m_non_od[, 3]
      a_tdp_4[, , i] <- a_tdp_4[, , i] * m_non_od[, 4]
      a_tdp_5[, , i] <- a_tdp_5[, , i] * m_non_od[, 5]
      a_tdp_6[, , i] <- a_tdp_6[, , i] * m_non_od[, 6]
      a_tdp_7[, , i] <- a_tdp_7[, , i] * m_non_od[, 7]
      a_tdp_8[, , i] <- a_tdp_8[, , i] * m_non_od[, 8]
      a_tdp_9[, , i] <- a_tdp_9[, , i] * m_non_od[, 9]
      a_tdp_10[, , i] <- a_tdp_10[, , i] * m_non_od[, 10]
      a_tdp_11[, , i] <- a_tdp_11[, , i] * m_non_od[, 11]
    }

    # Add overdose probabilities
    for (i in 1) {
      a_tdp_1[, odn, i] <- m_odn_first[, 1]
      a_tdp_2[, odn, i] <- m_odn_first[, 2]
      a_tdp_3[, odn, i] <- m_odn_first[, 3]
      a_tdp_4[, odn, i] <- m_odn_first[, 4]
      a_tdp_5[, odn, i] <- m_odn_first[, 5]
      a_tdp_6[, odn, i] <- m_odn_first[, 6]
      a_tdp_7[, odn, i] <- m_odn_first[, 7]
      a_tdp_8[, odn, i] <- m_odn_first[, 8]
      a_tdp_9[, odn, i] <- m_odn_first[, 9]
      a_tdp_10[, odn, i] <- m_odn_first[, 10]
      a_tdp_11[, odn, i] <- m_odn_first[, 11]

      a_tdp_1[, odf, i] <- m_odf_first[, 1]
      a_tdp_2[, odf, i] <- m_odf_first[, 2]
      a_tdp_3[, odf, i] <- m_odf_first[, 3]
      a_tdp_4[, odf, i] <- m_odf_first[, 4]
      a_tdp_5[, odf, i] <- m_odf_first[, 5]
      a_tdp_6[, odf, i] <- m_odf_first[, 6]
      a_tdp_7[, odf, i] <- m_odf_first[, 7]
      a_tdp_8[, odf, i] <- m_odf_first[, 8]
      a_tdp_9[, odf, i] <- m_odf_first[, 9]
      a_tdp_10[, odf, i] <- m_odf_first[, 10]
      a_tdp_11[, odf, i] <- m_odf_first[, 11]

      a_tdp_1[odf, odf, i] <- a_tdp_2[odf, odf, i] <- a_tdp_3[odf, odf, i] <- a_tdp_4[odf, odf, i] <- a_tdp_5[odf, odf, i] <- 1
      a_tdp_6[odf, odf, i] <- a_tdp_7[odf, odf, i] <- a_tdp_8[odf, odf, i] <- a_tdp_9[odf, odf, i] <- a_tdp_10[odf, odf, i] <- a_tdp_11[odf, odf, i] <- 1

      # For cohort rebalancing
      # Incident to prevalent
      a_tdp_1[ep1, ep2, i] <- a_tdp_1[ep2, ep2, i]
      a_tdp_2[ep1, ep2, i] <- a_tdp_2[ep2, ep2, i]
      a_tdp_3[ep1, ep2, i] <- a_tdp_3[ep2, ep2, i]
      a_tdp_4[ep1, ep2, i] <- a_tdp_4[ep2, ep2, i]
      a_tdp_5[ep1, ep2, i] <- a_tdp_5[ep2, ep2, i]
      a_tdp_6[ep1, ep2, i] <- a_tdp_6[ep2, ep2, i]
      a_tdp_7[ep1, ep2, i] <- a_tdp_7[ep2, ep2, i]
      a_tdp_8[ep1, ep2, i] <- a_tdp_8[ep2, ep2, i]
      a_tdp_9[ep1, ep2, i] <- a_tdp_9[ep2, ep2, i]
      a_tdp_10[ep1, ep2, i] <- a_tdp_10[ep2, ep2, i]
      a_tdp_11[ep1, ep2, i] <- a_tdp_11[ep2, ep2, i]
      # Prevalent to incident
      a_tdp_1[ep2, ep1, i] <- a_tdp_1[ep1, ep1, i]
      a_tdp_2[ep2, ep1, i] <- a_tdp_2[ep1, ep1, i]
      a_tdp_3[ep2, ep1, i] <- a_tdp_3[ep1, ep1, i]
      a_tdp_4[ep2, ep1, i] <- a_tdp_4[ep1, ep1, i]
      a_tdp_5[ep2, ep1, i] <- a_tdp_5[ep1, ep1, i]
      a_tdp_6[ep2, ep1, i] <- a_tdp_6[ep1, ep1, i]
      a_tdp_7[ep2, ep1, i] <- a_tdp_7[ep1, ep1, i]
      a_tdp_8[ep2, ep1, i] <- a_tdp_8[ep1, ep1, i]
      a_tdp_9[ep2, ep1, i] <- a_tdp_9[ep1, ep1, i]
      a_tdp_10[ep2, ep1, i] <- a_tdp_10[ep1, ep1, i]
      a_tdp_11[ep2, ep1, i] <- a_tdp_11[ep1, ep1, i]
    }
    ## All transitions 2+ months
    for (i in 2:n_t) {
      # Non-injection
      a_tdp_1[, odn, i] <- m_odn[, 1]
      a_tdp_2[, odn, i] <- m_odn[, 2]
      a_tdp_3[, odn, i] <- m_odn[, 3]
      a_tdp_4[, odn, i] <- m_odn[, 4]
      a_tdp_5[, odn, i] <- m_odn[, 5]
      a_tdp_6[, odn, i] <- m_odn[, 6]
      a_tdp_7[, odn, i] <- m_odn[, 7]
      a_tdp_8[, odn, i] <- m_odn[, 8]
      a_tdp_9[, odn, i] <- m_odn[, 9]
      a_tdp_10[, odn, i] <- m_odn[, 10]
      a_tdp_11[, odn, i] <- m_odn[, 11]

      # Non-injection
      a_tdp_1[, odf, i] <- m_odf[, 1]
      a_tdp_2[, odf, i] <- m_odf[, 2]
      a_tdp_3[, odf, i] <- m_odf[, 3]
      a_tdp_4[, odf, i] <- m_odf[, 4]
      a_tdp_5[, odf, i] <- m_odf[, 5]
      a_tdp_6[, odf, i] <- m_odf[, 6]
      a_tdp_7[, odf, i] <- m_odf[, 7]
      a_tdp_8[, odf, i] <- m_odf[, 8]
      a_tdp_9[, odf, i] <- m_odf[, 9]
      a_tdp_10[, odf, i] <- m_odf[, 10]
      a_tdp_11[, odf, i] <- m_odf[, 11]

      a_tdp_1[odf, odf, i] <- a_tdp_2[odf, odf, i] <- a_tdp_3[odf, odf, i] <- a_tdp_4[odf, odf, i] <- a_tdp_5[odf, odf, i] <- 1
      a_tdp_6[odf, odf, i] <- a_tdp_7[odf, odf, i] <- a_tdp_8[odf, odf, i] <- a_tdp_9[odf, odf, i] <- a_tdp_10[odf, odf, i] <- a_tdp_11[odf, odf, i] <- 1

      # For cohort rebalancing
      # Incident to prevalent
      a_tdp_1[ep1, ep2, i] <- a_tdp_1[ep2, ep2, i]
      a_tdp_2[ep1, ep2, i] <- a_tdp_2[ep2, ep2, i]
      a_tdp_3[ep1, ep2, i] <- a_tdp_3[ep2, ep2, i]
      a_tdp_4[ep1, ep2, i] <- a_tdp_4[ep2, ep2, i]
      a_tdp_5[ep1, ep2, i] <- a_tdp_5[ep2, ep2, i]
      a_tdp_6[ep1, ep2, i] <- a_tdp_6[ep2, ep2, i]
      a_tdp_7[ep1, ep2, i] <- a_tdp_7[ep2, ep2, i]
      a_tdp_8[ep1, ep2, i] <- a_tdp_8[ep2, ep2, i]
      a_tdp_9[ep1, ep2, i] <- a_tdp_9[ep2, ep2, i]
      a_tdp_10[ep1, ep2, i] <- a_tdp_10[ep2, ep2, i]
      a_tdp_11[ep1, ep2, i] <- a_tdp_11[ep2, ep2, i]
      # Prevalent to incident
      a_tdp_1[ep2, ep1, i] <- a_tdp_1[ep1, ep1, i]
      a_tdp_2[ep2, ep1, i] <- a_tdp_2[ep1, ep1, i]
      a_tdp_3[ep2, ep1, i] <- a_tdp_3[ep1, ep1, i]
      a_tdp_4[ep2, ep1, i] <- a_tdp_4[ep1, ep1, i]
      a_tdp_5[ep2, ep1, i] <- a_tdp_5[ep1, ep1, i]
      a_tdp_6[ep2, ep1, i] <- a_tdp_6[ep1, ep1, i]
      a_tdp_7[ep2, ep1, i] <- a_tdp_7[ep1, ep1, i]
      a_tdp_8[ep2, ep1, i] <- a_tdp_8[ep1, ep1, i]
      a_tdp_9[ep2, ep1, i] <- a_tdp_9[ep1, ep1, i]
      a_tdp_10[ep2, ep1, i] <- a_tdp_10[ep1, ep1, i]
      a_tdp_11[ep2, ep1, i] <- a_tdp_11[ep1, ep1, i]

      if (analytic_cohort == "bnx_only") {
        a_tdp_1[met, , i] <- a_tdp_2[met, , i] <- a_tdp_3[met, , i] <- a_tdp_4[met, , i] <- a_tdp_5[met, , i] <- a_tdp_6[met, , i] <- 0
        a_tdp_7[met, , i] <- a_tdp_8[met, , i] <- a_tdp_9[met, , i] <- a_tdp_10[met, , i] <- a_tdp_11[met, , i] <- 0
        a_tdp_1[oum, , i] <- a_tdp_2[oum, , i] <- a_tdp_3[oum, , i] <- a_tdp_4[oum, , i] <- a_tdp_5[oum, , i] <- a_tdp_6[oum, , i] <- 0
        a_tdp_7[oum, , i] <- a_tdp_8[oum, , i] <- a_tdp_9[oum, , i] <- a_tdp_10[oum, , i] <- a_tdp_11[oum, , i] <- 0
      } else if (analytic_cohort == "met_only") {
        a_tdp_1[bnx, , i] <- a_tdp_2[bnx, , i] <- a_tdp_3[bnx, , i] <- a_tdp_4[bnx, , i] <- a_tdp_5[bnx, , i] <- a_tdp_6[bnx, , i] <- 0
        a_tdp_7[bnx, , i] <- a_tdp_8[bnx, , i] <- a_tdp_9[bnx, , i] <- a_tdp_10[bnx, , i] <- a_tdp_11[bnx, , i] <- 0
        a_tdp_1[oub, , i] <- a_tdp_2[oub, , i] <- a_tdp_3[oub, , i] <- a_tdp_4[oub, , i] <- a_tdp_5[oub, , i] <- a_tdp_6[oub, , i] <- 0
        a_tdp_7[oub, , i] <- a_tdp_8[oub, , i] <- a_tdp_9[oub, , i] <- a_tdp_10[oub, , i] <- a_tdp_11[oub, , i] <- 0
      } else if (analytic_cohort == "cali") {
        a_tdp_1[bnx, , i] <- a_tdp_2[bnx, , i] <- a_tdp_3[bnx, , i] <- a_tdp_4[bnx, , i] <- a_tdp_5[bnx, , i] <- a_tdp_6[bnx, , i] <- 0
        a_tdp_7[bnx, , i] <- a_tdp_8[bnx, , i] <- a_tdp_9[bnx, , i] <- a_tdp_10[bnx, , i] <- a_tdp_11[bnx, , i] <- 0
        a_tdp_1[oub, , i] <- a_tdp_2[oub, , i] <- a_tdp_3[oub, , i] <- a_tdp_4[oub, , i] <- a_tdp_5[oub, , i] <- a_tdp_6[oub, , i] <- 0
        a_tdp_7[oub, , i] <- a_tdp_8[oub, , i] <- a_tdp_9[oub, , i] <- a_tdp_10[oub, , i] <- a_tdp_11[oub, , i] <- 0
      } else {
        print("No cohort selected")
      }
    }

    if (checks) {
      # State transitions post-overdose
      # First month (state-time)
      write.csv(a_tdp_1[, , 1], "checks/transitions/a_tdp_post_od_2010_wk1.csv")
      write.csv(a_tdp_7[, , 26], "checks/transitions/a_tdp_post_od_2016_wk26.csv")
      write.csv(a_tdp_11[, , 1], "checks/transitions/a_tdp_post_od_2020_wk1.csv")
    }

    #### Cohort rebalancing to mimic cohort dynamics over time
    # NOTE: Cohort re-balancing in both directions is only to capture years in which there are more incident than prevalent users
    # To re-balance in years where there were substantially more incidenct users (2011 & 2016)
    # This is a net outcome, so one of the transitions must be zero (either inc -> prev, or prev -> inc)
    # 2010
    # From incident to prevalent
    a_tdp_1[ep1, ep2, ] <- a_tdp_1[ep1, ep2, ] * p_inc_prev_2010 # move
    a_tdp_1[ep1, ep1, ] <- a_tdp_1[ep1, ep1, ] * (1 - p_inc_prev_2010) # stay
    # From prevalent to incident
    a_tdp_1[ep2, ep1, ] <- a_tdp_1[ep2, ep1, ] * p_prev_inc_2010 # move
    a_tdp_1[ep2, ep2, ] <- a_tdp_1[ep2, ep2, ] * (1 - p_prev_inc_2010) # stay

    # 2011
    # From incident to prevalent
    a_tdp_2[ep1, ep2, ] <- a_tdp_2[ep1, ep2, ] * p_inc_prev_2011 # move
    a_tdp_2[ep1, ep1, ] <- a_tdp_2[ep1, ep1, ] * (1 - p_inc_prev_2011) # stay
    # From prevalent to incident
    a_tdp_2[ep2, ep1, ] <- a_tdp_2[ep2, ep1, ] * p_prev_inc_2011 # move
    a_tdp_2[ep2, ep2, ] <- a_tdp_2[ep2, ep2, ] * (1 - p_prev_inc_2011) # stay

    # 2012
    # From incident to prevalent
    a_tdp_3[ep1, ep2, ] <- a_tdp_3[ep1, ep2, ] * p_inc_prev_2012 # move
    a_tdp_3[ep1, ep1, ] <- a_tdp_3[ep1, ep1, ] * (1 - p_inc_prev_2012) # stay
    # From prevalent to incident
    a_tdp_3[ep2, ep1, ] <- a_tdp_3[ep2, ep1, ] * p_prev_inc_2012 # move
    a_tdp_3[ep2, ep2, ] <- a_tdp_3[ep2, ep2, ] * (1 - p_prev_inc_2012) # stay

    # 2013
    # From incident to prevalent
    a_tdp_4[ep1, ep2, ] <- a_tdp_4[ep1, ep2, ] * p_inc_prev_2013 # move
    a_tdp_4[ep1, ep1, ] <- a_tdp_4[ep1, ep1, ] * (1 - p_inc_prev_2013) # stay
    # From prevalent to incident
    a_tdp_4[ep2, ep1, ] <- a_tdp_4[ep2, ep1, ] * p_prev_inc_2013 # move
    a_tdp_4[ep2, ep2, ] <- a_tdp_4[ep2, ep2, ] * (1 - p_prev_inc_2013) # stay

    # 2014
    # From incident to prevalent
    a_tdp_5[ep1, ep2, ] <- a_tdp_5[ep1, ep2, ] * p_inc_prev_2014 # move
    a_tdp_5[ep1, ep1, ] <- a_tdp_5[ep1, ep1, ] * (1 - p_inc_prev_2014) # stay
    # From prevalent to incident
    a_tdp_5[ep2, ep1, ] <- a_tdp_5[ep2, ep1, ] * p_prev_inc_2014 # move
    a_tdp_5[ep2, ep2, ] <- a_tdp_5[ep2, ep2, ] * (1 - p_prev_inc_2014) # stay

    # 2015
    # From incident to prevalent
    a_tdp_6[ep1, ep2, ] <- a_tdp_6[ep1, ep2, ] * p_inc_prev_2015 # move
    a_tdp_6[ep1, ep1, ] <- a_tdp_6[ep1, ep1, ] * (1 - p_inc_prev_2015) # stay
    # From prevalent to incident
    a_tdp_6[ep2, ep1, ] <- a_tdp_6[ep2, ep1, ] * p_prev_inc_2015 # move
    a_tdp_6[ep2, ep2, ] <- a_tdp_6[ep2, ep2, ] * (1 - p_prev_inc_2015) # stay

    # 2016
    # From incident to prevalent
    a_tdp_7[ep1, ep2, ] <- a_tdp_7[ep1, ep2, ] * p_inc_prev_2016 # move
    a_tdp_7[ep1, ep1, ] <- a_tdp_7[ep1, ep1, ] * (1 - p_inc_prev_2016) # stay
    # From prevalent to incident
    a_tdp_7[ep2, ep1, ] <- a_tdp_7[ep2, ep1, ] * p_prev_inc_2016 # move
    a_tdp_7[ep2, ep2, ] <- a_tdp_7[ep2, ep2, ] * (1 - p_prev_inc_2016) # stay

    # 2017
    # From incident to prevalent
    a_tdp_8[ep1, ep2, ] <- a_tdp_8[ep1, ep2, ] * p_inc_prev_2017 # move
    a_tdp_8[ep1, ep1, ] <- a_tdp_8[ep1, ep1, ] * (1 - p_inc_prev_2017) # stay
    # From prevalent to incident
    a_tdp_8[ep2, ep1, ] <- a_tdp_8[ep2, ep1, ] * p_prev_inc_2017 # move
    a_tdp_8[ep2, ep2, ] <- a_tdp_8[ep2, ep2, ] * (1 - p_prev_inc_2017) # stay

    # 2018
    # From incident to prevalent
    a_tdp_9[ep1, ep2, ] <- a_tdp_9[ep1, ep2, ] * p_inc_prev_2018 # move
    a_tdp_9[ep1, ep1, ] <- a_tdp_9[ep1, ep1, ] * (1 - p_inc_prev_2018) # stay
    # From prevalent to incident
    a_tdp_9[ep2, ep1, ] <- a_tdp_9[ep2, ep1, ] * p_prev_inc_2018 # move
    a_tdp_9[ep2, ep2, ] <- a_tdp_9[ep2, ep2, ] * (1 - p_prev_inc_2018) # stay

    # 2019
    # From incident to prevalent
    a_tdp_10[ep1, ep2, ] <- a_tdp_10[ep1, ep2, ] * p_inc_prev_2019 # move
    a_tdp_10[ep1, ep1, ] <- a_tdp_10[ep1, ep1, ] * (1 - p_inc_prev_2019) # stay
    # From prevalent to incident
    a_tdp_10[ep2, ep1, ] <- a_tdp_10[ep2, ep1, ] * p_prev_inc_2019 # move
    a_tdp_10[ep2, ep2, ] <- a_tdp_10[ep2, ep2, ] * (1 - p_prev_inc_2019) # stay

    # 2020
    # From incident to prevalent
    a_tdp_11[ep1, ep2, ] <- a_tdp_11[ep1, ep2, ] * p_inc_prev_2020 # move
    a_tdp_11[ep1, ep1, ] <- a_tdp_11[ep1, ep1, ] * (1 - p_inc_prev_2020) # stay
    # From prevalent to incident
    a_tdp_11[ep2, ep1, ] <- a_tdp_11[ep2, ep1, ] * p_prev_inc_2020 # move
    a_tdp_11[ep2, ep2, ] <- a_tdp_11[ep2, ep2, ] * (1 - p_prev_inc_2020) # stay

    if (checks) {
      # Full array at time = 1
      array_2010_1w <- a_tdp_1[, , 1]
      array_2016_1w <- a_tdp_7[, , 1]
      array_2020_1w <- a_tdp_11[, , 1]
      write.csv(array_2010_1w, "checks/full array/array_2010_1w.csv", row.names = TRUE)
      write.csv(array_2016_1w, "checks/full array/array_2016_1w.csv", row.names = TRUE)
      write.csv(array_2020_1w, "checks/full array/array_2020_1w.csv", row.names = TRUE)

      # Full array at time = max
      array_2010_52w <- a_tdp_1[, , 52]
      array_2016_52w <- a_tdp_7[, , 52]
      array_2020_52w <- a_tdp_11[, , 52]
      write.csv(array_2010_52w, "checks/full array/array_2010_52w.csv", row.names = TRUE)
      write.csv(array_2016_52w, "checks/full array/array_2016_52w.csv", row.names = TRUE)
      write.csv(array_2020_52w, "checks/full array/array_2020_52w.csv", row.names = TRUE)
    } else {

    }

    #### Check transition array ####
    # check_transition_probability(a_P = a_TDP, err_stop = err_stop, verbose = verbose) # check all probs [0, 1]
    # check_sum_of_transition_array(a_P = a_TDP, n_states = n_states, n_t = n_t, err_stop = err_stop, verbose = verbose) # check prob sums = 1

    #### Run Markov model ####
    # Create empty initial state vectors
    v_s_init <- rep(0, n_states)
    names(v_s_init) <- v_n_states

    #### Set initial state vector ####
    # v_s_init[bnx & ep1] <- v_init_dist$bnx_inc
    # v_s_init[met & ep1] <- v_init_dist$met_inc
    # v_s_init[oum & ep1] <- v_init_dist$oum_inc
    # v_s_init[oub & ep1] <- v_init_dist$oub_inc
    # v_s_init[ouo & ep1] <- v_init_dist$ouo_inc
    # v_s_init[odn & ep1] <- v_init_dist$odn_inc
    # v_s_init[odf & ep1] <- v_init_dist$odf_inc
    # v_s_init[abs & ep1] <- v_init_dist$abs_inc

    # v_s_init[bnx & ep2] <- v_init_dist$bnx_prev
    # v_s_init[met & ep2] <- v_init_dist$met_prev
    # v_s_init[oum & ep2] <- v_init_dist$oum_prev
    # v_s_init[oub & ep2] <- v_init_dist$oub_prev
    # v_s_init[ouo & ep2] <- v_init_dist$ouo_prev
    # v_s_init[odn & ep2] <- v_init_dist$odn_prev
    # v_s_init[odf & ep2] <- v_init_dist$odf_prev
    # v_s_init[abs & ep2] <- v_init_dist$abs_prev

    v_s_init[bnx & ep1] <- bnx_inc
    v_s_init[met & ep1] <- met_inc
    v_s_init[oum & ep1] <- oum_inc
    v_s_init[oub & ep1] <- oub_inc
    v_s_init[ouo & ep1] <- ouo_inc
    v_s_init[odn & ep1] <- odn_inc
    v_s_init[odf & ep1] <- odf_inc
    v_s_init[abs & ep1] <- abs_inc

    v_s_init[bnx & ep2] <- bnx_prev
    v_s_init[met & ep2] <- met_prev
    v_s_init[oum & ep2] <- oum_prev
    v_s_init[oub & ep2] <- oub_prev
    v_s_init[ouo & ep2] <- ouo_prev
    v_s_init[odn & ep2] <- odn_prev
    v_s_init[odf & ep2] <- odf_prev
    v_s_init[abs & ep2] <- abs_prev

    if (checks) {
      write.csv(v_s_init, "checks/initial/v_s_init.csv", row.names = TRUE)
    }

    # Create Markov Trace
    # Initialize population
    a_m_trace <- array(0,
      dim = c((n_t + 1), n_states, (n_t + 1)),
      dimnames = list(0:n_t, v_n_states, 0:n_t)
    )
    a_m_trace[1, , 1] <- v_s_init

    if (checks) {
      write.csv(a_m_trace[, , 1], "checks/initial/m_trace_init.csv")
    }

    # Calbrate 2012-2020 (to capture first year fentanyl)
    if (time_horizon == "cali") {
      # All model time periods
      # Year 1 (2012)
      for (i in 2:51) {
        # Time spent in given health state
        # First month (state-time)
        for (j in 1:(i - 1)) {
          # state-time-dependent transition probability (j) * age (model-time)-specific mortality (i) * model-time-specific overdose (track in separate matrix)

          # Conditional state transition matrix at time (j) ~ sojourn time is continuous time spent in health state
          # Remove deaths from all states (depends only on cohort age (i) and health state-specific SMR
          m_sojourn <- a_tdp_3[, , j] * m_alive[, i - 1] # state-time transition matrix at state-time j, re-weighted for model-time (age) varying mortality at each time point

          # Current health states from prior cycle
          # Row vector of individuals across all health states
          v_current_state <- as.vector(a_m_trace[i - 1, , j]) # all in current state

          # Transitions remaining in current states
          # Multiply current vector of state occupancy by diagonal of state-time-dependent array (i.e., the conditional probability of remaining in given health state)
          v_same_state <- as.vector(v_current_state * diag(m_sojourn)) # individuals remaining in state next period

          # Add remain probabilities to next period, increase sojourn time by 1
          # a_M_trace already initialized for i = 1 above
          # Add % individuals remaining in subsequent period only
          a_m_trace[i, , j + 1] <- v_same_state # add remain to next period

          # Set diagonal to zero as these have been counted
          # Reset only remain (diagonal) to zero and keep state transitions
          diag(m_sojourn) <- 0 # reset remain to 0 once counted

          # Transitions to new health states
          # Create vector of states following transition according to conditional state-transition matrix excluding remain (i.e., multiply same states by zero)
          v_new_state <- as.vector(v_current_state %*% m_sojourn) # populate new states post-transition (excluding remaining)

          # Add new state %'s to markov trace
          a_m_trace[i, , 1] <- v_new_state + a_m_trace[i, , 1] # add new state %'s to array
        }
      }
      # Year 2 (2013)
      for (i in 52:103) {
        for (j in 1:(i - 1)) {
          m_sojourn <- a_tdp_4[, , j] * m_alive[, i - 1]
          v_current_state <- as.vector(a_m_trace[i - 1, , j])
          v_same_state <- as.vector(v_current_state * diag(m_sojourn))
          a_m_trace[i, , j + 1] <- v_same_state
          diag(m_sojourn) <- 0
          v_new_state <- as.vector(v_current_state %*% m_sojourn)
          a_m_trace[i, , 1] <- v_new_state + a_m_trace[i, , 1]
        }
      }
      # Year 3 (2014)
      for (i in 104:155) {
        for (j in 1:(i - 1)) {
          m_sojourn <- a_tdp_5[, , j] * m_alive[, i - 1]
          v_current_state <- as.vector(a_m_trace[i - 1, , j]) # all in current state
          v_same_state <- as.vector(v_current_state * diag(m_sojourn)) # individuals remaining in state next period
          a_m_trace[i, , j + 1] <- v_same_state # add remain to next period
          diag(m_sojourn) <- 0 # reset remain to 0 once counted
          v_new_state <- as.vector(v_current_state %*% m_sojourn) # populate new states post-transition (excluding remaining)
          a_m_trace[i, , 1] <- v_new_state + a_m_trace[i, , 1] # add new state %'s to array
        }
      }
      # Year 4 (2015)
      for (i in 156:207) {
        # Time spent in given health state
        # First month (state-time)
        for (j in 1:(i - 1)) {
          m_sojourn <- a_tdp_6[, , j] * m_alive[, i - 1] # state-time transition matrix at state-time j, re-weighted for model-time (age) varying mortality at each time point, time-varying overdose (t=3)
          v_current_state <- as.vector(a_m_trace[i - 1, , j]) # all in current state
          v_same_state <- as.vector(v_current_state * diag(m_sojourn)) # individuals remaining in state next period
          a_m_trace[i, , j + 1] <- v_same_state # add remain to next period
          diag(m_sojourn) <- 0 # reset remain to 0 once counted
          v_new_state <- as.vector(v_current_state %*% m_sojourn) # populate new states post-transition (excluding remaining)
          a_m_trace[i, , 1] <- v_new_state + a_m_trace[i, , 1] # add new state %'s to array
        }
      }
      # Year 5 (2016)
      for (i in 208:259) {
        # Time spent in given health state
        # First month (state-time)
        for (j in 1:(i - 1)) {
          m_sojourn <- a_tdp_7[, , j] * m_alive[, i - 1] # state-time transition matrix at state-time j, re-weighted for model-time (age) varying mortality at each time point, time-varying overdose (t=3)
          v_current_state <- as.vector(a_m_trace[i - 1, , j]) # all in current state
          v_same_state <- as.vector(v_current_state * diag(m_sojourn)) # individuals remaining in state next period
          a_m_trace[i, , j + 1] <- v_same_state # add remain to next period
          diag(m_sojourn) <- 0 # reset remain to 0 once counted
          v_new_state <- as.vector(v_current_state %*% m_sojourn) # populate new states post-transition (excluding remaining)
          a_m_trace[i, , 1] <- v_new_state + a_m_trace[i, , 1] # add new state %'s to array
        }
      }
      # Year 6 (2017)
      for (i in 260:311) {
        # Time spent in given health state
        # First month (state-time)
        for (j in 1:(i - 1)) {
          m_sojourn <- a_tdp_8[, , j] * m_alive[, i - 1] # state-time transition matrix at state-time j, re-weighted for model-time (age) varying mortality at each time point, time-varying overdose (t=3)
          v_current_state <- as.vector(a_m_trace[i - 1, , j]) # all in current state
          v_same_state <- as.vector(v_current_state * diag(m_sojourn)) # individuals remaining in state next period
          a_m_trace[i, , j + 1] <- v_same_state # add remain to next period
          diag(m_sojourn) <- 0 # reset remain to 0 once counted
          v_new_state <- as.vector(v_current_state %*% m_sojourn) # populate new states post-transition (excluding remaining)
          a_m_trace[i, , 1] <- v_new_state + a_m_trace[i, , 1] # add new state %'s to array
        }
      }
      # Year 7 (2018)
      for (i in 312:363) {
        # Time spent in given health state
        # First month (state-time)
        for (j in 1:(i - 1)) {
          m_sojourn <- a_tdp_9[, , j] * m_alive[, i - 1] # state-time transition matrix at state-time j, re-weighted for model-time (age) varying mortality at each time point, time-varying overdose (t=3)
          v_current_state <- as.vector(a_m_trace[i - 1, , j]) # all in current state
          v_same_state <- as.vector(v_current_state * diag(m_sojourn)) # individuals remaining in state next period
          a_m_trace[i, , j + 1] <- v_same_state # add remain to next period
          diag(m_sojourn) <- 0 # reset remain to 0 once counted
          v_new_state <- as.vector(v_current_state %*% m_sojourn) # populate new states post-transition (excluding remaining)
          a_m_trace[i, , 1] <- v_new_state + a_m_trace[i, , 1] # add new state %'s to array
        }
      }
      # # Year 8 (2019)
      for (i in 364:415) {
        # Time spent in given health state
        # First month (state-time)
        for (j in 1:(i - 1)) {
          m_sojourn <- a_tdp_10[, , j] * m_alive[, i - 1] # state-time transition matrix at state-time j, re-weighted for model-time (age) varying mortality at each time point, time-varying overdose (t=3)
          v_current_state <- as.vector(a_m_trace[i - 1, , j]) # all in current state
          v_same_state <- as.vector(v_current_state * diag(m_sojourn)) # individuals remaining in state next period
          a_m_trace[i, , j + 1] <- v_same_state # add remain to next period
          diag(m_sojourn) <- 0 # reset remain to 0 once counted
          v_new_state <- as.vector(v_current_state %*% m_sojourn) # populate new states post-transition (excluding remaining)
          a_m_trace[i, , 1] <- v_new_state + a_m_trace[i, , 1] # add new state %'s to array
        }
      }
      # # Year 9 (2020)
      for (i in 416:(n_t)) {
        # Time spent in given health state
        # First month (state-time)
        for (j in 1:(i - 1)) {
          m_sojourn <- a_tdp_11[, , j] * m_alive[, i - 1] # state-time transition matrix at state-time j, re-weighted for model-time (age) varying mortality at each time point, time-varying overdose (t=3)
          v_current_state <- as.vector(a_m_trace[i - 1, , j]) # all in current state
          v_same_state <- as.vector(v_current_state * diag(m_sojourn)) # individuals remaining in state next period
          a_m_trace[i, , j + 1] <- v_same_state # add remain to next period
          diag(m_sojourn) <- 0 # reset remain to 0 once counted
          v_new_state <- as.vector(v_current_state %*% m_sojourn) # populate new states post-transition (excluding remaining)
          a_m_trace[i, , 1] <- v_new_state + a_m_trace[i, , 1] # add new state %'s to array
        }
      }
    } else if (time_horizon == "full") { # Non-calibration, run model for 2010-2020
      # 2010
      for (i in 2:51) {
        # Time spent in given health state
        # First month (state-time)
        for (j in 1:(i - 1)) {
          m_sojourn <- a_tdp_1[, , j] * m_alive[, i - 1] # state-time transition matrix at state-time j, re-weighted for model-time (age) varying mortality at each time point, time-varying overdose (t=3)
          v_current_state <- as.vector(a_m_trace[i - 1, , j]) # all in current state
          v_same_state <- as.vector(v_current_state * diag(m_sojourn)) # individuals remaining in state next period
          a_m_trace[i, , j + 1] <- v_same_state # add remain to next period
          diag(m_sojourn) <- 0 # reset remain to 0 once counted
          v_new_state <- as.vector(v_current_state %*% m_sojourn) # populate new states post-transition (excluding remaining)
          a_m_trace[i, , 1] <- v_new_state + a_m_trace[i, , 1] # add new state %'s to array
        }
      }
      # 2011
      for (i in 52:103) {
        # Time spent in given health state
        # First month (state-time)
        for (j in 1:(i - 1)) {
          m_sojourn <- a_tdp_2[, , j] * m_alive[, i - 1] # state-time transition matrix at state-time j, re-weighted for model-time (age) varying mortality at each time point, time-varying overdose (t=3)
          v_current_state <- as.vector(a_m_trace[i - 1, , j]) # all in current state
          v_same_state <- as.vector(v_current_state * diag(m_sojourn)) # individuals remaining in state next period
          a_m_trace[i, , j + 1] <- v_same_state # add remain to next period
          diag(m_sojourn) <- 0 # reset remain to 0 once counted
          v_new_state <- as.vector(v_current_state %*% m_sojourn) # populate new states post-transition (excluding remaining)
          a_m_trace[i, , 1] <- v_new_state + a_m_trace[i, , 1] # add new state %'s to array
        }
      }
      # 2012
      for (i in 104:155) {
        # Time spent in given health state
        # First month (state-time)
        for (j in 1:(i - 1)) {
          m_sojourn <- a_tdp_3[, , j] * m_alive[, i - 1] # state-time transition matrix at state-time j, re-weighted for model-time (age) varying mortality at each time point, time-varying overdose (t=3)
          v_current_state <- as.vector(a_m_trace[i - 1, , j]) # all in current state
          v_same_state <- as.vector(v_current_state * diag(m_sojourn)) # individuals remaining in state next period
          a_m_trace[i, , j + 1] <- v_same_state # add remain to next period
          diag(m_sojourn) <- 0 # reset remain to 0 once counted
          v_new_state <- as.vector(v_current_state %*% m_sojourn) # populate new states post-transition (excluding remaining)
          a_m_trace[i, , 1] <- v_new_state + a_m_trace[i, , 1] # add new state %'s to array
        }
      }
      # 2013
      for (i in 156:207) {
        # Time spent in given health state
        # First month (state-time)
        for (j in 1:(i - 1)) {
          m_sojourn <- a_tdp_4[, , j] * m_alive[, i - 1] # state-time transition matrix at state-time j, re-weighted for model-time (age) varying mortality at each time point, time-varying overdose (t=3)
          v_current_state <- as.vector(a_m_trace[i - 1, , j]) # all in current state
          v_same_state <- as.vector(v_current_state * diag(m_sojourn)) # individuals remaining in state next period
          a_m_trace[i, , j + 1] <- v_same_state # add remain to next period
          diag(m_sojourn) <- 0 # reset remain to 0 once counted
          v_new_state <- as.vector(v_current_state %*% m_sojourn) # populate new states post-transition (excluding remaining)
          a_m_trace[i, , 1] <- v_new_state + a_m_trace[i, , 1] # add new state %'s to array
        }
      }
      # 2014
      for (i in 208:259) {
        # Time spent in given health state
        # First month (state-time)
        for (j in 1:(i - 1)) {
          m_sojourn <- a_tdp_5[, , j] * m_alive[, i - 1] # state-time transition matrix at state-time j, re-weighted for model-time (age) varying mortality at each time point, time-varying overdose (t=3)
          v_current_state <- as.vector(a_m_trace[i - 1, , j]) # all in current state
          v_same_state <- as.vector(v_current_state * diag(m_sojourn)) # individuals remaining in state next period
          a_m_trace[i, , j + 1] <- v_same_state # add remain to next period
          diag(m_sojourn) <- 0 # reset remain to 0 once counted
          v_new_state <- as.vector(v_current_state %*% m_sojourn) # populate new states post-transition (excluding remaining)
          a_m_trace[i, , 1] <- v_new_state + a_m_trace[i, , 1] # add new state %'s to array
        }
      }
      # 2015
      for (i in 260:311) {
        # Time spent in given health state
        # First month (state-time)
        for (j in 1:(i - 1)) {
          m_sojourn <- a_tdp_6[, , j] * m_alive[, i - 1] # state-time transition matrix at state-time j, re-weighted for model-time (age) varying mortality at each time point, time-varying overdose (t=3)
          v_current_state <- as.vector(a_m_trace[i - 1, , j]) # all in current state
          v_same_state <- as.vector(v_current_state * diag(m_sojourn)) # individuals remaining in state next period
          a_m_trace[i, , j + 1] <- v_same_state # add remain to next period
          diag(m_sojourn) <- 0 # reset remain to 0 once counted
          v_new_state <- as.vector(v_current_state %*% m_sojourn) # populate new states post-transition (excluding remaining)
          a_m_trace[i, , 1] <- v_new_state + a_m_trace[i, , 1] # add new state %'s to array
        }
      }
      # 2016
      for (i in 312:363) {
        # Time spent in given health state
        # First month (state-time)
        for (j in 1:(i - 1)) {
          m_sojourn <- a_tdp_7[, , j] * m_alive[, i - 1] # state-time transition matrix at state-time j, re-weighted for model-time (age) varying mortality at each time point, time-varying overdose (t=3)
          v_current_state <- as.vector(a_m_trace[i - 1, , j]) # all in current state
          v_same_state <- as.vector(v_current_state * diag(m_sojourn)) # individuals remaining in state next period
          a_m_trace[i, , j + 1] <- v_same_state # add remain to next period
          diag(m_sojourn) <- 0 # reset remain to 0 once counted
          v_new_state <- as.vector(v_current_state %*% m_sojourn) # populate new states post-transition (excluding remaining)
          a_m_trace[i, , 1] <- v_new_state + a_m_trace[i, , 1] # add new state %'s to array
        }
      }
      # 2017
      for (i in 364:415) {
        # Time spent in given health state
        # First month (state-time)
        for (j in 1:(i - 1)) {
          m_sojourn <- a_tdp_8[, , j] * m_alive[, i - 1] # state-time transition matrix at state-time j, re-weighted for model-time (age) varying mortality at each time point, time-varying overdose (t=3)
          v_current_state <- as.vector(a_m_trace[i - 1, , j]) # all in current state
          v_same_state <- as.vector(v_current_state * diag(m_sojourn)) # individuals remaining in state next period
          a_m_trace[i, , j + 1] <- v_same_state # add remain to next period
          diag(m_sojourn) <- 0 # reset remain to 0 once counted
          v_new_state <- as.vector(v_current_state %*% m_sojourn) # populate new states post-transition (excluding remaining)
          a_m_trace[i, , 1] <- v_new_state + a_m_trace[i, , 1] # add new state %'s to array
        }
      }
      # 2018
      for (i in 416:467) {
        # Time spent in given health state
        # First month (state-time)
        for (j in 1:(i - 1)) {
          m_sojourn <- a_tdp_9[, , j] * m_alive[, i - 1] # state-time transition matrix at state-time j, re-weighted for model-time (age) varying mortality at each time point, time-varying overdose (t=3)
          v_current_state <- as.vector(a_m_trace[i - 1, , j]) # all in current state
          v_same_state <- as.vector(v_current_state * diag(m_sojourn)) # individuals remaining in state next period
          a_m_trace[i, , j + 1] <- v_same_state # add remain to next period
          diag(m_sojourn) <- 0 # reset remain to 0 once counted
          v_new_state <- as.vector(v_current_state %*% m_sojourn) # populate new states post-transition (excluding remaining)
          a_m_trace[i, , 1] <- v_new_state + a_m_trace[i, , 1] # add new state %'s to array
        }
      }
      # 2019
      for (i in 468:519) {
        # Time spent in given health state
        # First month (state-time)
        for (j in 1:(i - 1)) {
          m_sojourn <- a_tdp_10[, , j] * m_alive[, i - 1] # state-time transition matrix at state-time j, re-weighted for model-time (age) varying mortality at each time point, time-varying overdose (t=3)
          v_current_state <- as.vector(a_m_trace[i - 1, , j]) # all in current state
          v_same_state <- as.vector(v_current_state * diag(m_sojourn)) # individuals remaining in state next period
          a_m_trace[i, , j + 1] <- v_same_state # add remain to next period
          diag(m_sojourn) <- 0 # reset remain to 0 once counted
          v_new_state <- as.vector(v_current_state %*% m_sojourn) # populate new states post-transition (excluding remaining)
          a_m_trace[i, , 1] <- v_new_state + a_m_trace[i, , 1] # add new state %'s to array
        }
      }
      # 2020
      for (i in 520:(n_t)) {
        # Time spent in given health state
        # First month (state-time)
        for (j in 1:(i - 1)) {
          m_sojourn <- a_tdp_11[, , j] * m_alive[, i - 1] # state-time transition matrix at state-time j, re-weighted for model-time (age) varying mortality at each time point, time-varying overdose (t=3)
          v_current_state <- as.vector(a_m_trace[i - 1, , j]) # all in current state
          v_same_state <- as.vector(v_current_state * diag(m_sojourn)) # individuals remaining in state next period
          a_m_trace[i, , j + 1] <- v_same_state # add remain to next period
          diag(m_sojourn) <- 0 # reset remain to 0 once counted
          v_new_state <- as.vector(v_current_state %*% m_sojourn) # populate new states post-transition (excluding remaining)
          a_m_trace[i, , 1] <- v_new_state + a_m_trace[i, , 1] # add new state %'s to array
        }
      }
    } else if (time_horizon == "fent_era") { # Non-calibration, run model for 2012-2020
      # 2012
      for (i in 2:51) {
        # Time spent in given health state
        # First month (state-time)
        for (j in 1:(i - 1)) {
          m_sojourn <- a_tdp_3[, , j] * m_alive[, i - 1] # state-time transition matrix at state-time j, re-weighted for model-time (age) varying mortality at each time point, time-varying overdose (t=3)
          v_current_state <- as.vector(a_m_trace[i - 1, , j]) # all in current state
          v_same_state <- as.vector(v_current_state * diag(m_sojourn)) # individuals remaining in state next period
          a_m_trace[i, , j + 1] <- v_same_state # add remain to next period
          diag(m_sojourn) <- 0 # reset remain to 0 once counted
          v_new_state <- as.vector(v_current_state %*% m_sojourn) # populate new states post-transition (excluding remaining)
          a_m_trace[i, , 1] <- v_new_state + a_m_trace[i, , 1] # add new state %'s to array
        }
      }
      # 2013
      for (i in 52:103) {
        # Time spent in given health state
        # First month (state-time)
        for (j in 1:(i - 1)) {
          m_sojourn <- a_tdp_4[, , j] * m_alive[, i - 1] # state-time transition matrix at state-time j, re-weighted for model-time (age) varying mortality at each time point, time-varying overdose (t=3)
          v_current_state <- as.vector(a_m_trace[i - 1, , j]) # all in current state
          v_same_state <- as.vector(v_current_state * diag(m_sojourn)) # individuals remaining in state next period
          a_m_trace[i, , j + 1] <- v_same_state # add remain to next period
          diag(m_sojourn) <- 0 # reset remain to 0 once counted
          v_new_state <- as.vector(v_current_state %*% m_sojourn) # populate new states post-transition (excluding remaining)
          a_m_trace[i, , 1] <- v_new_state + a_m_trace[i, , 1] # add new state %'s to array
        }
      }
      # 2014
      for (i in 104:155) {
        # Time spent in given health state
        # First month (state-time)
        for (j in 1:(i - 1)) {
          m_sojourn <- a_tdp_5[, , j] * m_alive[, i - 1] # state-time transition matrix at state-time j, re-weighted for model-time (age) varying mortality at each time point, time-varying overdose (t=3)
          v_current_state <- as.vector(a_m_trace[i - 1, , j]) # all in current state
          v_same_state <- as.vector(v_current_state * diag(m_sojourn)) # individuals remaining in state next period
          a_m_trace[i, , j + 1] <- v_same_state # add remain to next period
          diag(m_sojourn) <- 0 # reset remain to 0 once counted
          v_new_state <- as.vector(v_current_state %*% m_sojourn) # populate new states post-transition (excluding remaining)
          a_m_trace[i, , 1] <- v_new_state + a_m_trace[i, , 1] # add new state %'s to array
        }
      }
      # 2015
      for (i in 156:207) {
        # Time spent in given health state
        # First month (state-time)
        for (j in 1:(i - 1)) {
          m_sojourn <- a_tdp_6[, , j] * m_alive[, i - 1] # state-time transition matrix at state-time j, re-weighted for model-time (age) varying mortality at each time point, time-varying overdose (t=3)
          v_current_state <- as.vector(a_m_trace[i - 1, , j]) # all in current state
          v_same_state <- as.vector(v_current_state * diag(m_sojourn)) # individuals remaining in state next period
          a_m_trace[i, , j + 1] <- v_same_state # add remain to next period
          diag(m_sojourn) <- 0 # reset remain to 0 once counted
          v_new_state <- as.vector(v_current_state %*% m_sojourn) # populate new states post-transition (excluding remaining)
          a_m_trace[i, , 1] <- v_new_state + a_m_trace[i, , 1] # add new state %'s to array
        }
      }
      # 2016
      for (i in 208:259) {
        # Time spent in given health state
        # First month (state-time)
        for (j in 1:(i - 1)) {
          m_sojourn <- a_tdp_7[, , j] * m_alive[, i - 1] # state-time transition matrix at state-time j, re-weighted for model-time (age) varying mortality at each time point, time-varying overdose (t=3)
          v_current_state <- as.vector(a_m_trace[i - 1, , j]) # all in current state
          v_same_state <- as.vector(v_current_state * diag(m_sojourn)) # individuals remaining in state next period
          a_m_trace[i, , j + 1] <- v_same_state # add remain to next period
          diag(m_sojourn) <- 0 # reset remain to 0 once counted
          v_new_state <- as.vector(v_current_state %*% m_sojourn) # populate new states post-transition (excluding remaining)
          a_m_trace[i, , 1] <- v_new_state + a_m_trace[i, , 1] # add new state %'s to array
        }
      }
      # 2017
      for (i in 260:311) {
        # Time spent in given health state
        # First month (state-time)
        for (j in 1:(i - 1)) {
          m_sojourn <- a_tdp_8[, , j] * m_alive[, i - 1] # state-time transition matrix at state-time j, re-weighted for model-time (age) varying mortality at each time point, time-varying overdose (t=3)
          v_current_state <- as.vector(a_m_trace[i - 1, , j]) # all in current state
          v_same_state <- as.vector(v_current_state * diag(m_sojourn)) # individuals remaining in state next period
          a_m_trace[i, , j + 1] <- v_same_state # add remain to next period
          diag(m_sojourn) <- 0 # reset remain to 0 once counted
          v_new_state <- as.vector(v_current_state %*% m_sojourn) # populate new states post-transition (excluding remaining)
          a_m_trace[i, , 1] <- v_new_state + a_m_trace[i, , 1] # add new state %'s to array
        }
      }
      # 2018
      for (i in 312:363) {
        # Time spent in given health state
        # First month (state-time)
        for (j in 1:(i - 1)) {
          m_sojourn <- a_tdp_9[, , j] * m_alive[, i - 1] # state-time transition matrix at state-time j, re-weighted for model-time (age) varying mortality at each time point, time-varying overdose (t=3)
          v_current_state <- as.vector(a_m_trace[i - 1, , j]) # all in current state
          v_same_state <- as.vector(v_current_state * diag(m_sojourn)) # individuals remaining in state next period
          a_m_trace[i, , j + 1] <- v_same_state # add remain to next period
          diag(m_sojourn) <- 0 # reset remain to 0 once counted
          v_new_state <- as.vector(v_current_state %*% m_sojourn) # populate new states post-transition (excluding remaining)
          a_m_trace[i, , 1] <- v_new_state + a_m_trace[i, , 1] # add new state %'s to array
        }
      }
      # 2019
      for (i in 364:415) {
        # Time spent in given health state
        # First month (state-time)
        for (j in 1:(i - 1)) {
          m_sojourn <- a_tdp_10[, , j] * m_alive[, i - 1] # state-time transition matrix at state-time j, re-weighted for model-time (age) varying mortality at each time point, time-varying overdose (t=3)
          v_current_state <- as.vector(a_m_trace[i - 1, , j]) # all in current state
          v_same_state <- as.vector(v_current_state * diag(m_sojourn)) # individuals remaining in state next period
          a_m_trace[i, , j + 1] <- v_same_state # add remain to next period
          diag(m_sojourn) <- 0 # reset remain to 0 once counted
          v_new_state <- as.vector(v_current_state %*% m_sojourn) # populate new states post-transition (excluding remaining)
          a_m_trace[i, , 1] <- v_new_state + a_m_trace[i, , 1] # add new state %'s to array
        }
      }
      # 2020
      for (i in 416:(n_t)) {
        # Time spent in given health state
        # First month (state-time)
        for (j in 1:(i - 1)) {
          m_sojourn <- a_tdp_11[, , j] * m_alive[, i - 1] # state-time transition matrix at state-time j, re-weighted for model-time (age) varying mortality at each time point, time-varying overdose (t=3)
          v_current_state <- as.vector(a_m_trace[i - 1, , j]) # all in current state
          v_same_state <- as.vector(v_current_state * diag(m_sojourn)) # individuals remaining in state next period
          a_m_trace[i, , j + 1] <- v_same_state # add remain to next period
          diag(m_sojourn) <- 0 # reset remain to 0 once counted
          v_new_state <- as.vector(v_current_state %*% m_sojourn) # populate new states post-transition (excluding remaining)
          a_m_trace[i, , 1] <- v_new_state + a_m_trace[i, , 1] # add new state %'s to array
        }
      }
    }

    # Collect trace for time-periods across all model states
    m_m_trace <- array(0,
      dim = c((n_t + 1), n_states),
      dimnames = list(0:n_t, v_n_states)
    )
    for (i in 1:n_t) {
      m_m_trace[i, ] <- rowSums(a_m_trace[i, , ])
    }

    v_agg_trace_states <- c("Alive", "Non-Overdose Death", "OOT", "ODN", "ODF", "OUM", "OUB", "OUO", "BNX", "MET", "ABS") # states to aggregate
    v_agg_trace_states_cali <- c("Non-Overdose Death", "ODN", "ODF", "OU", "OAT", "ABS") # states to aggregate
    v_cohort_balance_trace_states <- c("Incident user", "Experienced user") # states to aggregate

    n_agg_trace_states <- length(v_agg_trace_states)
    n_agg_trace_states_cali <- length(v_agg_trace_states_cali)
    n_cohort_balance_trace_states <- length(v_cohort_balance_trace_states)

    m_m_agg_trace <- array(0,
      dim = c((n_t + 1), n_agg_trace_states),
      dimnames = list(0:n_t, v_agg_trace_states)
    )
    m_m_agg_trace_cali <- array(0,
      dim = c((n_t + 1), n_agg_trace_states_cali),
      dimnames = list(0:n_t, v_agg_trace_states_cali)
    )
    m_m_cohort_balance_trace <- array(0,
      dim = c((n_t + 1), n_cohort_balance_trace_states),
      dimnames = list(0:n_t, v_cohort_balance_trace_states)
    )

    for (i in 1:n_t) {
      m_m_agg_trace[i, "BNX"] <- sum(m_m_trace[i, bnx])
      m_m_agg_trace[i, "MET"] <- sum(m_m_trace[i, met])
      m_m_agg_trace[i, "OUM"] <- sum(m_m_trace[i, oum])
      m_m_agg_trace[i, "OUB"] <- sum(m_m_trace[i, oub])
      m_m_agg_trace[i, "OUO"] <- sum(m_m_trace[i, ouo])
      m_m_agg_trace[i, "ABS"] <- sum(m_m_trace[i, abs])
      m_m_agg_trace[i, "ODN"] <- sum(m_m_trace[i, odn])
      m_m_agg_trace[i, "ODF"] <- sum(m_m_trace[i, odf])
      m_m_agg_trace[i, "OOT"] <- sum(m_m_trace[i, oot])
      m_m_agg_trace[i, "Non-Overdose Death"] <- 1 - sum(m_m_trace[i, ])
      m_m_agg_trace[i, "Alive"] <- sum(m_m_trace[i, !odf])
    }

    for (i in 1:n_t) {
      m_m_agg_trace_cali[i, "OAT"] <- sum(m_m_trace[i, oat])
      m_m_agg_trace_cali[i, "OU"] <- sum(m_m_trace[i, ou])
      m_m_agg_trace_cali[i, "ABS"] <- sum(m_m_trace[i, abs])
      m_m_agg_trace_cali[i, "ODN"] <- sum(m_m_trace[i, odn])
      m_m_agg_trace_cali[i, "ODF"] <- sum(m_m_trace[i, odf])
      m_m_agg_trace_cali[i, "Non-Overdose Death"] <- 1 - sum(m_m_trace[i, ])
    }

    for (i in 1:n_t) {
      m_m_cohort_balance_trace[i, "Incident user"] <- sum(m_m_trace[i, ep1])
      m_m_cohort_balance_trace[i, "Experienced user"] <- sum(m_m_trace[i, ep2])
      # m_m_cohort_balance_trace[i, "Incident (ODF)"] <- sum(m_m_trace[i, ep1 & odf])
      # m_m_cohort_balance_trace[i, "Prevalent (ODF)"] <- sum(m_m_trace[i, ep2 & odf])
      # m_m_cohort_balance_trace[i, "Non-Overdose Death"] <- 1 - sum(m_m_trace[i, ])
      # m_m_cohort_balance_trace[i, "Alive"] <- sum(m_m_trace[i, !odf])
    }

    return(list(
      l_index_s = l_index_s,
      m_m_trace = m_m_trace,
      m_m_agg_trace = m_m_agg_trace,
      m_m_agg_trace_cali = m_m_agg_trace_cali,
      m_m_cohort_balance_trace = m_m_cohort_balance_trace
    ))
  })
}

#' Check if transition array is valid
#'
#' \code{check_transition_probability} checks if individual transition probabilities are in range \[0, 1\].
#'
#' @param a_p A transition probability array.
#' @param err_stop Logical variable to stop model run if set up as TRUE. Default = FALSE.
#' @param verbose Logical variable to indicate print out of messages.
#' Default = FALSE
#'
#' @return
#' This function stops if transition probability array is not valid and shows which entries are invalid
#' @import utils
#' @export
check_transition_probability <- function(a_p,
                                         err_stop = FALSE,
                                         verbose = FALSE) {
  m_indices_notvalid <- arrayInd(
    which(a_p < 0 | a_p > 1),
    dim(a_p)
  )

  if (dim(m_indices_notvalid)[1] != 0) {
    v_rows_notval <- rownames(a_p)[m_indices_notvalid[, 1]]
    v_cols_notval <- colnames(a_p)[m_indices_notvalid[, 2]]
    v_cycles_notval <- dimnames(a_p)[[3]][m_indices_notvalid[, 3]]

    df_notvalid <- data.frame(
      `Transition probabilities not valid:` =
        matrix(paste0(
          paste(v_rows_notval, v_cols_notval, sep = "->"),
          "; at cycle ",
          v_cycles_notval
        ), ncol = 1),
      check.names = FALSE
    )

    if (err_stop) {
      stop(
        "Not valid transition probabilities\n",
        paste(capture.output(df_notvalid), collapse = "\n")
      )
    }

    if (verbose) {
      warning(
        "Not valid transition probabilities\n",
        paste(capture.output(df_notvalid), collapse = "\n")
      )
    }
  }
}

#' Check if the sum of transition probabilities from each state are equal to one.
#'
#' \code{check_sum_of_transition_array} checks if each of the rows of the
#' transition matrices sum to one.
#'
#' @param a_p Transition probability array.
#' @param n_states Number of health states.
#' @param n_t Number of time periods.
#' @param err_stop Logical variable to stop model run if set up as TRUE. Default = FALSE.
#' @param verbose Logical variable to indicate print out of messages. Default = FALSE.
#' @return
#' The transition probability array and the cohort trace matrix.
#' @import dplyr
#' @export
check_sum_of_transition_array <- function(a_p,
                                          n_states,
                                          n_t,
                                          err_stop = FALSE,
                                          verbose = FALSE) {
  valid <- (apply(a_p, 3, function(x) sum(rowSums(x))) == n_states)
  if (!isTRUE(all_equal(as.numeric(sum(valid)), as.numeric(n_t)))) {
    if (err_stop) {
      stop("This is not a valid transition Matrix")
    }

    if (verbose) {
      warning("This is not a valid transition Matrix")
    }
  }
}
