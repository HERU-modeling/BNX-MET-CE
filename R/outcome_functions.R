#' Outcomes
#'
#' \code{outcomes} implements functions to apply costs & QALYs to Markov trace.
#'
#' @param l_params_all List with all parameters
#' @param markov_model Run Markov model with l_params_all as inputs
#' @return
#' n_qalys_3mo:
#' n_qalys_6mo:
#' n_qalys_1yr:
#' n_qalys_2yr:
#' n_qalys_5yr:
#' n_qalys_10yr:
#' n_qalys_max:
#' @export
outcomes <- function(l_params_all,
                     v_params_calib = NULL,
                     v_params_dsa = NULL,
                     checks = FALSE,
                     psa = FALSE,
                     time_horizon = time_horizon,
                     ce_est = ce_est,
                     analytic_cohort = analytic_cohort) {
  # Substitute values of calibrated parameters in base-case with calibrated values
  # Combine with DSA ranges
  dsa_calib_check <- !is.null(v_params_calib)
  dsa_check <- !is.null(v_params_dsa)

  if (dsa_calib_check) {
    l_params_updated <- update_param_list(l_params_all = l_params_all, params_updated = v_params_calib)
  } else {
    print("Calibration parameters not provided")
  }

  if (dsa_check) {
    l_params_updated <- update_param_list(l_params_all = l_params_updated, params_updated = v_params_dsa)
  } else {

  }

  # l_params_all <- update_param_list(l_params_all = l_params_all, params_updated = v_params_updated)
  # PSA update already includes draws from posterior for cali parameters
  if (psa) {
    l_params_all <- l_params_all
  } else {
    l_params_all <- update_param_list(l_params_all = l_params_all, params_updated = l_params_updated)
  }

  # Substitute values of calibrated parameters in base-case with calibrated values
  # l_params_all <- update_param_list(l_params_all = l_params_all, params_updated = v_params_calib)

  # Run model
  l_out_markov <- markov_model(l_params_all = l_params_all, err_stop = FALSE, verbose = TRUE, checks = checks, time_horizon = time_horizon, ce_est = ce_est, analytic_cohort = analytic_cohort)

  # State indices
  l_index_s <- l_out_markov$l_index_s # Load state/strata indices to assign costs/QALYs
  bnx <- l_index_s$bnx
  met <- l_index_s$met
  oum <- l_index_s$oum
  oub <- l_index_s$oub
  ouo <- l_index_s$ouo
  odn <- l_index_s$odn
  odf <- l_index_s$odf
  abs <- l_index_s$abs

  #### Cost-effectiveness analysis ####
  with(as.list(l_params_all), {
    ### Matrices of individual costs/QALYs ###
    m_m_trace <- l_out_markov$m_m_trace # Load full markov trace output
    m_m_agg_trace <- l_out_markov$m_m_agg_trace # Load aggregated trace output
    m_m_cohort_balance_trace <- l_out_markov$m_m_cohort_balance_trace

    ### Cohort alive ###
    v_alive <- l_out_markov$m_m_agg_trace[, "Alive"] # cohort alive at time i

    # Max model periods
    n_t <- nrow(m_m_agg_trace) - 1

    ### Overdose deaths ###
    v_odf <- v_odf_adj <- v_odf_inc <- l_out_markov$m_m_agg_trace[, "ODF"] # cumulative overdose deaths at time i

    ### Non-fatal overdoses ###
    v_odn <- v_odn_adj <- l_out_markov$m_m_agg_trace[, "ODN"] # cumulative non-fatal overdoses at time i

    ### All-cause deaths ###
    v_acm <- v_acm_adj <- l_out_markov$m_m_agg_trace[, "ODF"] + l_out_markov$m_m_agg_trace[, "Non-Overdose Death"] # cumulative all-cause deaths at time i

    m_qalys <- m_qalys_adj <- m_m_trace

    # QALY weights set to 1 for all health states (no-weighting)
    m_qalys[, bnx] <- m_qalys[, bnx] * u_bnx
    m_qalys[, met] <- m_qalys[, met] * u_met
    m_qalys[, oum] <- m_qalys[, oum] * u_oum
    m_qalys[, oub] <- m_qalys[, oub] * u_oub
    m_qalys[, ouo] <- m_qalys[, ouo] * u_ouo
    m_qalys[, odn] <- m_qalys[, odn] * u_odn
    m_qalys[, odf] <- m_qalys[, odf] * u_odf
    m_qalys[, abs] <- m_qalys[, abs] * u_abs

    # Adjust yearly QALY weights to monthly
    m_qalys <- m_qalys / 52

    # Calulate deaths in each week (deaths at each time step are cumulative)
    # for (i in 1:n_t) {
    #   if (i != n_t) {
    #     v_odf_inc[i + 1] <- v_odf[i + 1] - v_odf[i]
    #     v_acm_inc[i + 1] <- v_acm[i + 1] - v_acm[i]
    #   } else {
    #     v_odf_inc[i + 1] <- 0
    #     v_acm_inc[i + 1] <- 0
    #   }
    # }

    v_qalys <- rowSums(m_qalys)

    ### Sum QALYs ###
    # Unadjusted and unscaled (cumulative)
    n_qalys_2010 <- sum(v_qalys[1:52])
    n_qalys_2011 <- sum(v_qalys[1:104])
    n_qalys_2012 <- sum(v_qalys[1:156])
    n_qalys_2013 <- sum(v_qalys[1:208])
    n_qalys_2014 <- sum(v_qalys[1:260])
    n_qalys_2015 <- sum(v_qalys[1:312])
    n_qalys_2016 <- sum(v_qalys[1:364])
    n_qalys_2017 <- sum(v_qalys[1:416])
    n_qalys_2018 <- sum(v_qalys[1:468])
    n_qalys_2019 <- sum(v_qalys[1:520])
    n_qalys_2020 <- sum(v_qalys[1:n_t])

    # Unadjusted annual (not cumulative)
    n_qalys_2010_ann <- sum(v_qalys[1:52])
    n_qalys_2011_ann <- sum(v_qalys[53:104])
    n_qalys_2012_ann <- sum(v_qalys[105:156])
    n_qalys_2013_ann <- sum(v_qalys[157:208])
    n_qalys_2014_ann <- sum(v_qalys[209:260])
    n_qalys_2015_ann <- sum(v_qalys[261:312])
    n_qalys_2016_ann <- sum(v_qalys[313:364])
    n_qalys_2017_ann <- sum(v_qalys[365:416])
    n_qalys_2018_ann <- sum(v_qalys[417:468])
    n_qalys_2019_ann <- sum(v_qalys[469:520])
    n_qalys_2020_ann <- sum(v_qalys[521:n_t])

    # Scaled annually
    n_qalys_2010_ann_scaled <- n_qalys_2010_ann * n_pop_scaling_2010
    n_qalys_2011_ann_scaled <- n_qalys_2011_ann * n_pop_scaling_2011
    n_qalys_2012_ann_scaled <- n_qalys_2012_ann * n_pop_scaling_2012
    n_qalys_2013_ann_scaled <- n_qalys_2013_ann * n_pop_scaling_2013
    n_qalys_2014_ann_scaled <- n_qalys_2014_ann * n_pop_scaling_2014
    n_qalys_2015_ann_scaled <- n_qalys_2015_ann * n_pop_scaling_2015
    n_qalys_2016_ann_scaled <- n_qalys_2016_ann * n_pop_scaling_2016
    n_qalys_2017_ann_scaled <- n_qalys_2017_ann * n_pop_scaling_2017
    n_qalys_2018_ann_scaled <- n_qalys_2018_ann * n_pop_scaling_2018
    n_qalys_2019_ann_scaled <- n_qalys_2019_ann * n_pop_scaling_2019
    n_qalys_2020_ann_scaled <- n_qalys_2020_ann * n_pop_scaling_2020

    ### Non-fatal overdose
    # Unadjusted
    # Non-fatal overdoses need to be summed across time points to generate cumulative estimates
    # Yearly non-fatal overdose
    n_odn_2010 <- sum(v_odn[1:52])
    n_odn_2011 <- sum(v_odn[1:104])
    n_odn_2012 <- sum(v_odn[1:156])
    n_odn_2013 <- sum(v_odn[1:208])
    n_odn_2014 <- sum(v_odn[1:260])
    n_odn_2015 <- sum(v_odn[1:312])
    n_odn_2016 <- sum(v_odn[1:364])
    n_odn_2017 <- sum(v_odn[1:416])
    n_odn_2018 <- sum(v_odn[1:468])
    n_odn_2019 <- sum(v_odn[1:520])
    n_odn_2020 <- sum(v_odn[1:n_t])

    # Unadjusted annual (not cumulative)
    n_odn_2010_ann <- sum(v_odn[1:52])
    n_odn_2011_ann <- sum(v_odn[53:104])
    n_odn_2012_ann <- sum(v_odn[105:156])
    n_odn_2013_ann <- sum(v_odn[157:208])
    n_odn_2014_ann <- sum(v_odn[209:260])
    n_odn_2015_ann <- sum(v_odn[261:312])
    n_odn_2016_ann <- sum(v_odn[313:364])
    n_odn_2017_ann <- sum(v_odn[365:416])
    n_odn_2018_ann <- sum(v_odn[417:468])
    n_odn_2019_ann <- sum(v_odn[469:520])
    n_odn_2020_ann <- sum(v_odn[521:n_t])

    # Scaled annually
    n_odn_2010_ann_scaled <- n_odn_2010_ann * n_pop_scaling_2010
    n_odn_2011_ann_scaled <- n_odn_2011_ann * n_pop_scaling_2011
    n_odn_2012_ann_scaled <- n_odn_2012_ann * n_pop_scaling_2012
    n_odn_2013_ann_scaled <- n_odn_2013_ann * n_pop_scaling_2013
    n_odn_2014_ann_scaled <- n_odn_2014_ann * n_pop_scaling_2014
    n_odn_2015_ann_scaled <- n_odn_2015_ann * n_pop_scaling_2015
    n_odn_2016_ann_scaled <- n_odn_2016_ann * n_pop_scaling_2016
    n_odn_2017_ann_scaled <- n_odn_2017_ann * n_pop_scaling_2017
    n_odn_2018_ann_scaled <- n_odn_2018_ann * n_pop_scaling_2018
    n_odn_2019_ann_scaled <- n_odn_2019_ann * n_pop_scaling_2019
    n_odn_2020_ann_scaled <- n_odn_2020_ann * n_pop_scaling_2020

    ### Overdose deaths ###
    # Unadjusted
    # Yearly fatal overdoses
    n_odf_2010 <- v_odf[52]
    n_odf_2011 <- v_odf[104]
    n_odf_2012 <- v_odf[156]
    n_odf_2013 <- v_odf[208]
    n_odf_2014 <- v_odf[260]
    n_odf_2015 <- v_odf[312]
    n_odf_2016 <- v_odf[364]
    n_odf_2017 <- v_odf[416]
    n_odf_2018 <- v_odf[468]
    n_odf_2019 <- v_odf[520]
    n_odf_2020 <- v_odf[n_t]

    # Unadjusted annual (not cumulative)
    n_odf_2010_ann <- v_odf[52]
    n_odf_2011_ann <- (v_odf[104] - v_odf[52])
    n_odf_2012_ann <- (v_odf[156] - v_odf[104])
    n_odf_2013_ann <- (v_odf[208] - v_odf[156])
    n_odf_2014_ann <- (v_odf[260] - v_odf[208])
    n_odf_2015_ann <- (v_odf[312] - v_odf[260])
    n_odf_2016_ann <- (v_odf[364] - v_odf[312])
    n_odf_2017_ann <- (v_odf[416] - v_odf[364])
    n_odf_2018_ann <- (v_odf[468] - v_odf[416])
    n_odf_2019_ann <- (v_odf[520] - v_odf[468])
    n_odf_2020_ann <- (v_odf[n_t] - v_odf[520])

    # Scaled annually
    n_odf_2010_ann_scaled <- n_odf_2010_ann * n_pop_scaling_2010
    n_odf_2011_ann_scaled <- n_odf_2011_ann * n_pop_scaling_2011
    n_odf_2012_ann_scaled <- n_odf_2012_ann * n_pop_scaling_2012
    n_odf_2013_ann_scaled <- n_odf_2013_ann * n_pop_scaling_2013
    n_odf_2014_ann_scaled <- n_odf_2014_ann * n_pop_scaling_2014
    n_odf_2015_ann_scaled <- n_odf_2015_ann * n_pop_scaling_2015
    n_odf_2016_ann_scaled <- n_odf_2016_ann * n_pop_scaling_2016
    n_odf_2017_ann_scaled <- n_odf_2017_ann * n_pop_scaling_2017
    n_odf_2018_ann_scaled <- n_odf_2018_ann * n_pop_scaling_2018
    n_odf_2019_ann_scaled <- n_odf_2019_ann * n_pop_scaling_2019
    n_odf_2020_ann_scaled <- n_odf_2020_ann * n_pop_scaling_2020

    # m_odf_comp <- cbind(v_odf, v_odf_adj, v_odf_inc)

    ### All-cause deaths ###
    # Unadjusted
    # Yearly all-cause deaths
    n_acm_2010 <- v_acm[52]
    n_acm_2011 <- v_acm[104]
    n_acm_2012 <- v_acm[156]
    n_acm_2013 <- v_acm[208]
    n_acm_2014 <- v_acm[260]
    n_acm_2015 <- v_acm[312]
    n_acm_2016 <- v_acm[364]
    n_acm_2017 <- v_acm[416]
    n_acm_2018 <- v_acm[468]
    n_acm_2019 <- v_acm[520]
    n_acm_2020 <- v_acm[n_t]

    # Unadjusted annual (not cumulative)
    n_acm_2010_ann <- v_acm[52]
    n_acm_2011_ann <- (v_acm[104] - v_acm[52])
    n_acm_2012_ann <- (v_acm[156] - v_acm[104])
    n_acm_2013_ann <- (v_acm[208] - v_acm[156])
    n_acm_2014_ann <- (v_acm[260] - v_acm[208])
    n_acm_2015_ann <- (v_acm[312] - v_acm[260])
    n_acm_2016_ann <- (v_acm[364] - v_acm[312])
    n_acm_2017_ann <- (v_acm[416] - v_acm[364])
    n_acm_2018_ann <- (v_acm[468] - v_acm[416])
    n_acm_2019_ann <- (v_acm[520] - v_acm[468])
    n_acm_2020_ann <- (v_acm[n_t] - v_acm[520])

    # Scaled annually
    n_acm_2010_ann_scaled <- n_acm_2010_ann * n_pop_scaling_2010
    n_acm_2011_ann_scaled <- n_acm_2011_ann * n_pop_scaling_2011
    n_acm_2012_ann_scaled <- n_acm_2012_ann * n_pop_scaling_2012
    n_acm_2013_ann_scaled <- n_acm_2013_ann * n_pop_scaling_2013
    n_acm_2014_ann_scaled <- n_acm_2014_ann * n_pop_scaling_2014
    n_acm_2015_ann_scaled <- n_acm_2015_ann * n_pop_scaling_2015
    n_acm_2016_ann_scaled <- n_acm_2016_ann * n_pop_scaling_2016
    n_acm_2017_ann_scaled <- n_acm_2017_ann * n_pop_scaling_2017
    n_acm_2018_ann_scaled <- n_acm_2018_ann * n_pop_scaling_2018
    n_acm_2019_ann_scaled <- n_acm_2019_ann * n_pop_scaling_2019
    n_acm_2020_ann_scaled <- n_acm_2020_ann * n_pop_scaling_2020

    ## Combined outcomes ##
    df_outcomes <- data.frame(
      # Cumulative life-years
      n_qalys_2010, n_qalys_2011, n_qalys_2012, n_qalys_2013, n_qalys_2014, n_qalys_2015,
      n_qalys_2016, n_qalys_2017, n_qalys_2018, n_qalys_2019, n_qalys_2020,
      # Annual life years (unscaled)
      n_qalys_2010_ann, n_qalys_2011_ann, n_qalys_2012_ann, n_qalys_2013_ann, n_qalys_2014_ann, n_qalys_2015_ann,
      n_qalys_2016_ann, n_qalys_2017_ann, n_qalys_2018_ann, n_qalys_2019_ann, n_qalys_2020_ann,
      # Annual life years (scaled)
      n_qalys_2010_ann_scaled, n_qalys_2011_ann_scaled, n_qalys_2012_ann_scaled, n_qalys_2013_ann_scaled, n_qalys_2014_ann_scaled, n_qalys_2015_ann_scaled,
      n_qalys_2016_ann_scaled, n_qalys_2017_ann_scaled, n_qalys_2018_ann_scaled, n_qalys_2019_ann_scaled, n_qalys_2020_ann_scaled,

      # Cumulative fatal overdoses
      n_odf_2010, n_odf_2011, n_odf_2012, n_odf_2013, n_odf_2014, n_odf_2015, n_odf_2016, n_odf_2017, n_odf_2018, n_odf_2019, n_odf_2020,
      # Annual fatal overdoses (unscaled)
      n_odf_2010_ann, n_odf_2011_ann, n_odf_2012_ann, n_odf_2013_ann, n_odf_2014_ann, n_odf_2015_ann, n_odf_2016_ann, n_odf_2017_ann, n_odf_2018_ann, n_odf_2019_ann, n_odf_2020_ann,
      # Annual fatal overdoses (scaled)
      n_odf_2010_ann_scaled, n_odf_2011_ann_scaled, n_odf_2012_ann_scaled, n_odf_2013_ann_scaled, n_odf_2014_ann_scaled, n_odf_2015_ann_scaled,
      n_odf_2016_ann_scaled, n_odf_2017_ann_scaled, n_odf_2018_ann_scaled, n_odf_2019_ann_scaled, n_odf_2020_ann_scaled,

      # Cumulative non-fatal overdoses
      n_odn_2010, n_odn_2011, n_odn_2012, n_odn_2013, n_odn_2014, n_odn_2015, n_odn_2016, n_odn_2017, n_odn_2018, n_odn_2019, n_odn_2020,
      # Annual fatal overdoses (unscaled)
      n_odn_2010_ann, n_odn_2011_ann, n_odn_2012_ann, n_odn_2013_ann, n_odn_2014_ann, n_odn_2015_ann, n_odn_2016_ann, n_odn_2017_ann, n_odn_2018_ann, n_odn_2019_ann, n_odn_2020_ann,
      # Annual fatal overdoses (scaled)
      n_odn_2010_ann_scaled, n_odn_2011_ann_scaled, n_odn_2012_ann_scaled, n_odn_2013_ann_scaled, n_odn_2014_ann_scaled, n_odn_2015_ann_scaled,
      n_odn_2016_ann_scaled, n_odn_2017_ann_scaled, n_odn_2018_ann_scaled, n_odn_2019_ann_scaled, n_odn_2020_ann_scaled,

      # Cumulative all-cause deaths
      n_acm_2010, n_acm_2011, n_acm_2012, n_acm_2013, n_acm_2014, n_acm_2015, n_acm_2016, n_acm_2017, n_acm_2018, n_acm_2019, n_acm_2020,
      # Annual all-cause deaths (unscaled)
      n_acm_2010_ann, n_acm_2011_ann, n_acm_2012_ann, n_acm_2013_ann, n_acm_2014_ann, n_acm_2015_ann, n_acm_2016_ann, n_acm_2017_ann, n_acm_2018_ann, n_acm_2019_ann, n_acm_2020_ann,
      # Annual all-cause deaths (scaled)
      n_acm_2010_ann_scaled, n_acm_2011_ann_scaled, n_acm_2012_ann_scaled, n_acm_2013_ann_scaled, n_acm_2014_ann_scaled, n_acm_2015_ann_scaled,
      n_acm_2016_ann_scaled, n_acm_2017_ann_scaled, n_acm_2018_ann_scaled, n_acm_2019_ann_scaled, n_acm_2020_ann_scaled
    )

    return(list(
      m_m_trace = m_m_trace,
      m_m_agg_trace = m_m_agg_trace,
      m_m_cohort_balance_trace = m_m_cohort_balance_trace,
      v_alive = v_alive,
      # Life-years
      n_qalys_2010 = n_qalys_2010,
      n_qalys_2011 = n_qalys_2011,
      n_qalys_2012 = n_qalys_2012,
      n_qalys_2013 = n_qalys_2013,
      n_qalys_2014 = n_qalys_2014,
      n_qalys_2015 = n_qalys_2015,
      n_qalys_2016 = n_qalys_2016,
      n_qalys_2017 = n_qalys_2017,
      n_qalys_2018 = n_qalys_2018,
      n_qalys_2019 = n_qalys_2019,
      n_qalys_2020 = n_qalys_2020,
      n_qalys_2010_ann = n_qalys_2010_ann,
      n_qalys_2011_ann = n_qalys_2011_ann,
      n_qalys_2012_ann = n_qalys_2012_ann,
      n_qalys_2013_ann = n_qalys_2013_ann,
      n_qalys_2014_ann = n_qalys_2014_ann,
      n_qalys_2015_ann = n_qalys_2015_ann,
      n_qalys_2016_ann = n_qalys_2016_ann,
      n_qalys_2017_ann = n_qalys_2017_ann,
      n_qalys_2018_ann = n_qalys_2018_ann,
      n_qalys_2019_ann = n_qalys_2019_ann,
      n_qalys_2020_ann = n_qalys_2020_ann,
      n_qalys_2010_ann_scaled = n_qalys_2010_ann_scaled,
      n_qalys_2011_ann_scaled = n_qalys_2011_ann_scaled,
      n_qalys_2012_ann_scaled = n_qalys_2012_ann_scaled,
      n_qalys_2013_ann_scaled = n_qalys_2013_ann_scaled,
      n_qalys_2014_ann_scaled = n_qalys_2014_ann_scaled,
      n_qalys_2015_ann_scaled = n_qalys_2015_ann_scaled,
      n_qalys_2016_ann_scaled = n_qalys_2016_ann_scaled,
      n_qalys_2017_ann_scaled = n_qalys_2017_ann_scaled,
      n_qalys_2018_ann_scaled = n_qalys_2018_ann_scaled,
      n_qalys_2019_ann_scaled = n_qalys_2019_ann_scaled,
      n_qalys_2020_ann_scaled = n_qalys_2020_ann_scaled,
      # Fatal overdoses
      n_odf_2010 = n_odf_2010,
      n_odf_2011 = n_odf_2011,
      n_odf_2012 = n_odf_2012,
      n_odf_2013 = n_odf_2013,
      n_odf_2014 = n_odf_2014,
      n_odf_2015 = n_odf_2015,
      n_odf_2016 = n_odf_2016,
      n_odf_2017 = n_odf_2017,
      n_odf_2018 = n_odf_2018,
      n_odf_2019 = n_odf_2019,
      n_odf_2020 = n_odf_2020,
      n_odf_2010_ann = n_odf_2010_ann,
      n_odf_2011_ann = n_odf_2011_ann,
      n_odf_2012_ann = n_odf_2012_ann,
      n_odf_2013_ann = n_odf_2013_ann,
      n_odf_2014_ann = n_odf_2014_ann,
      n_odf_2015_ann = n_odf_2015_ann,
      n_odf_2016_ann = n_odf_2016_ann,
      n_odf_2017_ann = n_odf_2017_ann,
      n_odf_2018_ann = n_odf_2018_ann,
      n_odf_2019_ann = n_odf_2019_ann,
      n_odf_2020_ann = n_odf_2020_ann,
      n_odf_2010_ann_scaled = n_odf_2010_ann_scaled,
      n_odf_2011_ann_scaled = n_odf_2011_ann_scaled,
      n_odf_2012_ann_scaled = n_odf_2012_ann_scaled,
      n_odf_2013_ann_scaled = n_odf_2013_ann_scaled,
      n_odf_2014_ann_scaled = n_odf_2014_ann_scaled,
      n_odf_2015_ann_scaled = n_odf_2015_ann_scaled,
      n_odf_2016_ann_scaled = n_odf_2016_ann_scaled,
      n_odf_2017_ann_scaled = n_odf_2017_ann_scaled,
      n_odf_2018_ann_scaled = n_odf_2018_ann_scaled,
      n_odf_2019_ann_scaled = n_odf_2019_ann_scaled,
      n_odf_2020_ann_scaled = n_odf_2020_ann_scaled,
      # Non-fatal overdoses
      n_odn_2010 = n_odn_2010,
      n_odn_2011 = n_odn_2011,
      n_odn_2012 = n_odn_2012,
      n_odn_2013 = n_odn_2013,
      n_odn_2014 = n_odn_2014,
      n_odn_2015 = n_odn_2015,
      n_odn_2016 = n_odn_2016,
      n_odn_2017 = n_odn_2017,
      n_odn_2018 = n_odn_2018,
      n_odn_2019 = n_odn_2019,
      n_odn_2020 = n_odn_2020,
      n_odn_2010_ann = n_odn_2010_ann,
      n_odn_2011_ann = n_odn_2011_ann,
      n_odn_2012_ann = n_odn_2012_ann,
      n_odn_2013_ann = n_odn_2013_ann,
      n_odn_2014_ann = n_odn_2014_ann,
      n_odn_2015_ann = n_odn_2015_ann,
      n_odn_2016_ann = n_odn_2016_ann,
      n_odn_2017_ann = n_odn_2017_ann,
      n_odn_2018_ann = n_odn_2018_ann,
      n_odn_2019_ann = n_odn_2019_ann,
      n_odn_2020_ann = n_odn_2020_ann,
      n_odn_2010_ann_scaled = n_odn_2010_ann_scaled,
      n_odn_2011_ann_scaled = n_odn_2011_ann_scaled,
      n_odn_2012_ann_scaled = n_odn_2012_ann_scaled,
      n_odn_2013_ann_scaled = n_odn_2013_ann_scaled,
      n_odn_2014_ann_scaled = n_odn_2014_ann_scaled,
      n_odn_2015_ann_scaled = n_odn_2015_ann_scaled,
      n_odn_2016_ann_scaled = n_odn_2016_ann_scaled,
      n_odn_2017_ann_scaled = n_odn_2017_ann_scaled,
      n_odn_2018_ann_scaled = n_odn_2018_ann_scaled,
      n_odn_2019_ann_scaled = n_odn_2019_ann_scaled,
      n_odn_2020_ann_scaled = n_odn_2020_ann_scaled,
      # All-cause deaths
      n_acm_2010 = n_acm_2010,
      n_acm_2011 = n_acm_2011,
      n_acm_2012 = n_acm_2012,
      n_acm_2013 = n_acm_2013,
      n_acm_2014 = n_acm_2014,
      n_acm_2015 = n_acm_2015,
      n_acm_2016 = n_acm_2016,
      n_acm_2017 = n_acm_2017,
      n_acm_2018 = n_acm_2018,
      n_acm_2019 = n_acm_2019,
      n_acm_2020 = n_acm_2020,
      n_acm_2010_ann = n_acm_2010_ann,
      n_acm_2011_ann = n_acm_2011_ann,
      n_acm_2012_ann = n_acm_2012_ann,
      n_acm_2013_ann = n_acm_2013_ann,
      n_acm_2014_ann = n_acm_2014_ann,
      n_acm_2015_ann = n_acm_2015_ann,
      n_acm_2016_ann = n_acm_2016_ann,
      n_acm_2017_ann = n_acm_2017_ann,
      n_acm_2018_ann = n_acm_2018_ann,
      n_acm_2019_ann = n_acm_2019_ann,
      n_acm_2020_ann = n_acm_2020_ann,
      n_acm_2010_ann_scaled = n_acm_2010_ann_scaled,
      n_acm_2011_ann_scaled = n_acm_2011_ann_scaled,
      n_acm_2012_ann_scaled = n_acm_2012_ann_scaled,
      n_acm_2013_ann_scaled = n_acm_2013_ann_scaled,
      n_acm_2014_ann_scaled = n_acm_2014_ann_scaled,
      n_acm_2015_ann_scaled = n_acm_2015_ann_scaled,
      n_acm_2016_ann_scaled = n_acm_2016_ann_scaled,
      n_acm_2017_ann_scaled = n_acm_2017_ann_scaled,
      n_acm_2018_ann_scaled = n_acm_2018_ann_scaled,
      n_acm_2019_ann_scaled = n_acm_2019_ann_scaled,
      n_acm_2020_ann_scaled = n_acm_2020_ann_scaled,
      df_outcomes = df_outcomes
    ))
  })
}

#' inc_outcomes
#'
#' \code{inc_outcomes} calculate incremental health differences between scenarios
#'
#' @param outcomes_comp List with outputs from comparator
#' @param outcomes_int List with outputs from intervention
#' @return
#' @export
inc_outcomes <- function(outcomes_comp, outcomes_int) {
  # Max model periods
  n_t <- nrow(outcomes_comp$m_m_agg_trace) - 1

  # Calculate average proportion alive for adjustment
  n_alive_int_mean_2010 <- mean(outcomes_int$v_alive[1:52])
  n_alive_comp_mean_2010 <- mean(outcomes_comp$v_alive[1:52])
  n_alive_int_mean_2012 <- mean(outcomes_int$v_alive[1:156])
  n_alive_comp_mean_2012 <- mean(outcomes_comp$v_alive[1:156])
  n_alive_int_mean_2014 <- mean(outcomes_int$v_alive[1:260])
  n_alive_comp_mean_2014 <- mean(outcomes_comp$v_alive[1:260])
  n_alive_int_mean_2016 <- mean(outcomes_int$v_alive[1:364])
  n_alive_comp_mean_2016 <- mean(outcomes_comp$v_alive[1:364])
  n_alive_int_mean_2018 <- mean(outcomes_int$v_alive[1:468])
  n_alive_comp_mean_2018 <- mean(outcomes_comp$v_alive[1:468])
  n_alive_int_mean_2020 <- mean(outcomes_int$v_alive[1:n_t])
  n_alive_comp_mean_2020 <- mean(outcomes_comp$v_alive[1:n_t])

  #################################
  ### Total outcomes (unscaled) ###
  #################################
  # 2010
  # QALYs
  n_qalys_2010_int <- outcomes_int$n_qalys_2010_ann
  n_qalys_2010_comp <- outcomes_comp$n_qalys_2010_ann
  # Fatal overdoses
  n_odf_2010_int <- outcomes_int$n_odf_2010_ann
  n_odf_2010_comp <- outcomes_comp$n_odf_2010_ann
  # Non-fatal overdoses
  n_odn_2010_int <- outcomes_int$n_odn_2010_ann
  n_odn_2010_comp <- outcomes_comp$n_odn_2010_ann
  # All-cause deaths
  n_acm_2010_int <- outcomes_int$n_acm_2010_ann
  n_acm_2010_comp <- outcomes_comp$n_acm_2010_ann

  # 2012
  # QALYs
  n_qalys_2012_int <- sum(outcomes_int$n_qalys_2010_ann, outcomes_int$n_qalys_2011_ann, outcomes_int$n_qalys_2012_ann)
  n_qalys_2012_comp <- sum(outcomes_comp$n_qalys_2010_ann, outcomes_comp$n_qalys_2011_ann, outcomes_comp$n_qalys_2012_ann)
  # Fatal overdoses
  n_odf_2012_int <- sum(outcomes_int$n_odf_2010_ann, outcomes_int$n_odf_2011_ann, outcomes_int$n_odf_2012_ann)
  n_odf_2012_comp <- sum(outcomes_comp$n_odf_2010_ann, outcomes_comp$n_odf_2011_ann, outcomes_comp$n_odf_2012_ann)
  # Non-fatal overdoses
  n_odn_2012_int <- sum(outcomes_int$n_odn_2010_ann, outcomes_int$n_odn_2011_ann, outcomes_int$n_odn_2012_ann)
  n_odn_2012_comp <- sum(outcomes_comp$n_odn_2010_ann, outcomes_comp$n_odn_2011_ann, outcomes_comp$n_odn_2012_ann)
  # All-cause deaths
  n_acm_2012_int <- sum(outcomes_int$n_acm_2010_ann, outcomes_int$n_acm_2011_ann, outcomes_int$n_acm_2012_ann)
  n_acm_2012_comp <- sum(outcomes_comp$n_acm_2010_ann, outcomes_comp$n_acm_2011_ann, outcomes_comp$n_acm_2012_ann)

  # 2014
  # QALYs
  n_qalys_2014_int <- sum(outcomes_int$n_qalys_2010_ann, outcomes_int$n_qalys_2011_ann, outcomes_int$n_qalys_2012_ann, outcomes_int$n_qalys_2013_ann, outcomes_int$n_qalys_2014_ann)
  n_qalys_2014_comp <- sum(outcomes_comp$n_qalys_2010_ann, outcomes_comp$n_qalys_2011_ann, outcomes_comp$n_qalys_2012_ann, outcomes_comp$n_qalys_2013_ann, outcomes_comp$n_qalys_2014_ann)
  # Fatal overdoses
  n_odf_2014_int <- sum(outcomes_int$n_odf_2010_ann, outcomes_int$n_odf_2011_ann, outcomes_int$n_odf_2012_ann, outcomes_int$n_odf_2013_ann, outcomes_int$n_odf_2014_ann)
  n_odf_2014_comp <- sum(outcomes_comp$n_odf_2010_ann, outcomes_comp$n_odf_2011_ann, outcomes_comp$n_odf_2012_ann, outcomes_comp$n_odf_2013_ann, outcomes_comp$n_odf_2014_ann)
  # Non-fatal overdoses
  n_odn_2014_int <- sum(outcomes_int$n_odn_2010_ann, outcomes_int$n_odn_2011_ann, outcomes_int$n_odn_2012_ann, outcomes_int$n_odn_2013_ann, outcomes_int$n_odn_2014_ann)
  n_odn_2014_comp <- sum(outcomes_comp$n_odn_2010_ann, outcomes_comp$n_odn_2011_ann, outcomes_comp$n_odn_2012_ann, outcomes_comp$n_odn_2013_ann, outcomes_comp$n_odn_2014_ann)
  # All-cause deaths
  n_acm_2014_int <- sum(outcomes_int$n_acm_2010_ann, outcomes_int$n_acm_2011_ann, outcomes_int$n_acm_2012_ann, outcomes_int$n_acm_2013_ann, outcomes_int$n_acm_2014_ann)
  n_acm_2014_comp <- sum(outcomes_comp$n_acm_2010_ann, outcomes_comp$n_acm_2011_ann, outcomes_comp$n_acm_2012_ann, outcomes_comp$n_acm_2013_ann, outcomes_comp$n_acm_2014_ann)

  # 2016
  # QALYs
  n_qalys_2016_int <- sum(
    outcomes_int$n_qalys_2010_ann, outcomes_int$n_qalys_2011_ann, outcomes_int$n_qalys_2012_ann,
    outcomes_int$n_qalys_2013_ann, outcomes_int$n_qalys_2014_ann, outcomes_int$n_qalys_2015_ann, outcomes_int$n_qalys_2016_ann
  )
  n_qalys_2016_comp <- sum(
    outcomes_comp$n_qalys_2010_ann, outcomes_comp$n_qalys_2011_ann, outcomes_comp$n_qalys_2012_ann,
    outcomes_comp$n_qalys_2013_ann, outcomes_comp$n_qalys_2014_ann, outcomes_comp$n_qalys_2015_ann, outcomes_comp$n_qalys_2016_ann
  )
  # Fatal overdoses
  n_odf_2016_int <- sum(
    outcomes_int$n_odf_2010_ann, outcomes_int$n_odf_2011_ann, outcomes_int$n_odf_2012_ann,
    outcomes_int$n_odf_2013_ann, outcomes_int$n_odf_2014_ann, outcomes_int$n_odf_2015_ann, outcomes_int$n_odf_2016_ann
  )
  n_odf_2016_comp <- sum(
    outcomes_comp$n_odf_2010_ann, outcomes_comp$n_odf_2011_ann, outcomes_comp$n_odf_2012_ann,
    outcomes_comp$n_odf_2013_ann, outcomes_comp$n_odf_2014_ann, outcomes_comp$n_odf_2015_ann, outcomes_comp$n_odf_2016_ann
  )
  # Non-fatal overdoses
  n_odn_2016_int <- sum(
    outcomes_int$n_odn_2010_ann, outcomes_int$n_odn_2011_ann, outcomes_int$n_odn_2012_ann,
    outcomes_int$n_odn_2013_ann, outcomes_int$n_odn_2014_ann, outcomes_int$n_odn_2015_ann, outcomes_int$n_odn_2016_ann
  )
  n_odn_2016_comp <- sum(
    outcomes_comp$n_odn_2010_ann, outcomes_comp$n_odn_2011_ann, outcomes_comp$n_odn_2012_ann,
    outcomes_comp$n_odn_2013_ann, outcomes_comp$n_odn_2014_ann, outcomes_comp$n_odn_2015_ann, outcomes_comp$n_odn_2016_ann
  )
  # All-cause deaths
  n_acm_2016_int <- sum(
    outcomes_int$n_acm_2010_ann, outcomes_int$n_acm_2011_ann, outcomes_int$n_acm_2012_ann,
    outcomes_int$n_acm_2013_ann, outcomes_int$n_acm_2014_ann, outcomes_int$n_acm_2015_ann, outcomes_int$n_acm_2016_ann
  )
  n_acm_2016_comp <- sum(
    outcomes_comp$n_acm_2010_ann, outcomes_comp$n_acm_2011_ann, outcomes_comp$n_acm_2012_ann,
    outcomes_comp$n_acm_2013_ann, outcomes_comp$n_acm_2014_ann, outcomes_comp$n_acm_2015_ann, outcomes_comp$n_acm_2016_ann
  )

  # 2018
  # QALYs
  n_qalys_2018_int <- sum(
    outcomes_int$n_qalys_2010_ann, outcomes_int$n_qalys_2011_ann, outcomes_int$n_qalys_2012_ann, outcomes_int$n_qalys_2013_ann,
    outcomes_int$n_qalys_2014_ann, outcomes_int$n_qalys_2015_ann, outcomes_int$n_qalys_2016_ann, outcomes_int$n_qalys_2017_ann, outcomes_int$n_qalys_2018_ann
  )
  n_qalys_2018_comp <- sum(
    outcomes_comp$n_qalys_2010_ann, outcomes_comp$n_qalys_2011_ann, outcomes_comp$n_qalys_2012_ann, outcomes_comp$n_qalys_2013_ann,
    outcomes_comp$n_qalys_2014_ann, outcomes_comp$n_qalys_2015_ann, outcomes_comp$n_qalys_2016_ann, outcomes_comp$n_qalys_2017_ann, outcomes_comp$n_qalys_2018_ann
  )
  # Fatal overdoses
  n_odf_2018_int <- sum(
    outcomes_int$n_odf_2010_ann, outcomes_int$n_odf_2011_ann, outcomes_int$n_odf_2012_ann, outcomes_int$n_odf_2013_ann,
    outcomes_int$n_odf_2014_ann, outcomes_int$n_odf_2015_ann, outcomes_int$n_odf_2016_ann, outcomes_int$n_odf_2017_ann, outcomes_int$n_odf_2018_ann
  )
  n_odf_2018_comp <- sum(
    outcomes_comp$n_odf_2010_ann, outcomes_comp$n_odf_2011_ann, outcomes_comp$n_odf_2012_ann, outcomes_comp$n_odf_2013_ann,
    outcomes_comp$n_odf_2014_ann, outcomes_comp$n_odf_2015_ann, outcomes_comp$n_odf_2016_ann, outcomes_comp$n_odf_2017_ann, outcomes_comp$n_odf_2018_ann
  )
  # Non-fatal overdoses
  n_odn_2018_int <- sum(
    outcomes_int$n_odn_2010_ann, outcomes_int$n_odn_2011_ann, outcomes_int$n_odn_2012_ann, outcomes_int$n_odn_2013_ann,
    outcomes_int$n_odn_2014_ann, outcomes_int$n_odn_2015_ann, outcomes_int$n_odn_2016_ann, outcomes_int$n_odn_2017_ann, outcomes_int$n_odn_2018_ann
  )
  n_odn_2018_comp <- sum(
    outcomes_comp$n_odn_2010_ann, outcomes_comp$n_odn_2011_ann, outcomes_comp$n_odn_2012_ann, outcomes_comp$n_odn_2013_ann,
    outcomes_comp$n_odn_2014_ann, outcomes_comp$n_odn_2015_ann, outcomes_comp$n_odn_2016_ann, outcomes_comp$n_odn_2017_ann, outcomes_comp$n_odn_2018_ann
  )
  # All-cause deaths
  n_acm_2018_int <- sum(
    outcomes_int$n_acm_2010_ann, outcomes_int$n_acm_2011_ann, outcomes_int$n_acm_2012_ann, outcomes_int$n_acm_2013_ann,
    outcomes_int$n_acm_2014_ann, outcomes_int$n_acm_2015_ann, outcomes_int$n_acm_2016_ann, outcomes_int$n_acm_2017_ann, outcomes_int$n_acm_2018_ann
  )
  n_acm_2018_comp <- sum(
    outcomes_comp$n_acm_2010_ann, outcomes_comp$n_acm_2011_ann, outcomes_comp$n_acm_2012_ann, outcomes_comp$n_acm_2013_ann,
    outcomes_comp$n_acm_2014_ann, outcomes_comp$n_acm_2015_ann, outcomes_comp$n_acm_2016_ann, outcomes_comp$n_acm_2017_ann, outcomes_comp$n_acm_2018_ann
  )

  # 2020
  # QALYs
  n_qalys_2020_int <- sum(
    outcomes_int$n_qalys_2010_ann, outcomes_int$n_qalys_2011_ann, outcomes_int$n_qalys_2012_ann, outcomes_int$n_qalys_2013_ann, outcomes_int$n_qalys_2014_ann, outcomes_int$n_qalys_2015_ann,
    outcomes_int$n_qalys_2016_ann, outcomes_int$n_qalys_2017_ann, outcomes_int$n_qalys_2018_ann, outcomes_int$n_qalys_2019_ann, outcomes_int$n_qalys_2020_ann
  )
  n_qalys_2020_comp <- sum(
    outcomes_comp$n_qalys_2010_ann, outcomes_comp$n_qalys_2011_ann, outcomes_comp$n_qalys_2012_ann, outcomes_comp$n_qalys_2013_ann, outcomes_comp$n_qalys_2014_ann, outcomes_comp$n_qalys_2015_ann,
    outcomes_comp$n_qalys_2016_ann, outcomes_comp$n_qalys_2017_ann, outcomes_comp$n_qalys_2018_ann, outcomes_comp$n_qalys_2019_ann, outcomes_comp$n_qalys_2020_ann
  )
  # Fatal overdoses
  n_odf_2020_int <- sum(
    outcomes_int$n_odf_2010_ann, outcomes_int$n_odf_2011_ann, outcomes_int$n_odf_2012_ann, outcomes_int$n_odf_2013_ann, outcomes_int$n_odf_2014_ann, outcomes_int$n_odf_2015_ann,
    outcomes_int$n_odf_2016_ann, outcomes_int$n_odf_2017_ann, outcomes_int$n_odf_2018_ann, outcomes_int$n_odf_2019_ann, outcomes_int$n_odf_2020_ann
  )
  n_odf_2020_comp <- sum(
    outcomes_comp$n_odf_2010_ann, outcomes_comp$n_odf_2011_ann, outcomes_comp$n_odf_2012_ann, outcomes_comp$n_odf_2013_ann, outcomes_comp$n_odf_2014_ann, outcomes_comp$n_odf_2015_ann,
    outcomes_comp$n_odf_2016_ann, outcomes_comp$n_odf_2017_ann, outcomes_comp$n_odf_2018_ann, outcomes_comp$n_odf_2019_ann, outcomes_comp$n_odf_2020_ann
  )
  # Non-fatal overdoses
  n_odn_2020_int <- sum(
    outcomes_int$n_odn_2010_ann, outcomes_int$n_odn_2011_ann, outcomes_int$n_odn_2012_ann, outcomes_int$n_odn_2013_ann, outcomes_int$n_odn_2014_ann, outcomes_int$n_odn_2015_ann,
    outcomes_int$n_odn_2016_ann, outcomes_int$n_odn_2017_ann, outcomes_int$n_odn_2018_ann, outcomes_int$n_odn_2019_ann, outcomes_int$n_odn_2020_ann
  )
  n_odn_2020_comp <- sum(
    outcomes_comp$n_odn_2010_ann, outcomes_comp$n_odn_2011_ann, outcomes_comp$n_odn_2012_ann, outcomes_comp$n_odn_2013_ann, outcomes_comp$n_odn_2014_ann, outcomes_comp$n_odn_2015_ann,
    outcomes_comp$n_odn_2016_ann, outcomes_comp$n_odn_2017_ann, outcomes_comp$n_odn_2018_ann, outcomes_comp$n_odn_2019_ann, outcomes_comp$n_odn_2020_ann
  )
  # All-cause deaths
  n_acm_2020_int <- sum(
    outcomes_int$n_acm_2010_ann, outcomes_int$n_acm_2011_ann, outcomes_int$n_acm_2012_ann, outcomes_int$n_acm_2013_ann, outcomes_int$n_acm_2014_ann, outcomes_int$n_acm_2015_ann,
    outcomes_int$n_acm_2016_ann, outcomes_int$n_acm_2017_ann, outcomes_int$n_acm_2018_ann, outcomes_int$n_acm_2019_ann, outcomes_int$n_acm_2020_ann
  )
  n_acm_2020_comp <- sum(
    outcomes_comp$n_acm_2010_ann, outcomes_comp$n_acm_2011_ann, outcomes_comp$n_acm_2012_ann, outcomes_comp$n_acm_2013_ann, outcomes_comp$n_acm_2014_ann, outcomes_comp$n_acm_2015_ann,
    outcomes_comp$n_acm_2016_ann, outcomes_comp$n_acm_2017_ann, outcomes_comp$n_acm_2018_ann, outcomes_comp$n_acm_2019_ann, outcomes_comp$n_acm_2020_ann
  )

  ###############################
  ### Total outcomes (scaled) ###
  ###############################
  # 2010
  # QALYs
  n_qalys_2010_int_scaled <- outcomes_int$n_qalys_2010_ann_scaled
  n_qalys_2010_comp_scaled <- outcomes_comp$n_qalys_2010_ann_scaled
  # Fatal overdoses
  n_odf_2010_int_scaled <- outcomes_int$n_odf_2010_ann_scaled
  n_odf_2010_comp_scaled <- outcomes_comp$n_odf_2010_ann_scaled
  # Fatal overdoses
  n_odn_2010_int_scaled <- outcomes_int$n_odn_2010_ann_scaled
  n_odn_2010_comp_scaled <- outcomes_comp$n_odn_2010_ann_scaled
  # All-cause deaths
  n_acm_2010_int_scaled <- outcomes_int$n_acm_2010_ann_scaled
  n_acm_2010_comp_scaled <- outcomes_comp$n_acm_2010_ann_scaled

  # 2012
  # QALYs
  n_qalys_2012_int_scaled <- sum(outcomes_int$n_qalys_2010_ann_scaled, outcomes_int$n_qalys_2011_ann_scaled, outcomes_int$n_qalys_2012_ann_scaled)
  n_qalys_2012_comp_scaled <- sum(outcomes_comp$n_qalys_2010_ann_scaled, outcomes_comp$n_qalys_2011_ann_scaled, outcomes_comp$n_qalys_2012_ann_scaled)
  # Fatal overdoses
  n_odf_2012_int_scaled <- sum(outcomes_int$n_odf_2010_ann_scaled, outcomes_int$n_odf_2011_ann_scaled, outcomes_int$n_odf_2012_ann_scaled)
  n_odf_2012_comp_scaled <- sum(outcomes_comp$n_odf_2010_ann_scaled, outcomes_comp$n_odf_2011_ann_scaled, outcomes_comp$n_odf_2012_ann_scaled)
  # Non-fatal overdoses
  n_odn_2012_int_scaled <- sum(outcomes_int$n_odn_2010_ann_scaled, outcomes_int$n_odn_2011_ann_scaled, outcomes_int$n_odn_2012_ann_scaled)
  n_odn_2012_comp_scaled <- sum(outcomes_comp$n_odn_2010_ann_scaled, outcomes_comp$n_odn_2011_ann_scaled, outcomes_comp$n_odn_2012_ann_scaled)
  # All-cause deaths
  n_acm_2012_int_scaled <- sum(outcomes_int$n_acm_2010_ann_scaled, outcomes_int$n_acm_2011_ann_scaled, outcomes_int$n_acm_2012_ann_scaled)
  n_acm_2012_comp_scaled <- sum(outcomes_comp$n_acm_2010_ann_scaled, outcomes_comp$n_acm_2011_ann_scaled, outcomes_comp$n_acm_2012_ann_scaled)

  # 2014
  # QALYs
  n_qalys_2014_int_scaled <- sum(
    outcomes_int$n_qalys_2010_ann_scaled, outcomes_int$n_qalys_2011_ann_scaled, outcomes_int$n_qalys_2012_ann_scaled,
    outcomes_int$n_qalys_2013_ann_scaled, outcomes_int$n_qalys_2014_ann_scaled
  )
  n_qalys_2014_comp_scaled <- sum(
    outcomes_comp$n_qalys_2010_ann_scaled, outcomes_comp$n_qalys_2011_ann_scaled, outcomes_comp$n_qalys_2012_ann_scaled,
    outcomes_comp$n_qalys_2013_ann_scaled, outcomes_comp$n_qalys_2014_ann_scaled
  )
  # Fatal overdoses
  n_odf_2014_int_scaled <- sum(
    outcomes_int$n_odf_2010_ann_scaled, outcomes_int$n_odf_2011_ann_scaled, outcomes_int$n_odf_2012_ann_scaled,
    outcomes_int$n_odf_2013_ann_scaled, outcomes_int$n_odf_2014_ann_scaled
  )
  n_odf_2014_comp_scaled <- sum(
    outcomes_comp$n_odf_2010_ann_scaled, outcomes_comp$n_odf_2011_ann_scaled, outcomes_comp$n_odf_2012_ann_scaled,
    outcomes_comp$n_odf_2013_ann_scaled, outcomes_comp$n_odf_2014_ann_scaled
  )
  # Non-fatal overdoses
  n_odn_2014_int_scaled <- sum(
    outcomes_int$n_odn_2010_ann_scaled, outcomes_int$n_odn_2011_ann_scaled, outcomes_int$n_odn_2012_ann_scaled,
    outcomes_int$n_odn_2013_ann_scaled, outcomes_int$n_odn_2014_ann_scaled
  )
  n_odn_2014_comp_scaled <- sum(
    outcomes_comp$n_odn_2010_ann_scaled, outcomes_comp$n_odn_2011_ann_scaled, outcomes_comp$n_odn_2012_ann_scaled,
    outcomes_comp$n_odn_2013_ann_scaled, outcomes_comp$n_odn_2014_ann_scaled
  )
  # All-cause deaths
  n_acm_2014_int_scaled <- sum(
    outcomes_int$n_acm_2010_ann_scaled, outcomes_int$n_acm_2011_ann_scaled, outcomes_int$n_acm_2012_ann_scaled,
    outcomes_int$n_acm_2013_ann_scaled, outcomes_int$n_acm_2014_ann_scaled
  )
  n_acm_2014_comp_scaled <- sum(
    outcomes_comp$n_acm_2010_ann_scaled, outcomes_comp$n_acm_2011_ann_scaled, outcomes_comp$n_acm_2012_ann_scaled,
    outcomes_comp$n_acm_2013_ann_scaled, outcomes_comp$n_acm_2014_ann_scaled
  )

  # 2016
  # QALYs
  n_qalys_2016_int_scaled <- sum(
    outcomes_int$n_qalys_2010_ann_scaled, outcomes_int$n_qalys_2011_ann_scaled, outcomes_int$n_qalys_2012_ann_scaled,
    outcomes_int$n_qalys_2013_ann_scaled, outcomes_int$n_qalys_2014_ann_scaled, outcomes_int$n_qalys_2015_ann_scaled, outcomes_int$n_qalys_2016_ann_scaled
  )
  n_qalys_2016_comp_scaled <- sum(
    outcomes_comp$n_qalys_2010_ann_scaled, outcomes_comp$n_qalys_2011_ann_scaled, outcomes_comp$n_qalys_2012_ann_scaled,
    outcomes_comp$n_qalys_2013_ann_scaled, outcomes_comp$n_qalys_2014_ann_scaled, outcomes_comp$n_qalys_2015_ann_scaled, outcomes_comp$n_qalys_2016_ann_scaled
  )
  # Fatal overdoses
  n_odf_2016_int_scaled <- sum(
    outcomes_int$n_odf_2010_ann_scaled, outcomes_int$n_odf_2011_ann_scaled, outcomes_int$n_odf_2012_ann_scaled,
    outcomes_int$n_odf_2013_ann_scaled, outcomes_int$n_odf_2014_ann_scaled, outcomes_int$n_odf_2015_ann_scaled, outcomes_int$n_odf_2016_ann_scaled
  )
  n_odf_2016_comp_scaled <- sum(
    outcomes_comp$n_odf_2010_ann_scaled, outcomes_comp$n_odf_2011_ann_scaled, outcomes_comp$n_odf_2012_ann_scaled,
    outcomes_comp$n_odf_2013_ann_scaled, outcomes_comp$n_odf_2014_ann_scaled, outcomes_comp$n_odf_2015_ann_scaled, outcomes_comp$n_odf_2016_ann_scaled
  )
  # Non-fatal overdoses
  n_odn_2016_int_scaled <- sum(
    outcomes_int$n_odn_2010_ann_scaled, outcomes_int$n_odn_2011_ann_scaled, outcomes_int$n_odn_2012_ann_scaled,
    outcomes_int$n_odn_2013_ann_scaled, outcomes_int$n_odn_2014_ann_scaled, outcomes_int$n_odn_2015_ann_scaled, outcomes_int$n_odn_2016_ann_scaled
  )
  n_odn_2016_comp_scaled <- sum(
    outcomes_comp$n_odn_2010_ann_scaled, outcomes_comp$n_odn_2011_ann_scaled, outcomes_comp$n_odn_2012_ann_scaled,
    outcomes_comp$n_odn_2013_ann_scaled, outcomes_comp$n_odn_2014_ann_scaled, outcomes_comp$n_odn_2015_ann_scaled, outcomes_comp$n_odn_2016_ann_scaled
  )
  # All-cause deaths
  n_acm_2016_int_scaled <- sum(
    outcomes_int$n_acm_2010_ann_scaled, outcomes_int$n_acm_2011_ann_scaled, outcomes_int$n_acm_2012_ann_scaled,
    outcomes_int$n_acm_2013_ann_scaled, outcomes_int$n_acm_2014_ann_scaled, outcomes_int$n_acm_2015_ann_scaled, outcomes_int$n_acm_2016_ann_scaled
  )
  n_acm_2016_comp_scaled <- sum(
    outcomes_comp$n_acm_2010_ann_scaled, outcomes_comp$n_acm_2011_ann_scaled, outcomes_comp$n_acm_2012_ann_scaled,
    outcomes_comp$n_acm_2013_ann_scaled, outcomes_comp$n_acm_2014_ann_scaled, outcomes_comp$n_acm_2015_ann_scaled, outcomes_comp$n_acm_2016_ann_scaled
  )

  # 2018
  # QALYs
  n_qalys_2018_int_scaled <- sum(
    outcomes_int$n_qalys_2010_ann_scaled, outcomes_int$n_qalys_2011_ann_scaled, outcomes_int$n_qalys_2012_ann_scaled, outcomes_int$n_qalys_2013_ann_scaled,
    outcomes_int$n_qalys_2014_ann_scaled, outcomes_int$n_qalys_2015_ann_scaled, outcomes_int$n_qalys_2016_ann_scaled, outcomes_int$n_qalys_2017_ann_scaled, outcomes_int$n_qalys_2018_ann_scaled
  )
  n_qalys_2018_comp_scaled <- sum(
    outcomes_comp$n_qalys_2010_ann_scaled, outcomes_comp$n_qalys_2011_ann_scaled, outcomes_comp$n_qalys_2012_ann_scaled, outcomes_comp$n_qalys_2013_ann_scaled,
    outcomes_comp$n_qalys_2014_ann_scaled, outcomes_comp$n_qalys_2015_ann_scaled, outcomes_comp$n_qalys_2016_ann_scaled, outcomes_comp$n_qalys_2017_ann_scaled, outcomes_comp$n_qalys_2018_ann_scaled
  )
  # Fatal overdoses
  n_odf_2018_int_scaled <- sum(
    outcomes_int$n_odf_2010_ann_scaled, outcomes_int$n_odf_2011_ann_scaled, outcomes_int$n_odf_2012_ann_scaled, outcomes_int$n_odf_2013_ann_scaled,
    outcomes_int$n_odf_2014_ann_scaled, outcomes_int$n_odf_2015_ann_scaled, outcomes_int$n_odf_2016_ann_scaled, outcomes_int$n_odf_2017_ann_scaled, outcomes_int$n_odf_2018_ann_scaled
  )
  n_odf_2018_comp_scaled <- sum(
    outcomes_comp$n_odf_2010_ann_scaled, outcomes_comp$n_odf_2011_ann_scaled, outcomes_comp$n_odf_2012_ann_scaled, outcomes_comp$n_odf_2013_ann_scaled,
    outcomes_comp$n_odf_2014_ann_scaled, outcomes_comp$n_odf_2015_ann_scaled, outcomes_comp$n_odf_2016_ann_scaled, outcomes_comp$n_odf_2017_ann_scaled, outcomes_comp$n_odf_2018_ann_scaled
  )
  # Non-fatal overdoses
  n_odn_2018_int_scaled <- sum(
    outcomes_int$n_odn_2010_ann_scaled, outcomes_int$n_odn_2011_ann_scaled, outcomes_int$n_odn_2012_ann_scaled, outcomes_int$n_odn_2013_ann_scaled,
    outcomes_int$n_odn_2014_ann_scaled, outcomes_int$n_odn_2015_ann_scaled, outcomes_int$n_odn_2016_ann_scaled, outcomes_int$n_odn_2017_ann_scaled, outcomes_int$n_odn_2018_ann_scaled
  )
  n_odn_2018_comp_scaled <- sum(
    outcomes_comp$n_odn_2010_ann_scaled, outcomes_comp$n_odn_2011_ann_scaled, outcomes_comp$n_odn_2012_ann_scaled, outcomes_comp$n_odn_2013_ann_scaled,
    outcomes_comp$n_odn_2014_ann_scaled, outcomes_comp$n_odn_2015_ann_scaled, outcomes_comp$n_odn_2016_ann_scaled, outcomes_comp$n_odn_2017_ann_scaled, outcomes_comp$n_odn_2018_ann_scaled
  )
  # All-cause deaths
  n_acm_2018_int_scaled <- sum(
    outcomes_int$n_acm_2010_ann_scaled, outcomes_int$n_acm_2011_ann_scaled, outcomes_int$n_acm_2012_ann_scaled, outcomes_int$n_acm_2013_ann_scaled,
    outcomes_int$n_acm_2014_ann_scaled, outcomes_int$n_acm_2015_ann_scaled, outcomes_int$n_acm_2016_ann_scaled, outcomes_int$n_acm_2017_ann_scaled, outcomes_int$n_acm_2018_ann_scaled
  )
  n_acm_2018_comp_scaled <- sum(
    outcomes_comp$n_acm_2010_ann_scaled, outcomes_comp$n_acm_2011_ann_scaled, outcomes_comp$n_acm_2012_ann_scaled, outcomes_comp$n_acm_2013_ann_scaled,
    outcomes_comp$n_acm_2014_ann_scaled, outcomes_comp$n_acm_2015_ann_scaled, outcomes_comp$n_acm_2016_ann_scaled, outcomes_comp$n_acm_2017_ann_scaled, outcomes_comp$n_acm_2018_ann_scaled
  )

  # 2020
  # QALYs
  n_qalys_2020_int_scaled <- sum(
    outcomes_int$n_qalys_2010_ann_scaled, outcomes_int$n_qalys_2011_ann_scaled, outcomes_int$n_qalys_2012_ann_scaled, outcomes_int$n_qalys_2013_ann_scaled, outcomes_int$n_qalys_2014_ann_scaled, outcomes_int$n_qalys_2015_ann_scaled,
    outcomes_int$n_qalys_2016_ann_scaled, outcomes_int$n_qalys_2017_ann_scaled, outcomes_int$n_qalys_2018_ann_scaled, outcomes_int$n_qalys_2019_ann_scaled, outcomes_int$n_qalys_2020_ann_scaled
  )
  n_qalys_2020_comp_scaled <- sum(
    outcomes_comp$n_qalys_2010_ann_scaled, outcomes_comp$n_qalys_2011_ann_scaled, outcomes_comp$n_qalys_2012_ann_scaled, outcomes_comp$n_qalys_2013_ann_scaled, outcomes_comp$n_qalys_2014_ann_scaled, outcomes_comp$n_qalys_2015_ann_scaled,
    outcomes_comp$n_qalys_2016_ann_scaled, outcomes_comp$n_qalys_2017_ann_scaled, outcomes_comp$n_qalys_2018_ann_scaled, outcomes_comp$n_qalys_2019_ann_scaled, outcomes_comp$n_qalys_2020_ann_scaled
  )
  # Fatal overdoses
  n_odf_2020_int_scaled <- sum(
    outcomes_int$n_odf_2010_ann_scaled, outcomes_int$n_odf_2011_ann_scaled, outcomes_int$n_odf_2012_ann_scaled, outcomes_int$n_odf_2013_ann_scaled, outcomes_int$n_odf_2014_ann_scaled, outcomes_int$n_odf_2015_ann_scaled,
    outcomes_int$n_odf_2016_ann_scaled, outcomes_int$n_odf_2017_ann_scaled, outcomes_int$n_odf_2018_ann_scaled, outcomes_int$n_odf_2019_ann_scaled, outcomes_int$n_odf_2020_ann_scaled
  )
  n_odf_2020_comp_scaled <- sum(
    outcomes_comp$n_odf_2010_ann_scaled, outcomes_comp$n_odf_2011_ann_scaled, outcomes_comp$n_odf_2012_ann_scaled, outcomes_comp$n_odf_2013_ann_scaled, outcomes_comp$n_odf_2014_ann_scaled, outcomes_comp$n_odf_2015_ann_scaled,
    outcomes_comp$n_odf_2016_ann_scaled, outcomes_comp$n_odf_2017_ann_scaled, outcomes_comp$n_odf_2018_ann_scaled, outcomes_comp$n_odf_2019_ann_scaled, outcomes_comp$n_odf_2020_ann_scaled
  )
  # Non-fatal overdoses
  n_odn_2020_int_scaled <- sum(
    outcomes_int$n_odn_2010_ann_scaled, outcomes_int$n_odn_2011_ann_scaled, outcomes_int$n_odn_2012_ann_scaled, outcomes_int$n_odn_2013_ann_scaled, outcomes_int$n_odn_2014_ann_scaled, outcomes_int$n_odn_2015_ann_scaled,
    outcomes_int$n_odn_2016_ann_scaled, outcomes_int$n_odn_2017_ann_scaled, outcomes_int$n_odn_2018_ann_scaled, outcomes_int$n_odn_2019_ann_scaled, outcomes_int$n_odn_2020_ann_scaled
  )
  n_odn_2020_comp_scaled <- sum(
    outcomes_comp$n_odn_2010_ann_scaled, outcomes_comp$n_odn_2011_ann_scaled, outcomes_comp$n_odn_2012_ann_scaled, outcomes_comp$n_odn_2013_ann_scaled, outcomes_comp$n_odn_2014_ann_scaled, outcomes_comp$n_odn_2015_ann_scaled,
    outcomes_comp$n_odn_2016_ann_scaled, outcomes_comp$n_odn_2017_ann_scaled, outcomes_comp$n_odn_2018_ann_scaled, outcomes_comp$n_odn_2019_ann_scaled, outcomes_comp$n_odn_2020_ann_scaled
  )
  # All-cause deaths
  n_acm_2020_int_scaled <- sum(
    outcomes_int$n_acm_2010_ann_scaled, outcomes_int$n_acm_2011_ann_scaled, outcomes_int$n_acm_2012_ann_scaled, outcomes_int$n_acm_2013_ann_scaled, outcomes_int$n_acm_2014_ann_scaled, outcomes_int$n_acm_2015_ann_scaled,
    outcomes_int$n_acm_2016_ann_scaled, outcomes_int$n_acm_2017_ann_scaled, outcomes_int$n_acm_2018_ann_scaled, outcomes_int$n_acm_2019_ann_scaled, outcomes_int$n_acm_2020_ann_scaled
  )
  n_acm_2020_comp_scaled <- sum(
    outcomes_comp$n_acm_2010_ann_scaled, outcomes_comp$n_acm_2011_ann_scaled, outcomes_comp$n_acm_2012_ann_scaled, outcomes_comp$n_acm_2013_ann_scaled, outcomes_comp$n_acm_2014_ann_scaled, outcomes_comp$n_acm_2015_ann_scaled,
    outcomes_comp$n_acm_2016_ann_scaled, outcomes_comp$n_acm_2017_ann_scaled, outcomes_comp$n_acm_2018_ann_scaled, outcomes_comp$n_acm_2019_ann_scaled, outcomes_comp$n_acm_2020_ann_scaled
  )

  #######################################
  ### Incremental outcomes (unscaled) ###
  #######################################
  # 2010 (unadjusted)
  n_inc_qalys_2010 <- n_qalys_2010_int - n_qalys_2010_comp
  n_inc_odf_2010 <- n_odf_2010_int - n_odf_2010_comp
  n_inc_odn_2010 <- n_odn_2010_int - n_odn_2010_comp
  n_inc_acm_2010 <- n_acm_2010_int - n_acm_2010_comp
  # 2010 (adjusted)
  n_inc_qalys_adj_2010 <- n_inc_qalys_2010 / n_alive_int_mean_2010
  n_inc_odf_adj_2010 <- n_inc_odf_2010 / n_alive_int_mean_2010
  n_inc_odn_adj_2010 <- n_inc_odn_2010 / n_alive_int_mean_2010
  n_inc_acm_adj_2010 <- n_inc_acm_2010 / n_alive_int_mean_2010

  # 2012 (unadjusted)
  n_inc_qalys_2012 <- n_qalys_2012_int - n_qalys_2012_comp
  n_inc_odf_2012 <- n_odf_2012_int - n_odf_2012_comp
  n_inc_odn_2012 <- n_odn_2012_int - n_odn_2012_comp
  n_inc_acm_2012 <- n_acm_2012_int - n_acm_2012_comp
  # 2012 (adjusted)
  n_inc_qalys_adj_2012 <- n_inc_qalys_2012 / n_alive_int_mean_2012
  n_inc_odf_adj_2012 <- n_inc_odf_2012 / n_alive_int_mean_2012
  n_inc_odn_adj_2012 <- n_inc_odn_2012 / n_alive_int_mean_2012
  n_inc_acm_adj_2012 <- n_inc_acm_2012 / n_alive_int_mean_2012

  # 2014 (unadjusted)
  n_inc_qalys_2014 <- n_qalys_2014_int - n_qalys_2014_comp
  n_inc_odf_2014 <- n_odf_2014_int - n_odf_2014_comp
  n_inc_odn_2014 <- n_odn_2014_int - n_odn_2014_comp
  n_inc_acm_2014 <- n_acm_2014_int - n_acm_2014_comp
  # 2014 (adjusted)
  n_inc_qalys_adj_2014 <- n_inc_qalys_2014 / n_alive_int_mean_2014
  n_inc_odf_adj_2014 <- n_inc_odf_2014 / n_alive_int_mean_2014
  n_inc_odn_adj_2014 <- n_inc_odn_2014 / n_alive_int_mean_2014
  n_inc_acm_adj_2014 <- n_inc_acm_2014 / n_alive_int_mean_2014

  # 2016 (unadjusted)
  n_inc_qalys_2016 <- n_qalys_2016_int - n_qalys_2016_comp
  n_inc_odf_2016 <- n_odf_2016_int - n_odf_2016_comp
  n_inc_odn_2016 <- n_odn_2016_int - n_odn_2016_comp
  n_inc_acm_2016 <- n_acm_2016_int - n_acm_2016_comp
  # 2016 (adjusted)
  n_inc_qalys_adj_2016 <- n_inc_qalys_2016 / n_alive_int_mean_2016
  n_inc_odf_adj_2016 <- n_inc_odf_2016 / n_alive_int_mean_2016
  n_inc_odn_adj_2016 <- n_inc_odn_2016 / n_alive_int_mean_2016
  n_inc_acm_adj_2016 <- n_inc_acm_2016 / n_alive_int_mean_2016

  # 2018 (unadjusted)
  n_inc_qalys_2018 <- n_qalys_2018_int - n_qalys_2018_comp
  n_inc_odf_2018 <- n_odf_2018_int - n_odf_2018_comp
  n_inc_odn_2018 <- n_odn_2018_int - n_odn_2018_comp
  n_inc_acm_2018 <- n_acm_2018_int - n_acm_2018_comp
  # 2018 (adjusted)
  n_inc_qalys_adj_2018 <- n_inc_qalys_2018 / n_alive_int_mean_2018
  n_inc_odf_adj_2018 <- n_inc_odf_2018 / n_alive_int_mean_2018
  n_inc_odn_adj_2018 <- n_inc_odn_2018 / n_alive_int_mean_2018
  n_inc_acm_adj_2018 <- n_inc_acm_2018 / n_alive_int_mean_2018

  # 2020 (unadjusted)
  n_inc_qalys_2020 <- n_qalys_2020_int - n_qalys_2020_comp
  n_inc_odf_2020 <- n_odf_2020_int - n_odf_2020_comp
  n_inc_odn_2020 <- n_odn_2020_int - n_odn_2020_comp
  n_inc_acm_2020 <- n_acm_2020_int - n_acm_2020_comp
  # 2020 (adjusted)
  n_inc_qalys_adj_2020 <- n_inc_qalys_2020 / n_alive_int_mean_2020
  n_inc_odf_adj_2020 <- n_inc_odf_2020 / n_alive_int_mean_2020
  n_inc_odn_adj_2020 <- n_inc_odn_2020 / n_alive_int_mean_2020
  n_inc_acm_adj_2020 <- n_inc_acm_2020 / n_alive_int_mean_2020

  #####################################
  ### Incremental outcomes (scaled) ###
  #####################################
  # 2010 (unadjusted)
  n_inc_qalys_2010_scaled <- n_qalys_2010_int_scaled - n_qalys_2010_comp_scaled
  n_inc_odf_2010_scaled <- n_odf_2010_int_scaled - n_odf_2010_comp_scaled
  n_inc_odn_2010_scaled <- n_odn_2010_int_scaled - n_odn_2010_comp_scaled
  n_inc_acm_2010_scaled <- n_acm_2010_int_scaled - n_acm_2010_comp_scaled
  # 2010 (adjusted)
  n_inc_qalys_adj_2010_scaled <- n_inc_qalys_2010_scaled / n_alive_int_mean_2010
  n_inc_odf_adj_2010_scaled <- n_inc_odf_2010_scaled / n_alive_int_mean_2010
  n_inc_odn_adj_2010_scaled <- n_inc_odn_2010_scaled / n_alive_int_mean_2010
  n_inc_acm_adj_2010_scaled <- n_inc_acm_2010_scaled / n_alive_int_mean_2010

  # 2012 (unadjusted)
  n_inc_qalys_2012_scaled <- n_qalys_2012_int_scaled - n_qalys_2012_comp_scaled
  n_inc_odf_2012_scaled <- n_odf_2012_int_scaled - n_odf_2012_comp_scaled
  n_inc_odn_2012_scaled <- n_odn_2012_int_scaled - n_odn_2012_comp_scaled
  n_inc_acm_2012_scaled <- n_acm_2012_int_scaled - n_acm_2012_comp_scaled
  # 2012 (adjusted)
  n_inc_qalys_adj_2012_scaled <- n_inc_qalys_2012_scaled / n_alive_int_mean_2012
  n_inc_odf_adj_2012_scaled <- n_inc_odf_2012_scaled / n_alive_int_mean_2012
  n_inc_odn_adj_2012_scaled <- n_inc_odn_2012_scaled / n_alive_int_mean_2012
  n_inc_acm_adj_2012_scaled <- n_inc_acm_2012_scaled / n_alive_int_mean_2012

  # 2014 (unadjusted)
  n_inc_qalys_2014_scaled <- n_qalys_2014_int_scaled - n_qalys_2014_comp_scaled
  n_inc_odf_2014_scaled <- n_odf_2014_int_scaled - n_odf_2014_comp_scaled
  n_inc_odn_2014_scaled <- n_odn_2014_int_scaled - n_odn_2014_comp_scaled
  n_inc_acm_2014_scaled <- n_acm_2014_int_scaled - n_acm_2014_comp_scaled
  # 2014 (adjusted)
  n_inc_qalys_adj_2014_scaled <- n_inc_qalys_2014_scaled / n_alive_int_mean_2014
  n_inc_odf_adj_2014_scaled <- n_inc_odf_2014_scaled / n_alive_int_mean_2014
  n_inc_odn_adj_2014_scaled <- n_inc_odn_2014_scaled / n_alive_int_mean_2014
  n_inc_acm_adj_2014_scaled <- n_inc_acm_2014_scaled / n_alive_int_mean_2014

  # 2016 (unadjusted)
  n_inc_qalys_2016_scaled <- n_qalys_2016_int_scaled - n_qalys_2016_comp_scaled
  n_inc_odf_2016_scaled <- n_odf_2016_int_scaled - n_odf_2016_comp_scaled
  n_inc_odn_2016_scaled <- n_odn_2016_int_scaled - n_odn_2016_comp_scaled
  n_inc_acm_2016_scaled <- n_acm_2016_int_scaled - n_acm_2016_comp_scaled
  # 2016 (adjusted)
  n_inc_qalys_adj_2016_scaled <- n_inc_qalys_2016_scaled / n_alive_int_mean_2016
  n_inc_odf_adj_2016_scaled <- n_inc_odf_2016_scaled / n_alive_int_mean_2016
  n_inc_odn_adj_2016_scaled <- n_inc_odn_2016_scaled / n_alive_int_mean_2016
  n_inc_acm_adj_2016_scaled <- n_inc_acm_2016_scaled / n_alive_int_mean_2016

  # 2018 (unadjusted)
  n_inc_qalys_2018_scaled <- n_qalys_2018_int_scaled - n_qalys_2018_comp_scaled
  n_inc_odf_2018_scaled <- n_odf_2018_int_scaled - n_odf_2018_comp_scaled
  n_inc_odn_2018_scaled <- n_odn_2018_int_scaled - n_odn_2018_comp_scaled
  n_inc_acm_2018_scaled <- n_acm_2018_int_scaled - n_acm_2018_comp_scaled
  # 2018 (adjusted)
  n_inc_qalys_adj_2018_scaled <- n_inc_qalys_2018_scaled / n_alive_int_mean_2018
  n_inc_odf_adj_2018_scaled <- n_inc_odf_2018_scaled / n_alive_int_mean_2018
  n_inc_odn_adj_2018_scaled <- n_inc_odn_2018_scaled / n_alive_int_mean_2018
  n_inc_acm_adj_2018_scaled <- n_inc_acm_2018_scaled / n_alive_int_mean_2018

  # 2020 (unadjusted)
  n_inc_qalys_2020_scaled <- n_qalys_2020_int_scaled - n_qalys_2020_comp_scaled
  n_inc_odf_2020_scaled <- n_odf_2020_int_scaled - n_odf_2020_comp_scaled
  n_inc_odn_2020_scaled <- n_odn_2020_int_scaled - n_odn_2020_comp_scaled
  n_inc_acm_2020_scaled <- n_acm_2020_int_scaled - n_acm_2020_comp_scaled
  # 2020 (adjusted)
  n_inc_qalys_adj_2020_scaled <- n_inc_qalys_2020_scaled / n_alive_int_mean_2020
  n_inc_odf_adj_2020_scaled <- n_inc_odf_2020_scaled / n_alive_int_mean_2020
  n_inc_odn_adj_2020_scaled <- n_inc_odn_2020_scaled / n_alive_int_mean_2020
  n_inc_acm_adj_2020_scaled <- n_inc_acm_2020_scaled / n_alive_int_mean_2020

  ## Combined dataframes ##
  df_incremental <- data.frame(
    # QALYs
    n_inc_qalys_2010, n_inc_qalys_2012, n_inc_qalys_2014, n_inc_qalys_2016, n_inc_qalys_2018, n_inc_qalys_2020,
    n_inc_qalys_adj_2010, n_inc_qalys_adj_2012, n_inc_qalys_adj_2014, n_inc_qalys_adj_2016, n_inc_qalys_adj_2018, n_inc_qalys_adj_2020,
    # Fatal overdoses
    n_inc_odf_2010, n_inc_odf_2012, n_inc_odf_2014, n_inc_odf_2016, n_inc_odf_2018, n_inc_odf_2020,
    n_inc_odf_adj_2010, n_inc_odf_adj_2012, n_inc_odf_adj_2014, n_inc_odf_adj_2016, n_inc_odf_adj_2018, n_inc_odf_adj_2020,
    # Non-fatal overdoses
    n_inc_odn_2010, n_inc_odn_2012, n_inc_odn_2014, n_inc_odn_2016, n_inc_odn_2018, n_inc_odn_2020,
    n_inc_odn_adj_2010, n_inc_odn_adj_2012, n_inc_odn_adj_2014, n_inc_odn_adj_2016, n_inc_odn_adj_2018, n_inc_odn_adj_2020,
    # All-cause deaths
    n_inc_acm_2010, n_inc_acm_2012, n_inc_acm_2014, n_inc_acm_2016, n_inc_acm_2018, n_inc_acm_2020,
    n_inc_acm_adj_2010, n_inc_acm_adj_2012, n_inc_acm_adj_2014, n_inc_acm_adj_2016, n_inc_acm_adj_2018, n_inc_acm_adj_2020
  )
  df_incremental_scaled <- data.frame(
    # QALYs
    n_inc_qalys_2010_scaled, n_inc_qalys_2012_scaled, n_inc_qalys_2014_scaled, n_inc_qalys_2016_scaled, n_inc_qalys_2018_scaled, n_inc_qalys_2020_scaled,
    n_inc_qalys_adj_2010_scaled, n_inc_qalys_adj_2012_scaled, n_inc_qalys_adj_2014_scaled, n_inc_qalys_adj_2016_scaled, n_inc_qalys_adj_2018_scaled, n_inc_qalys_adj_2020_scaled,
    # Fatal overdoses
    n_inc_odf_2010_scaled, n_inc_odf_2012_scaled, n_inc_odf_2014_scaled, n_inc_odf_2016_scaled, n_inc_odf_2018_scaled, n_inc_odf_2020_scaled,
    n_inc_odf_adj_2010_scaled, n_inc_odf_adj_2012_scaled, n_inc_odf_adj_2014_scaled, n_inc_odf_adj_2016_scaled, n_inc_odf_adj_2018_scaled, n_inc_odf_adj_2020_scaled,
    # Non-fatal overdoses
    n_inc_odn_2010_scaled, n_inc_odn_2012_scaled, n_inc_odn_2014_scaled, n_inc_odn_2016_scaled, n_inc_odn_2018_scaled, n_inc_odn_2020_scaled,
    n_inc_odn_adj_2010_scaled, n_inc_odn_adj_2012_scaled, n_inc_odn_adj_2014_scaled, n_inc_odn_adj_2016_scaled, n_inc_odn_adj_2018_scaled, n_inc_odn_adj_2020_scaled,
    # All-cause deaths
    n_inc_acm_2010_scaled, n_inc_acm_2012_scaled, n_inc_acm_2014_scaled, n_inc_acm_2016_scaled, n_inc_acm_2018_scaled, n_inc_acm_2020_scaled,
    n_inc_acm_adj_2010_scaled, n_inc_acm_adj_2012_scaled, n_inc_acm_adj_2014_scaled, n_inc_acm_adj_2016_scaled, n_inc_acm_adj_2018_scaled, n_inc_acm_adj_2020_scaled
  )

  return(list(
    df_incremental = df_incremental,
    df_incremental_scaled = df_incremental_scaled
  ))
}

# Function to combine parallel outputs
combine_custom_twsa <- function(LL1, LL2) {
  df_outcomes_met_twsa <- rbind(LL1$df_outcomes_met_twsa, LL2$df_outcomes_met_twsa)
  df_outcomes_bnx_twsa <- rbind(LL1$df_outcomes_bnx_twsa, LL2$df_outcomes_bnx_twsa)
  df_incremental_twsa <- rbind(LL1$df_incremental_twsa, LL2$df_incremental_twsa)
  df_incremental_twsa_scaled <- rbind(LL1$df_incremental_twsa_scaled, LL2$df_incremental_twsa_scaled)

  return(list(
    df_outcomes_met_twsa = df_outcomes_met_twsa,
    df_outcomes_bnx_twsa = df_outcomes_bnx_twsa,
    df_incremental_twsa = df_incremental_twsa,
    df_incremental_twsa_scaled = df_incremental_twsa_scaled
  ))
}
