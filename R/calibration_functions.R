#' Generate model outputs for calibration from a parameter set
#'
#' \code{calibration_out} computes model outputs to be used for calibration
#' routines.
#'
#' @param v_params_calib is a vector of parameters that need to be calibrated.
#' @param l_params_all is a list with all parameters of the decision model.
#' @return
#' A list with all cause deaths, and non-fatal overdoses.
#' @export
calibration_out <- function(v_params_updated,
                            l_params_all,
                            time_horizon = "cali",
                            ce_est = "cali",
                            analytic_cohort = "cali") {
  # Substitute values of calibrated parameters in base-case with calibrated values
  l_params_all <- update_param_list(l_params_all = l_params_all, params_updated = v_params_updated)

  # Run model with updated calibrated parameters
  l_out_markov <- markov_model(l_params_all = l_params_all, time_horizon = time_horizon, ce_est = ce_est, analytic_cohort = analytic_cohort)

  # TODO: REMOVE
  # m_trace <- l_out_markov$m_m_agg_trace

  #### Model outputs ####
  ### Overdose deaths ###
  v_odf <- l_out_markov$m_m_agg_trace[, "ODF"] # cumulative deaths at time i

  ### Non-overdose deaths ###
  v_dno <- l_out_markov$m_m_agg_trace[, "Non-Overdose Death"] # cumulative deaths at time i

  ### Non-fatal overdoses ###
  v_odn <- l_out_markov$m_m_agg_trace[, "ODN"] # cumulative non-fatal overdoses at time i

  ### Time spent out of treatment ###
  v_oot <- l_out_markov$m_m_agg_trace[, "OOT"]

  ### Time spent in long-term abstinence ###
  v_abs <- l_out_markov$m_m_agg_trace[, "ABS"]

  ### Cohort alive ###
  v_alive <- l_out_markov$m_m_agg_trace[, "Alive"] # cohort alive at time i

  ### Select time-points ###
  ### Overdose deaths (2011-2020) ###
  n_per_2012 <- l_params_all$n_per_2012
  n_per_2013 <- l_params_all$n_per_2013
  n_per_2014 <- l_params_all$n_per_2014
  n_per_2015 <- l_params_all$n_per_2015
  n_per_2016 <- l_params_all$n_per_2016
  n_per_2017 <- l_params_all$n_per_2017
  n_per_2018 <- l_params_all$n_per_2018
  n_per_2019 <- l_params_all$n_per_2019
  n_per_2020 <- l_params_all$n_per_2020

  ### Cohort size (2010-2020) ###
  ### Note: cohort size (number alive) calculated at beginning of year versus overdoses calculated at end of year
  # Adjustment starts with year 2
  n_alive_t1 <- v_alive[n_per_2012]
  n_alive_t2 <- v_alive[n_per_2013]
  n_alive_t3 <- v_alive[n_per_2014]
  n_alive_t4 <- v_alive[n_per_2015]
  n_alive_t5 <- v_alive[n_per_2016]
  n_alive_t6 <- v_alive[n_per_2017]
  n_alive_t7 <- v_alive[n_per_2018]
  n_alive_t8 <- v_alive[n_per_2019]
  n_alive_t9 <- v_alive[n_per_2020]

  ### Subset output by time-points ###
  ### Overdose deaths ###
  # Yearly fatal overdoses
  n_dno1 <- v_dno[n_per_2012]
  n_dno2 <- (v_dno[n_per_2013] - v_dno[n_per_2012]) / n_alive_t2
  n_dno3 <- (v_dno[n_per_2014] - v_dno[n_per_2013]) / n_alive_t3
  n_dno4 <- (v_dno[n_per_2015] - v_dno[n_per_2014]) / n_alive_t4
  n_dno5 <- (v_dno[n_per_2016] - v_dno[n_per_2015]) / n_alive_t5
  n_dno6 <- (v_dno[n_per_2017] - v_dno[n_per_2016]) / n_alive_t6
  n_dno7 <- (v_dno[n_per_2018] - v_dno[n_per_2017]) / n_alive_t7
  n_dno8 <- (v_dno[n_per_2019] - v_dno[n_per_2018]) / n_alive_t8
  n_dno9 <- (v_dno[n_per_2020] - v_dno[n_per_2019]) / n_alive_t9

  ### Non-overdose deaths ###
  # Yearly
  n_odf1 <- v_odf[n_per_2012]
  n_odf2 <- (v_odf[n_per_2013] - v_odf[n_per_2012]) / n_alive_t2
  n_odf3 <- (v_odf[n_per_2014] - v_odf[n_per_2013]) / n_alive_t3
  n_odf4 <- (v_odf[n_per_2015] - v_odf[n_per_2014]) / n_alive_t4
  n_odf5 <- (v_odf[n_per_2016] - v_odf[n_per_2015]) / n_alive_t5
  n_odf6 <- (v_odf[n_per_2017] - v_odf[n_per_2016]) / n_alive_t6
  n_odf7 <- (v_odf[n_per_2018] - v_odf[n_per_2017]) / n_alive_t7
  n_odf8 <- (v_odf[n_per_2019] - v_odf[n_per_2018]) / n_alive_t8
  n_odf9 <- (v_odf[n_per_2020] - v_odf[n_per_2019]) / n_alive_t9

  ### Non-fatal overdose ###
  # Non-fatal overdoses need to be summed across time points to generate cumulative estimates
  # Yearly non-fatal overdose
  n_odn1 <- sum(v_odn[c(1:n_per_2012)])
  n_odn2 <- sum(v_odn[c((n_per_2012 + 1):n_per_2013)]) / n_alive_t2
  n_odn3 <- sum(v_odn[c((n_per_2013 + 1):n_per_2014)]) / n_alive_t3
  n_odn4 <- sum(v_odn[c((n_per_2014 + 1):n_per_2015)]) / n_alive_t4
  n_odn5 <- sum(v_odn[c((n_per_2015 + 1):n_per_2016)]) / n_alive_t5
  n_odn6 <- sum(v_odn[c((n_per_2016 + 1):n_per_2017)]) / n_alive_t6
  n_odn7 <- sum(v_odn[c((n_per_2017 + 1):n_per_2018)]) / n_alive_t7
  n_odn8 <- sum(v_odn[c((n_per_2018 + 1):n_per_2019)]) / n_alive_t8
  n_odn9 <- sum(v_odn[c((n_per_2019 + 1):n_per_2020)]) / n_alive_t9

  ### Average yearly time spent out of treatment ###
  n_oot1 <- mean(v_oot[c(1:n_per_2012)])
  n_oot2 <- mean(v_oot[c((n_per_2012 + 1):n_per_2013)]) / n_alive_t2
  n_oot3 <- mean(v_oot[c((n_per_2013 + 1):n_per_2014)]) / n_alive_t3
  n_oot4 <- mean(v_oot[c((n_per_2014 + 1):n_per_2015)]) / n_alive_t4
  n_oot5 <- mean(v_oot[c((n_per_2015 + 1):n_per_2016)]) / n_alive_t5
  n_oot6 <- mean(v_oot[c((n_per_2016 + 1):n_per_2017)]) / n_alive_t6
  n_oot7 <- mean(v_oot[c((n_per_2017 + 1):n_per_2018)]) / n_alive_t7
  n_oot8 <- mean(v_oot[c((n_per_2018 + 1):n_per_2019)]) / n_alive_t8
  n_oot9 <- mean(v_oot[c((n_per_2019 + 1):n_per_2020)]) / n_alive_t9

  ### Average yearly time spent in abstinence ###
  n_abs1 <- mean(v_abs[c(1:n_per_2012)])
  n_abs2 <- mean(v_abs[c((n_per_2012 + 1):n_per_2013)]) / n_alive_t2
  n_abs3 <- mean(v_abs[c((n_per_2013 + 1):n_per_2014)]) / n_alive_t3
  n_abs4 <- mean(v_abs[c((n_per_2014 + 1):n_per_2015)]) / n_alive_t4
  n_abs5 <- mean(v_abs[c((n_per_2015 + 1):n_per_2016)]) / n_alive_t5
  n_abs6 <- mean(v_abs[c((n_per_2016 + 1):n_per_2017)]) / n_alive_t6
  n_abs7 <- mean(v_abs[c((n_per_2017 + 1):n_per_2018)]) / n_alive_t7
  n_abs8 <- mean(v_abs[c((n_per_2018 + 1):n_per_2019)]) / n_alive_t8
  n_abs9 <- mean(v_abs[c((n_per_2019 + 1):n_per_2020)]) / n_alive_t9

  #### Return Output ####
  l_out <- list(
    fatal_overdose = c(n_odf1, n_odf2, n_odf3, n_odf4, n_odf5, n_odf6, n_odf7, n_odf8, n_odf9),
    death_non_od = c(n_dno1, n_dno2, n_dno3, n_dno4, n_dno5, n_dno6, n_dno7, n_dno8, n_dno9),
    overdose = c(n_odn1, n_odn2, n_odn3, n_odn4, n_odn5, n_odn6, n_odn7, n_odn8, n_odn9),
    time_oot = c(n_oot1, n_oot2, n_oot3, n_oot4, n_oot5, n_oot6, n_oot7, n_oot8, n_oot9),
    time_abs = c(n_abs1, n_abs2, n_abs3, n_abs4, n_abs5, n_abs6, n_abs7, n_abs8, n_abs9) # ,
    # TODO: REMOVE m_trace = m_trace
  )
  return(l_out)
}

## Sample prior distribution ##
sample.prior <- function(n_samp,
                         v_param_names = v_cali_param_names,
                         v_alpha = v_par1,
                         v_beta = v_par2) {
  n_param <- length(v_param_names)
  # random latin hypercube sampling
  m_lhs_unit <- lhs::randomLHS(n = n_samp, k = n_param)
  m_param_samp <- matrix(nrow = n_samp, ncol = n_param)
  colnames(m_param_samp) <- v_param_names

  # draw parameters
  draws <- data.frame(
    n_oat_od = qgamma(m_lhs_unit[, 1], shape = v_alpha[1], scale = v_beta[1]), # baseline OD rate in all OAT
    n_fent_prev_od_mult = qunif(m_lhs_unit[, 2], min = v_alpha[2], max = v_beta[2]), # fentanyl multiplier for OD rate
    n_fent_delta_od_mult = qunif(m_lhs_unit[, 3], min = v_alpha[3], max = v_beta[3]), # fentanyl multiplier for OD rate
    n_fatal_od_oat = qgamma(m_lhs_unit[, 4], shape = v_alpha[4], scale = v_beta[4]), # fatal OD rate
    hr_oat = qgamma(m_lhs_unit[, 5], shape = v_alpha[5], scale = v_beta[5]), # non overdose mortality mult in OAT
    p_weibull_scale_ou = qgamma(m_lhs_unit[, 6], shape = v_alpha[6], scale = v_beta[6]), # weibull scale on remain out of treatment
    p_weibull_shape_ou = qgamma(m_lhs_unit[, 7], shape = v_alpha[7], scale = v_beta[7]), # weibull shape on remain out of treatment
    p_weibull_scale_abs = qgamma(m_lhs_unit[, 8], shape = v_alpha[8], scale = v_beta[8]), # weibull scale on remain out of treatment
    p_weibull_shape_abs = qgamma(m_lhs_unit[, 9], shape = v_alpha[9], scale = v_beta[9]) # weibull shape on remain out of treatment
  )

  return(as.matrix(draws))
}

#### Log prior ####
log_prior <- function(v_params,
                      v_param_names = v_cali_param_names,
                      v_alpha = v_par1,
                      v_beta = v_par2) {
  if (is.null(dim(v_params))) { # If vector, change to matrix
    v_params <- t(v_params)
  }
  n_param <- length(v_param_names)
  n_samp <- nrow(v_params)
  colnames(v_params) <- v_param_names
  lprior <- rep(0, n_samp)

  lprior <- lprior + dgamma(v_params[, 1], shape = v_alpha[1], scale = v_beta[1], log = TRUE)
  lprior <- lprior + dunif(v_params[, 2], min = v_alpha[2], max = v_beta[2], log = TRUE)
  lprior <- lprior + dunif(v_params[, 3], min = v_alpha[3], max = v_beta[3], log = TRUE)
  lprior <- lprior + dgamma(v_params[, 4], shape = v_alpha[4], scale = v_beta[4], log = TRUE)
  lprior <- lprior + dgamma(v_params[, 5], shape = v_alpha[5], scale = v_beta[5], log = TRUE)
  lprior <- lprior + dgamma(v_params[, 6], shape = v_alpha[6], scale = v_beta[6], log = TRUE)
  lprior <- lprior + dgamma(v_params[, 7], shape = v_alpha[7], scale = v_beta[7], log = TRUE)
  lprior <- lprior + dgamma(v_params[, 8], shape = v_alpha[8], scale = v_beta[8], log = TRUE)
  lprior <- lprior + dgamma(v_params[, 9], shape = v_alpha[9], scale = v_beta[9], log = TRUE)

  return(lprior)
}

#' Evaluate prior of calibrated parameters
prior <- function(v_params) {
  v_prior <- exp(log_prior(v_params))
  return(v_prior)
}

#' Log-likelihood function for a parameter set
log_lik <- function(v_params) { # User defined
  if (is.null(dim(v_params))) { # If vector, change to matrix
    v_params <- t(v_params)
  }

  n_samp <- nrow(v_params)
  v_target_names <- c("Fatal overdoses", "Non-overdose deaths", "Non-fatal overdoses", "Time out of treatment", "Time in long-term abstinence") # "Non-Overdose Deaths"
  n_target <- length(v_target_names)
  v_llik <- matrix(0, nrow = n_samp, ncol = n_target)
  colnames(v_llik) <- v_target_names
  v_llik_overall <- numeric(n_samp)
  for (j in 1:n_samp) { # j=1
    jj <- tryCatch(
      {
        ###   Run model for parameter set "v_params" ###
        l_model_res <- calibration_out(
          v_params_calib = v_params[j, ],
          l_params_all = l_params_all
        )

        ###  Calculate log-likelihood of model outputs to targets  ###
        ## Uses calibration weights from input file for each year (set all to 1 for equal weight)
        ## TARGET 1: Fatal overdoses ("fatal_overdose")
        ## Normal log-likelihood
        v_llik[j, "Fatal overdoses"] <- sum(dnorm(
          x = l_cali_targets$df_odf$pe,
          mean = l_model_res$fatal_overdose,
          sd = l_cali_targets$df_odf$se,
          log = TRUE
        ) * l_cali_targets$df_odf$weight)
        ## TARGET 2: Non-overdose deaths ("death_non_od")
        ## Normal log-likelihood
        v_llik[j, "Non-overdose deaths"] <- sum(dnorm(
          x = l_cali_targets$df_dno$pe,
          mean = l_model_res$death_non_od,
          sd = l_cali_targets$df_dno$se,
          log = TRUE
        ) * l_cali_targets$df_dno$weight)
        ## TARGET 3: Non-fatal overdoses ("overdose")
        ## Normal log-likelihood
        v_llik[j, "Non-fatal overdoses"] <- sum(dnorm(
          x = l_cali_targets$df_odn$pe,
          mean = l_model_res$overdose,
          sd = l_cali_targets$df_odn$se,
          log = TRUE
        ) * l_cali_targets$df_odn$weight)
        ## TARGET 4: Time spent out-of-treatment ("time_oot")
        ## Normal log-likelihood
        v_llik[j, "Time out of treatment"] <- sum(dnorm(
          x = l_cali_targets$df_oot$pe,
          mean = l_model_res$time_oot,
          sd = l_cali_targets$df_oot$se,
          log = TRUE
        ) * l_cali_targets$df_oot$weight)
        ## TARGET 5: Time in long-term abstinence
        # Note: this target is only one point-estimate average over the entire period
        # Removing the sum over all time-periods should reduce to a single pe
        # Can model as a vector if necessary to match other targets (same across time points)
        ## Normal log-likelihood
        v_llik[j, "Time in long-term abstinence"] <- sum(dnorm(
          x = l_cali_targets$df_abs$pe,
          mean = l_model_res$time_abs,
          sd = l_cali_targets$df_abs$se,
          log = TRUE
        ) * l_cali_targets$df_abs$weight)

        ## targets different weights
        v_weights <- c(1, 1, 0.5, 1, 1) # 100% fatal overdose; 50% non-fatal overdose
        # v_weights <- rep(1, n_target) # set to 1 for equal weight
        ## weighted sum
        v_llik_overall[j] <- v_llik[j, ] %*% v_weights
      },
      error = function(e) NA
    )
    if (is.na(jj)) {
      v_llik_overall <- -Inf
    }
  } ## End loop over sampled parameter sets
  ## return GOF
  return(v_llik_overall)
}

#' Parallel evaluation of log-likelihood function for a sets of parameters
#'
#' \code{log_lik_par} computes a log-likelihood value for one (or multiple)
#' parameter set(s) using parallel computation.
#'
#' @param v_params Vector (or matrix) of model parameters.
#' @return
#' A scalar (or vector) with log-likelihood values.
log_lik_par <- function(v_params,
                        ...) {
  if (is.null(dim(v_params))) { # If vector, change to matrix
    v_params <- t(v_params)
  }
  n_samp <- nrow(v_params)

  ### Get OS
  os <- get_os()
  no_cores <- parallel::detectCores() - 1
  print(paste0("Parallelized Likelihood calculations on ", os, " using ", no_cores, " cores"))
  n_time_init_likpar <- Sys.time()

  if (os == "macosx") {
    # Initialize cluster object
    cl <- parallel::makeForkCluster(no_cores)
    doParallel::registerDoParallel(cl)
    v_llk <- foreach::foreach(i = 1:n_samp, .combine = c) %dopar% {
      log_lik(v_params[i, ]) # i = 1
    }
    n_time_end_likpar <- Sys.time()
  }
  if (os == "windows") {
    # Initialize cluster object
    cl <- parallel::makeCluster(no_cores)
    doParallel::registerDoParallel(cl)
    opts <- list(attachExportEnv = TRUE)
    v_llk <- foreach::foreach(
      i = 1:n_samp, .combine = c,
      .export = ls(globalenv()),
      .packages = c(),
      .options.snow = opts
    ) %dopar% {
      log_lik(v_params[i, ])
    }
    n_time_end_likpar <- Sys.time()
  }
  if (os == "linux") {
    # Initialize cluster object
    cl <- parallel::makeCluster(no_cores)
    doMC::registerDoMC(cl)
    v_llk <- foreach::foreach(i = 1:n_samp, .combine = c) %dopar% {
      log_lik(v_params[i, ])
    }
    n_time_end_likpar <- Sys.time()
  }
  parallel::stopCluster(cl)
  n_time_total_likpar <- difftime(n_time_end_likpar, n_time_init_likpar,
    units = "mins"
  )
  print(paste0("Runtime: ", round(n_time_total_likpar, 2), " mins."))
  #-# Try this: # PO
  rm(cl) # PO
  gc() # PO
  #-#           # PO
  return(v_llk)
}

#' Likelihood
# likelihood <- function(v_params) {
#   v_like <- exp(log_lik(v_params))
#   return(v_like)
# }
# Calculate likelihood using parallel
likelihood <- function(v_params) {
  v_like <- exp(log_lik_par(v_params))
  return(v_like)
}

#' Evaluate log-posterior of calibrated parameters
log_post <- function(v_params) {
  v_lpost <- log_prior(v_params) + log_lik_par(v_params)
  return(v_lpost)
}

#' Evaluate posterior of calibrated parameters
posterior <- function(v_params) {
  v_posterior <- exp(log_post(v_params))
  return(v_posterior)
}

# Function to combine parallel outputs
combine_custom_cali <- function(LL1, LL2) {
  m_model_targets_odf <- rbind(LL1$m_model_targets_odf, LL2$m_model_targets_odf)
  m_model_targets_dno <- rbind(LL1$m_model_targets_dno, LL2$m_model_targets_dno)
  m_model_targets_odn <- rbind(LL1$m_model_targets_odn, LL2$m_model_targets_odn)
  m_model_targets_oot <- rbind(LL1$m_model_targets_oot, LL2$m_model_targets_oot)
  m_model_targets_abs <- rbind(LL1$m_model_targets_abs, LL2$m_model_targets_abs)

  return(list(
    m_model_targets_odf = m_model_targets_odf,
    m_model_targets_dno = m_model_targets_dno,
    m_model_targets_odn = m_model_targets_odn,
    m_model_targets_oot = m_model_targets_oot,
    m_model_targets_abs = m_model_targets_abs
  ))
}

#' Get operating system
#'
#' @return
#' A string with the operating system.
#' @export
get_os <- function() {
  sysinf <- Sys.info()
  if (!is.null(sysinf)) {
    os <- sysinf["sysname"]
    if (os == "Darwin") {
      os <- "MacOSX"
    }
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os)) {
      os <- "osx"
    }
    if (grepl("linux-gnu", R.version$os)) {
      os <- "linux"
    }
  }
  tolower(os)
}
