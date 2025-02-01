#' Generate PSA dataset of CEA parameters
#'
#' \code{generate_psa_params} generates PSA input dataset by sampling decision
#' model parameters from their distributions. The sample of the calibrated
#' parameters is a draw from their posterior distribution obtained with the
#' IMIS algorithm.
#'
#' Parameters that are not sampled in PSA do not need to be defined here, they will
#' default to their original input value for each simulation.
#' @param n_sim Number of PSA samples.
#' @param seed Seed for reproducibility of Monte Carlo sampling.
#' @param n_pop Sample size to determine dirichlet distribution variance.
#' @param scenario Choose analysis scenario (options: 'tx_switch'; 'bnx_only'; 'met_only')
#' @return
#' A data frame with \code{n_sim} rows and {n_states} columns of parameters for PSA.
#' Each row is a parameter set sampled from distributions that characterize
#' their uncertainty
#' @examples
#' generate_psa_params()
#' @export
generate_psa_params <- function(n_sim = n_sim, seed = seed, n_pop = n_pop, # scenario = scenario,
                                file.death_hr = NULL,
                                file.weibull = NULL,
                                file.ce_tx = NULL,
                                file.ce_death = NULL,
                                file.overdose = NULL,
                                file.fentanyl = NULL,
                                file.naloxone = NULL,
                                file.imis_output = NULL) {
  # Load files with parameter distribution values
  df_death_hr <- read.csv(file = file.death_hr, row.names = 1, header = TRUE) # Mortality hazard ratios
  df_weibull <- read.csv(file = file.weibull, row.names = 1, header = TRUE) # Weibull shape and scale
  df_ce_tx <- read.csv(file = file.ce_tx, row.names = 1, header = TRUE) # Comparative effectiveness parameters (treatment retention)
  df_ce_death <- read.csv(file = file.ce_death, row.names = 1, header = TRUE) # Comparative effectiveness parameters (mortality)
  df_overdose <- read.csv(file = file.overdose, row.names = 1, header = TRUE) # Overdose params
  df_fentanyl <- read.csv(file = file.fentanyl, row.names = 1, header = TRUE) # Fentanyl params
  df_naloxone <- read.csv(file = file.naloxone, row.names = 1, header = TRUE) # Time-varying naloxone parameters for calibration

  ## Load calibrated parameters

  load(file = file.imis_output)
  df_calib_post <- as.data.frame(m_calib_post)

  # Number of simulations
  n_sim <- n_sim


  if (n_sim != nrow(df_calib_post)) {
    warning("Number of PSA simulations and posterior draws not equal")
    df_calib_post <- df_calib_post[1:n_sim, ] # Truncate to match simulations
  }

  # Set seed for random number generator
  seed <- seed
  if (!missing(seed)) {
    set.seed(seed)
  }

  # Set sample-size for trial-based parameter uncertainty
  n_pop <- n_pop
  # Function to generate lognormal parameter
  location <- function(m = m, s = s) {
    log(m^2 / sqrt(s^2 + m^2))
  }
  shape <- function(m = m, s = s) {
    shape <- sqrt(log(1 + (s^2 / m^2)))
  }

  #### First week multiplier for overdose in opioid use ####
  n_ou_od_mult <- rgamma(n_sim, shape = df_overdose["shape", "ou_od_mult"], scale = df_overdose["scale", "ou_od_mult"])
  #### s
  hr_death_pp_inc <- rlnorm(n_sim, location(m = df_ce_death["pe", "pp_inc"], s = df_ce_death["sd", "pp_inc"]), shape(m = df_ce_death["pe", "pp_inc"], s = df_ce_death["sd", "pp_inc"]))
  hr_death_pp_prev <- rlnorm(n_sim, location(m = df_ce_death["pe", "pp_prev"], s = df_ce_death["sd", "pp_prev"]), shape(m = df_ce_death["pe", "pp_prev"], s = df_ce_death["sd", "pp_prev"]))

  df_psa_params <- data.frame(
    ### Calibrated parameters
    df_calib_post,

    # Non-overdose death
    hr_dno_ou_oat = rlnorm(n_sim, location(m = df_death_hr["pe", "dno_ou_oat"], s = df_death_hr["sd", "dno_ou_oat"]), shape(m = df_death_hr["pe", "dno_ou_oat"], s = df_death_hr["sd", "dno_ou_oat"])),

    # Comparative effectiveness
    # Log-normal distribution
    # Incidence
    hr_tx_itt_ps_inc = rlnorm(n_sim, location(m = df_ce_tx["pe", "itt_ps_inc"], s = df_ce_tx["sd", "itt_ps_inc"]), shape(m = df_ce_tx["pe", "itt_ps_inc"], s = df_ce_tx["sd", "itt_ps_inc"])),
    hr_tx_pp_inc = rlnorm(n_sim, location(m = df_ce_tx["pe", "pp_inc"], s = df_ce_tx["sd", "pp_inc"]), shape(m = df_ce_tx["pe", "pp_inc"], s = df_ce_tx["sd", "pp_inc"])),
    hr_tx_pp_hd_inc = rlnorm(n_sim, location(m = df_ce_tx["pe", "pp_hd_inc"], s = df_ce_tx["sd", "pp_hd_inc"]), shape(m = df_ce_tx["pe", "pp_hd_inc"], s = df_ce_tx["sd", "pp_hd_inc"])),
    # Prevalent
    hr_tx_itt_ps_prev = rlnorm(n_sim, location(m = df_ce_tx["pe", "itt_ps_prev"], s = df_ce_tx["sd", "itt_ps_prev"]), shape(m = df_ce_tx["pe", "itt_ps_prev"], s = df_ce_tx["sd", "itt_ps_prev"])),
    hr_tx_pp_prev = rlnorm(n_sim, location(m = df_ce_tx["pe", "pp_prev"], s = df_ce_tx["sd", "pp_prev"]), shape(m = df_ce_tx["pe", "pp_prev"], s = df_ce_tx["sd", "pp_prev"])),
    hr_tx_pp_hd_prev = rlnorm(n_sim, location(m = df_ce_tx["pe", "pp_hd_prev"], s = df_ce_tx["sd", "pp_hd_prev"]), shape(m = df_ce_tx["pe", "pp_hd_prev"], s = df_ce_tx["sd", "pp_hd_prev"])),
    # Incident
    hr_death_itt_ps_inc = hr_death_pp_inc,
    hr_death_pp_inc = hr_death_pp_inc,
    hr_death_pp_hd_inc = rlnorm(n_sim, location(m = df_ce_death["pe", "pp_hd_inc"], s = df_ce_death["sd", "pp_hd_inc"]), shape(m = df_ce_death["pe", "pp_hd_inc"], s = df_ce_death["sd", "pp_hd_inc"])),
    # Prevalent
    hr_death_itt_ps_prev = hr_death_pp_prev,
    hr_death_pp_prev = hr_death_pp_prev,
    hr_death_pp_hd_prev = rlnorm(n_sim, location(m = df_ce_death["pe", "pp_hd_prev"], s = df_ce_death["sd", "pp_hd_prev"]), shape(m = df_ce_death["pe", "pp_hd_prev"], s = df_ce_death["sd", "pp_hd_prev"])),

    # Weibull
    # Shape
    # Incident user (EP1)
    p_weibull_shape_met_inc = rnorm(n_sim, mean = df_weibull["pe", "met_shape_inc"], sd = df_weibull["sd", "met_shape_inc"]),
    p_weibull_shape_oat_inc = rnorm(n_sim, mean = df_weibull["pe", "oat_shape_inc"], sd = df_weibull["sd", "oat_shape_inc"]),

    # Prevalent user (EP2+)
    p_weibull_shape_met_prev = rnorm(n_sim, mean = df_weibull["pe", "met_shape_prev"], sd = df_weibull["sd", "met_shape_prev"]),
    p_weibull_shape_oat_prev = rnorm(n_sim, mean = df_weibull["pe", "oat_shape_prev"], sd = df_weibull["sd", "oat_shape_prev"]),

    # Scale
    # Incident user (EP1)
    p_weibull_scale_met_inc = rnorm(n_sim, mean = df_weibull["pe", "met_scale_inc"], sd = df_weibull["sd", "met_scale_inc"]),
    p_weibull_scale_oat_inc = rnorm(n_sim, mean = df_weibull["pe", "oat_scale_inc"], sd = df_weibull["sd", "oat_scale_inc"]),

    # Prevalent user (EP2+)
    p_weibull_scale_met_prev = rnorm(n_sim, mean = df_weibull["pe", "met_scale_prev"], sd = df_weibull["sd", "met_scale_prev"]),
    p_weibull_scale_oat_prev = rnorm(n_sim, mean = df_weibull["pe", "oat_scale_prev"], sd = df_weibull["sd", "oat_scale_prev"]),

    ### Overdose ###
    n_oum_od_mult = n_ou_od_mult,
    n_oub_od_mult = n_ou_od_mult,
    n_ouo_od_mult = n_ou_od_mult,
    n_ou_oat_odn_mult = rgamma(n_sim, shape = df_overdose["shape", "ou_oat_odn_mult"], scale = df_overdose["scale", "ou_oat_odn_mult"]),
    n_ou_oat_odf_mult = rgamma(n_sim, shape = df_overdose["shape", "ou_oat_odf_mult"], scale = df_overdose["scale", "ou_oat_odf_mult"])
  )
  return(df_psa_params)
}
