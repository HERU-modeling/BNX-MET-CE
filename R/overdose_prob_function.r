#### Overdose probability ####
#' Probability of non-fatal and fatal overdose
#'
#' \code{overdose} is used to calculate overdose probabilities from health states. This function also requires additional overdose/fentanyl/naloxone parameters included in `l_params_all`
#'
#' @param l_params_all Import parameter list from model
#' @param rate Baseline overdose rate for health states
#' @param rate_fatal Fatal overdose rate
#' @param rate_fent Fentanyl overdose rate
#' @param multiplier Multiplier for elevated overdose in first month of health state
#' @param fent_mult Multiplier for overdose rate in health states when exposed to fentanyl
#' @param first_week Logical parameter to switch between week 1 and week 2+ for parameter estimation
#' @param fatal Logical parameter to switch between fatal/non-fatal overdose
#' @param injection Logical parameter to adjust rate calculation for injection/non-injection use
#' @param time Time period for time-varying parameters
#'
#' @return
#' `p_OD` monthly probability of fatal or non-fatal overdose from a given health state
#' @export
overdose <- function(
    l_params_all,
    rate,
    rate_fatal,
    multiplier,
    fent_mult,
    fent_delta_mult,
    time,
    ce_fatal_od_mult,
    first_week = FALSE,
    fatal = FALSE,
    bnx = FALSE) {
  with(as.list(l_params_all), {
    # Probability of naloxone use
    p_nx_rev <- v_nx_rev[time]

    # Modifier for CE estimate for BNX fatal OD rate
    # Applied to individuals in treatment only and the week following treatment dropout
    if (!is.null(ce_fatal_od_mult)) {
      p_fatal_od <- 1 - exp(-(rate_fatal * ce_fatal_od_mult))
    } else {
      # Probability of death following all OD
      p_fatal_od <- 1 - exp(-(rate_fatal))
    }
    p_fent_exp <- v_fent_prev[time]
    p_fent_delta <- v_fent_delta[time]

    # Adjustment for change in fentanyl prevalence
    # Prevalence and change are incorporated at this stage
    fent_delta_mult_adj <- exp(fent_delta_mult * p_fent_delta)
    fent_prev_mult_adj <- exp(fent_mult * p_fent_exp)
    # fent_mult_adj <- fent_mult * fent_delta_mult_yr

    # Convert input monthly rates to monthly probabilities - multiply rates by first month multiplier before converting
    if (first_week == TRUE) {
      p_base_od <- 1 - exp(-(rate * multiplier * fent_prev_mult_adj * fent_delta_mult_adj))
    } else if (first_week == FALSE) {
      p_base_od <- 1 - exp(-(rate * fent_prev_mult_adj * fent_delta_mult_adj))
    }

    # Naloxone effect on fatal overdose
    p_fatal_od_nx <- p_fatal_od * (1 - p_nx_rev)

    # Calculate fatal and non-fatal overdose probabilities
    if (fatal == TRUE) {
      p_od <- p_base_od * p_fatal_od_nx
    } else {
      p_od <- p_base_od * (1 - p_fatal_od_nx)
    }
    return(p_od)
  })
}
