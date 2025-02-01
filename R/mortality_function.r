#### Mortality ####
#' Mortality probability estimates
#'
#' \code{mort} is used to populate mortality probability vectors.
#'
#' @param l_params_all Import parameter list from model
#' @param hr
#' @param per
#' @return
#' Mortality vectors for each age applied to model periods (months or weeks), includes state-specific hr.
#' Overdose deaths tracked as "ODF"
#' @export
mort <- function(
    l_params_all,
    hr = hr,
    per = per) {
  with(as.list(l_params_all), {
    v_mort <- rep((1 - exp(-v_r_mort_by_age[n_age_init:(n_age_max_full - 1), ] * (1 / per) * hr)), each = per) # weekly
    return(v_mort)
  })
}
