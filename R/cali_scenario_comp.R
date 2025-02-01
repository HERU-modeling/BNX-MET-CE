#' cali_output_comp
#'
#' Description: This function calculates the calibration outcomes for different scenarios.
#'
#' @param scenario The scenario to calculate calibration outcomes for. Choose one of: "cali_itt", "cali_pp", "met_itt", "met_pp", "bnx_itt", "bnx_pp".
#' @return A list containing the calibration outcomes for the specified scenario.
#' @export
#'
cali_output_comp <- function(
    scenario,
    l_cali_targets) {
  # Initialize lists
  l_cali_outcomes_odf <- list()
  l_cali_outcomes_dno <- list()
  l_cali_outcomes_odn <- list()

  # Import calibration targets
  l_cali_targets <- l_cali_targets

  # Output #1: Calibration (ITT)
  if (scenario == "cali_itt") {
    l_params_all <- l_params_cali_itt
    m_params_updated <- m_calib_post_itt
    ce_est <- "cali"
    analytic_cohort <- "cali"
    # Output #2: Calibration (PP)
  } else if (scenario == "cali_pp") {
    l_params_all <- l_params_cali_pp
    m_params_updated <- m_calib_post_pp
    ce_est <- "cali"
    analytic_cohort <- "cali"
    # Output #3: MET (ITT)
  } else if (scenario == "met_itt") {
    l_params_all <- l_params_met_itt
    m_params_updated <- as.matrix(df_psa_params_itt)
    ce_est <- "itt_ps"
    analytic_cohort <- "met_only"
    # Output #4: MET (PP)
  } else if (scenario == "met_pp") {
    l_params_all <- l_params_met_pp
    m_params_updated <- as.matrix(df_psa_params_pp)
    ce_est <- "pp"
    analytic_cohort <- "met_only"
    # Output #5: BNX (ITT)
  } else if (scenario == "bnx_itt") {
    l_params_all <- l_params_bnx_itt
    m_params_updated <- as.matrix(df_psa_params_itt)
    ce_est <- "itt_ps"
    analytic_cohort <- "bnx_only"
    # Output #6: BNX (PP)
  } else if (scenario == "bnx_pp") {
    l_params_all <- l_params_bnx_pp
    m_params_updated <- as.matrix(df_psa_params_pp)
    ce_est <- "pp"
    analytic_cohort <- "bnx_only"
  } else {
    stop("Must select valid scenario")
  }

  for (j in (0:(n_blocks - 1))) {
    l_cali_target_fit <- foreach(i = (n_start + 1 + j * n_block_size):(n_start + (j + 1) * n_block_size), .combine = combine_custom_cali_output, .packages = "tidyr") %dopar% {
      l_model_target_fit <- calibration_out(
        v_params_updated = m_params_updated[i, ],
        l_params_all = l_params_all,
        ce_est = ce_est,
        analytic_cohort = analytic_cohort
      )
      m_model_targets_odf <- l_model_target_fit$fatal_overdose
      m_model_targets_dno <- l_model_target_fit$death_non_od
      m_model_targets_odn <- l_model_target_fit$overdose

      return(list(
        m_model_targets_odf = m_model_targets_odf,
        m_model_targets_dno = m_model_targets_dno,
        m_model_targets_odn = m_model_targets_odn
      ))
    }

    m_cali_outcomes_odf <- l_cali_target_fit$m_model_targets_odf
    m_cali_outcomes_dno <- l_cali_target_fit$m_model_targets_dno
    m_cali_outcomes_odn <- l_cali_target_fit$m_model_targets_odn

    l_cali_outcomes_odf[[j + 1]] <- m_cali_outcomes_odf
    l_cali_outcomes_dno[[j + 1]] <- m_cali_outcomes_dno
    l_cali_outcomes_odn[[j + 1]] <- m_cali_outcomes_odn

    out <- paste0("Block ", (j + 1), " of ", n_blocks, " complete.") # Status
    print(out)
  }

  m_outcomes_odf_comb <- m_outcomes_dno_comb <- m_outcomes_odn_comb <- matrix(0, nrow = 0, ncol = 9)

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

  return(list(
    m_outcomes_odf_comb = m_outcomes_odf_comb,
    m_outcomes_dno_comb = m_outcomes_dno_comb,
    m_outcomes_odn_comb = m_outcomes_odn_comb
  ))
}

# Function to combine parallel outputs
combine_custom_cali_output <- function(LL1, LL2) {
  m_model_targets_odf <- rbind(LL1$m_model_targets_odf, LL2$m_model_targets_odf)
  m_model_targets_dno <- rbind(LL1$m_model_targets_dno, LL2$m_model_targets_dno)
  m_model_targets_odn <- rbind(LL1$m_model_targets_odn, LL2$m_model_targets_odn)

  return(list(
    m_model_targets_odf = m_model_targets_odf,
    m_model_targets_dno = m_model_targets_dno,
    m_model_targets_odn = m_model_targets_odn
  ))
}
