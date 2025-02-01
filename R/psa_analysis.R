#' PSA analysis
#'
#' \code{psa_analysis} implements PSA analysis for each analytic scenario
#'
#' @param scenario analytic scenarios ("ITT_PS"; "ITT_IV"; "PP")
#' @return
#' df_outcomes_MET_PSA_comb:
#' df_outcomes_BNX_PSA_comb:
#' df_incremental_PSA_comb:
#' @export
psa_analysis <- function(ce_est) { # choose one of three analytic scenarios ("ITT_PS"; "ITT_IV"; "PP")
  # Initialize lists
  l_outcomes_met_psa <- list()
  l_outcomes_bnx_psa <- list()
  l_incremental_psa <- list()
  l_incremental_psa_scaled <- list()

  if (ce_est == "itt_ps") {
    l_params_met <- l_params_met_itt
    l_params_bnx <- l_params_bnx_itt
    df_psa_params_met <- df_psa_params_itt
    df_psa_params_bnx <- df_psa_params_itt
  } else if (ce_est == "pp_hd") {
    l_params_met <- l_params_met_pp
    l_params_bnx <- l_params_bnx_pp
    df_psa_params_met <- df_psa_params_pp
    df_psa_params_bnx <- df_psa_params_pp
  } else if (ce_est == "pp") {
    l_params_met <- l_params_met_pp
    l_params_bnx <- l_params_bnx_pp
    df_psa_params_met <- df_psa_params_pp
    df_psa_params_bnx <- df_psa_params_pp
  } else {
    stop("Must select valid scenario")
  }

  for (j in (0:(n_blocks - 1))) {
    l_psa <- foreach(i = (n_start + 1 + j * n_block_size):(n_start + (j + 1) * n_block_size), .combine = combine_custom, .packages = "tidyr") %dopar% {
      # Update parameter set for each scenario with next set of PSA drawn parameters (i rows in PSA parameter data frame)
      l_psa_input_met <- update_param_list(l_params_all = l_params_met, params_updated = df_psa_params_met[i, ])
      l_psa_input_bnx <- update_param_list(l_params_all = l_params_bnx, params_updated = df_psa_params_bnx[i, ])

      # Run model and generate outputs
      l_outcomes_met <- outcomes(l_params_all = l_psa_input_met, v_params_calib = NULL, psa = TRUE, time_horizon = "full", ce_est = ce_est, analytic_cohort = "met_only")
      l_outcomes_bnx <- outcomes(l_params_all = l_psa_input_bnx, v_params_calib = NULL, psa = TRUE, time_horizon = "full", ce_est = ce_est, analytic_cohort = "bnx_only")
      df_outcomes_met_psa <- l_outcomes_met$df_outcomes
      df_outcomes_bnx_psa <- l_outcomes_bnx$df_outcomes

      # Calculate incremental outcomes
      l_inc_outcomes <- inc_outcomes(outcomes_comp = l_outcomes_met, outcomes_int = l_outcomes_bnx)
      df_incremental_psa <- l_inc_outcomes$df_incremental
      df_incremental_psa_scaled <- l_inc_outcomes$df_incremental_scaled

      return(list(
        df_outcomes_met_psa = df_outcomes_met_psa,
        df_outcomes_bnx_psa = df_outcomes_bnx_psa,
        df_incremental_psa = df_incremental_psa,
        df_incremental_psa_scaled = df_incremental_psa_scaled
      ))
    }

    df_outcomes_met_psa <- l_psa$df_outcomes_met_psa
    df_outcomes_bnx_psa <- l_psa$df_outcomes_bnx_psa
    df_incremental_psa <- l_psa$df_incremental_psa
    df_incremental_psa_scaled <- l_psa$df_incremental_psa_scaled

    l_outcomes_met_psa[[j + 1]] <- df_outcomes_met_psa
    l_outcomes_bnx_psa[[j + 1]] <- df_outcomes_bnx_psa
    l_incremental_psa[[j + 1]] <- df_incremental_psa
    l_incremental_psa_scaled[[j + 1]] <- df_incremental_psa_scaled

    out <- paste0("Block ", (j + 1), " of ", n_blocks, " complete.") # Status
    print(out)
  }

  # combine output data sets
  # initialize empty data frames
  df_outcomes_met_psa_comb <- df_outcomes_bnx_psa_comb <- df_incremental_psa_comb <- df_incremental_psa_scaled_comb <- data.frame()

  for (i in 1:n_blocks) {
    df_temp <- l_outcomes_met_psa[[i]]
    df_outcomes_met_psa_comb <- rbind(df_outcomes_met_psa_comb, df_temp)
  }
  for (i in 1:n_blocks) {
    df_temp <- l_outcomes_bnx_psa[[i]]
    df_outcomes_bnx_psa_comb <- rbind(df_outcomes_bnx_psa_comb, df_temp)
  }
  for (i in 1:n_blocks) {
    df_temp <- l_incremental_psa[[i]]
    df_incremental_psa_comb <- rbind(df_incremental_psa_comb, df_temp)
  }
  for (i in 1:n_blocks) {
    df_temp <- l_incremental_psa_scaled[[i]]
    df_incremental_psa_scaled_comb <- rbind(df_incremental_psa_scaled_comb, df_temp)
  }


  return(list(
    df_outcomes_met_psa_comb = df_outcomes_met_psa_comb,
    df_outcomes_bnx_psa_comb = df_outcomes_bnx_psa_comb,
    df_incremental_psa_comb = df_incremental_psa_comb,
    df_incremental_psa_scaled_comb = df_incremental_psa_scaled_comb
  ))
}

# Function to combine parallel outputs
combine_custom <- function(LL1, LL2) {
  df_outcomes_met_psa <- rbind(LL1$df_outcomes_met_psa, LL2$df_outcomes_met_psa)
  df_outcomes_bnx_psa <- rbind(LL1$df_outcomes_bnx_psa, LL2$df_outcomes_bnx_psa)
  df_incremental_psa <- rbind(LL1$df_incremental_psa, LL2$df_incremental_psa)
  df_incremental_psa_scaled <- rbind(LL1$df_incremental_psa_scaled, LL2$df_incremental_psa_scaled)

  return(list(
    df_outcomes_met_psa = df_outcomes_met_psa,
    df_outcomes_bnx_psa = df_outcomes_bnx_psa,
    df_incremental_psa = df_incremental_psa,
    df_incremental_psa_scaled = df_incremental_psa_scaled
  ))
}
