#' Trace plots
#'
#' \code{trace_plots} implements function to plot Markov trace.
#'
#' @param outcomes List with model outcomes
#' @return
#' main_states_trace_plot
#' sero_states_trace_plot
#' main_states_time
#' plots
#' layout
#' @export

trace_plots <- function(
    outcomes,
    analysis) {
  ## Prepare data
  if (analysis == "cali") {
    ## Calibration
    df_m_agg_trace_cali <- as.data.frame(outcomes$m_m_agg_trace_cali)
    df_m_agg_trace_cali$time <- as.numeric(rownames(df_m_agg_trace_cali))
    df_m_agg_trace_cali$year <- df_m_agg_trace_cali$time / 52
    df_m_agg_trace_plot_cali <- df_m_agg_trace_cali %>% gather(state, proportion, "Non-Overdose Death", "ODN", "ODF", "OU", "OAT", "ABS") # health states to plot
    df_m_agg_state_time_cali <- df_m_agg_trace_cali %>% gather(state, proportion, "ODN", "OU", "OAT", "ABS") # alive health states to plot
    df_m_agg_state_time_cali <- df_m_agg_state_time_cali %>%
      group_by(state) %>%
      dplyr::summarise_at(.vars = vars(proportion), .funs = list(sum)) %>%
      mutate(
        percentage = round((proportion / sum(proportion)) * 100, 1),
        time_alive = sum(proportion)
      )
    # Preserve order for plotting
    state_order_trace_cali <- factor(df_m_agg_trace_plot_cali$state, levels = c("Non-Overdose Death", "ODF", "ODN", "OU", "OAT", "ABS"))
    state_order_time_cali <- factor(df_m_agg_state_time_cali$state, levels = c("ODN", "OU", "OAT", "ABS"))
  } else if (analysis == "full") {
    ## Full model
    df_m_agg_trace <- as.data.frame(outcomes$m_m_agg_trace)
    df_m_agg_trace$time <- as.numeric(rownames(df_m_agg_trace))
    df_m_agg_trace$year <- df_m_agg_trace$time / 52
    df_m_agg_trace_plot <- df_m_agg_trace %>% gather(state, proportion, "Non-Overdose Death", "ODN", "ODF", "OUM", "OUB", "OUO", "BNX", "MET", "ABS") # health states to plot
    df_m_agg_state_time <- df_m_agg_trace %>% gather(state, proportion, "ODN", "OUM", "OUB", "OUO", "BNX", "MET", "ABS") # alive health states to plot
    df_m_agg_state_time <- df_m_agg_state_time %>%
      group_by(state) %>%
      dplyr::summarise_at(.vars = vars(proportion), .funs = list(sum)) %>%
      mutate(
        percentage = round((proportion / sum(proportion)) * 100, 1),
        time_alive = sum(proportion)
      )
    # Preserve order for plotting
    state_order_trace <- factor(df_m_agg_trace_plot$state, levels = c("Non-Overdose Death", "ODF", "ODN", "OUM", "OUB", "OUO", "BNX", "MET", "ABS"))
    state_order_time <- factor(df_m_agg_state_time$state, levels = c("ODN", "OUM", "OUB", "OUO", "BNX", "MET", "ABS"))
  } else if (analysis == "cohort_balance") {
    ## Cohort balance
    # FIXME: need to adjust proportions for alive
    # TODO: turn into line plot with dots for observed data
    df_m_agg_cohort_balance_trace <- as.data.frame(outcomes$m_m_cohort_balance_trace)
    df_m_agg_cohort_balance_trace$time <- as.numeric(rownames(df_m_agg_cohort_balance_trace))
    df_m_agg_cohort_balance_trace$year <- df_m_agg_cohort_balance_trace$time / 52
    df_m_agg_cohort_balance_trace_plot <- df_m_agg_cohort_balance_trace %>% gather(state, proportion, "Incident user", "Experienced user") # health states to plot
    # Preserve order for plotting
    state_order_trace_cohort_balance <- factor(df_m_agg_cohort_balance_trace_plot$state, levels = c("Incident user", "Experienced user"))
  } else {
    print("No analysis selected")
  }

  # Palette #3                Death       ODF        ODN                OUM          OUB           OUO         BNX            MET         ABS
  state_colours_trace3 <- c("grey50", "#d73027", "#fdae61", "#74add1", "#ffffb3", "#a8ddb5", "#FFB612", "#003087", "#bebada") # colour pallette 3
  #                                 Death         ODF          ODN          OU           OAT          ABS
  state_colours_trace3_cali <- c("grey50", "#d73027", "#fdae61", "#a8ddb5", "#046A38", "#bebada") # colour pallette 3
  #                            ODN           OUM          OUB           OUO         BNX            MET         ABS
  state_colours_time3 <- c("#fdae61", "#74add1", "#ffffb3", "#a8ddb5", "#FFB612", "#003087", "#bebada") # colour pallette 3
  #                            ODN            OU          OAT          ABS
  state_colours_time3_cali <- c("#fdae61", "#a8ddb5", "#046A38", "#bebada") # colour pallette 3

  #                             Incident user   Experienced user
  state_colours_cohort_balance <- c("#d73027", "#4575b4") # colour pallette 3

  # Labels
  # Health state occupancy
  state_trace_label <- c("Non-Overdose Death", "Fatal overdose", "Non-fatal overdose", "Out-of-tx (MET)", "Out-of-tx (BNX)", "Out-of-tx", "BNX", "MET", "Abstinence")
  state_trace_label_cali <- c("Non-Overdose Death", "Fatal overdose", "Non-fatal overdose", "Out-of-tx", "OAT", "Abstinence")
  state_trace_label_cohort_balance <- c("Incident user", "Experienced user")
  # Time spent in health states
  state_time_label <- c("ODN", "OOT-M", "OOT-B", "OOT", "BNX", "MET", "ABS")
  state_time_label_cali <- c("ODN", "OOT", "OAT", "ABS")

  ### Markov trace plots ###
  # Base model states
  if (analysis == "cali") {
    main_states_trace_plot <- ggplot(df_m_agg_trace_plot_cali, aes(x = year, y = proportion, fill = state_order_trace_cali)) +
      theme_bw() +
      theme(legend.position = "bottom") +
      xlab("Year") +
      ylab("Proportion in state") +
      geom_area() +
      scale_fill_manual(
        name = "Health States",
        labels = state_trace_label_cali,
        values = state_colours_trace3_cali
      ) +
      scale_x_continuous(
        breaks = seq(from = 0, to = 8, by = 1),
        labels = c("2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020*")
      )
  } else if (analysis == "full") {
    main_states_trace_plot <- ggplot(df_m_agg_trace_plot, aes(x = year, y = proportion, fill = state_order_trace)) +
      theme_bw() +
      theme(legend.position = "bottom") +
      xlab("Year") +
      ylab("Proportion in state") +
      geom_area() +
      scale_fill_manual(
        name = "Health States",
        labels = state_trace_label,
        values = state_colours_trace3
      ) +
      scale_x_continuous(
        breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
        labels = c("2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020*")
      )
  } else if (analysis == "cohort_balance") {
    main_states_trace_plot <- ggplot(df_m_agg_cohort_balance_trace_plot, aes(x = year, y = proportion, fill = state_order_trace_cohort_balance)) +
      theme_bw() +
      theme(legend.position = "bottom") +
      xlab("Year") +
      ylab("Proportion in state") +
      geom_area() +
      scale_fill_manual(
        name = "Health States",
        labels = state_trace_label_cohort_balance,
        values = state_colours_cohort_balance
      ) +
      scale_x_continuous(
        breaks = seq(from = 0, to = 10, by = 1),
        labels = c("2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020*")
      )
  } else {
    print("No analysis selected")
  }

  ### Time spent in health states ###
  if (analysis == "cali") {
    main_states_time <- ggplot(df_m_agg_state_time_cali, aes(x = state_order_time_cali, y = proportion, fill = state_order_time_cali)) +
      theme_bw() +
      theme(legend.position = "none") +
      xlab("Health State") +
      ylab("Time") +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = state_colours_time3_cali) +
      scale_x_discrete(labels = state_time_label_cali) +
      geom_text(aes(label = paste0(round(proportion, 1), " (", percentage, "%)")), hjust = -0.25, size = 3.5) +
      annotate("text", x = 1.25, y = 375, label = paste0(round((df_m_agg_state_time_cali$time_alive) / 52, 2), " years alive"), size = 3.5) +
      coord_flip(ylim = c(0, 400))
  } else if (analysis == "full") {
    main_states_time <- ggplot(df_m_agg_state_time, aes(x = state_order_time, y = proportion, fill = state_order_time)) +
      theme_bw() +
      theme(legend.position = "none") +
      xlab("Health State") +
      ylab("Time") +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = state_colours_time3) +
      scale_x_discrete(labels = state_time_label) +
      geom_text(aes(label = paste0(round(proportion, 1), " (", percentage, "%)")), hjust = -0.25, size = 3.5) +
      annotate("text", x = 1.25, y = 375, label = paste0(round((df_m_agg_state_time$time_alive) / 52, 2), " years alive"), size = 3.5) +
      coord_flip(ylim = c(0, 400))
  } else if (analysis == "cohort_balance") {
    main_states_time <- NULL
  } else {
    print("No analysis selected")
  }


  ### Combined plot ###
  # if (exists("main_states_time")) {
  #   plots <- list()
  #   plots[[1]] <- main_states_trace_plot
  #   plots[[2]] <- main_states_time
  #   layout <- matrix(c(1, 1, 2), nrow = 3, byrow = TRUE)
  # } else {
  #   plots <- list()
  #   plots[[1]] <- main_states_trace_plot
  #   layout <- matrix(c(1, 1), nrow = 2, byrow = TRUE)
  # }

  # plots <- list()
  # plots[[1]] <- main_states_trace_plot
  # plots[[2]] <- main_states_time
  # layout <- matrix(c(1, 1, 2), nrow = 3, byrow = TRUE)

  return(list(
    main_states_trace_plot,
    main_states_time # ,
    # plots,
    # layout
  ))
}
