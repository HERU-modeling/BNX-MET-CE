# BNX-MET-CE

[`BNX-MET-CE`](https://github.com/HERU-modeling/BNX-MET-CE) provides the necessary code and data to run a semi-Markov cohort model
and reproduce results for the comparative effectiveness of BNX versus methadone, presented in: "A decision-analytic model-based 
evaluation of the long-term comparative effectiveness of buprenorphine/naloxone versus methadone in British Columbia, Canada."

# R code modules

The R scripts to run the model and conduct the analysis are organized
into two sections:

1.  [R](https://github.com/HERU-modeling/BNX-MET-CE/tree/main/R)
    contains R functions for core modules.
    - You should not need to modify or run these scripts directly.

2.  [Analysis](https://github.com/HERU-modeling/BNX-MET-CE/tree/main/analysis)
    contains R scripts to run all analyses. Note that some scripts such as the calibration modules or PSA may take a long time to run (~12-24 hours or more), depending on computer speed and number of processors available.
    - Scripts to produce manuscript results
      - [07_run_cali_scenario_comp.R](https://github.com/HERU-modeling/BNX-MET-CE/blob/main/analysis/07_run_cali_scenario_comp.R) will produce the results for Figure 2.
      - [04_psa_run_model.R](https://github.com/HERU-modeling/BNX-MET-CE/blob/main/analysis/04_psa_run_model.R) and [04_psa_plot_results.R](https://github.com/HERU-modeling/BNX-MET-CE/blob/main/analysis/04_psa_plot_results.R) will produce the results for Figure 3.
      - [06_twsa_run_model.R](https://github.com/HERU-modeling/BNX-MET-CE/blob/main/analysis/06_twsa_run_model.R) and [06_twsa_plot_results.R](https://github.com/HERU-modeling/BNX-MET-CE/blob/main/analysis/06_twsa_plot_results.R) will produce the results for Figure 4.
    - Scripts to produce supplementary results
      - [01_calibration_itt.R](https://github.com/HERU-modeling/BNX-MET-CE/blob/main/analysis/01_calibration_itt.R) will produce the results for Figures S1, S10, and S11.
      - [01_calibration_pp.R](https://github.com/HERU-modeling/BNX-MET-CE/blob/main/analysis/01_calibration_pp.R) will produce the results for Figures S2, S12, and S13.
      - [02_trace_plots_itt.R](https://github.com/HERU-modeling/BNX-MET-CE/blob/main/analysis/02_trace_plots_itt.R) will produce the results for Figures S4, and S5.
      - [02_trace_plots_pp.R](https://github.com/HERU-modeling/BNX-MET-CE/blob/main/analysis/02_trace_plots_pp.R) will produce the results for Figures S6, and S7.
      - [05_owsa_run_model.R](https://github.com/HERU-modeling/BNX-MET-CE/blob/main/analysis/05_owsa_run_model.R) and [05_owsa_tornado.R](https://github.com/HERU-modeling/BNX-MET-CE/blob/main/analysis/05_owsa_tornado.R) will produce the results for Figures S8 and S9.
      - [08_cohort_balance.R](https://github.com/HERU-modeling/BNX-MET-CE/blob/main/analysis/08_cohort_balance.R) will produce the results for Figure S3.