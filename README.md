# MERS-CoV-mortality-paper-frequentist

# Flexible Frequentist Risk Modelling of MERS-CoV Mortality (South Korea, 2015)

This repository contains reproducible **R code** for the paper:

> **Flexible Frequentist Risk Modelling of Case-Fatality in the 2015 South Korean MERS-CoV Outbreak:  
> A Spline-Based Logistic Framework with Internal Validation and Simulation Benchmarking**  
> Debashis Chatterjee, Department of Statistics, Visva–Bharati University.

The project implements a **fully frequentist**, spline-based logistic regression framework for modelling case-fatality risk, with:

- A **Monte Carlo simulation study** that mimics a MERS-like outbreak and compares:
  - a *proposed* flexible spline logistic model, vs.  
  - a *naïve* logit-linear model in age and sex.

- A **real-data analysis** of the 2015 South Korean MERS-CoV outbreak using the
  `mers_korea_2015` dataset from the `outbreaks` R package.

Both components reproduce the tables and figures used in the manuscript.

---

## Repository structure

A suggested (and typical) organisation is:

```text
.
├── R/
│   ├── 01_simulation.R                # Simulation DGM, model fitting, metrics, plots
│   └── 02_mers_real_data_analysis.R   # Real MERS-CoV analysis, figures, tables, zip output
├── results_sim/
│   ├── plots/                         # Simulation plots (AUC boxplot, Brier boxplot, calibration plot, etc.)
│   └── tables/                        # Simulation summary tables, raw metrics, paired differences
├── mers_real_results/
│   ├── plots/                         # Real-data plots (epi curve, age/delay distributions, ROC, calibration, partial effects)
│   ├── tables/                        # Real-data tables (outcome, age, delay, model performance, coefficients, calibration deciles)
│   └── models/                        # Saved fitted model objects (.rds)
├── mers_real_results.zip              # Zipped bundle of real-data output (created by the script)
├── README.md
└── (optionally) LICENSE
