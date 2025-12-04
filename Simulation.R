############################################################
## Simulation study for MERS mortality risk modelling
## - Proposed method: logistic regression with splines
## - Comparator: misspecified classical logistic regression
## Author: Dr. Debashis Chatterjee
############################################################

## 0. Setup -------------------------------------------------

## Load required packages
library(tidyverse)
library(splines)      # for spline bases
library(pROC)         # AUC
# library(rms)        # no longer needed for simulation; using splines::ns instead
library(survival)     # (if later you want Cox models)
library(glmnet)       # if you want to add penalised models later

## Set seed for reproducibility
set.seed(12345)

## Create output directories
out_base   <- "results_sim"
plots_dir  <- file.path(out_base, "plots")
tables_dir <- file.path(out_base, "tables")

dir.create(out_base,   showWarnings = FALSE)
dir.create(plots_dir,  showWarnings = FALSE)
dir.create(tables_dir, showWarnings = FALSE)


############################################################
## 1. Data-generating mechanism (DGM)
############################################################

## Function to generate one MERS-like dataset
simulate_mers_like <- function(n = 200) {
  ## Age: roughly 25 to 85, skewed towards older
  age <- pmin(pmax(rnorm(n, mean = 60, sd = 15), 25), 85)
  
  ## Sex: 0 = female, 1 = male
  sex <- rbinom(n, size = 1, prob = 0.65)
  
  ## Place of infection: hospital vs community vs abroad
  place_levels <- c("hospital", "community", "abroad")
  place <- sample(place_levels, size = n, replace = TRUE,
                  prob = c(0.6, 0.3, 0.1))
  place <- factor(place, levels = place_levels)
  
  ## Onset-to-diagnosis delay (days): right-skewed, e.g. gamma
  delay <- rgamma(n, shape = 3, rate = 0.5)  # mean ~6 days
  delay <- pmin(delay, 20)                   # truncate extreme tails
  
  ## Network-like covariates: in-degree and out-degree of infector
  ## For simplicity we generate them as Poisson with small means.
  kin  <- rpois(n, lambda = 0.3)             # in-degree (number of known infectors)
  kout <- rpois(n, lambda = 1.0)             # out-degree of the infector (super-spreading)
  
  ## TRUE log-odds of death with non-linearities and interactions
  ## These choices create a scenario where flexible model is "correct"
  ## and the naive linear model is misspecified.
  
  ## Centre age & delay for nicer coefficients
  age_c   <- (age - 60) / 10        # age in decades centred at 60
  delay_c <- (delay - 6) / 4        # delay in units of 4 days centred at 6
  
  ## Non-linear effects
  f_age   <- 0.5 * age_c + 0.3 * (age_c^2)
  f_delay <- 0.4 * delay_c + 0.6 * pmax(delay_c - 0.5, 0)^2
  
  ## Effect of sex (male higher risk)
  beta_sex <- 0.6
  
  ## Effect of place of infection
  ## baseline: hospital; community and abroad have different risks
  beta_place <- c(
    hospital  = 0.0,   # reference
    community = -0.3,
    abroad    = 0.5
  )
  
  ## Network effects
  beta_kin  <- 0.25    # more known infectors = higher risk
  beta_kout <- 0.15    # infector’s out-degree (as proxy for super-spreaders)
  
  ## Intercept tuned to give ~20–30% mortality
  beta0 <- -1.0
  
  eta_true <- beta0 +
    f_age +
    f_delay +
    beta_sex * sex +
    beta_place[as.character(place)] +
    beta_kin  * kin +
    beta_kout * kout
  
  p_death <- 1 / (1 + exp(-eta_true))
  
  ## Outcome
  y <- rbinom(n, size = 1, prob = p_death)
  
  ## Return data frame with true probabilities (for reference if needed)
  tibble(
    y        = y,
    age      = age,
    sex      = factor(sex, levels = c(0, 1), labels = c("F", "M")),
    place    = place,
    delay    = delay,
    kin      = kin,
    kout     = kout,
    p_true   = p_death,
    eta_true = eta_true
  )
}


############################################################
## 2. Models
############################################################

## 2.1 Proposed model: Logistic regression with flexible splines
##     + key covariates, fitted via glm + splines::ns
fit_proposed_model <- function(dat) {
  ## Natural cubic splines with df = 4 for age and delay
  fit <- glm(
    y ~ ns(age, df = 4) + ns(delay, df = 4) +
      sex + place + kin + kout,
    data   = dat,
    family = binomial()
  )
  
  ## Predicted probabilities
  p_hat <- as.numeric(fitted(fit))
  
  list(
    fit   = fit,
    p_hat = p_hat
  )
}

## 2.2 Comparator (not-so-good) model:
##     Naive logistic regression with only linear age and sex
##     (no delay non-linearity, no place, no network)
fit_naive_model <- function(dat) {
  fit <- glm(
    y ~ age + sex,
    data   = dat,
    family = binomial()
  )
  p_hat <- as.numeric(fitted(fit))
  
  list(
    fit   = fit,
    p_hat = p_hat
  )
}


############################################################
## 3. Performance metrics
############################################################

compute_metrics <- function(y, p_hat) {
  ## AUC
  roc_obj <- pROC::roc(response = y, predictor = p_hat, quiet = TRUE)
  auc_val <- as.numeric(pROC::auc(roc_obj))
  
  ## Brier score
  brier <- mean((y - p_hat)^2)
  
  ## Calibration intercept & slope:
  ## Fit y ~ logit(p_hat) by logistic regression
  eps      <- 1e-6
  logit_ph <- log(p_hat + eps) - log(1 - p_hat + eps)
  calib_fit <- glm(y ~ logit_ph, family = binomial())
  
  calib_intercept <- coef(calib_fit)[1]
  calib_slope     <- coef(calib_fit)[2]
  
  tibble(
    auc              = auc_val,
    brier            = brier,
    calib_intercept  = calib_intercept,
    calib_slope      = calib_slope
  )
}


############################################################
## 4. One-shot demonstration (single dataset)
############################################################

## Generate one dataset and show model fits + calibration plot
n_demo   <- 200
dat_demo <- simulate_mers_like(n = n_demo)

## Fit proposed and naive models
mod_prop_demo  <- fit_proposed_model(dat_demo)
mod_naive_demo <- fit_naive_model(dat_demo)

## Compute metrics for demo
metrics_prop_demo  <- compute_metrics(dat_demo$y, mod_prop_demo$p_hat)
metrics_naive_demo <- compute_metrics(dat_demo$y, mod_naive_demo$p_hat)

print("=== Demo dataset: performance metrics ===")
print("Proposed model:")
print(metrics_prop_demo)

print("Naive model:")
print(metrics_naive_demo)

## Calibration plot for demo dataset
## Group patients into deciles of predicted risk (proposed model)
calib_data_demo <- tibble(
  y       = dat_demo$y,
  p_prop  = mod_prop_demo$p_hat,
  p_naive = mod_naive_demo$p_hat
) %>%
  mutate(
    decile_prop = ntile(p_prop, 10)
  ) %>%
  group_by(decile_prop) %>%
  summarise(
    p_hat_mean_prop = mean(p_prop),
    y_mean_prop     = mean(y),
    .groups = "drop"
  ) %>%
  rename(decile = decile_prop)

## Plot calibration for proposed model
p_calib_demo <- ggplot(calib_data_demo,
                       aes(x = p_hat_mean_prop, y = y_mean_prop)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(size = 2) +
  geom_line() +
  labs(
    x = "Mean predicted risk (proposed, deciles)",
    y = "Observed mortality proportion",
    title = "Calibration plot (demo dataset, proposed model)"
  ) +
  theme_minimal()

## Save and show the plot
calib_demo_file <- file.path(plots_dir, "calibration_demo_proposed.pdf")
ggsave(calib_demo_file, p_calib_demo, width = 6, height = 4)

print(p_calib_demo)  ## show.plt


############################################################
## 5. Simulation loop
############################################################

nsim      <- 500        # number of simulated datasets
n_per_sim <- 200        # sample size per dataset

## Storage for metrics
sim_results <- tibble(
  sim   = integer(),
  model = character(),
  auc   = double(),
  brier = double(),
  calib_intercept = double(),
  calib_slope     = double()
)

for (s in 1:nsim) {
  if (s %% 50 == 0) {
    cat("Simulation", s, "of", nsim, "\n")
  }
  
  dat <- simulate_mers_like(n = n_per_sim)
  
  ## Proposed model
  prop_res <- fit_proposed_model(dat)
  met_prop <- compute_metrics(dat$y, prop_res$p_hat) %>%
    mutate(sim = s, model = "Proposed")
  
  ## Naive model
  naive_res <- fit_naive_model(dat)
  met_naive <- compute_metrics(dat$y, naive_res$p_hat) %>%
    mutate(sim = s, model = "Naive")
  
  sim_results <- bind_rows(sim_results, met_prop, met_naive)
}

## Save raw simulation results
sim_results_file <- file.path(tables_dir, "simulation_results_raw.csv")
write.csv(sim_results, sim_results_file, row.names = FALSE)

print("=== Head of raw simulation results ===")
print(head(sim_results))


############################################################
## 6. Summarise simulation results (tables)
############################################################

## Summary by model
summary_table <- sim_results %>%
  group_by(model) %>%
  summarise(
    auc_mean   = mean(auc),
    auc_sd     = sd(auc),
    brier_mean = mean(brier),
    brier_sd   = sd(brier),
    calib_intercept_mean = mean(calib_intercept),
    calib_intercept_sd   = sd(calib_intercept),
    calib_slope_mean     = mean(calib_slope),
    calib_slope_sd       = sd(calib_slope),
    .groups = "drop"
  )

summary_table_file <- file.path(tables_dir, "simulation_summary_table.csv")
write.csv(summary_table, summary_table_file, row.names = FALSE)

print("=== Simulation summary table (to copy into paper) ===")
print(summary_table)


## Pairwise differences: Proposed - Naive
diff_table <- sim_results %>%
  select(sim, model, auc, brier, calib_slope) %>%
  pivot_wider(
    id_cols = sim,
    names_from = model,
    values_from = c(auc, brier, calib_slope)
  ) %>%
  mutate(
    diff_auc         = auc_Proposed - auc_Naive,
    diff_brier       = brier_Naive - brier_Proposed,   # positive => Proposed better
    diff_calib_slope = abs(calib_slope_Proposed - 1) - abs(calib_slope_Naive - 1)
  )

diff_summary <- diff_table %>%
  summarise(
    mean_diff_auc        = mean(diff_auc),
    sd_diff_auc          = sd(diff_auc),
    prop_auc_better      = mean(diff_auc > 0),
    mean_diff_brier      = mean(diff_brier),
    sd_diff_brier        = sd(diff_brier),
    prop_brier_better    = mean(diff_brier > 0),
    mean_diff_calib_slope = mean(diff_calib_slope),
    sd_diff_calib_slope   = sd(diff_calib_slope),
    prop_calib_slope_better = mean(diff_calib_slope < 0)
  )

diff_summary_file <- file.path(tables_dir, "simulation_diff_summary.csv")
write.csv(diff_summary, diff_summary_file, row.names = FALSE)

print("=== Differences (Proposed minus Naive) summary ===")
print(diff_summary)


############################################################
## 7. Plots comparing methods
############################################################

## Boxplots for AUC
p_auc_box <- sim_results %>%
  ggplot(aes(x = model, y = auc, fill = model)) +
  geom_boxplot(alpha = 0.7) +
  coord_cartesian(ylim = c(0.5, 1.0)) +
  labs(
    title = "Simulation: AUC distribution by model",
    x = "",
    y = "AUC"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

auc_plot_file <- file.path(plots_dir, "simulation_auc_boxplot.pdf")
ggsave(auc_plot_file, p_auc_box, width = 6, height = 4)

print(p_auc_box)   ## show.plt


## Boxplots for Brier score
p_brier_box <- sim_results %>%
  ggplot(aes(x = model, y = brier, fill = model)) +
  geom_boxplot(alpha = 0.7) +
  labs(
    title = "Simulation: Brier score distribution by model",
    x = "",
    y = "Brier score"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

brier_plot_file <- file.path(plots_dir, "simulation_brier_boxplot.pdf")
ggsave(brier_plot_file, p_brier_box, width = 6, height = 4)

print(p_brier_box)  ## show.plt


## Calibration slope distribution
p_calib_slope_box <- sim_results %>%
  ggplot(aes(x = model, y = calib_slope, fill = model)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_boxplot(alpha = 0.7) +
  labs(
    title = "Simulation: Calibration slope by model",
    x = "",
    y = "Calibration slope"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

calib_slope_plot_file <- file.path(plots_dir, "simulation_calibration_slope_boxplot.pdf")
ggsave(calib_slope_plot_file, p_calib_slope_box, width = 6, height = 4)

print(p_calib_slope_box)  ## show.plt


############################################################
## 8. Notes for manuscript
############################################################
## - Use `summary_table` in a main simulation results table.
## - Use `diff_summary` to argue that the proposed method
##   systematically improves AUC and Brier score, and yields
##   calibration slopes closer to 1, compared to the naive model.
## - Plots in `results_sim/plots` can be used as figures.
############################################################
