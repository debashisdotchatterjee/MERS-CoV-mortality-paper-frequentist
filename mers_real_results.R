############################################################
## Real-data analysis: MERS Korea 2015
## - Apply proposed spline logistic model vs naive baseline
## - Generate descriptive tables and colorful plots
## - Save outputs + zip folder
############################################################

## 0. Setup -------------------------------------------------

## Load required packages
library(tidyverse)
library(outbreaks)
library(splines)      # for spline bases
library(pROC)         # for AUC / ROC
library(broom)        # for tidy model summaries
# library(janitor)    # optional for nice tables
# library(patchwork)  # optional if you want plot grids

## Set seed (not strictly needed, but good practice when sampling)
set.seed(12345)

## Output directories
out_base    <- "mers_real_results"
plots_dir   <- file.path(out_base, "plots")
tables_dir  <- file.path(out_base, "tables")
models_dir  <- file.path(out_base, "models")

dir.create(out_base,   showWarnings = FALSE)
dir.create(plots_dir,  showWarnings = FALSE)
dir.create(tables_dir, showWarnings = FALSE)
dir.create(models_dir, showWarnings = FALSE)


############################################################
## 1. Load data and prepare analysis dataset
############################################################

## Load dataset
data("mers_korea_2015")
linelist <- mers_korea_2015$linelist
contacts <- mers_korea_2015$contacts

## Quick check (optional)
# str(linelist)
# str(contacts)

## Outcome: Death vs Alive
## outcome is factor with levels "Alive", "Dead"
linelist <- linelist %>%
  mutate(
    y        = if_else(outcome == "Dead", 1L, 0L),
    y_factor = factor(outcome, levels = c("Alive", "Dead"))
  )

## Delays (in days)
linelist <- linelist %>%
  mutate(
    delay_diag   = as.numeric(dt_diag   - dt_onset),
    delay_report = as.numeric(dt_report - dt_onset)
  )

## Remove clearly invalid delays (e.g. negative)
linelist <- linelist %>%
  mutate(
    delay_diag   = if_else(delay_diag   < 0, NA_real_, delay_diag),
    delay_report = if_else(delay_report < 0, NA_real_, delay_report)
  )

## 1.1 Network covariates from contacts ---------------------

## contacts typically has columns: from, to, exposure, diff_dt_onset
# colnames(contacts)

## In-degree: number of times ID appears as 'to'
deg_in <- contacts %>%
  count(to, name = "kin")

## Out-degree: number of times ID appears as 'from'
deg_out <- contacts %>%
  count(from, name = "kout")

## Join network features back to line list via id
linelist_net <- linelist %>%
  left_join(deg_in,  by = c("id" = "to")) %>%
  left_join(deg_out, by = c("id" = "from")) %>%
  mutate(
    kin  = replace_na(kin,  0L),
    kout = replace_na(kout, 0L)
  )

## 1.2 Final analysis dataset -------------------------------

dat_full <- linelist_net %>%
  mutate(
    sex          = sex,
    place_infect = place_infect,
    hospital     = loc_hosp
  ) %>%
  select(
    id,
    y, y_factor,
    age,
    sex,
    place_infect,
    hospital,
    dt_onset,
    dt_report,
    week_report,
    delay_diag,
    delay_report,
    kin, kout
  )

## Keep complete cases for the proposed model’s key covariates
dat_cc <- dat_full %>%
  filter(
    !is.na(y),
    !is.na(age),
    !is.na(sex),
    !is.na(place_infect),
    !is.na(delay_diag)
  )

## Quick check sample size
nrow(dat_cc)


############################################################
## 2. Descriptive epidemiology
############################################################

## 2.1 Outcome distribution
table_outcome <- dat_cc %>%
  count(y_factor, name = "n") %>%
  mutate(prop = n / sum(n))

## Save + show
write.csv(table_outcome,
          file.path(tables_dir, "table_outcome_distribution.csv"),
          row.names = FALSE)
print("=== Outcome distribution (Alive vs Dead) ===")
print(table_outcome)


## 2.2 Age distribution (overall and by outcome)
table_age_summary <- dat_cc %>%
  group_by(y_factor) %>%
  summarise(
    n         = n(),
    age_mean  = mean(age, na.rm = TRUE),
    age_sd    = sd(age,   na.rm = TRUE),
    age_median = median(age, na.rm = TRUE),
    age_min   = min(age,  na.rm = TRUE),
    age_max   = max(age,  na.rm = TRUE),
    .groups   = "drop"
  )

write.csv(table_age_summary,
          file.path(tables_dir, "table_age_by_outcome.csv"),
          row.names = FALSE)
print("=== Age summary by outcome ===")
print(table_age_summary)


## 2.3 Delay to diagnosis by outcome
table_delay_summary <- dat_cc %>%
  group_by(y_factor) %>%
  summarise(
    n             = n(),
    delay_mean    = mean(delay_diag, na.rm = TRUE),
    delay_sd      = sd(delay_diag,   na.rm = TRUE),
    delay_median  = median(delay_diag, na.rm = TRUE),
    delay_min     = min(delay_diag,  na.rm = TRUE),
    delay_max     = max(delay_diag,  na.rm = TRUE),
    .groups       = "drop"
  )

write.csv(table_delay_summary,
          file.path(tables_dir, "table_delay_by_outcome.csv"),
          row.names = FALSE)
print("=== Onset-to-diagnosis delay summary by outcome ===")
print(table_delay_summary)


## 2.4 Place of infection by outcome
table_place_outcome <- dat_cc %>%
  count(place_infect, y_factor, name = "n") %>%
  group_by(y_factor) %>%
  mutate(prop_within_outcome = n / sum(n)) %>%
  ungroup()

write.csv(table_place_outcome,
          file.path(tables_dir, "table_place_by_outcome.csv"),
          row.names = FALSE)
print("=== Place of infection by outcome ===")
print(table_place_outcome)


## 2.5 Epi curve by onset date, colored by outcome ----------
p_epi_curve <- ggplot(dat_cc, aes(x = dt_onset, fill = y_factor)) +
  geom_histogram(binwidth = 1, colour = "black", alpha = 0.8) +
  scale_fill_brewer(palette = "Set1", name = "Outcome") +
  labs(
    x = "Date of symptom onset",
    y = "Number of cases",
    title = "Epidemic curve by onset date and outcome"
  ) +
  theme_minimal()

epi_curve_file <- file.path(plots_dir, "epi_curve_onset_by_outcome.pdf")
ggsave(epi_curve_file, p_epi_curve, width = 7, height = 4.5)
print(p_epi_curve)  ## show.plt


## 2.6 Age distribution by outcome (hist + density) --------
p_age_hist <- ggplot(dat_cc, aes(x = age, fill = y_factor)) +
  geom_histogram(alpha = 0.7, position = "identity", bins = 20, colour = "black") +
  scale_fill_brewer(palette = "Dark2", name = "Outcome") +
  labs(
    x = "Age (years)",
    y = "Count",
    title = "Age distribution by outcome"
  ) +
  theme_minimal()

age_hist_file <- file.path(plots_dir, "age_hist_by_outcome.pdf")
ggsave(age_hist_file, p_age_hist, width = 6, height = 4)
print(p_age_hist)  ## show.plt

p_age_density <- ggplot(dat_cc, aes(x = age, colour = y_factor, fill = y_factor)) +
  geom_density(alpha = 0.3) +
  scale_colour_brewer(palette = "Set1", name = "Outcome") +
  scale_fill_brewer(palette = "Set1", name = "Outcome") +
  labs(
    x = "Age (years)",
    y = "Density",
    title = "Age density by outcome"
  ) +
  theme_minimal()

age_density_file <- file.path(plots_dir, "age_density_by_outcome.pdf")
ggsave(age_density_file, p_age_density, width = 6, height = 4)
print(p_age_density)  ## show.plt


## 2.7 Delay distribution by outcome ------------------------
p_delay_hist <- ggplot(dat_cc, aes(x = delay_diag, fill = y_factor)) +
  geom_histogram(alpha = 0.7, position = "identity", bins = 15, colour = "black") +
  scale_fill_brewer(palette = "Paired", name = "Outcome") +
  labs(
    x = "Onset-to-diagnosis delay (days)",
    y = "Count",
    title = "Delay to diagnosis by outcome"
  ) +
  theme_minimal()

delay_hist_file <- file.path(plots_dir, "delay_diag_hist_by_outcome.pdf")
ggsave(delay_hist_file, p_delay_hist, width = 6, height = 4)
print(p_delay_hist)  ## show.plt

p_delay_box <- ggplot(dat_cc, aes(x = y_factor, y = delay_diag, fill = y_factor)) +
  geom_boxplot(alpha = 0.8) +
  scale_fill_brewer(palette = "Set2", name = "Outcome") +
  labs(
    x = "Outcome",
    y = "Onset-to-diagnosis delay (days)",
    title = "Delay to diagnosis by outcome"
  ) +
  theme_minimal()

delay_box_file <- file.path(plots_dir, "delay_diag_box_by_outcome.pdf")
ggsave(delay_box_file, p_delay_box, width = 5, height = 4)
print(p_delay_box)  ## show.plt


## 2.8 Age vs delay scatter, colored by outcome -------------
p_age_delay <- ggplot(dat_cc,
                      aes(x = age, y = delay_diag, colour = y_factor)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_colour_brewer(palette = "Set1", name = "Outcome") +
  labs(
    x = "Age (years)",
    y = "Onset-to-diagnosis delay (days)",
    title = "Age vs delay to diagnosis, by outcome"
  ) +
  theme_minimal()

age_delay_file <- file.path(plots_dir, "age_vs_delay_by_outcome.pdf")
ggsave(age_delay_file, p_age_delay, width = 6, height = 4)
print(p_age_delay)  ## show.plt


############################################################
## 3. Fit proposed and naive logistic regression models
############################################################

## Proposed model: splines + full covariates
mod_prop <- glm(
  y ~ ns(age, df = 4) + ns(delay_diag, df = 4) +
    sex + place_infect + kin + kout,
  data   = dat_cc,
  family = binomial()
)

## Naive model: linear age + sex only
mod_naive <- glm(
  y ~ age + sex,
  data   = dat_cc,
  family = binomial()
)

## Save models
saveRDS(mod_prop,  file.path(models_dir, "mod_proposed_logit.rds"))
saveRDS(mod_naive, file.path(models_dir, "mod_naive_logit.rds"))


############################################################
## 4. Model summaries (tables)
############################################################

## 4.1 Coefficients (tidy)
tidy_prop  <- broom::tidy(mod_prop) %>%
  mutate(OR = exp(estimate),
         OR_low = exp(estimate - 1.96 * std.error),
         OR_high = exp(estimate + 1.96 * std.error))

tidy_naive <- broom::tidy(mod_naive) %>%
  mutate(OR = exp(estimate),
         OR_low = exp(estimate - 1.96 * std.error),
         OR_high = exp(estimate + 1.96 * std.error))

write.csv(tidy_prop,
          file.path(tables_dir, "model_proposed_coefficients.csv"),
          row.names = FALSE)
write.csv(tidy_naive,
          file.path(tables_dir, "model_naive_coefficients.csv"),
          row.names = FALSE)

print("=== Proposed model coefficients (logit scale + OR) ===")
print(tidy_prop)

print("=== Naive model coefficients (logit scale + OR) ===")
print(tidy_naive)


############################################################
## 5. Model performance: AUC, Brier, calibration
############################################################

## Add fitted probabilities to dataset
dat_cc <- dat_cc %>%
  mutate(
    p_prop  = predict(mod_prop,  type = "response"),
    p_naive = predict(mod_naive, type = "response")
  )

## Helper: compute metrics
compute_metrics <- function(y, p_hat) {
  roc_obj <- pROC::roc(response = y, predictor = p_hat, quiet = TRUE)
  auc_val <- as.numeric(pROC::auc(roc_obj))
  brier   <- mean((y - p_hat)^2)
  eps      <- 1e-6
  logit_ph <- log(p_hat + eps) - log(1 - p_hat + eps)
  calib_fit <- glm(y ~ logit_ph, family = binomial())
  calib_intercept <- coef(calib_fit)[1]
  calib_slope     <- coef(calib_fit)[2]
  tibble(
    auc             = auc_val,
    brier           = brier,
    calib_intercept = calib_intercept,
    calib_slope     = calib_slope
  )
}

metrics_prop  <- compute_metrics(dat_cc$y, dat_cc$p_prop)  %>%
  mutate(model = "Proposed")
metrics_naive <- compute_metrics(dat_cc$y, dat_cc$p_naive) %>%
  mutate(model = "Naive")

model_perf <- bind_rows(metrics_prop, metrics_naive) %>%
  select(model, everything())

write.csv(model_perf,
          file.path(tables_dir, "model_performance_real_data.csv"),
          row.names = FALSE)

print("=== Model performance on real MERS data ===")
print(model_perf)


############################################################
## 6. Calibration plots (deciles)
############################################################

## Decile calibration for both models
calib_deciles <- dat_cc %>%
  transmute(
    y,
    p_prop,
    p_naive
  ) %>%
  mutate(
    decile_prop  = ntile(p_prop,  10),
    decile_naive = ntile(p_naive, 10)
  )

calib_prop <- calib_deciles %>%
  group_by(decile_prop) %>%
  summarise(
    mean_p = mean(p_prop),
    obs    = mean(y),
    .groups = "drop"
  ) %>%
  mutate(model = "Proposed")

calib_naive <- calib_deciles %>%
  group_by(decile_naive) %>%
  summarise(
    mean_p = mean(p_naive),
    obs    = mean(y),
    .groups = "drop"
  ) %>%
  mutate(model = "Naive")

calib_both <- bind_rows(
  calib_prop %>% rename(decile = decile_prop),
  calib_naive %>% rename(decile = decile_naive)
)

write.csv(calib_both,
          file.path(tables_dir, "calibration_deciles_real_data.csv"),
          row.names = FALSE)

## Plot calibration curves
p_calib_both <- ggplot(calib_both,
                       aes(x = mean_p, y = obs, colour = model)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(size = 2) +
  geom_line() +
  scale_colour_brewer(palette = "Set1", name = "Model") +
  labs(
    x = "Mean predicted risk (by decile)",
    y = "Observed mortality proportion",
    title = "Calibration curves (deciles) for proposed vs naive models"
  ) +
  theme_minimal()

calib_both_file <- file.path(plots_dir, "calibration_deciles_proposed_vs_naive.pdf")
ggsave(calib_both_file, p_calib_both, width = 6, height = 4)
print(p_calib_both)  ## show.plt


############################################################
## 7. ROC curves plot
############################################################

roc_prop  <- pROC::roc(response = dat_cc$y, predictor = dat_cc$p_prop,  quiet = TRUE)
roc_naive <- pROC::roc(response = dat_cc$y, predictor = dat_cc$p_naive, quiet = TRUE)

roc_prop_df <- tibble(
  fpr = 1 - roc_prop$specificities,
  tpr = roc_prop$sensitivities,
  model = "Proposed"
)
roc_naive_df <- tibble(
  fpr = 1 - roc_naive$specificities,
  tpr = roc_naive$sensitivities,
  model = "Naive"
)

roc_both_df <- bind_rows(roc_prop_df, roc_naive_df)

p_roc_both <- ggplot(roc_both_df,
                     aes(x = fpr, y = tpr, colour = model)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_line(size = 1.2) +
  scale_colour_brewer(palette = "Set1", name = "Model") +
  coord_equal() +
  labs(
    x = "1 - Specificity (FPR)",
    y = "Sensitivity (TPR)",
    title = "ROC curves: proposed vs naive"
  ) +
  theme_minimal()

roc_both_file <- file.path(plots_dir, "roc_proposed_vs_naive.pdf")
ggsave(roc_both_file, p_roc_both, width = 6, height = 5)
print(p_roc_both)  ## show.plt


############################################################
## 8. Partial-effect plots (age and delay)
############################################################

## Create a grid for age at fixed sex / place / delay
age_seq <- seq(min(dat_cc$age, na.rm = TRUE),
               max(dat_cc$age, na.rm = TRUE),
               length.out = 100)

## Choose reference levels: sex = most frequent, place_infect = most frequent
ref_sex   <- dat_cc %>% count(sex)          %>% slice_max(n, n = 1) %>% pull(sex)
ref_place <- dat_cc %>% count(place_infect) %>% slice_max(n, n = 1) %>% pull(place_infect)

## Use median delay, kin, kout as typical values
ref_delay <- median(dat_cc$delay_diag, na.rm = TRUE)
ref_kin   <- median(dat_cc$kin,         na.rm = TRUE)
ref_kout  <- median(dat_cc$kout,        na.rm = TRUE)

## Data frame for prediction over age
pred_age_df <- tibble(
  age          = age_seq,
  delay_diag   = ref_delay,
  sex          = ref_sex,
  place_infect = ref_place,
  kin          = ref_kin,
  kout         = ref_kout
)

pred_age_df <- pred_age_df %>%
  mutate(
    p_prop  = predict(mod_prop,  newdata = pred_age_df, type = "response"),
    p_naive = predict(mod_naive, newdata = pred_age_df, type = "response")
  ) %>%
  pivot_longer(cols = c(p_prop, p_naive),
               names_to = "model",
               values_to = "p") %>%
  mutate(
    model = recode(model,
                   p_prop  = "Proposed",
                   p_naive = "Naive")
  )

p_partial_age <- ggplot(pred_age_df, aes(x = age, y = p, colour = model)) +
  geom_line(size = 1.2) +
  scale_colour_brewer(palette = "Set1", name = "Model") +
  labs(
    x = "Age (years)",
    y = "Predicted probability of death",
    title = "Predicted mortality risk vs age\n(holding delay, sex, place, kin, kout fixed)"
  ) +
  theme_minimal()

partial_age_file <- file.path(plots_dir, "partial_effect_age.pdf")
ggsave(partial_age_file, p_partial_age, width = 6, height = 4.5)
print(p_partial_age)  ## show.plt


## Now partial effect of delay at fixed age
delay_seq <- seq(
  min(dat_cc$delay_diag, na.rm = TRUE),
  max(dat_cc$delay_diag, na.rm = TRUE),
  length.out = 100
)

ref_age <- median(dat_cc$age, na.rm = TRUE)

pred_delay_df <- tibble(
  age          = ref_age,
  delay_diag   = delay_seq,
  sex          = ref_sex,
  place_infect = ref_place,
  kin          = ref_kin,
  kout         = ref_kout
)

pred_delay_df <- pred_delay_df %>%
  mutate(
    p_prop  = predict(mod_prop,  newdata = pred_delay_df, type = "response"),
    p_naive = predict(mod_naive, newdata = pred_delay_df, type = "response")
  ) %>%
  pivot_longer(cols = c(p_prop, p_naive),
               names_to = "model",
               values_to = "p") %>%
  mutate(
    model = recode(model,
                   p_prop  = "Proposed",
                   p_naive = "Naive")
  )

p_partial_delay <- ggplot(pred_delay_df,
                          aes(x = delay_diag, y = p, colour = model)) +
  geom_line(size = 1.2) +
  scale_colour_brewer(palette = "Set1", name = "Model") +
  labs(
    x = "Onset-to-diagnosis delay (days)",
    y = "Predicted probability of death",
    title = "Predicted mortality risk vs delay\n(holding age, sex, place, kin, kout fixed)"
  ) +
  theme_minimal()

partial_delay_file <- file.path(plots_dir, "partial_effect_delay.pdf")
ggsave(partial_delay_file, p_partial_delay, width = 6, height = 4.5)
print(p_partial_delay)  ## show.plt


############################################################
## 9. Zip all results for convenient download
############################################################

## This uses base R's zip; requires a 'zip' utility to be installed.
## It creates mers_real_results.zip in the working directory.

all_files_to_zip <- list.files(out_base, recursive = TRUE, full.names = TRUE)

zip_file <- "mers_real_results.zip"
utils::zip(zipfile = zip_file, files = all_files_to_zip)

print(paste("Zipped results written to:", zip_file))
############################################################
## End of script
############################################################
