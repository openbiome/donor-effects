#!/usr/bin/env Rscript

library(tidyverse)
library(exact2x2)

# Define functions -----------------------------------------------------

logit <- function(p) log(p / (1 - p))
invlogit <- function(x) exp(x) / (1 + exp(x))

# Chi-square test should give a p-value of 1, not NA, in cases
# like: chisq.test(matrix(c(1, 0, 0, 0), nrow = 2))
chisq_test <- function(x) {
  test <- chisq.test(x)
  if (is.na(test$p.value)) test$p.value <- 1
  test
}

# Fisher-Freeman-Halton test if possible, and chi-square if not
test_f <- function(x) {
  tryCatch(
    fisher.test(x),
    error = function(y) chisq_test(x)
  )
}

# Simulate outcomes from a single trial given the donor's efficacies
# and the number of patients per donor
simulate_outcomes <- function(donor_p, patients_per_donor) {
  rbinom(donor_p, patients_per_donor, donor_p) %>%
    { matrix(c(., patients_per_donor - .), nrow = 2, byrow = TRUE) }
}

# Simulate many clinicial trials
#   mean_log_odds : mean efficacy across donors
#   placebo_log_odds : rough lower limit of donor efficacy
#   n_steps : number of slices in SD-log-odds space, for estimating intercept with 80% power
#   n_iter : number of clinical trials per step
#   serialize : cache results of the simulations?
simulate_trials <- function(mean_log_odds, placebo_log_odds, n_patients, n_donors, n_steps = 10, n_iter = 1000, serialize = TRUE) {
  stopifnot(mean_log_odds > placebo_log_odds)

  # Store the parameters, which uniquely identify the simulation
  control <- list(
    mean_log_odds = mean_log_odds,
    placebo_log_odds = placebo_log_odds,
    n_donors = n_donors,
    n_patients = n_patients,
    n_steps = n_steps,
    n_iter = n_iter
  )

  # Check if this simulation has been performed before; if so, load it
  hash <- digest::digest(control)
  fn <- str_glue("cache/{hash}.rds")
  if (serialize && file.exists(fn)) return(readRDS(fn)$results)

  # Run n_iter simulations for the n_steps slices between 0 and mean - placebo
  # For each simulation, draw random donor efficacies, simulate outcomes, and
  # compute the statistical test
  results <- crossing(
    sd_log_odds = seq(1e-6, mean_log_odds - placebo_log_odds, length.out = n_steps),
    iter = 1:n_iter
  ) %>%
    mutate(
      patients_per_donor = n_patients / n_donors,
      donor_log_odds = pmap(list(n_donors, mean_log_odds, sd_log_odds), rnorm),
      donor_p = map(donor_log_odds, invlogit),
      outcomes = map2(donor_p, patients_per_donor, simulate_outcomes),
      test = map(outcomes, test_f),
      p_value = map_dbl(test, ~ .$p.value),
      sig = p_value < 0.05
    )

  if (serialize) {
    list(control = control, results = results) %>%
      saveRDS(fn)
  }

  results
}

# Run simulations -----------------------------------------------------

simulations <- crossing(
  mean_log_odds = logit(0.30), # mean 30% efficacy
  placebo_log_odds = logit(0.05), # 5% placebo rate
  n_patients = c(24, 48, 96, 192, 384),
  n_donors = c(2, 3, 4, 6, 12, 24)
) %>%
  mutate(simulations = pmap(list(mean_log_odds, placebo_log_odds, n_patients, n_donors), simulate_trials))


# Compute powers ------------------------------------------------------

# Given some (x, y) coordinates, find the value x0 where y = y0 using
# linear interpolation between the y values that bracket y0
approxroot <- function(x, y, y0 = 0.8) {
  if (max(y) < y0) return(list(status = 'always under power', root = NA))
  if (all(y > y0)) return(list(status = 'always over power', root = NA))

  left <- max(which(y < y0))
  right <- min(which(y > y0))
  idx <- c(left, right)

  f <- approxfun(x[idx], y[idx])
  obj <- function(x) f(x) - y0

  tryCatch(
    list(status = 'success', root = uniroot(obj, x[idx])$root),
    error = function(y) list(status = 'uniroot failed', root = NA)
  )
}

# Mean difference in treatment efficacy (MDTE)
mdte <- function(mean_log_odds, sd_log_odds, n = 1e6) {
  if (is.na(sd_log_odds)) return(NA)
  p <- invlogit(rnorm(n * 2, mean_log_odds, sd_log_odds))
  median(abs(p[c(TRUE, FALSE)] - p[c(FALSE, TRUE)]))
}

# For each set of simulations, find the log odds SD (or, equivalently, the
# MDTE) for 80% power
powers <- simulations %>%
  unnest() %>%
  group_by(mean_log_odds, placebo_log_odds, n_patients, n_donors, sd_log_odds) %>%
  summarize(power = mean(sig)) %>%
  nest() %>%
  mutate(
    power_root = map(data, ~ approxroot(.$sd_log_odds, .$power)),
    status = map_chr(power_root, ~ .$status),
    sd_log_odds_power80 = map_dbl(power_root, ~ .$root),
    effect_size = map2_dbl(mean_log_odds, sd_log_odds_power80, mdte)
  )

# check that all the root findings were successful
stopifnot(all(powers$status %in% c('success', 'always under power')))

# Save the results
powers %>%
  select(n_patients, n_donors, effect_size) %>%
  write_tsv('powers.tsv')
