#!/usr/bin/env Rscript --vanilla

# Analysis of Rossen et al.

library(tidyverse)
library(exact2x2)

# Load data -----------------------------------------------------------

# Parameters from the whole trial
n_patients <- 23
n_success <- 7

# Data about the 3 donors mentioned in the paper
data <- tibble(
  donor = LETTERS[1:3],
  total = c(8, 3, 2),
  success = c(4, 1, 0),
  fail = total - success
)

# 2x2 tests comparing the 3 donors against all others -----------------

results <- data %>%
  mutate(
    other_success = n_success - success,
    other_fail = n_patients - n_success - fail,
    outcomes = pmap(
      list(success, fail, other_success, other_fail),
      ~ matrix(c(..1, ..2, ..3, ..4), nrow = 2)
    ),
    test = map(outcomes, partial(exact2x2, midp = TRUE)),
    midp = map_dbl(test, ~ .$p.value),
    odds_ratio = map_dbl(test, ~ .$estimate),
    or_lci = map_dbl(test, ~ .$conf.int[1]),
    or_uci = map_dbl(test, ~ .$conf.int[2])
  )

cat("Outcomes from Rossen et al. --------------------------------------\n")
results %>% select(donor, success, fail, other_success, other_fail)

cat("\nFisher's exact tests for each donor ------------------------------\n")
results %>% select(donor, odds_ratio, or_lci, or_uci, midp)
