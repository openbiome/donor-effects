#!/usr/bin/env Rscript

library(tidyverse)
library(lme4)
source("../utils.R")

# Donor "D" in the paper is actually 7 different donors, one per patient
data <- tribble(
  ~donor, ~n_remission, ~n_response, ~n_total,
  "A",    6,            8,           12,
  "B",    3,            4,           7,
  "C",    1,            3,           4,
  "D1",   1,            1,           1,
  "D2",   1,            1,           1,
  "D3",   1,            1,           1,
  "D4",   0,            1,           1,
  "D5",   0,            1,           1,
  "D6",   0,            1,           1,
  "D7",   0,            0,           1
)

# Run FFH test using remission
data %>%
  mutate(n_fail = n_total - n_remission) %>%
  select(n_remission, n_fail) %>%
  { fisher.test(as.matrix(.)) }

# Run FFH test using response
data %>%
  mutate(n_fail = n_total - n_response) %>%
  select(n_response, n_fail) %>%
  { fisher.test(as.matrix(.)) }


# Run a mixed logistic regression using remission
long_data <- data %>%
  mutate(n_success = n_remission, n_fail = n_total - n_success) %>%
  mutate(outcome = map2(n_success, n_fail, ~ c(rep(1, .x), rep(0, .y)))) %>%
  select(donor, outcome) %>%
  unnest()

model <- glmer(outcome ~ (1 | donor), family = "binomial", data = long_data)

model
print(effect_sizes(model, by = "donor"))
