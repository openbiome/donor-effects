#!/usr/bin/env Rscript

library(tidyverse)
library(lme4)
source("../utils.R")

data <- tribble(
  ~patient_id, ~outcome, ~donor,
  "a", 1, "4",
  "b", 0, "2",
  "c", 0, "3",
  "d", 1, "1",
  "e", 1, "4",
  "f", 1, "3",
  "g", 1, "1",
  "h", 1, "1",
  "i", 0, "3",
  "j", 1, "2",
  "k", 1, "2",
  "l", 1, "4",
  "m", 1, "1",
  "n", 0, "2",
  "o", 0, "1",
  "p", 0, "4",
  "q", 1, "3",
  "r", 0, "2"
)

# FFH test for outcome by donor
data %>%
  count(donor, outcome) %>%
  pivot_wider(names_from = outcome, values_from = n) %>%
  select(`0`, `1`) %>%
  { fisher.test(as.matrix(.)) }

model <- glmer(outcome ~ (1 | donor), family = "binomial", data = data)

model
print(effect_sizes(model, by = "donor"))
