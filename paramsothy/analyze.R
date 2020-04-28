#!/usr/bin/env Rscript

library(exact2x2)
library(lme4)
library(tidyverse)

source("../utils.R")

# Re-do the analysis on donor D54 as reported in the original paper
telegraph("D54 comparison from the original paper")

paper_data <- c(14, 38 - 14, 7, 40 - 7) %>%
  matrix(nrow = 2)

exact2x2(paper_data, midp = TRUE)

# Read in the full data
patient_data <- read_tsv("patient-data.tsv") %>%
  filter(treatment != "no_rescue")

total_patients <- nrow(patient_data)
total_success <- sum(patient_data$outcome)
total_fail <- total_patients - total_success

# Data by donor

donor_data <- patient_data %>%
  gather("donor", "present", A:M) %>%
  filter(present == 1) %>%
  group_by(donor) %>%
  summarize(
    n = n(),
    success = sum(outcome),
    fail = n - success,
    other_success = total_success - success,
    other_fail = total_fail - fail
  ) %>%
  mutate(
    x = pmap(list(success, fail, other_success, other_fail), c),
    contingency_table = map(x, ~ matrix(., nrow = 2)),
    test = map(contingency_table, ~ exact2x2(., midp = TRUE)),
    p_value = map_dbl(test, ~ .$p.value),
    fdr = p.adjust(p_value, "BH"),
    bonferroni_p = p.adjust(p_value)
  )

donor_data %>%
  select(donor, success, fail, p_value, fdr, bonferroni_p)

# Data by pool
# Combine columns A through M into a single column "pool"

pool_data <- patient_data %>%
  unite("pool", A:M)

telegraph("Omnibus test of variance in efficacy by pool")

pool_data %>%
  group_by(pool) %>%
  summarize(
    success = sum(outcome),
    fail = n() - success
  ) %>%
  select(success, fail) %>%
  as.matrix() %>%
  fisher.test()

telegraph("Mixed model of efficacy with random pool effect")

model <- glmer(outcome ~ (1 | pool), family = "binomial", data = pool_data)
summary(model)

telegraph("Typical deviations")
print(effect_sizes(model, by = "pool"))
