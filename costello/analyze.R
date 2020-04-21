#!/usr/bin/env Rscript

library(lme4)
library(exact2x2)
library(tidyverse)
source("../utils.R")

# Include only data from patients that received FMT
data <- read_tsv("patient-data.tsv") %>%
  filter(fmt == 1)

# For each pool, count (i) the number of patients who got this pool and
# had a positive outcome, (ii) a negative outcome, (iii) who got other
# pools and had a positive outcome, and (iv) a negative outcome
pool_table <- data %>%
  group_by(pool) %>%
  summarize(
    success = sum(outcome),
    fail = n() - success
  ) %>%
  mutate(
    other_success = sum(success) - success,
    other_fail = sum(success) - fail
  )

# 2 x P test ----------------------------------------------------------
pool_matrix <- pool_table %>%
  select(success, fail) %>%
  as.matrix()

telegraph("Fisher-Freeman-Halton test on 2xP table of outcomes by pool")
fisher.test(pool_matrix)


# 2x2 tests by donor --------------------------------------------------
telegraph("Fisher's exact test on 2x2 tables of outcomes by donor")

data %>%
  gather("donor", "present", starts_with("donor")) %>%
  group_by(donor, present) %>%
  summarize(
    success = sum(outcome),
    fail = n() - success
  ) %>%
  group_by(donor) %>%
  group_nest() %>%
  mutate(
    donor_matrix = map(data, ~ matrix(c(.$success, .$fail), nrow = 2)),
    test = map(donor_matrix, partial(exact2x2, midp = TRUE)),
    midp = map_dbl(test, ~ .$p.value),
    estimate = map_dbl(test, ~ .$estimate),
    cil = map_dbl(test, ~ .$conf.int[1]),
    ciu = map_dbl(test, ~ .$conf.int[2])
  ) %>%
  select(donor, midp, estimate, cil, ciu)


# Logistic regression of outcomes by pool -----------------------------

model <- glmer(outcome ~ (1 | pool), family = "binomial", data = data)

telegraph("Logistic regression of outcomes by pool")
summary(model)

telegraph("Typical deviations")
print(effect_sizes(model))
