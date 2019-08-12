#!/usr/bin/env Rscript

library(lme4)
library(tidyr)
library(exact2x2)
import::from(scales, percent)

data <- read_tsv('data.tsv') %>%
  filter(fmt == 1)

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

cat("\nFisher-Freeman-Halton test on 2xP table of outcomes by pool -----------\n")
fisher.test(pool_matrix)


# 2x2 tests by donor --------------------------------------------------
cat("\nFisher's exact test on 2x2 tables of outcomes by donor ----------------\n")

data %>%
  gather('donor', 'present', starts_with('donor')) %>%
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

model <- glmer(outcome ~ (1 | pool), family = 'binomial', data = data)

cat("\nLogistic regression of outcomes by pool -----------------------------\n")
summary(model)

cat("\nTypical deviations --------------------------------------------------\n")

invlogit <- function(x) exp(x) / (1 + exp(x))
fix_est <- fixef(model)
ran_sd_ci <- confint(model, parm = '.sig01')

typical <- function(mean, sd, n = 1e6) {
  p1 <- invlogit(rnorm(n, mean, sd))
  p2 <- invlogit(rnorm(n, mean, sd))
  median(abs(p1 - p2))
}

cat(str_glue("Baseline efficacy: {fix_est} -> {percent(invlogit(fix_est))}\n\n"))
cat(str_glue("Typical deviation (lower 95% CI): {ran_sd_ci[1]} -> {percent(typical(fix_est, ran_sd_ci[1]))}\n\n"))
cat(str_glue("Typical deviation (upper 95% CI): {ran_sd_ci[2]} -> {percent(typical(fix_est, ran_sd_ci[2]))}\n\n"))
