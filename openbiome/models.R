#!/usr/bin/env Rscript --vanilla

library(tidyverse)
library(lme4)
library(scales)

data <- read_tsv('data/data.tsv') %>%
  mutate(route = factor(route, levels = c('lower', 'upper', 'capsule', 'both')))

cat("\nMultivariate mixed model ----------------------------------------------\n")
model <- glmer(outcome ~ cdi_type + route + (1 | partner) + (1 | donor), family = 'binomial', data = data)
summary(model)

cat("\nLikelihood ratio test -------------------------------------------------\n")
null_model <- update(model, outcome ~ cdi_type + route + (1 | partner))
lmtest::lrtest(null_model, model)

cat("\nTypical deviates ------------------------------------------------------\n")

invlogit <- function(x) exp(x) / (1 + exp(x))
fix_est <- model %>% fixef %>% `[`('(Intercept)')
ran_sd <- model %>% VarCorr %$% donor %>% attr('stddev') %>% unname
ran_sd_ci <- model %>% confint(oldNames = FALSE, parm = 'sd_(Intercept)|donor')

typical <- function(mean, sd, n = 1e6) {
  p1 <- invlogit(rnorm(n, mean, sd))
  p2 <- invlogit(rnorm(n, mean, sd))
  median(abs(p1 - p2))
}

cat(str_glue("Baseline efficacy: {fix_est} -> {percent(invlogit(fix_est))}\n\n"))
cat(str_glue("Typical deviation: {ran_sd} -> {percent(typical(fix_est, ran_sd))}\n\n"))
cat(str_glue("Typical deviation (lower 95% CI): {ran_sd_ci[1]} -> {percent(typical(fix_est, ran_sd_ci[1]))}\n\n"))
cat(str_glue("Typical deviation (upper 95% CI): {ran_sd_ci[2]} -> {percent(typical(fix_est, ran_sd_ci[2]))}\n\n"))
