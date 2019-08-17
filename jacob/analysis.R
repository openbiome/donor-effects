#!/usr/bin/env Rscript --vanilla

library(tidyverse)
library(exact2x2)
library(lme4)
library(scales)

data <- read_tsv('data/data.tsv')


cat(" 2 x P test of pools ------------------------------ \n\n")

pool_test <- data %>%
  group_by(pool) %>%
  summarize(
    total = n(),
    success = sum(outcome),
    fail = total - success
  ) %>%
  select(success, fail) %>%
  as.matrix() %T>%
  # check that table counts add up to number of patients
  { stopifnot(sum(.) == nrow(data)) } %>%
  fisher.test()

pool_test


cat(" 2x2 tests of donors ------------------------------ \n\n")

total_success <- sum(data$outcome)
total_fail <- sum(!data$outcome)

donorwise <- data %>%
  tidyr::gather('donor', 'donor_present', starts_with('donor')) %>%
  filter(donor_present == 1) %>%
  group_by(donor) %>%
  summarize(
    total = n(),
    success = sum(outcome),
    fail = total - success
  ) %>%
  mutate(
    other_success = total_success - success,
    other_fail = total_fail - fail,
    test_matrix = pmap(list(success, fail, other_success, other_fail), ~ matrix(c(..1, ..2, ..3, ..4), nrow = 2)),
    n_patients = map_dbl(test_matrix, sum),
    test = map(test_matrix, ~ exact2x2(., midp = TRUE)),
    test_p = map_dbl(test, ~ .$p.value),
    sig_adj = p.adjust(test_p, 'BH') < 0.05
  )

# check that each 2x2 table has the right number of counts
stopifnot(all(donorwise$n_patients == nrow(data)))

donorwise %>%
  select(donor, test_p)


cat("\n\n Mixed pool model ----------------------------------\n\n")

model <- glmer(outcome ~ (1 | pool), family = 'binomial', data = data)
summary(model)

cat("\n\nTypical deviations ---------------------------------\n")

invlogit <- function(lo) exp(lo) / (1 + exp(lo))

estimate <- model %>% summary %>% coef %>% { .['(Intercept)', 'Estimate'] }
ran_sd <- model %>% VarCorr %$% pool %>% attr('stddev') %>% unname
ran_sd_ci <- confint(model, oldNames = FALSE, parm = 'sd_(Intercept)|pool')

typical <- function(mean, sd, n = 1e6) {
  p1 <- invlogit(rnorm(n, mean, sd))
  p2 <- invlogit(rnorm(n, mean, sd))
  median(abs(p1 - p2))
}

cat(str_glue("Baseline efficacy: {estimate} -> {percent(invlogit(estimate))}\n\n"))
cat(str_glue("Typical deviation: {ran_sd} -> {percent(typical(estimate, ran_sd))}\n\n"))
cat(str_glue("Typical deviation (lower 95% CI): {ran_sd_ci[1]} -> {percent(typical(estimate, ran_sd_ci[1]))}\n\n"))
cat(str_glue("Typical deviation (upper 95% CI): {ran_sd_ci[2]} -> {percent(typical(estimate, ran_sd_ci[2]))}\n\n"))
