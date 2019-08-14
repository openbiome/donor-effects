#!/usr/bin/env Rscript

library(exact2x2, quietly = TRUE)
library(lme4, quietly = TRUE)
library(scales, quietly = TRUE)

cat("\n\nD54 comparison from the original paper -----------------------\n")

paper_data <- c(14, 38 - 14, 7, 40 - 7) %>%
  matrix(nrow = 2)

exact2x2(paper_data, midp = TRUE)

data <- read_tsv('data/data.tsv', col_types = cols(patient = 'c', treatment = 'c', .default = 'i')) %>%
  filter(treatment != 'no_rescue')

n_patients <- nrow(data)
n_success <- data %>% filter(outcome == 1) %>% nrow

# get all columns, starting at A until the end
donors <- names(data) %>% { .[match('A', .):length(.)] }


cat("\n\nOmnibus test of variance in efficacy by pool -----------------\n")

pools <- data %>%
  select(donors) %>%
  apply(1, function(x) str_c(x, collapse = '')) %>%
  { LETTERS[match(., unique(.))] }

data %<>% mutate(pool = pools)

data %>%
  group_by(pool) %>%
  summarize(
    success = sum(outcome),
    fail = n() - success
  ) %>%
  select(success, fail) %>%
  as.matrix() %>%
  fisher.test()


cat("\n\nMixed model of efficacy with random pool effect --------------\n")

model <- glmer(outcome ~ (1 | pool), family = 'binomial', data = data)

summary(model)

cat("\n\nTypical deviations -------------------------------------------\n")

invlogit <- function(x) exp(x) / (1 + exp(x))
fix_est <- fixef(model)
ran_sd_ci <- confint(model, parm = 'sd_(Intercept)|pool', oldNames = FALSE)

typical <- function(mean, sd, n = 1e6) {
  p1 <- invlogit(rnorm(n, mean, sd))
  p2 <- invlogit(rnorm(n, mean, sd))
  median(abs(p1 - p2))
}

cat("\n")
cat(str_glue("Baseline efficacy: {fix_est} -> {percent(invlogit(fix_est))}\n\n"))
cat(str_glue("Typical deviation (lower 95% CI): {ran_sd_ci[1]} -> {percent(typical(fix_est, ran_sd_ci[1]))}\n\n"))
cat(str_glue("Typical deviation (upper 95% CI): {ran_sd_ci[2]} -> {percent(typical(fix_est, ran_sd_ci[2]))}\n\n"))
