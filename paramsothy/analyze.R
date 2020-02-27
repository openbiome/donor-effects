#!/usr/bin/env Rscript

library(exact2x2)
library(lme4)
library(tidyverse)

source('../utils/utils.R')

# Re-do the analysis on donor D54 as reported in the original paper
cat("\n\nD54 comparison from the original paper -----------------------\n")

paper_data <- c(14, 38 - 14, 7, 40 - 7) %>%
  matrix(nrow = 2)

exact2x2(paper_data, midp = TRUE)

# Read in the full data
raw_data <- read_tsv('data.tsv', col_types = cols(patient = 'c', treatment = 'c', .default = 'i'))

# Check that the last donor is "M"
stopifnot(last(names(raw_data)) == "M")

data <- raw_data %>%
  filter(treatment != 'no_rescue') %>%
  # combine columns A through M into a single column "pool"
  unite("pool", A:M)

n_patients <- nrow(data)
n_success <- sum(data$outcome)

cat("\n\nOmnibus test of variance in efficacy by pool -----------------\n")

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

print(effect_sizes(model, by = 'pool'))
