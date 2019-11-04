#!/usr/bin/env Rscript --vanilla

library(lme4)
library(lmtest)
library(tidyverse)
source('../utils/utils.R')

data <- read_tsv('private/data.tsv')

cat(str_glue("There are {nrow(data)} patients, treated with {length(unique(data$donor))} donors at {length(unique(data$partner))} partners\n"))

cat("\nMultivariate mixed model ----------------------------------------------\n")
model <- glmer(outcome ~ cdi_type + route + (1 | partner) + (1 | donor), family = 'binomial', data = data)
summary(model)

cat("\nLikelihood ratio tests ------------------------------------------------\n")
null_model <- update(model, outcome ~ cdi_type + route + (1 | partner))
lrtest(null_model, model)

null_model <- update(model, outcome ~ cdi_type + route + (1 | donor))
lrtest(null_model, model)

cat("\nDonor effect effect size ----------------------------------------------\n")
print(effect_sizes(model, by = 'donor', accuracy = 1))

cat("\nVariation between clinical partners -----------------------------------\n")
print(effect_sizes(model, by = 'partner', accuracy = 1))
