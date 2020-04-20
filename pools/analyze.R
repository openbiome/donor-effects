#!/usr/bin/env Rscript

library(lme4)
library(tidyverse)

source("../utils.R")

# Paramsothy ----------------------------------------------------------

paramsothy <- read_tsv(
  "../paramsothy/patient-data.tsv",
  col_types = cols(patient = "c", treatment = "c", .default = "i")
) %>%
  filter(treatment != "no_rescue")

# get all columns, starting at A until the end
donors <- names(paramsothy) %>% { .[match("A", .):length(.)] }

# Add pool identifiers
pools <- paramsothy %>%
  select(donors) %>%
  apply(1, function(x) str_c(x, collapse = "")) %>%
  { LETTERS[match(., unique(.))] }

paramsothy <- paramsothy %>% mutate(pool = pools)

# Costello and Jacob --------------------------------------------------

costello <- read_tsv("../costello/patient-data.tsv") %>%
  filter(fmt == 1)

jacob <- read_tsv("../jacob/patient-data.tsv")

# Makre sure pools in different studies have the same names
data <- bind_rows(paramsothy, costello, jacob, .id = "study") %>%
  rename(pool_in_study = pool) %>%
  mutate(pool = str_c(study, "_", pool_in_study))


# Model ---------------------------------------------------------------

model <- glmer(outcome ~ study + (1 | pool), family = "binomial", data = data)
summary(model)
print(effect_sizes(model))
