#!/usr/bin/env Rscript

library(lme4)
library(tidyverse)

source('../utils/utils.R')

# Paramsothy ----------------------------------------------------------

paramsothy <- read_tsv(
  '../paramsothy/data.tsv',
  col_types = cols(patient = 'c', treatment = 'c', .default = 'i')
) %>%
  filter(treatment != 'no_rescue')

# get all columns, starting at A until the end
donors <- names(paramsothy) %>% { .[match('A', .):length(.)] }

# Add pool identifiers
pools <- paramsothy %>%
  select(donors) %>%
  apply(1, function(x) str_c(x, collapse = '')) %>%
  { LETTERS[match(., unique(.))] }

paramsothy <- paramsothy %>% mutate(pool = pools)

# Costello and Jacob --------------------------------------------------

costello <- read_tsv('../costello/data.tsv') %>%
  filter(fmt == 1)

jacob <- read_tsv('../jacob/data/data.tsv')

# Function to ensure that pools with the same name in different studies
# get a different global name
name_pool <- function(df, nm) {
  df %>%
    mutate(
      pool = str_c(nm, '_', pool),
      study = nm
    ) %>%
    select(study, pool, outcome)
}

data <- bind_rows(
  name_pool(paramsothy, 'paramsothy'),
  name_pool(costello, 'costello'),
  name_pool(jacob, 'jacob')
)


# Model ---------------------------------------------------------------

model <- glmer(outcome ~ study + (1 | pool), family = 'binomial', data = data)
summary(model)
print(effect_sizes(model))
