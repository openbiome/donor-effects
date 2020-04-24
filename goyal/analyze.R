#!/usr/bin/env Rscript

library(tidyverse)
source("../utils.R")

n_patients <- 21

# Read the diversity file produced by Qiime
diversity <- read_tsv("diversity-data/alpha-diversity.tsv", skip = 1, col_names = c("sample_id", "diversity")) %>%
  mutate(patient_id = as.numeric(str_replace(sample_id, "^patient", "")))

stopifnot(nrow(diversity) == n_patients)

# Load patient data
patients <- read_tsv("patient-data.tsv")
stopifnot(nrow(patients) == n_patients)

# Join and export data
data <- patients %>%
  left_join(diversity, by = "patient_id")

telegraph("Summary of diversity")
data %>%
  group_by(outcome) %>%
  summarize(
    n_patients = n(),
    mean_diversity = mean(diversity),
    sem_diversity = sd(diversity) / sqrt(n_patients)
  )

telegraph("Wilcoxon of patients' diversity by outcome")
wilcox.test(diversity ~ outcome, data = data, conf.int = TRUE)

telegraph("Model of patients by outcome")
model <- glm(outcome ~ diversity, data = data, family = "binomial")
summary(model)

telegraph("Model predictions for extreme diversities")

tibble(
  extreme = c("min observed diversity", "max observed diversity"),
  diversity = c(min(data$diversity), max(data$diversity)),
  pred = map(diversity, ~ predict(model, newdata = tibble(diversity = .), se.fit = TRUE)),
  fit = map_dbl(pred, ~ .$fit),
  se = map_dbl(pred, ~ .$se),
  cil = fit - 1.96 * se,
  ciu = fit + 1.96 * se
) %>%
  mutate_at(c("fit", "cil", "ciu"), invlogit) %>%
  select(extreme, diversity, fit, cil, ciu)
