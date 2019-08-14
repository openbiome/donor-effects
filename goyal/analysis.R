#!/usr/bin/env Rscript

data <- read_tsv('data.tsv')

cat("\n\nSummary of diversity -----------------------------------------------\n")
data %>%
  group_by(outcome) %>%
  summarize(
    n_patients = n(),
    mean_diversity = mean(diversity),
    sem_diversity = sd(diversity) / sqrt(n_patients)
  )

cat("\n\nWilcoxon of patients' diversity by outcome -----------------------------\n")
wilcox.test(diversity ~ outcome, data = data, conf.int = TRUE)

cat("\n\nModel of patients by outcome -------------------------------------------\n")

model <- glm(outcome ~ diversity, data = data, family = 'binomial')
summary(model)

cat("\nModel predictions for extreme diversities ------------------------------\n")

invlogit <- function(lo) exp(lo) / (1 + exp(lo))

tibble(
  extreme = c('min observed diversity', 'max observed diversity'),
  diversity = c(min(data$diversity), max(data$diversity)),
  pred = map(diversity, ~ predict(model, newdata = tibble(diversity = .), se.fit = TRUE)),
  fit = map_dbl(pred, ~ .$fit),
  se = map_dbl(pred, ~ .$se),
  cil = fit - 1.96 * se,
  ciu = fit + 1.96 * se
) %>%
  mutate_at(c('fit', 'cil', 'ciu'), invlogit) %>%
  select(extreme, diversity, fit, cil, ciu)
