#!/usr/bin/env Rscript --vanilla

library(tidyverse)

patients <- read_tsv('data/patients.tsv')
donations <- read_tsv('data/donations.tsv')

patients
donations

# Make sure we have the same counts as the original publication -------
expected_counts <- tibble(
  response = c('no_response', 'remission'),
  n = c(12L, 16L)
)

actual_counts <- patients %>%
  filter(response != 'partial_response') %>%
  left_join(donations, by = 'patient') %>%
  count(response)

stopifnot(all_equal(expected_counts, actual_counts))


cat("\nCounts of donors and patients ---------------------------------------")

cat(str_glue("There are {length(unique(donations$donor))} donors\n\n"))
cat("Counts of patients who received FMT, by outcome:\n")
count(patients, response)


cat("\n\nTest of response by donation diversity ----------------------------------\n")

donors <- donations %>%
  group_by(donor, patient) %>%
  summarize_at('diversity', mean)

data <- patients %>%
  left_join(donors, by = 'patient')

cat("NR-RE\n")
data %>%
  filter(response %in% c('no_response', 'remission')) %>%
  wilcox.test(diversity ~ response, data = ., conf.int = TRUE)

cat("NR-PR\n")
data %>%
  filter(response %in% c('no_response', 'partial_response')) %>%
  wilcox.test(diversity ~ response, data = ., conf.int = TRUE)

cat("PR-RE\n")
data %>%
  filter(response %in% c('partial_response', 'remission')) %>%
  wilcox.test(diversity ~ response, data = ., conf.int = TRUE)


cat("\n\nModel of response, predicting by diversity ----------------------\n")

model_data <- data %>%
  filter(response %in% c('remission', 'no_response')) %>%
  mutate(y = recode(response, no_response = 0, remission = 1))

model <- glm(y ~ diversity, data = model_data, family = 'binomial')

summary(model)

invlogit <- function(lo) exp(lo) / (1 + exp(lo))

tibble(
  diversity = range(model_data$diversity),
  pred = map(diversity, ~ predict(model, newdata = tibble(diversity = .), se.fit = TRUE)),
  fit = map_dbl(pred, ~ .$fit),
  se = map_dbl(pred, ~.$se),
  cil = fit - 1.96 * se,
  ciu = fit + 1.96 * se
) %>%
  mutate_at(c('fit', 'cil', 'ciu'), invlogit)


# Plot of diversities by outcome --------------------------------------

plot <- data %>%
  ggplot(aes(response, diversity)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_point(shape = 1, size = 3) +
  scale_y_continuous(
    "Diversity (number of ASVs)",
    limits = c(0, 1050),
    expand = c(0, 0)
  ) +
  scale_x_discrete(
    "",
    breaks = c('no_response', 'partial_response', 'remission'),
    labels = c('No\nresponse', 'Partial\nresponse', 'Remission')
  ) +
  theme_classic() +
  theme(
    panel.grid.major.x = element_line(color = 'gray'),
    axis.text.x = element_text(size = 10, color = 'black'),
    axis.text.y = element_text(size = 8, color = 'black')
  )

ggsave(
  'diversity-by-response.pdf', plot = plot,
  height = 3, width = 4, units = 'in'
)
