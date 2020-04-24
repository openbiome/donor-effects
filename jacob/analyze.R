#!/usr/bin/env Rscript

library(exact2x2)
library(lme4)
library(tidyverse)
library(scales)
library(patchwork)
source("../utils.R")

n_patients <- 20

patient_data <- read_tsv("patient-data.tsv")

# check that the number of patients is what's reported in the paper
stopifnot(nrow(patient_data) == n_patients)

# Table tests and GLM -------------------------------------------------

telegraph("2xP test of pools")
pool_test <- patient_data %>%
  group_by(pool) %>%
  summarize(
    total = n(),
    success = sum(outcome),
    fail = total - success
  ) %>%
  select(success, fail) %>%
  as.matrix() %T>%
  # check that table counts add up to number of patients
  { stopifnot(sum(.) == nrow(patient_data)) } %>%
  fisher.test()

pool_test

telegraph("2x2 tests of donors")

total_success <- sum(patient_data$outcome)
total_fail <- sum(!patient_data$outcome)

donorwise <- patient_data %>%
  tidyr::gather("donor", "donor_present", starts_with("donor")) %>%
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
    sig_adj = p.adjust(test_p, "BH") < 0.05
  )

# check that each 2x2 table has the right number of counts
stopifnot(all(donorwise$n_patients == n_patients))

donorwise %>%
  select(donor, test_p)

telegraph("Mixed model by pool")

model <- glmer(outcome ~ (1 | pool), family = "binomial", data = patient_data)

summary(model)
print(effect_sizes(model))

# Diversity model -----------------------------------------------------

raw_diversity_data <- read_tsv("diversity-data/alpha-diversity.tsv") %>%
  set_names(c("patient", "diversity")) %>%
  mutate_at("patient", ~ as.integer(str_replace(., "^patient", "")))

stopifnot(nrow(raw_diversity_data) == n_patients)

diversity_data <- patient_data %>%
  left_join(raw_diversity_data, by = "patient") %>%
  filter(!is.na(diversity))

print(diversity_data)

telegraph("Distribution of diversities by patient outcomes")
diversity_data %>%
  group_by(outcome) %>%
  summarize(
    n = n(),
    mean_diversity = mean(diversity),
    sem_diversity = sd(diversity) / sqrt(n)
  )

telegraph("Wilcox test on diversity data")
wilcox.test(diversity ~ outcome, data = diversity_data, conf.int = TRUE)

telegraph("Logisitic regression")
invlogit <- function(lo) exp(lo) / (1 + exp(lo))

model <- glm(outcome ~ diversity, family = "binomial", data = diversity_data)
summary(model)

tibble(
  diversity = range(diversity_data$diversity),
  pred = map(diversity, ~ predict(model, newdata = tibble(diversity = .), se.fit = TRUE)),
  fit = map_dbl(pred, ~ .$fit),
  se = map_dbl(pred, ~.$se),
  cil = fit - 1.96 * se,
  ciu = fit + 1.96 * se
) %>%
  mutate_at(c("fit", "cil", "ciu"), invlogit)

# Make the figure -----------------------------------------------------

alpha_plot <- diversity_data %>%
  ggplot(aes(factor(outcome), diversity)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_point(shape = 1, size = 3) +
  scale_x_discrete(
    "",
    breaks = c(0, 1),
    labels = c("No response", "Response")
  ) +
  scale_y_continuous(
    "Diversity (Shannon)",
    limits = c(3.5, 6),
    expand = c(0, 0)
  ) +
  cowplot::theme_half_open()

telegraph("PERMANOVA")

distance_matrix <- read_tsv("diversity-data/distance-matrix.tsv")
sample_ids <- sprintf("patient%02i", patient_data$patient)

results <- permanova(distance_matrix, sample_ids, patient_data$outcome)

results$test

beta_plot <- results$mds %>%
  ggplot(aes(coord1, coord2, color = factor(outcome))) +
  geom_point(shape = 1, size = 3, stroke = 1) +
  scale_color_manual(
    "Outcome",
    breaks = c(0, 1),
    labels = c("No response", "Response"),
    values = c("#d7191c", "#2c7bb6")
  ) +
  labs(x = "coordinate 1", y = "coordinate 2") +
  cowplot::theme_half_open() +
  theme(
    legend.position = "top"
  )

# Assemble the two-part plot ------------------------------------------

plot <- alpha_plot + beta_plot +
  plot_layout(widths = c(1, 1.5)) +
  plot_annotation(tag_level = "a")

ggsave("plot.pdf", width = 19, units = "cm")
