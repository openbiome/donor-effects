#!/usr/bin/env Rscript

library(exact2x2)
library(lme4)
library(tidyverse)
library(scales)
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
  set_names(c("patient_id", "diversity")) %>%
  mutate_at("patient_id", ~ as.integer(str_match(., "patient(\\d+)")[, 2]))

stopifnot(nrow(raw_diversity_data) == n_patients)

diversity_data <- patient_data %>%
  left_join(raw_diversity_data, by = "patient_id") %>%
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

plot <- diversity_data %>%
  ggplot(aes(factor(outcome), diversity)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_point(shape = 1, size = 3) +
  scale_x_discrete(
    "",
    breaks = c(0, 1),
    labels = c("No remission", "Remission")
  ) +
  scale_y_continuous(
    "Diversity (number of ASVs)",
    limits = c(0, 500),
    expand = c(0, 0)
  ) +
  theme_classic() +
  theme(
    panel.grid.major.x = element_line(color = "gray"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 8, color = "black")
  )

ggsave("plot.pdf", height = 3, width = 4, units = "in")
