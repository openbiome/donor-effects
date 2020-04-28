#!/usr/bin/env Rscript

library(tidyverse)
library(patchwork)
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

alpha_plot <- data %>%
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
    limits = c(5, 8),
    expand = c(0, 0)
  ) +
  cowplot::theme_half_open()

telegraph("PERMANOVA")

distance_matrix <- read_tsv("diversity-data/distance-matrix.tsv")
sample_ids <- sprintf("patient%02i", patients$patient_id)

results <- permanova(distance_matrix, sample_ids, patients$outcome)

results$test

beta_plot <- results$mds %>%
  ggplot(aes(coord1, coord2, color = factor(outcome))) +
  geom_point(shape = 1, size = 3, stroke = 1) +
  scale_color_manual(
    "Outcome",
    breaks = c(0, 1),
    labels = c("No remission", "Remission"),
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
