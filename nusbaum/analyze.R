#!/usr/bin/env Rscript

library(tidyverse)
library(scales)
library(patchwork)
source("../utils.R")

id_to_outcome <- function(x) {
  case_when(
    str_detect(x, "^Responder") ~ 1,
    str_detect(x, "^NonResponder") ~ 0
  )
}

# Diversity model -----------------------------------------------------
telegraph("Wilcox test on diversity data")

diversity_data <- read_tsv("diversity-data/alpha-diversity.tsv") %>%
  set_names(c("patient", "diversity")) %>%
  mutate(outcome = id_to_outcome(patient))

diversity_data

wilcox.test(diversity ~ outcome, data = diversity_data, conf.int = TRUE)

# Make the figure -----------------------------------------------------

alpha_plot <- diversity_data %>%
  mutate_at("outcome", factor) %>%
  ggplot(aes(outcome, diversity, color = outcome)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_point(shape = 1, size = 3) +
  scale_x_discrete(
    "",
    breaks = c(0, 1),
    labels = c("No response", "Response"),
  ) +
  scale_y_continuous(
    "Diversity (Shannon)",
    limits = c(4.0, 6.5),
    expand = c(0, 0)
  ) +
  scale_color_manual(values = c(red, blue)) +
  cowplot::theme_half_open() +
  theme(legend.position = "none")

# PERMANOVA -----------------------------------------------------------
telegraph("PERMANOVA")

distance_matrix <- read_tsv("diversity-data/distance-matrix.tsv")
sample_ids <- names(distance_matrix)[-1]
outcomes <- id_to_outcome(sample_ids)

results <- permanova(distance_matrix, sample_ids, outcomes)

results$test

beta_plot <- results$mds %>%
  ggplot(aes(coord1, coord2, color = factor(outcome))) +
  geom_point(shape = 1, size = 3, stroke = 1) +
  scale_color_manual(
    "Outcome",
    breaks = c(0, 1),
    labels = c("No response", "Response"),
    values = c(red, blue)
  ) +
  labs(x = "coordinate 1", y = "coordinate 2") +
  coord_fixed() +
  cowplot::theme_half_open() +
  theme(legend.position = "top")

# Assemble the two-part plot ------------------------------------------

plot <- alpha_plot + beta_plot +
  plot_layout(widths = c(1, 1.5)) +
  plot_annotation(tag_level = "a")

ggsave("plot.pdf", width = 19, height = 10, units = "cm")
