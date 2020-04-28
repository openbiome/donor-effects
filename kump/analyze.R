#!/usr/bin/env Rscript

library(tidyverse)
library(patchwork)
source("../utils.R")

# Clean data ----------------------------------------------------------

meta <- read_tsv("kump2018.metadata.tsv") %>%
  filter(DonorID != "Control", Matter == "Donorstool") %>%
  select(
    sample = `#SampleID`,
    patient = PatientID,
    response = Response,
    donor = DonorID,
    timepoint = Sampling_day
  ) %>%
  mutate(
    donor = str_replace(donor, "^D", "donor"),
    patient = str_replace(patient, "^P", "patient"),
    timepoint = str_replace(timepoint, "^d", "day"),
    response = recode(response, NR = "no_response", PR = "partial_response", RE = "remission")
  )

patients <- meta %>%
  select(patient, response) %>%
  distinct() %>%
  arrange(patient)

cat(str_glue("\n\nNumber of patients with donor samples: {nrow(patients)}\n\n"))

raw_diversity_data <- read_tsv("diversity-data/alpha-diversity.tsv") %>%
  set_names(c("sample_id", "diversity")) %>%
  separate(sample_id, c("patient", "donor", "timepoint"))

donations <- meta %>%
  left_join(raw_diversity_data, by = c("patient", "donor", "timepoint")) %>%
  select(patient, donor, timepoint, diversity) %>%
  arrange(patient)


# Make sure we have the same counts as the original publication --------

expected_counts <- tibble(
  response = c("no_response", "remission"),
  n = c(12L, 16L)
)

actual_counts <- patients %>%
  filter(response != "partial_response") %>%
  left_join(donations, by = "patient") %>%
  count(response)

stopifnot(all_equal(expected_counts, actual_counts))


telegraph("Counts of donors and patients")
cat(str_glue("There are {length(unique(donations$donor))} donors\n\n"))
cat("Counts of patients who received FMT, by outcome:\n")
count(patients, response)

telegraph("Test of response by donation diversity")

donors <- donations %>%
  group_by(donor, patient) %>%
  summarize(
    n_samples = n(),
    diversity = mean(diversity)
  )

data <- patients %>%
  left_join(donors, by = "patient")

cat("Omnibus test\n")
with(data, { kruskal.test(diversity, response) })

cat("NR-RE\n")
data %>%
  filter(response %in% c("no_response", "remission")) %>%
  wilcox.test(diversity ~ response, data = ., conf.int = TRUE)

cat("NR-PR\n")
data %>%
  filter(response %in% c("no_response", "partial_response")) %>%
  wilcox.test(diversity ~ response, data = ., conf.int = TRUE)

cat("PR-RE\n")
data %>%
  filter(response %in% c("partial_response", "remission")) %>%
  wilcox.test(diversity ~ response, data = ., conf.int = TRUE)

# Plot of diversities by outcome --------------------------------------

alpha_plot <- data %>%
  ggplot(aes(response, diversity)) +
  geom_point(shape = 1, size = 3) +
  scale_y_continuous(
    "Diversity (Shannon)",
    limits = c(5, 7),
    expand = c(0, 0)
  ) +
  scale_x_discrete(
    "",
    breaks = c("no_response", "partial_response", "remission"),
    labels = c("No\nresponse", "Partial\nresponse", "Remission")
  ) +
  cowplot::theme_half_open()

# PERMANOVA -----------------------------------------------------------

telegraph("PERMANOVA")

distance_matrix <- read_tsv("diversity-data/distance-matrix.tsv")

permanova_data <- patients %>%
  filter(response != "partial_response") %>%
  mutate(outcome = as.numeric(response == "remission")) %>%
  left_join(donations, by = "patient") %>%
  group_by(patient) %>%
  filter(timepoint == min(timepoint)) %>%
  ungroup() %>%
  mutate(sample_id = str_c(patient, donor, timepoint, sep = "-"))

results <- permanova(distance_matrix, permanova_data$sample_id, permanova_data$outcome)

results$test

beta_plot <- results$mds %>%
  ggplot(aes(coord1, coord2, color = factor(outcome))) +
  geom_point(shape = 1, size = 3, stroke = 1) +
  scale_color_manual(
    "Outcome",
    breaks = c(0, 1),
    labels = c("No response", "Remission"),
    values = c("#d7191c", "#2c7bb6")
  ) +
  labs(x = "coordinate 1", y = "coordinate 2") +
  cowplot::theme_half_open() +
  theme(
    legend.position = "top"
  )

plot <- alpha_plot + beta_plot +
  plot_layout(widths = c(1, 1.25)) +
  plot_annotation(tag_level = "a")

ggsave("plot.pdf", width = 19, units = "cm")
