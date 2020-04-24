#!/usr/bin/env Rscript

library(tidyverse)
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
  select(donor, patient, timepoint, diversity)

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
  summarize_at("diversity", mean)

data <- patients %>%
  left_join(donors, by = "patient")

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


telegraph("Model of response, predicting by diversity")

model_data <- data %>%
  filter(response %in% c("remission", "no_response")) %>%
  mutate(y = recode(response, no_response = 0, remission = 1))

model <- glm(y ~ diversity, data = model_data, family = "binomial")

summary(model)

tibble(
  diversity = range(model_data$diversity),
  pred = map(diversity, ~ predict(model, newdata = tibble(diversity = .), se.fit = TRUE)),
  fit = map_dbl(pred, ~ .$fit),
  se = map_dbl(pred, ~.$se),
  cil = fit - 1.96 * se,
  ciu = fit + 1.96 * se
) %>%
  mutate_at(c("fit", "cil", "ciu"), invlogit)


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
    breaks = c("no_response", "partial_response", "remission"),
    labels = c("No\nresponse", "Partial\nresponse", "Remission")
  ) +
  theme_classic() +
  theme(
    panel.grid.major.x = element_line(color = "gray"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 8, color = "black")
  )

ggsave("plot.pdf", height = 3, width = 4, units = "in")
