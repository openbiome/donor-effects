#!/usr/bin/env Rscript

# Read the diversity file produced by Qiime
raw <- read_tsv('raw/alpha-diversity.tsv', skip = 1, col_names = c('sample_id', 'diversity'))

# Clean up the diversity table
diversity <- raw %>%
  mutate(patient_id = as.numeric(str_match(sample_id, 'FMT\\.0\\d\\.0*(\\d+)\\.D\\d*')[, 2])) %>%
  group_by(patient_id) %>%
  summarize_at('diversity', mean)

# Load patient data
patients <- read_tsv('patients.tsv')

# Join and export data
data <- patients %>%
  left_join(diversity, by = 'patient_id')

write_tsv(data, 'data.tsv')
