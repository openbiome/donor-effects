#!/usr/bin/env Rscript --vanilla

# Get list of patients: unique combinations of patient / donor / response
# and list of donor samples: donor / sample / timepoint

library(tidyverse)

meta <- read_tsv('raw/kump2018.metadata.tsv') %>%
  filter(DonorID != 'Control', Matter == 'Donorstool') %>%
  select(
    sample = `#SampleID`,
    patient = PatientID,
    response = Response,
    donor = DonorID,
    timepoint = Sampling_day
  ) %>%
  # change donors from "DB" to just "B", patients from "PA" to "A"
  mutate(
    donor = str_replace(donor, '^D', 'donor'),
    patient = str_replace(patient, '^P', 'patient'),
    timepoint = str_replace(timepoint, '^d', 'day'),
    response = recode(response, NR = 'no_response', PR = 'partial_response', RE = 'remission')
  )

patients <- meta %>%
  select(patient, response) %>%
  distinct() %>%
  arrange(patient)

write_tsv(patients, 'patients.tsv')

diversity <- read_tsv('raw/alpha-diversity.tsv') %>%
  separate(X1, c('patient', 'donor', 'timepoint')) %>%
  rename(diversity = observed_otus)

donations <- meta %>%
  left_join(diversity, by = c('patient', 'donor', 'timepoint')) %>%
  select(donor, patient, timepoint, diversity)

write_tsv(patients, 'patients.tsv')
write_tsv(donations, 'donations.tsv')
