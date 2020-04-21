#!/usr/bin/env Rscript

library(optparse)
library(tidyverse)

option_list <- c(
  make_option("--metadata"),
  make_option("--filereport"),
  make_option("--manifest"),
  make_option("--samples")
)

opts <- parse_args(OptionParser(option_list = option_list))

# Get the "sample ID" (called "sample_alias" in ENA) for donor samples
meta <- read_tsv(opts$metadata) %>%
  filter(DonorID != "Control", Matter == "Donorstool") %>%
  mutate(
    patient = str_replace(PatientID, "^P", "patient"),
    donor = str_replace(DonorID, "^D", "donor"),
    day = str_replace(Sampling_day, "^d", "day"),
    sample_id = str_c(patient, donor, day, sep = "-")
  ) %>%
  select(sample_alias = `#SampleID`, sample_id)

raw_filereport <- read_tsv(options$filereport)

filereport <- meta %>%
  left_join(raw_filereport, by = "sample_alias") %>%
  select(sample_id, fastq_ftp, fastq_md5) %>%
  # split the URLs and MD5s into separate rows
  mutate_at(c("fastq_ftp", "fastq_md5"), ~ str_split(., ";")) %>%
  mutate(direction = map(fastq_ftp, ~ c("forward", "reverse"))) %>%
  unnest() %>%
  mutate(
    url = str_c("ftp://", fastq_ftp),
    filepath = sprintf("fastq/%s-%s.fastq.gz", sample_id, direction),
    absolute_filepath = str_c("$PWD/", filepath),
  )

manifest <- filereport %>%
  select(
    `sample-id` = sample_id,
    `absolute-filepath` = absolute_filepath,
    direction
  )

samples <- filereport %>%
  select(filepath, url, md5 = fastq_md5)

write_csv(manifest, opts$manifest)
write_csv(samples, opts$samples)
