#!/usr/bin/env Rscript

library(optparse)
library(tidyverse)

option_list <- list(
  make_option("--metadata"),
  make_option("--filereport"),
  make_option("--manifest"),
  make_option("--samples")
)

opts <- parse_args(OptionParser(option_list = option_list))

# donor samples of the form FMT.01.123.D, where "123" is the patient ID
patient_regex <- regex("^FMT\\.0[13]\\.(\\d{3})\\.D1?$")

filereport <- read_tsv(opts$filereport) %>%
  # keep only donor samples
  filter(str_detect(library_name, patient_regex)) %>%
  mutate(patient_id = as.numeric(str_match(library_name, patient_regex)[, 2])) %>%
  mutate(
    sample_id = sprintf("patient%02i", patient_id),
    url = str_c("ftp://", fastq_ftp),
    filepath = str_c("fastq/", run_accession, ".fastq.gz"),
    absolute_filepath = str_c("$PWD/", filepath)
  )

manifest <- filereport %>%
  mutate(direction = "forward") %>%
  select(
    `sample-id` = sample_id,
    `absolute-filepath` = absolute_filepath,
    direction
  )

samples <- filereport %>%
  select(filepath, url, md5 = fastq_md5)

write_csv(manifest, opts$manifest)
write_csv(samples, opts$samples)
