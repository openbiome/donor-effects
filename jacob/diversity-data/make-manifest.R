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

patients <- read_tsv(opts$metadata)

filereport <- read_tsv(opts$filereport) %>%
  mutate(patient_id = as.numeric(str_match(library_name, 'Longman\\.FMT\\.(\\d+)\\.Donor')[, 2])) %>%
  # keep only donor samples
  filter(
    str_detect(library_name, 'Donor'),
    patient_id %in% patients$patient_id
  ) %>%
  select(patient_id, fastq_ftp, fastq_md5) %>%
  # split the URLs and MD5s into separate rows
  mutate_at(c('fastq_ftp', 'fastq_md5'), ~ str_split(., ';')) %>%
  mutate(direction = map(fastq_ftp, ~ c('forward', 'reverse'))) %>%
  unnest(cols = c(fastq_ftp, fastq_md5, direction)) %>%
  mutate(
    sample_id = sprintf('patient%02i', patient_id),
    url = str_c('ftp://', fastq_ftp),
    filepath = sprintf('fastq/%s_%s.fastq.gz', sample_id, direction),
    absolute_filepath = str_c('$PWD/', filepath),
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
