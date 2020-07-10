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

raw_filereport <- read_tsv(opts$filereport)

filereport <- raw_filereport %>%
  select(library_name, fastq_ftp, fastq_md5) %>%
  separate(library_name, c("sequencing_type", "outcome", "id", "sample_type")) %>%
  filter(sequencing_type == "16S", sample_type == "Donor") %>%
  mutate(sample_id = str_c(outcome, id)) %>%
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
