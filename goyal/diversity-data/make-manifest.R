#!/usr/bin/env Rscript --vanilla

library(optparse)
library(tidyverse)

option_list <- list(
  make_option("--metadata"),
  make_option("--filereport"),
  make_option("--manifest"),
  make_option("--samples")
)

opts <- parse_args(OptionParser(option_list = option_list))

filereport <- read_tsv(opts$filereport) %>%
  # keep only donor samples
  filter(str_detect(library_name, "\\.D\\d?$")) %>%
  mutate(
    url = str_c("ftp://", fastq_ftp),
    filepath = str_c("fastq/", run_accession, ".fastq.gz"),
    absolute_filepath = str_c("$PWD/", filepath)
  )

manifest <- filereport %>%
  mutate(direction = "forward") %>%
  select(
    `sample-id` = library_name,
    `absolute-filepath` = absolute_filepath,
    direction
  )

samples <- filereport %>%
  select(filepath, url, md5 = fastq_md5)

write_csv(manifest, opts$manifest)
write_csv(samples, opts$samples)
