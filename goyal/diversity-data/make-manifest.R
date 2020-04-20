#!/usr/bin/env Rscript --vanilla

library(tidyverse)

filereport <- read_tsv(snakemake@input$filereport) %>%
  # keep only donor samples
  filter(str_detect(library_name, '\\.D\\d?$')) %>%
  mutate(
    url = str_c('ftp://', fastq_ftp),
    filepath = str_c('fastq/', run_accession, '.fastq.gz'),
    absolute_filepath = str_c('$PWD/', filepath)
  )

manifest <- filereport %>%
  mutate(direction = 'forward') %>%
  select(
    `sample-id` = library_name,
    `absolute-filepath` = absolute_filepath,
    direction
  )

samples <- filereport %>%
  select(filepath, url, fastq_md5)

write_csv(manifest, snakemake@output$manifest)
write_csv(samples, snakemake@output$samples)
