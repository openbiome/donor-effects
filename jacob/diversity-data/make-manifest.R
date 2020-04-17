#!/usr/bin/env Rscript --vanilla

library(tidyverse)

patients <- read_tsv(snakemake@input$patient_data)

filereport <- read_tsv(snakemake@input$filereport) %>%
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
  unnest() %>%
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

write_csv(manifest, snakemake@output$manifest)
write_csv(samples, snakemake@output$samples)
