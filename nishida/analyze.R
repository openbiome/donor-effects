#!/usr/bin/env Rscript

library(tidyverse)

# Read in the data from Nishida et al. Table 4
table4 <- read_tsv("table4.tsv")

table4 %>%
  mutate(MHC_significant = p.adjust(table4$p_value, "BH") < 0.05)
