#!/usr/bin/env Rscript

library(exact2x2, quietly = TRUE)

# Donor B's performance, compared to all other donors
data <- c(7, 18 - 7, 2, 20 - 2) %>%
  matrix(nrow = 2)

# Fisher-Freeman-Halton test with mid-p value
exact2x2(data, midp = TRUE)
