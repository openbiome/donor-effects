#!/usr/bin/env Rscript --vanilla

# Adaptive trial simulation

library(tidyverse)
library(exact2x2)

n_donors <- 6
n_explore_flips <- 4
n_total_flips <- 38
flip_p <- (7 + 2) / (18 + 20)
alpha <- 0.05

max_n_flips <- n_total_flips - (n_donors - 1) * n_explore_flips
stopifnot(max_n_flips > 0)

simulate <- function() {
  flips <- rbinom(n_donors * max_n_flips, 1, flip_p) %>%
    matrix(nrow = n_donors)

  exploration <- flips[, 1:n_explore_flips] %>% rowSums()
  best_donor_i <- which.max(exploration)

  best_donor_flips <- flips[best_donor_i, ]
  other_donor_flips <- flips[-best_donor_i, 1:n_explore_flips]

  good_success <- sum(best_donor_flips == 1)
  good_fail <- sum(best_donor_flips == 0)
  bad_success <- sum(other_donor_flips == 1)
  bad_fail <- sum(other_donor_flips == 0)

  matrix(c(good_success, good_fail, bad_success, bad_fail), nrow = 2)
}

results <- tibble(iter = 1:1e4) %>%
  mutate(
    data = rerun(n(), simulate()),
    midp = map_dbl(data, ~ exact2x2(., midp = TRUE)$p.value),
    significant = midp < 0.05
  )

cat("Proportion of simulated trials with midp < 0.05 ------------------")
results %>% { binom.test(sum(.$significant), nrow(.)) }
