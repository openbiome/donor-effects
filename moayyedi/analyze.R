#!/usr/bin/env Rscript --vanilla

library(tidyverse)
library(exact2x2)
source("../utils.R", chdir=TRUE)

# Table tests ---------------------------------------------------------
telegraph("Table tests")

# Donor B's performance, compared to all other donors
data <- c(7, 18 - 7, 2, 20 - 2) %>%
  matrix(nrow = 2)

# Fisher-Freeman-Halton test with mid-p value
exact2x2(data, midp = TRUE)

# Adaptive trial simulation
telegraph("Adaptive trial simulation")

# Parameters
n_donors <- 6
n_explore_flips <- 4
n_total_flips <- 38
flip_p <- (7 + 2) / (18 + 20) # probability of "heads" = positive patient outcome
alpha <- 0.05
n_iter <- 1e5

max_n_flips <- n_total_flips - (n_donors - 1) * n_explore_flips
stopifnot(max_n_flips > 0)

simulate <- function() {
  # Simulate all the flips, exploration & exploitation
  flips <- rbinom(n_donors * max_n_flips, 1, flip_p) %>%
    matrix(nrow = n_donors)

  # Summarize donors' performance in exploration round
  exploration <- flips[, 1:n_explore_flips] %>% rowSums()
  best_donor_i <- which.max(exploration)

  # Compare the exploration-best donor with other donors
  best_donor_flips <- flips[best_donor_i, ]
  other_donor_flips <- flips[-best_donor_i, 1:n_explore_flips]

  good_success <- sum(best_donor_flips == 1)
  good_fail <- sum(best_donor_flips == 0)
  bad_success <- sum(other_donor_flips == 1)
  bad_fail <- sum(other_donor_flips == 0)

  matrix(c(good_success, good_fail, bad_success, bad_fail), nrow = 2)
}

results <- tibble(iter = 1:n_iter) %>%
  mutate(
    data = map(iter, ~ simulate()),
    midp = map_dbl(data, ~ exact2x2(., midp = TRUE)$p.value),
    significant = midp < 0.05
  )

cat("Proportion of simulated trials with midp < 0.05 ------------------")
results %>% { binom.test(sum(.$significant), nrow(.)) }
