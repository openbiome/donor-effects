#!/usr/bin/env Rscript --vanilla

library(tidyverse)
library(scales)

powers <- read_tsv('powers.tsv')

# Bounds on the colors
lower_bound <- 0.10
upper_bound <- 0.30
stopifnot(all(between(powers$effect_size, lower_bound, upper_bound), na.rm = TRUE))

plot <- powers %>%
  mutate(label = if_else(
    is.na(effect_size),
    '',
    as.character(round(effect_size * 100))
  )) %>%
  ggplot(aes(factor(n_donors), factor(n_patients))) +
  geom_tile(aes(fill = effect_size)) +
  geom_text(aes(label = label), size = 8 * 25.4 / 72) +
  scale_fill_gradient(
    name = 'Effect size (p.p.)',
    limits = c(lower_bound, upper_bound + 1e-6),
    labels = function(x) round(x * 100),
    low = 'white', high = 'red',
    na.value = 'gray75'
  ) +
  coord_fixed() +
  xlab('Number of donors') +
  ylab('Number of FMT-treated patients') +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    axis.text = element_text(size = 8, color = 'black'),
    legend.text = element_text(size = 8, color = 'black'),
    legend.title = element_text(size = 8, color = 'black')
  )

ggsave('power-landscape.pdf', plot = plot, width = 4, height = 2.8, units = 'in')
