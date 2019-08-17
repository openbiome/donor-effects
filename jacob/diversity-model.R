#!/usr/bin/env Rscript --vanilla

library(lme4)
library(scales)
library(tidyverse)

patients <- read_tsv('data/data.tsv')
diversity <- read_tsv('data/raw/alpha-diversity.tsv') %>%
  mutate(patient_id = as.integer(str_match(X1, 'patient(\\d+)')[, 2])) %>%
  select(patient_id, diversity = observed_otus)

data <- patients %>%
  left_join(diversity, by = 'patient_id') %>%
  filter(!is.na(diversity))

data %>%
  group_by(outcome) %>%
  summarize(
    n = n(),
    mean_diversity = mean(diversity),
    sem_diversity = sd(diversity) / sqrt(n)
  )

wilcox.test(diversity ~ outcome, data = data, conf.int = TRUE)

model <- glm(outcome ~ diversity, family = 'binomial', data = data)
summary(model)

invlogit <- function(lo) exp(lo) / (1 + exp(lo))

tibble(
  diversity = range(data$diversity),
  pred = map(diversity, ~ predict(model, newdata = tibble(diversity = .), se.fit = TRUE)),
  fit = map_dbl(pred, ~ .$fit),
  se = map_dbl(pred, ~.$se),
  cil = fit - 1.96 * se,
  ciu = fit + 1.96 * se
) %>%
  mutate_at(c('fit', 'cil', 'ciu'), invlogit)

# figure

plot <- data %>%
  ggplot(aes(factor(outcome), diversity)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_point(shape = 1, size = 3) +
  scale_x_discrete(
    "",
    breaks = c(0, 1),
    labels = c("No remission", "Remission")
  ) +
  scale_y_continuous(
    "Diversity (number of ASVs)",
    limits = c(0, 500),
    expand = c(0, 0)
  ) +
  theme_classic() +
  theme(
    panel.grid.major.x = element_line(color = 'gray'),
    axis.text.x = element_text(size = 10, color = 'black'),
    axis.text.y = element_text(size = 8, color = 'black')
  )

ggsave(
  'diversity-by-outcome.pdf', plot = plot,
  height = 3, width = 4, units = 'in'
)
