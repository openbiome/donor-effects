invlogit <- function(x) exp(x) / (1 + exp(x))

median_difference <- function(mean, sd, n = 1e6) {
  p1 <- invlogit(rnorm(n, mean, sd))
  p2 <- invlogit(rnorm(n, mean, sd))
  median(abs(p1 - p2))
}

effect_sizes <- function(model, by = 'pool', accuracy = 0) {
  estimate <- model %>% summary %>% coef %>% { .['(Intercept)', 'Estimate'] }
  ran_sd <- model %>% VarCorr %>% `[[`(by) %>% attr('stddev') %>% unname
  parm <- str_glue('sd_(Intercept)|{by}')
  ran_sd_ci <- confint(model, oldNames = FALSE, parm = parm)

  baseline_efficacy <- invlogit(estimate)
  effect_size <- median_difference(estimate, ran_sd)
  effect_size_cil <- median_difference(estimate, ran_sd_ci[1])
  effect_size_ciu <- median_difference(estimate, ran_sd_ci[2])

  efficacy_format <- str_glue('%.{accuracy}f%%')
  effect_size_format <- str_glue('%.{accuracy}f p.p')

  show <- function(x) sprintf(effect_size_format, x * 100)

  list(
    baseline_efficacy = sprintf(efficacy_format, baseline_efficacy * 100),
    effect_size = show(effect_size),
    effect_size_cil = show(effect_size_cil),
    effect_size_ciu = show(effect_size_ciu)
  )
}
