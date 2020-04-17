# Interface conveniences ----------------------------------------------

telegraph <- function(message, width = 80, left_buffer = 2) {
  dashes <- function(x) str_c(rep("—", x), collapse = "")

  right_buffer <- width - str_length(message) - left_buffer - 2
  if (right_buffer <= 0) stop("Message too long")

  cat(str_c("\n", dashes(left_buffer), " ", message, " ", dashes(right_buffer), "\n", sep = ""))
}

# Analysis tools ------------------------------------------------------

invlogit <- function(x) exp(x) / (1 + exp(x))

# Given a log odds mean ("estimate") and standard deviation, what is the difference in
# the probabilities (i.e., invlogit of log odds) between the alpha-th and (1-alpha)-th
# percentiles?
effect_size_helper <- function(estimate, sd, alpha = 0.1) {
  invlogit(estimate + sd * qnorm(1 - alpha)) - invlogit(estimate + sd * qnorm(alpha))
}

effect_sizes <- function(model, by = "pool", accuracy = 0, alpha = 0.1) {
  estimate <- coef(summary(model))["(Intercept)", "Estimate"]
  baseline_efficacy <- invlogit(estimate)

  ran_sd <- model %>% VarCorr %>% `[[`(by) %>% attr("stddev") %>% unname
  parm <- str_glue("sd_(Intercept)|{by}")
  ran_sd_ci <- confint(model, oldNames = FALSE, parm = parm)

  effect_size <- effect_size_helper(estimate, ran_sd)
  effect_size_cil <- effect_size_helper(estimate, ran_sd_ci[1])
  effect_size_ciu <- effect_size_helper(estimate, ran_sd_ci[2])

  efficacy_format <- str_glue("%.{accuracy}f%%")
  effect_size_format <- str_glue("%.{accuracy}f p.p")

  show <- function(x) sprintf(effect_size_format, x * 100)

  list(
    baseline_efficacy = sprintf(efficacy_format, baseline_efficacy * 100),
    effect_size = show(effect_size),
    effect_size_cil = show(effect_size_cil),
    effect_size_ciu = show(effect_size_ciu)
  )
}
