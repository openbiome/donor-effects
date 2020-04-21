#!/usr/bin/env Rscript

options(repos = c(CRAN = "https://cran.rstudio.com/"))

needed_packages <- c("tidyverse", "scales", "optparse", "lmtest", "lme4")
existing_packages <- installed.packages()[, 1]
missing_packages <- setdiff(needed_packages, existing_packages)
install.packages(missing_packages)
