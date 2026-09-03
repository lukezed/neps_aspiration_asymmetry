# script/00_packages.R
# Sourced by every other script. Libraries and global options only.

library(here)        # project-root-relative paths
library(haven)       # NEPS ships SPSS (.sav) and Stata (.dta)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(loo)

theme_set(bayesplot::theme_default())

options(mc.cores = 4)                     # 4 chains on the performance cores
options(cmdstanr_write_stan_file_dir = here("stan"))

neps_na <- function(x) {
  x <- haven::zap_labels(x)
  x[!is.na(x) & x < 0] <- NA
  x
}

points_to_note <- function(p) pmin(pmax((17 - p) / 3, 1), 6)

if (!dir.exists(here("models"))) dir.create(here("models"))

