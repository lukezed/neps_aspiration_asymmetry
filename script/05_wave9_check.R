
source(here::here("script", "00_packages.R"))

# Every cell, on whichever build the analysis is carrying. AUM_SUFFIX="" reruns
# the check against the four-covariate build for comparison.
SUFFIX <- Sys.getenv("AUM_SUFFIX", "_cov")
CELLS <- paste0(c("aum_math_realistic", "aum_read_realistic",
                  "aum_math_idealistic", "aum_read_idealistic"), SUFFIX)
model <- cmdstan_model(here("stan", "aum.stan"))

MATRICES <- c("observed", "y", "v_obs", "obs_a", "obs_sc", "obs_x",
              "obs_news", "s_a", "s_sc", "s_x", "news_a")

drop_wave <- function(d, wave) {
  j <- which(d$waves == wave)
  d[MATRICES] <- lapply(d[MATRICES], function(m) m[, -j, drop = FALSE])
  d$fielded_a <- d$fielded_a[-j]
  d$fielded_sc <- d$fielded_sc[-j]
  d$fielded_x <- d$fielded_x[-j]
  d$waves <- d$waves[-j]
  d$n_wave <- length(d$waves)
  d
}

for (cell in CELLS) {
  d <- drop_wave(readRDS(here("data", "processed", paste0(cell, ".rds"))), 9)
  nm <- paste0(cell, "_m5_drop9")
  message("fitting ", nm, ": ", d$n_person, " students, ", d$n_wave, " waves")

  fit <- model$sample(
    data = c(d[c("n_person", "n_wave", "n_cov", "observed", "y", "v_obs",
                 "fielded_a", "fielded_sc", "fielded_x",
                 "obs_a", "obs_sc", "obs_x",
                 "s_a", "s_sc", "s_x", "news_a", "obs_news",
                 "z", "n_mod", "w")],
             list(model_id = 5L)),
    chains = 4, parallel_chains = 4,
    iter_warmup = 1000, iter_sampling = 1000, refresh = 0, seed = 1
  )
  fit$save_object(here("models", paste0(nm, ".rds")))

  # From the draws, not fit$summary(): its quantile column names fail silently.
  dr <- as_draws_df(fit$draws(c("beta_a", "b_pos", "b_neg")))
  keep <- which(d$fielded_a == 1)[-1]
  cat("\n", nm, ": mark coefficient by wave\n", sep = "")
  for (t in keep) {
    x <- dr[[sprintf("beta_a[%d]", t)]]
    cat(sprintf("  W%-2d  %.3f [%.3f, %.3f]\n", d$waves[t],
                mean(x), quantile(x, 0.025), quantile(x, 0.975)))
  }
  first <- keep[1]; last <- keep[length(keep)]
  rise <- dr[[sprintf("beta_a[%d]", last)]] - dr[[sprintf("beta_a[%d]", first)]]
  cat(sprintf("  rise from W%d to W%d: %.3f [%.3f, %.3f], P(rise) = %.3f\n",
              d$waves[first], d$waves[last], mean(rise),
              quantile(rise, .025), quantile(rise, .975), mean(rise > 0)))

  stopifnot("wave 9 still carries a self-concept measurement" =
              sum(d$fielded_sc) == 1)

  asym <- dr[["b_pos[1]"]] - dr[["b_neg[1]"]]
  cat(sprintf("  b+ minus b-: %.3f [%.3f, %.3f], Pr = %.3f\n",
              mean(asym), quantile(asym, .025), quantile(asym, .975),
              mean(asym > 0)))
}

# -----------------------------------------------------------------------------
# RESULT (re-run 2026-08-23; four cells, model 5 on waves 4 to 8).
# 0 divergent, 0 max-treedepth, 16 chains, no rejected proposals in sampling.
#
# b_pos - b_neg with wave 9 absent, against the full-data value:
#   math_realistic   +.074 [+.045, +.104]   Pr 1.000   (full data +.080)
#   read_realistic   +.117 [+.091, +.143]   Pr 1.000   (full data +.135)
#   math_idealistic  +.034 [+.005, +.062]   Pr  .992   (full data +.037)
#   read_idealistic  +.103 [+.078, +.127]   Pr 1.000   (full data +.104)
#
# Change in the mark coefficient, W5 to W8, with wave 9 absent:
#   math_realistic  +.009 [-.034, +.053]  Pr .663
#   read_realistic  -.004 [-.043, +.034]  Pr .424
#   math_idealistic +.024 [-.018, +.067]  Pr .873
#   read_idealistic +.000 [-.037, +.038]  Pr .512
#   Covers zero in all four, without the wave whose instrument changed.
#
# Per-wave mark coefficients are printed by the script and tabulated by
# 04_compare.R as the tbl-drop9 block.
# -----------------------------------------------------------------------------

