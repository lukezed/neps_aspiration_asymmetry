
source(here::here("script", "00_packages.R"))

# Every cell, on whichever build the analysis is carrying. AUM_SUFFIX="" reruns
# the check against the four-covariate build for comparison.
SUFFIX <- Sys.getenv("AUM_SUFFIX", "_cov")
CELLS <- paste0(c("aum_math_realistic", "aum_read_realistic",
                  "aum_math_idealistic", "aum_read_idealistic"), SUFFIX)
model <- cmdstan_model(here("stan", "aum.stan"))

lag_mark <- function(d) {
  shift <- function(m) cbind(0, m[, -ncol(m), drop = FALSE])
  keep <- d$fielded_a                                  # same transitions as before
  d$s_a <- shift(d$s_a)
  d$obs_a <- shift(d$obs_a) * rep(keep, each = nrow(d$obs_a))
  d$s_a[d$obs_a == 0] <- 0
  # The change has to be lagged with the level it belongs to, otherwise the
  # transition would take last year's mark alongside this year's change.
  d$news_a <- shift(d$news_a)
  d$obs_news <- shift(d$obs_news) * rep(keep, each = nrow(d$obs_news))
  d$news_a[d$obs_news == 0] <- 0
  d
}

for (cell in CELLS) {
  d <- readRDS(here("data", "processed", paste0(cell, ".rds")))
  lagged <- lag_mark(d)
  message("\n=== ", cell, " ===")
  stopifnot(
    "lagging changed the set of achievement transitions" =
      identical(lagged$fielded_a, d$fielded_a),
    "the change score was not lagged with the level" =
      all(lagged$news_a[lagged$obs_news == 0] == 0) &&
      sum(lagged$obs_news) < sum(d$obs_news),
    "a lagged cell carries a value where it is not observed" =
      all(lagged$s_a[lagged$obs_a == 0] == 0),
    "lagging lost most of the signal" =
      sum(lagged$obs_a) > 0.6 * sum(d$obs_a)
  )
  message("marks observed: ", sum(d$obs_a), " concurrent, ", sum(lagged$obs_a), " lagged")

  fit <- model$sample(
    data = c(lagged[c("n_person", "n_wave", "n_cov", "observed", "y", "v_obs",
                      "fielded_a", "fielded_sc", "fielded_x",
                      "obs_a", "obs_sc", "obs_x",
                      "s_a", "s_sc", "s_x", "news_a", "obs_news",
                   "z", "n_mod", "w")],
             list(model_id = 5L)),
    chains = 4, parallel_chains = 4,
    iter_warmup = 1000, iter_sampling = 1000, refresh = 0, seed = 1
  )
  fit$save_object(here("models", paste0(cell, "_m5_lagged.rds")))

  dr <- as_draws_df(fit$draws(c("beta_a", "beta_sc", "b_pos", "b_neg")))
  ci <- function(x) sprintf("%+.3f [%+.3f, %+.3f]",
                            mean(x), quantile(x, .025), quantile(x, .975))
  keep <- which(d$fielded_a == 1)[-1]
  cat("\nmark coefficient, using the previous wave's mark\n")
  for (t in keep) {
    cat(sprintf("  W%-2d  %s\n", d$waves[t], ci(dr[[sprintf("beta_a[%d]", t)]])))
  }
  rise <- dr[[sprintf("beta_a[%d]", keep[length(keep)])]] - dr[[sprintf("beta_a[%d]", keep[1])]]
  cat(sprintf("  rise: %s, Pr(rise) = %.3f\n", ci(rise), mean(rise > 0)))

}

# The contrast the manuscript reports is recomputed from the saved fits by
# 04_compare.R against the concurrent ones, so it is not repeated here.

# -----------------------------------------------------------------------------
# RESULT (2026-08-23; four cells, model 5 refitted with achievement at t-1).
# 0 divergent, 0 max-treedepth, 16 chains, no rejected proposals in sampling.
#
# The three quantities the manuscript reports are printed by 04_compare.R as
# the tbl-lagged block. Under the lag, at Grade 12:
#   improved stays clear of zero in all four cells
#   slipped continues to cover zero in all four
#   the change in the gap stays clear of zero in mathematics, covers it in reading
#
# Mark coefficient by wave under the lag, read_idealistic:
#   W5 +.046  W6 +.031  W7 +.036  W8 +.030  W9 +.043
#   rise W5 to W9: -0.003 [-0.041, +0.037], Pr = 0.440
#
# Marks retained: 25,031 lagged against 28,166 concurrent (read_idealistic).
# -----------------------------------------------------------------------------

