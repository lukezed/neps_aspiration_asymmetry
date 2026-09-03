
source(here::here("script", "00_packages.R"))

# One cell of the design per call, so a run can be resumed or parallelised.
DATASET <- Sys.getenv("AUM_DATASET", "aum_math_realistic")
MODELS <- as.integer(strsplit(Sys.getenv("AUM_MODELS", "0,1,2"), ",")[[1]])

dat <- readRDS(here("data", "processed", paste0(DATASET, ".rds")))
model <- cmdstan_model(here("stan", "aum.stan"))

stan_input <- function(model_id) {
  c(dat[c("n_person", "n_wave", "n_cov", "observed", "y", "v_obs",
          "fielded_a", "fielded_sc", "fielded_x", "obs_a", "obs_sc", "obs_x",
          "s_a", "s_sc", "s_x", "news_a", "obs_news",
          "z", "n_mod", "w")],
    list(model_id = model_id))
}

for (model_id in MODELS) {
  label <- paste0(DATASET, "_m", model_id)
  message("fitting ", label)

  fit <- model$sample(
    data = stan_input(model_id),
    chains = 4, parallel_chains = 4,
    iter_warmup = 1000, iter_sampling = 1000,
    refresh = 100, seed = 1
  )

  # save_object() writes the draws to models/, which is gitignored. Without it
  # the draws live only in the session's temp directory and vanish with it.
  fit$save_object(here("models", paste0(label, ".rds")))

  diag <- fit$diagnostic_summary(quiet = TRUE)
  summ <- fit$summary(c("sigma_eps", "sigma_eta"))
  message(sprintf(
    "  %s: %d divergent, %d max-treedepth, max rhat %.3f",
    label, sum(diag$num_divergent), sum(diag$num_max_treedepth),
    max(summ$rhat, na.rm = TRUE)
  ))
}

cat("done; draws in models/\n")

# -----------------------------------------------------------------------------
# RESULT (2026-08-03)
#
# This script produces draws rather than numbers, so its result is the state of
# the fits it wrote. The convergence sweep that reads them lives at the foot of
# 04_compare.R, which is where the manuscript's one sentence on convergence
# comes from; running this script alone verifies nothing.
#
# 26 reported fits, four cells by six specifications plus the two diagnostics:
# 0 divergent transitions, 0 iterations at maximum treedepth, worst R-hat 1.008,
# smallest bulk ESS 2,119. The AR-only fits from 07_ar_only.R are the one
# non-convergence in models/ and are not reported.
#
# Note: a cell takes roughly 13 minutes per model on four chains. The full grid
# is run one cell per invocation through AUM_DATASET and AUM_MODELS so that a
# failure costs one fit rather than the grid.
# -----------------------------------------------------------------------------

