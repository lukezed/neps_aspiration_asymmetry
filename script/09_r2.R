
source(here::here("script", "00_packages.R"))

SUFFIX <- Sys.getenv("AUM_SUFFIX", "_cov")
CELLS <- paste0(c("aum_math_realistic", "aum_read_realistic",
                  "aum_math_idealistic", "aum_read_idealistic"), SUFFIX)
MODELS <- c("m4", "m5")

model <- cmdstan_model(here("stan", "aum.stan"))

stan_input <- function(dat, model_id) {
  c(dat[c("n_person", "n_wave", "n_cov", "observed", "y", "v_obs",
          "fielded_a", "fielded_sc", "fielded_x", "obs_a", "obs_sc", "obs_x",
          "s_a", "s_sc", "s_x", "news_a", "obs_news",
          "z", "n_mod", "w")],
    list(model_id = model_id))
}

for (cell in CELLS) {
  dat <- readRDS(here("data", "processed", paste0(cell, ".rds")))
  for (m in MODELS) {
    label <- paste0(cell, "_", m)
    fit <- readRDS(here("models", paste0(label, ".rds")))
    gq <- model$generate_quantities(
      fitted_params = fit$draws(),
      data = stan_input(dat, as.integer(sub("m", "", m))),
      seed = 1)
    r2 <- as.numeric(gq$draws("bayes_r2"))
    saveRDS(r2, here("models", paste0(label, "_r2.rds")))
    message(sprintf("%-32s R2 = %.3f [%.3f, %.3f]", label,
                    mean(r2), quantile(r2, .025), quantile(r2, .975)))
  }
}

# -----------------------------------------------------------------------------
# RESULT (2026-08-23). Bayesian R-squared, model 5, on the one-step-ahead
# predictive mean. Identical to three decimals under model 4.
#
#   math_realistic   .195 [.187, .202]
#   read_realistic   .301 [.294, .309]
#   math_idealistic  .197 [.190, .203]
#   read_idealistic  .247 [.240, .254]
# -----------------------------------------------------------------------------

