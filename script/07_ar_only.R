
source(here::here("script", "00_packages.R"))

model <- cmdstan_model(here("stan", "aum.stan"))
CELLS <- c("aum_math_realistic", "aum_read_realistic",
           "aum_math_idealistic", "aum_read_idealistic")

for (cl in CELLS) {
  d <- readRDS(here("data", "processed", paste0(cl, ".rds")))
  d$fielded_a <- rep(0L, d$n_wave)
  d$fielded_sc <- rep(0L, d$n_wave)
  d$fielded_x <- rep(0L, d$n_wave)

  message("fitting ", cl, "_ar")
  fit <- model$sample(
    data = c(d[c("n_person", "n_wave", "n_cov", "observed", "y", "v_obs",
                 "fielded_a", "fielded_sc", "fielded_x",
                 "obs_a", "obs_sc", "obs_x",
                 "s_a", "s_sc", "s_x", "z", "n_mod", "w")],
             list(model_id = 4L)),
    chains = 4, parallel_chains = 4,
    iter_warmup = 1000, iter_sampling = 1000, refresh = 0, seed = 1
  )
  fit$save_object(here("models", paste0(cl, "_ar.rds")))
  dg <- fit$diagnostic_summary(quiet = TRUE)
  message(sprintf("  %s_ar: %d divergent, max rhat %.3f", cl,
                  sum(dg$num_divergent),
                  max(fit$summary(c("phi", "sigma_eps"))$rhat, na.rm = TRUE)))
}

# -----------------------------------------------------------------------------
# RESULT (2026-08-03; realistic aspiration, both domains, 4 chains x 1000+1000)
#
# DOES NOT CONVERGE. NOT REPORTED. Recorded so the absence is on the record.
#
#   autoregressive coefficient, mathematics    rhat   bulk ESS
#     phi[1]  +1.01 [+0.22, +1.82]             1.00      9,033
#     phi[2]  +0.84 [+0.70, +1.16]             1.53          7
#     phi[3]  +0.87 [+0.80, +0.99]             1.52          7
#     phi[4]  +0.80 [+0.73, +0.95]             1.53          7
#     phi[5]  +0.64 [-0.11, +0.94]             1.53          7
#     phi[6]  -0.85 [-5.77, +0.75]             1.53          7
#     phi[7]  +0.97 [+0.81, +1.33]             1.53          7
#     phi[8]  +0.91 [+0.84, +1.02]             1.46          8
#
#   reading: phi[6] = -1.05 [-6.90, +0.86], rhat 1.53
#
# Eight free phi with no signals leaves the wave-9 transition unidentified;
# phi[6] trades against the innovation variance. Bulk ESS 7 of 4,000.
#
# The fits stay in models/ and 04_compare.R excludes them by name.
# -----------------------------------------------------------------------------

