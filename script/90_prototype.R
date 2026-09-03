
source(here::here("script", "00_packages.R"))

n_person <- 6864
n_wave <- 8
fielded_a <- c(1, 1, 1, 1, 1, 1, 0, 0)
fielded_sc <- c(0, 1, 0, 0, 0, 1, 0, 0)

# Share of students whose signal is observed where it was fielded. The real
# panel runs roughly 0.4 to 0.8 depending on wave.
p_obs_signal <- 0.6
p_obs_aspiration <- 0.65

# True values. Achievement permeability declines across the six signal waves,
# which is the developmental hypothesis the design has to be able to detect.
true <- list(
  lambda_a = c(0.28, 0.24, 0.20, 0.16, 0.12, 0.10, 0, 0),
  lambda_sc = c(0, 0.18, 0, 0, 0, 0.20, 0, 0),
  sigma_nu = 0.9,
  sigma_eps = 0.5,
  sigma_eta = 0.3,
  sigma_a = 1.0,
  gamma_0 = 0
)

set.seed(20260802)

# Moderators of the updating rate. Each needs a directional prediction to earn
# its place: gender, socioeconomic status, migration background.
w <- cbind(
  female = rbinom(n_person, 1, 0.5),
  ses = rnorm(n_person),
  migrant = rbinom(n_person, 1, 0.25)
)
z <- matrix(rnorm(n_person), ncol = 1)                 # one initial-state covariate

# Signals depend on the covariate, which is what makes the auxiliary model do
# real work: the imputed mean differs across students rather than being one number.
s_a <- 0.5 * z[, 1] + matrix(rnorm(n_person * n_wave, 0, true$sigma_a), n_person, n_wave)
SC_R <- 0.65
s_sc <- 0.5 * z[, 1] + SC_R * (s_a - 0.5 * z[, 1]) +
  sqrt(1 - SC_R^2) * matrix(rnorm(n_person * n_wave, 0, true$sigma_a), n_person, n_wave)

CROSS_R <- 0.42
s_x <- CROSS_R * s_a +
  sqrt(1 - CROSS_R^2) * matrix(rnorm(n_person * n_wave, 0, true$sigma_a), n_person, n_wave)

b_star <- matrix(NA_real_, n_person, n_wave)
b_star[, 1] <- true$gamma_0 + 0.4 * z[, 1] + rnorm(n_person, 0, true$sigma_nu)
for (t in 2:n_wave) {
  la <- true$lambda_a[t] * fielded_a[t]
  ls <- true$lambda_sc[t] * fielded_sc[t]
  b_star[, t] <- (1 - la - ls) * b_star[, t - 1] +
    la * s_a[, t] + ls * s_sc[, t] +
    rnorm(n_person, 0, true$sigma_eta)
}

y <- b_star + matrix(rnorm(n_person * n_wave, 0, true$sigma_eps), n_person, n_wave)

# Missingness is intermittent rather than monotone, as in the real panel: 3,797
# of 7,146 students have interior gaps.
observed <- matrix(rbinom(n_person * n_wave, 1, p_obs_aspiration), n_person, n_wave)
obs_a <- matrix(rbinom(n_person * n_wave, 1, p_obs_signal), n_person, n_wave) *
  matrix(rep(fielded_a, each = n_person), n_person, n_wave)
obs_sc <- matrix(rbinom(n_person * n_wave, 1, 0.887), n_person, n_wave) * obs_a *
  matrix(rep(fielded_sc, each = n_person), n_person, n_wave)
obs_x <- obs_a

obs_news <- obs_a * cbind(0L, obs_a[, -n_wave, drop = FALSE])
news_a <- (s_a - cbind(0, s_a[, -n_wave, drop = FALSE])) * obs_news

# The generating values stay in b_star; what the model sees is blanked out.
y[observed == 0] <- 0
s_a_obs <- s_a; s_a_obs[obs_a == 0] <- 0
s_sc_obs <- s_sc; s_sc_obs[obs_sc == 0] <- 0
s_x_obs <- s_x; s_x_obs[obs_x == 0] <- 0

stan_data <- list(
  n_person = n_person,
  n_wave = n_wave,
  n_cov = ncol(z),
  observed = observed,
  y = y,
  # Zero throughout, which is what the KldB outcome supplies: one occupation
  # code maps to one O*NET occupation, so no averaging error is carried.
  v_obs = matrix(0, n_person, n_wave),
  fielded_a = fielded_a,
  fielded_sc = fielded_sc,
  fielded_x = fielded_a,
  obs_a = obs_a,
  obs_sc = obs_sc,
  obs_x = obs_x,
  s_a = s_a_obs,
  s_sc = s_sc_obs,
  s_x = s_x_obs,
  news_a = news_a,
  obs_news = obs_news,
  z = z,
  n_mod = ncol(w),
  w = w,
  model_id = 2L
)

model <- cmdstan_model(here("stan", "aum.stan"))

fit <- model$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 500,
  refresh = 250,
  seed = 1
)

# --- self-check --------------------------------------------------------------
# The prototype is only useful if it recovers what generated the data. Fail loudly.

est <- fit$summary(c("lambda_a", "lambda_sc", "sigma_eps", "sigma_eta",
                     "sigma_a", "delta"))
print(est, n = 40)

lam_a_hat <- est$mean[grepl("^lambda_a", est$variable)]
lam_sc_hat <- est$mean[grepl("^lambda_sc", est$variable)]
delta_hat <- est$mean[grepl("^delta", est$variable)]
signal_waves <- which(fielded_a == 1)[-1]              # wave 1 has no transition
sc_waves <- which(fielded_sc == 1)

stopifnot(
  "lambda_a not recovered" =
    max(abs(lam_a_hat[signal_waves] - true$lambda_a[signal_waves])) < 0.06,
  "lambda_sc not recovered" =
    max(abs(lam_sc_hat[sc_waves] - true$lambda_sc[sc_waves])) < 0.06,
  "declining trend not recovered" =
    cor(signal_waves, lam_a_hat[signal_waves]) < -0.8,
  # Data are generated with no moderator effects, so delta must recover zero.
  "spurious moderator effect" = max(abs(delta_hat)) < 0.15,
  "signal dispersion not recovered" =
    abs(est$mean[est$variable == "sigma_a"] - true$sigma_a) < 0.1,
  "sampler unhealthy" = max(est$rhat, na.rm = TRUE) < 1.01
)

cat("\nrecovery OK with", round(100 * (1 - p_obs_signal)), "% of fielded signals missing;",
    "max |bias| in lambda_a =",
    round(max(abs(lam_a_hat[signal_waves] - true$lambda_a[signal_waves])), 3), "\n")

# Per-student marginal log likelihood is what leave-one-student-out CV runs on.
loo_out <- loo(fit$draws("log_lik"))
print(loo_out)

for (id in c(0L, 1L)) {
  smoke <- model$sample(
    data = modifyList(stan_data, list(model_id = id)),
    chains = 2, parallel_chains = 2,
    iter_warmup = 200, iter_sampling = 200,
    refresh = 0, seed = 1, show_messages = FALSE
  )
  ll <- smoke$summary("log_lik")
  stopifnot(
    "model returned non-finite log_lik" = all(is.finite(ll$mean)),
    "model failed to converge in smoke test" = max(ll$rhat, na.rm = TRUE) < 1.05
  )
  cat("model_id", id, "OK; elpd_loo =",
      round(loo(smoke$draws("log_lik"))$estimates["elpd_loo", "Estimate"], 1), "\n")
}

cross_fit <- function(data, label) {
  f <- model$sample(
    data = data, chains = 4, parallel_chains = 4,
    iter_warmup = 500, iter_sampling = 500, refresh = 0, seed = 1
  )
  dr <- as_draws_df(f$draws("beta_x"))
  ts <- which(fielded_a == 1)[-1]
  m <- sapply(ts, function(t) mean(dr[[sprintf("beta_x[%d]", t)]]))
  lo <- sapply(ts, function(t) quantile(dr[[sprintf("beta_x[%d]", t)]], 0.025))
  hi <- sapply(ts, function(t) quantile(dr[[sprintf("beta_x[%d]", t)]], 0.975))
  cat("\n", label, ": beta_x at the six signal waves\n", sep = "")
  for (k in seq_along(ts)) {
    cat(sprintf("  wave %d  %+.3f [%+.3f, %+.3f]\n", ts[k], m[k], lo[k], hi[k]))
  }
  stopifnot("beta_x summary is empty" = length(m) == length(ts) && all(is.finite(m)))
  list(fit = f, est = m, lo = lo, hi = hi)
}

null_case <- cross_fit(modifyList(stan_data, list(model_id = 3L)), "no true effect")

stopifnot(
  "cross-domain effect invented where none exists" =
    max(abs(null_case$est)) < 0.05,
  "cross-domain interval excludes zero where the effect is zero" =
    all(null_case$lo < 0 & null_case$hi > 0)
)

cat("\nmodel 3 OK: no false positive at r =", CROSS_R, "\n")

PHI <- 0.75
BETA <- 0.06
B_SYM <- 0.08

b5 <- matrix(NA_real_, n_person, n_wave)
b5[, 1] <- b_star[, 1]
for (t in 2:n_wave) {
  b5[, t] <- PHI * b5[, t - 1] +
    fielded_a[t] * (BETA * s_a[, t] + B_SYM * news_a[, t]) +
    fielded_sc[t] * BETA * s_sc[, t] +
    rnorm(n_person, 0, true$sigma_eta)
}
y5 <- b5 + matrix(rnorm(n_person * n_wave, 0, true$sigma_eps), n_person, n_wave)
y5[observed == 0] <- 0

fit5 <- model$sample(
  data = modifyList(stan_data, list(model_id = 5L, y = y5)),
  chains = 4, parallel_chains = 4,
  iter_warmup = 500, iter_sampling = 500, refresh = 0, seed = 1
)

d5 <- as_draws_df(fit5$draws(c("b_pos", "b_neg")))
gap <- d5[["b_pos[1]"]] - d5[["b_neg[1]"]]
cat(sprintf("\nsymmetric truth (b = %.2f both directions): b_pos - b_neg = %+.3f [%+.3f, %+.3f]\n",
            B_SYM, mean(gap), quantile(gap, .025), quantile(gap, .975)))
cat(sprintf("  b_pos %+.3f   b_neg %+.3f\n",
            mean(d5[["b_pos[1]"]]), mean(d5[["b_neg[1]"]])))

stopifnot(
  "asymmetry invented where the response is symmetric" = abs(mean(gap)) < 0.03,
  "asymmetry interval excludes zero where there is none" =
    quantile(gap, .025) < 0 && quantile(gap, .975) > 0,
  "the symmetric coefficient is not recovered" =
    abs(mean(d5[["b_pos[1]"]]) - B_SYM) < 0.03,
  "the sampler did not converge" =
    max(fit5$summary(c("phi", "b_pos", "b_neg"))$rhat, na.rm = TRUE) < 1.01
)
cat("\nmodel 5 OK: no asymmetry invented at a true coefficient of", B_SYM, "\n")

real <- readRDS(here("data", "processed", "aum_math_realistic.rds"))
grid <- sort(unique(real$y[real$observed == 1]))
real_y <- real$y; real_y[real$observed == 0] <- NA
real_zero <- mean(diff(t(real_y)) == 0, na.rm = TRUE)
STICKY <- 0.365            # share of pairs naming the identical occupation

snap <- function(v) grid[findInterval(v, (grid[-1] + grid[-length(grid)]) / 2) + 1L]

y5r <- matrix(snap(as.vector(y5)), n_person, n_wave)
for (t in 2:n_wave) {
  same <- runif(n_person) < STICKY
  y5r[same, t] <- y5r[same, t - 1]
}
zero_change <- mean(diff(t(y5r)) == 0)
cat(sprintf("\nsimulated outcome: %d attainable values, %.1f%% of adjacent pairs unchanged\n",
            length(grid), 100 * zero_change))
cat(sprintf("  real outcome:    %d attainable values, %.1f%% unchanged (same basis)\n",
            length(grid), 100 * real_zero))
y5r[observed == 0] <- 0

fit5r <- model$sample(
  data = modifyList(stan_data, list(model_id = 5L, y = y5r)),
  chains = 4, parallel_chains = 4,
  iter_warmup = 500, iter_sampling = 500, refresh = 0, seed = 1
)

d5r <- as_draws_df(fit5r$draws(c("b_pos", "b_neg")))
gap_r <- d5r[["b_pos[1]"]] - d5r[["b_neg[1]"]]
cat(sprintf("\nsymmetric truth on a bounded lumpy outcome: b_pos - b_neg = %+.3f [%+.3f, %+.3f]\n",
            mean(gap_r), quantile(gap_r, .025), quantile(gap_r, .975)))
cat(sprintf("  b_pos %+.3f   b_neg %+.3f\n",
            mean(d5r[["b_pos[1]"]]), mean(d5r[["b_neg[1]"]])))
cat(sprintf("  for comparison, the asymmetry reported in the real data is +.083 [+.057, +.110]\n"))

stopifnot(
  # The simulated outcome has to be at least as lumpy as the real one, or the
  # check is easier than the objection it answers.
  "the simulated outcome is not lumpy enough to be the test it claims to be" =
    zero_change > real_zero - 0.05,
  "asymmetry invented by a bounded lumpy outcome" = abs(mean(gap_r)) < 0.03,
  "asymmetry interval excludes zero on a symmetric truth" =
    quantile(gap_r, .025) < 0 && quantile(gap_r, .975) > 0,
  "the sampler did not converge" =
    max(fit5r$summary(c("phi", "b_pos", "b_neg"))$rhat, na.rm = TRUE) < 1.01
)
cat("\nmodel 5 OK: a bounded lumpy outcome does not manufacture the asymmetry\n")

# -----------------------------------------------------------------------------
# RESULT (re-run 2026-08-23; n = 6,864 simulated students, seed 1, reproduces
# the 2026-08-03 run exactly). Supports Appendix B, Recovery on Simulated Data.
#
# Model 2, parameter recovery
#   max |bias| in lambda_a  0.032  against a tolerance of 0.06, with 40% of
#   fielded signals missing; lambda_sc recovered; delta recovers zero
#
# Model 3, cross-domain coefficient, generated with no cross-domain effect and
# the two marks correlated 0.42:
#   +0.005 [-0.017,+0.027]   -0.010 [-0.030,+0.010]   -0.019 [-0.039,+0.003]
#   -0.004 [-0.026,+0.018]   +0.020 [-0.000,+0.040]      every interval covers 0
#
# Model 5, symmetric truth, b = 0.080 in both directions:
#   b_pos - b_neg = -0.007 [-0.021, +0.006]
#   b_pos +0.070   b_neg +0.078
#
# Model 5, same truth on an outcome snapped to the 32 attainable scores:
#   b_pos - b_neg = -0.009 [-0.023, +0.004]
#   b_pos +0.023   b_neg +0.033
#   simulated 48.4% of adjacent pairs unchanged against a real 51.2%
#
# Reading: the interval covers zero and its upper end, +.004, is an order of
# magnitude below the +.083 [+.057, +.110] reported in the real data. A bounded
# outcome that moves in steps of about three points, with nearly half of
# consecutive reports not moving at all, does not manufacture a directional
# response. This is the check the experimental belief-updating literature was
# faulted for lacking, and it is the one a reviewer will ask for, because the
# argument that undid that literature was precisely that its paradigm produced
# apparent asymmetry from a symmetric process.
#
# Note: discretisation attenuates. Both coefficients come back at roughly a
# third of the value that generated them, 0.023 and 0.033 against 0.080, because
# a response smaller than one step of the scale cannot be expressed. That is
# worth knowing when reading the size of any coefficient in this paper. It does
# not license a claim that the true asymmetry is larger than the estimated one:
# the truth here is a difference of zero, so this run measures the attenuation
# of the levels and says nothing about the attenuation of a non-zero difference.
# -----------------------------------------------------------------------------

