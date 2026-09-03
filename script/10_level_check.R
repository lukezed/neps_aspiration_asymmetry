
source(here::here("script", "00_packages.R"))

SUFFIX <- Sys.getenv("AUM_SUFFIX", "_cov")

CELLS <- sapply(c(math_realistic  = "aum_math_realistic",
                  read_realistic  = "aum_read_realistic",
                  math_idealistic = "aum_math_idealistic",
                  read_idealistic = "aum_read_idealistic"),
                function(x) paste0(x, SUFFIX))

ONET_SD <- sapply(CELLS, function(f)
  readRDS(here("data", "processed", paste0(f, ".rds")))$onet_sd)

ci <- function(x, d = 3) sprintf(paste0("%+.", d, "f [%+.", d, "f, %+.", d, "f]"),
                                 mean(x), quantile(x, .025), quantile(x, .975))

draws_of <- function(cell, model, pars)
  as_draws_df(readRDS(here("models", paste0(CELLS[[cell]], "_", model, ".rds")))$draws(pars))

# --- 1. does the asymmetry survive the interaction? ---------------------------

cat("\n=== b_pos - b_neg, with and without the level interaction ===\n")
cat(sprintf("%-16s %-26s %-26s %s\n", "cell", "m5 (reported)", "m6 (level freed)",
            "change in the gap"))
gaps <- list()
for (k in names(CELLS)) {
  d5 <- draws_of(k, "m5", c("b_pos", "b_neg"))
  d6 <- draws_of(k, "m6", c("b_pos", "b_neg"))
  g5 <- d5[["b_pos[1]"]] - d5[["b_neg[1]"]]
  g6 <- d6[["b_pos[1]"]] - d6[["b_neg[1]"]]
  gaps[[k]] <- list(m5 = g5, m6 = g6)
  cat(sprintf("%-16s %-26s %-26s %+.1f%%\n", k, ci(g5), ci(g6),
              100 * (mean(g6) - mean(g5)) / mean(g5)))
}

# --- 2. the interaction itself ------------------------------------------------

cat("\n=== b_lvl: does the response to a change depend on the level? ===\n")
cat(sprintf("%-16s %-26s %s\n", "cell", "b_lvl", "Pr(!= 0 in the reported sign)"))
for (k in names(CELLS)) {
  d <- draws_of(k, "m6", "b_lvl")
  b <- d[["b_lvl[1]"]]
  cat(sprintf("%-16s %-26s %.3f\n", k, ci(b), max(mean(b > 0), mean(b < 0))))
}

# --- 3. the two directions separately, in relatedness points ------------------

G9 <- 2L; G12 <- 6L
cat("\n=== Movement after a one-grade change, in relatedness points, model 6 ===\n")
cat(sprintf("%-16s %-6s %-26s %s\n", "cell", "grade", "improved", "slipped"))
for (k in names(CELLS)) {
  d <- draws_of(k, "m6", c("beta_a", "b_pos", "b_neg"))
  mv <- function(t, s) ONET_SD[[k]] * s *
    (d[[sprintf("beta_a[%d]", t)]] + d[[if (s > 0) "b_pos[1]" else "b_neg[1]"]])
  for (tt in c(G9, G12))
    cat(sprintf("%-16s %-6s %-26s %s\n", if (tt == G9) k else "",
                if (tt == G9) "9" else "12", ci(mv(tt, +1), 2), ci(mv(tt, -1), 2)))
}

# --- 4. convergence of the new fits -------------------------------------------

cat("\n=== Convergence, model 6 ===\n")
for (k in names(CELLS)) {
  f <- readRDS(here("models", paste0(CELLS[[k]], "_m6.rds")))
  dg <- f$diagnostic_summary(quiet = TRUE)
  s <- f$summary(c("b_pos", "b_neg", "b_lvl", "sigma_eps", "sigma_eta"))
  cat(sprintf("%-16s %d divergent, %d treedepth, max rhat %.3f, min ess %.0f\n",
              k, sum(dg$num_divergent), sum(dg$num_max_treedepth),
              max(s$rhat, na.rm = TRUE), min(s$ess_bulk, na.rm = TRUE)))
}

cat("\n")

# -----------------------------------------------------------------------------
# RESULT (re-run 2026-08-23; reproduces the 2026-08-13 run).
# Model 6: 0 divergent, 0 max-treedepth, max rhat 1.003, min bulk ESS 2222.
#
# b_pos - b_neg, level held constant (m5) and freed (m6):
#   math_realistic   +.080 [+.055, +.106]   +.091 [+.053, +.129]
#   read_realistic   +.135 [+.113, +.157]   +.157 [+.125, +.189]
#   math_idealistic  +.037 [+.012, +.063]   +.050 [+.014, +.086]
#   read_idealistic  +.104 [+.083, +.125]   +.134 [+.102, +.165]
#
# b_lvl, the change interacted with the mark it changed from:
#   math_realistic   +.007 [-.010, +.024]   Pr .774
#   read_realistic   +.014 [-.001, +.029]   Pr .969
#   math_idealistic  +.008 [-.009, +.024]   Pr .818
#   read_idealistic  +.019 [+.004, +.033]   Pr .991   <- the one that excludes zero
#
# The contrast survives in all four cells and grows once the level is freed.
# -----------------------------------------------------------------------------

