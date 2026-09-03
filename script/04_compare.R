
source(here::here("script", "00_packages.R"))

# The analysis carries the covariate version. AUM_SUFFIX="" reruns everything on
# the four-covariate build for comparison; nothing else in this script changes.
SUFFIX <- Sys.getenv("AUM_SUFFIX", "_cov")

CELLS <- sapply(c(math_realistic  = "aum_math_realistic",
                  read_realistic  = "aum_read_realistic",
                  math_idealistic = "aum_math_idealistic",
                  read_idealistic = "aum_read_idealistic"),
                function(x) paste0(x, SUFFIX))

# Wave-4 SD of the outcome, the divisor 02_build.R standardised by. Read, not typed.
ONET_SD <- sapply(CELLS, function(f)
  readRDS(here("data", "processed", paste0(f, ".rds")))$onet_sd)
stopifnot("outcome SD outside the 0-100 O*NET scale" = all(ONET_SD > 5 & ONET_SD < 30))

MARK_SD <- sapply(CELLS, function(f)
  readRDS(here("data", "processed", paste0(f, ".rds")))$anchoring["sd", "s_a"])
stopifnot("mark SD outside the 1-6 grade scale" = all(MARK_SD > 0.3 & MARK_SD < 2))

WAVES <- c(4:9, 11, 13)
G9 <- 2L                                    # wave index of Grade 9
G12 <- 6L                                   # wave index of Grade 12

fmt3 <- function(x) sub("^(-?)0\\.", "\\1.", sprintf("%.3f", x))

csv_row <- function(name, vals) paste0('"', name, '",',
                                     paste0('"', vals, '"', collapse = ","))

ci <- function(x, d = 3) sprintf(paste0("%+.", d, "f [%+.", d, "f, %+.", d, "f]"),
                                 mean(x), quantile(x, .025), quantile(x, .975))
cit <- function(x) sprintf("%s [%s, %s]", fmt3(mean(x)),
                           fmt3(quantile(x, .025)), fmt3(quantile(x, .975)))
draws_of <- function(cell, model, pars)
  as_draws_df(readRDS(here("models", paste0(CELLS[[cell]], "_", model, ".rds")))$draws(pars))

# --- 1. the handover -----------------------------------------------------------

cat("\n=== Coefficients on the mark and on self-concept, Grades 9 and 12 ===\n")
cat(sprintf("%-14s %-24s %-24s %-24s %s\n", "cell",
            "mark G9", "mark G12", "self-concept G9", "self-concept G12"))
for (k in names(CELLS)) {
  d <- draws_of(k, "m5", c("beta_a", "beta_sc"))
  g <- function(n, t) d[[sprintf("%s[%d]", n, t)]]
  cat(sprintf("%-14s %-24s %-24s %-24s %s\n", k,
              ci(g("beta_a", G9)), ci(g("beta_a", G12)),
              ci(g("beta_sc", G9)), ci(g("beta_sc", G12))))
}

cat("\n=== Change in each coefficient, Grade 9 to Grade 12 ===\n")
cat(sprintf("%-14s %-26s %-12s %-26s %s\n", "cell",
            "mark", "Pr(rise)", "self-concept", "Pr(fall)"))
for (k in names(CELLS)) {
  d <- draws_of(k, "m5", c("beta_a", "beta_sc"))
  ch <- function(n) d[[sprintf("%s[%d]", n, G12)]] - d[[sprintf("%s[%d]", n, G9)]]
  cat(sprintf("%-14s %-26s %-12.3f %-26s %.3f\n", k,
              ci(ch("beta_a")), mean(ch("beta_a") > 0),
              ci(ch("beta_sc")), mean(ch("beta_sc") < 0)))
}

cat("\n=== The handover: mark minus self-concept, and how it moves ===\n")
cat(sprintf("%-14s %-24s %-24s %-26s %s\n", "cell",
            "gap at G9", "gap at G12", "change in the gap", "Pr(gap grows)"))
for (k in names(CELLS)) {
  d <- draws_of(k, "m5", c("beta_a", "beta_sc"))
  gap <- function(t) d[[sprintf("beta_a[%d]", t)]] - d[[sprintf("beta_sc[%d]", t)]]
  ch <- gap(G12) - gap(G9)
  cat(sprintf("%-14s %-24s %-24s %-26s %.3f\n", k,
              ci(gap(G9)), ci(gap(G12)), ci(ch), mean(ch > 0)))
}

# --- 2. the asymmetry -----------------------------------------------------------

# Movement, not coefficient. Coefficient on a change is (beta_a + b); a slip is a
# change of -1, so its movement is that times -1. Printed unsigned it inverts.
cat("\n=== Movement of the named occupation, in O*NET points, for a mark",
    "one grade better and one grade worse than last year ===\n")
cat(sprintf("%-16s %-6s %-24s %-24s %s\n", "cell", "grade",
            "one grade better", "one grade worse", "b+ minus b-"))
for (k in names(CELLS)) {
  d <- draws_of(k, "m5", c("beta_a", "b_pos", "b_neg"))
  gap <- d[["b_pos[1]"]] - d[["b_neg[1]"]]
  for (t in c(G9, G12)) {
    ba <- d[[sprintf("beta_a[%d]", t)]]
    cat(sprintf("%-16s %-6s %-24s %-24s %s\n",
                if (t == G9) k else "", if (t == G9) "9" else "12",
                ci((ba + d[["b_pos[1]"]]) * ONET_SD[[k]] / MARK_SD[[k]] * (+1), 2),
                ci((ba + d[["b_neg[1]"]]) * ONET_SD[[k]] / MARK_SD[[k]] * (-1), 2),
                if (t == G9) sprintf("%s Pr = %.3f", ci(gap), mean(gap > 0)) else ""))
  }
}

# Three distinct quantities: (beta_a + b), its conversion to points, and b itself.

cat("\n=== Coefficient on the change in the mark, beta_a + b ===\n")
cat(sprintf("%-20s %-6s %-26s %-26s %s\n", "cell", "grade",
            "mark improved", "mark slipped", "slipped covers 0"))
for (k in names(CELLS)) {
  d <- draws_of(k, "m5", c("beta_a", "b_pos", "b_neg"))
  for (t in c(G9, G12)) {
    ba <- d[[sprintf("beta_a[%d]", t)]]
    dn <- ba + d[["b_neg[1]"]]
    cat(sprintf("%-20s %-6s %-26s %-26s %s\n",
                if (t == G9) k else "", if (t == G9) "9" else "12",
                ci(ba + d[["b_pos[1]"]]), ci(dn),
                if (quantile(dn, .025) < 0 && quantile(dn, .975) > 0) "yes" else "NO"))
  }
}

# --- 3. model comparison --------------------------------------------------------

# Signed asymmetric minus symmetric.
cat("\n=== Table, model comparison. Paste into paper.qmd's tbl-loso chunk ===\n")
cat('"","$\\Delta$elpd","$SE$","$\\Delta$elpd / $SE$"\n')
CELLNAME <- c(math_realistic = "Mathematics, realistic",
              read_realistic = "Reading, realistic",
              math_idealistic = "Mathematics, idealistic",
              read_idealistic = "Reading, idealistic")
kmax <- 0
for (k in names(CELLS)) {
  loos <- lapply(c("m7", "m5"), function(m)
    loo(readRDS(here("models", paste0(CELLS[[k]], "_", m, ".rds")))$draws("log_lik")))
  kmax <- max(kmax, sapply(loos, function(l) max(pareto_k_values(l))))
  cmp <- loo_compare(setNames(loos, c("m7", "m5")))
  d <- abs(cmp["m7", "elpd_diff"]); se <- cmp["m7", "se_diff"]
  cat(csv_row(CELLNAME[[k]], c(sprintf("%.1f", d), sprintf("%.1f", se),
                               sprintf("%.1f", d / se))), "\n")
}
cat(sprintf("\n  worst Pareto k across all eight fits: %.2f\n", kmax))

cat("\n=== Leave-one-student-out comparison ===\n")
for (k in names(CELLS)) {
  fits <- lapply(c("m4", "m7", "m5"), function(m)
    readRDS(here("models", paste0(CELLS[[k]], "_", m, ".rds"))))
  names(fits) <- c("no_change_term", "symmetric", "asymmetric")
  loos <- lapply(fits, function(f) loo(f$draws("log_lik")))
  cmp <- loo_compare(loos)
  cat("\n", k, "  (Pareto k > 0.7: ",
      paste(sapply(loos, function(l) sum(pareto_k_values(l) > 0.7)), collapse = "/"),
      ")\n", sep = "")
  print(round(cmp[, c("elpd_diff", "se_diff")], 1))
}

# --- 4. persistence -------------------------------------------------------------

cat("\n=== Autoregressive coefficient, Grades 9 and 12 ===\n")
for (k in names(CELLS)) {
  d <- draws_of(k, "m5", "phi")
  x <- d[[sprintf("phi[%d]", G12)]] - d[[sprintf("phi[%d]", G9)]]
  cat(sprintf("  %-14s G9 %.3f  G12 %.3f  change %s Pr = %.2f\n", k,
              mean(d[[sprintf("phi[%d]", G9)]]), mean(d[[sprintf("phi[%d]", G12)]]),
              ci(x), mean(x > 0)))
}

# --- 5. convergence -------------------------------------------------------------

cat("\n=== Convergence across every fit the paper reports ===\n")
# Reported fits only. _r2.rds holds a numeric vector, not a fit; excluded by pattern.
FITS <- grep("_cov(_m[456])(_drop9|_lagged)?\\.rds$",
             list.files(here("models"), "\\.rds$"), value = TRUE)
stopifnot("the reported-fit sweep lost fits" = length(FITS) == 20)
worst <- do.call(rbind, lapply(FITS, function(f) {
  fit <- readRDS(here("models", f))
  s <- fit$summary(setdiff(fit$metadata()$stan_variables, c("lp__", "log_lik")))
  d <- fit$diagnostic_summary(quiet = TRUE)
  data.frame(fit = sub("\\.rds$", "", f), div = sum(d$num_divergent),
             tree = sum(d$num_max_treedepth),
             rhat = max(s$rhat, na.rm = TRUE), ess = min(s$ess_bulk, na.rm = TRUE))
}))
reported <- worst
body <- reported[!grepl("_drop9|_lagged|_m6", reported$fit), ]
print(worst[order(-worst$rhat), ], row.names = FALSE, digits = 4)
cat(sprintf("\n  reported fits: %d (%d in the body, %d refits in Appendix C),\n",
            nrow(reported), nrow(body), nrow(reported) - nrow(body)))
cat(sprintf("  all with %d divergent and %d max-treedepth;\n",
            sum(reported$div), sum(reported$tree)))
cat(sprintf("  worst rhat %.3f, smallest bulk ESS %.0f (%.0f in the body)\n",
            max(reported$rhat), min(reported$ess), min(body$ess)))
stopifnot(
  "a reported fit failed a convergence diagnostic" =
    all(reported$div == 0) && all(reported$tree == 0) &&
    max(reported$rhat) < 1.01 && min(reported$ess) > 400
)

# -----------------------------------------------------------------------------
# RESULT (2026-08-23; _cov build, 4 chains x 1000+1000).
#
# 1. SIGNAL COEFFICIENTS, Grades 9 and 12.
#                     mark G9                mark G12               sc G9                  sc G12
#   math_realistic    +.051 [+.020, +.081]   +.080 [+.048, +.111]   +.080 [+.052, +.107]   -.014 [-.041, +.013]
#   read_realistic    +.095 [+.070, +.122]   +.093 [+.062, +.124]   +.057 [+.029, +.084]   -.011 [-.040, +.018]
#   math_idealistic   +.015 [-.012, +.043]   +.037 [+.007, +.068]   +.084 [+.060, +.109]   -.011 [-.038, +.015]
#   read_idealistic   +.072 [+.046, +.097]   +.061 [+.029, +.092]   +.033 [+.007, +.059]   -.002 [-.032, +.028]
#
# 2. CHANGE IN EACH COEFFICIENT, G9 to G12.
#                     mark                      Pr(rise)   self-concept              Pr(fall)
#   math_realistic    +.029 [-.015, +.072]      .906       -.094 [-.132, -.054]      1.000
#   read_realistic    -.003 [-.042, +.037]      .455       -.067 [-.106, -.026]       .999
#   math_idealistic   +.022 [-.018, +.062]      .851       -.096 [-.131, -.059]      1.000
#   read_idealistic   -.011 [-.052, +.028]      .302       -.036 [-.076, +.006]       .954
#
# 3. THE HANDOVER, mark minus self-concept.
#                     gap G9                   gap G12                  change                   Pr
#   math_realistic    -.028 [-.080, +.023]     +.094 [+.044, +.144]     +.122 [+.051, +.193]     1.000
#   read_realistic    +.039 [-.007, +.085]     +.104 [+.051, +.156]     +.065 [-.006, +.134]      .965
#   math_idealistic   -.069 [-.116, -.022]     +.049 [-.000, +.098]     +.118 [+.051, +.184]     1.000
#   read_idealistic   +.039 [-.007, +.083]     +.063 [+.009, +.117]     +.025 [-.047, +.096]      .749
#
# 4. MOVEMENT IN O*NET POINTS PER ONE GRADE, Grade 12 (corrected 2026-09-03;
#    previously divided only by the outcome SD, which reported the response per
#    standard deviation of the mark and called it a grade. The mark's Grade 8 SD
#    is .98 grades in maths and .83 in German, so reading was understated by a
#    fifth).
#                     one grade better        one grade worse         b+ minus b-
#   math_realistic    +1.60 [+1.09, +2.10]    -0.40 [-0.93, +0.14]    +.080 [+.055, +.106]  Pr 1.000
#   read_realistic    +2.39 [+1.87, +2.91]    -0.17 [-0.71, +0.38]    +.135 [+.113, +.157]  Pr 1.000
#   math_idealistic   +0.98 [+0.47, +1.50]    -0.39 [-0.94, +0.15]    +.037 [+.012, +.063]  Pr  .998
#   read_idealistic   +1.76 [+1.23, +2.27]    -0.10 [-0.63, +0.45]    +.104 [+.083, +.125]  Pr 1.000
#
# 5. MODEL COMPARISON, m7 (symmetric) against m5 (asymmetric).
#   math_realistic    17.4 (SE  6.2)  ratio 2.8
#   read_realistic    73.2 (SE 12.8)  ratio 5.7
#   math_idealistic    3.3 (SE  2.9)  ratio 1.1
#   read_idealistic   44.1 (SE  9.6)  ratio 4.6
#   Worst Pareto k: 0.24 over these eight fits, 0.26 including m4.
#
# 6. PERSISTENCE, phi at Grades 9 and 12.
#   math_realistic    .698 -> .690   -.008 [-.074, +.055]  Pr .41
#   read_realistic    .745 -> .780   +.035 [-.022, +.092]  Pr .89
#   math_idealistic   .738 -> .704   -.033 [-.097, +.030]  Pr .15
#   read_idealistic   .785 -> .821   +.037 [-.024, +.099]  Pr .88
#
# 7. CONVERGENCE. 20 reported fits, 8 in the body and 12 refits in Appendix C.
#   0 divergent and 0 max-treedepth everywhere.
#   Worst rhat 1.009. Smallest bulk ESS 1783, and 2654 in the body.
#
# --- 6. the coefficient table ---------------------------------------------------
# Model 5's population-level parameters, four cells, in the CSV form paper.qmd
# embeds. Unestimated rows are dropped: mark at waves 5-9, self-concept at 5 and 9,
# and wave 4 is the initial state so beta_a[1] and phi[1] only sample their priors.

GRADE <- c("Grade 8", "Grade 9", "Grade 9", "Grade 10", "Grade 11", "Grade 12",
           "2019", "2021")
lab <- function(t) sprintf("%s (W%d)", GRADE[t], WAVES[t])

rows <- c(
  lapply(2:6, function(t) list(lab(t), "beta_a", t)),
  lapply(c(2, 6), function(t) list(lab(t), "beta_sc", t)),
  list(list("Grade 9 (W5)",  "gap", 2),
       list("Grade 12 (W9)", "gap", 6)),
  # Movement, not coefficient, and not b_pos/b_neg, which invert if read alone.
  lapply(2:6, function(t) list(paste0(lab(t), ", improved"), "chg", c(t, +1))),
  lapply(2:6, function(t) list(paste0(lab(t), ", slipped"),  "chg", c(t, -1))),
  lapply(2:8, function(t) list(lab(t), "phi", t)),
  # sigma_nu belongs with the Grade 8 predictors below, not here.
  list(list("Aspiration report, $\\sigma_{\\varepsilon}$", "sigma_eps", NA),
       list("Latent state, wave to wave, $\\sigma_{\\eta}$", "sigma_eta", NA)),
  list(list("Intercept, $\\gamma_{0}$",                  "gamma_0", NA),
       list("Mathematics achievement, Grade 5 (W1)",           "gamma", 1),
       list("Self-concept, Grade 5 (W1)",               "gamma", 2),
       list("Domain interest, Grade 6 (W2)",            "gamma", 3),
       list("Parental status (HISEI)",                 "gamma", 4),
       list("Female",                                  "gamma", 5),
       list("Migration background",                    "gamma", 6),
       list("Gymnasium",                               "gamma", 7),
       list("Latent state at Grade 8, $\\sigma_{\\nu}$", "sigma_nu", NA))
)

cell_draws <- lapply(names(CELLS), function(k)
  draws_of(k, "m5", c("beta_a", "beta_sc", "phi", "b_pos", "b_neg",
                      "sigma_eps", "sigma_eta", "sigma_nu", "gamma_0", "gamma")))
names(cell_draws) <- names(CELLS)

est <- function(k, par, idx) {
  d <- cell_draws[[k]]
  v <- if (identical(par, "lit")) return(idx[[k]])
       else if (identical(par, "gap"))
         d[[sprintf("beta_a[%d]", idx)]] - d[[sprintf("beta_sc[%d]", idx)]]
       else if (identical(par, "chg"))
         idx[2] * (d[[sprintf("beta_a[%d]", idx[1])]] +
                   d[[if (idx[2] > 0) "b_pos[1]" else "b_neg[1]"]])
       else if (length(idx) == 1 && is.na(idx)) d[[par]]
       else d[[sprintf("%s[%d]", par, idx)]]
  sprintf("%s [%s, %s]", fmt3(mean(v)),
          fmt3(quantile(v, .025)), fmt3(quantile(v, .975)))
}

r2row <- sapply(names(CELLS), function(k) {
  f <- here("models", paste0(CELLS[[k]], "_m5_r2.rds"))
  if (!file.exists(f)) return("--")
  v <- readRDS(f)
  sprintf("%s [%s, %s]", fmt3(mean(v)), fmt3(quantile(v, .025)), fmt3(quantile(v, .975)))
})
rows <- c(rows, list(list("Bayesian $R^{2}$", "lit", r2row)))

cat("\n=== Table 3, coefficients. Paste into paper.qmd's tbl-coef chunk ===\n")
cat('"","Maths, realistic","Reading, realistic","Maths, idealistic","Reading, idealistic"\n')
for (r in rows)
  cat(csv_row(r[[1]], sapply(names(CELLS), function(k) est(k, r[[2]], r[[3]]))), "\n")

cat("\n  pack_rows: mark 1-5, self-concept 6-7, gap 8-9, movement 10-19,",
    "autoregression 20-26, residual SDs 27-28, Grade 8 starting point 29-37,",
    "model fit 38\n")

# --- 7. the two sensitivity refits ---------------------------------------------
# Fits come from 05_wave9_check.R and 06_timing_check.R.

cat("\n=== Table, wave-9 check. Paste into paper.qmd's tbl-drop9 chunk ===\n")
cat('"","Grade 9 (W5)","Grade 9 (W6)","Grade 10 (W7)","Grade 11 (W8)"\n')
CELLNAME <- c(math_realistic = "Mathematics, realistic",
              read_realistic = "Reading, realistic",
              math_idealistic = "Mathematics, idealistic",
              read_idealistic = "Reading, idealistic")
for (k in names(CELLS)) {
  full <- draws_of(k, "m5", "beta_a")
  d9 <- as_draws_df(readRDS(here("models",
          paste0(CELLS[[k]], "_m5_drop9.rds")))$draws("beta_a"))
  cat(csv_row(CELLNAME[[k]], rep("", 4)), "\n")
  cat(csv_row("All waves",
      sapply(2:5, function(t) cit(full[[sprintf("beta_a[%d]", t)]]))), "\n")
  cat(csv_row("Wave 9 dropped",
      sapply(2:5, function(t) cit(d9[[sprintf("beta_a[%d]", t)]]))), "\n")
}

cat("\n=== Table, lagged check. Paste into paper.qmd's tbl-lagged chunk ===\n")
cat('"","Mark improved","Mark slipped","Change in the gap"\n')
for (k in names(CELLS)) {
  lg <- function(par) as_draws_df(readRDS(here("models",
          paste0(CELLS[[k]], "_m5_lagged.rds")))$draws(par))
  gap <- function(d) (d[["beta_a[6]"]] - d[["beta_sc[6]"]]) -
                     (d[["beta_a[2]"]] - d[["beta_sc[2]"]])
  row <- function(get) {
    a <- get(draws_of(k, "m5", c("beta_a", "beta_sc", "b_pos", "b_neg")))
    a
  }
  ac <- draws_of(k, "m5", c("beta_a", "beta_sc", "b_pos", "b_neg"))
  al <- lg(c("beta_a", "beta_sc", "b_pos", "b_neg"))
  cat(csv_row(CELLNAME[[k]], rep("", 3)), "\n")
  # Movement, so the slipped column is negated, matching Table 3 and Figure 3.
  cat(csv_row("Concurrent achievement", c(
      cit(  ac[["beta_a[6]"]] + ac[["b_pos[1]"]]),
      cit(-(ac[["beta_a[6]"]] + ac[["b_neg[1]"]])),
      cit(gap(ac)))), "\n")
  cat(csv_row("Previous wave's achievement", c(
      cit(  al[["beta_a[6]"]] + al[["b_pos[1]"]]),
      cit(-(al[["beta_a[6]"]] + al[["b_neg[1]"]])),
      cit(gap(al)))), "\n")
}

cat("\n=== wave-9 check: achievement coefficient, interior waves ===\n")
GRADE_AT <- c("G8", "G9", "G9", "G10", "G11", "G12")
cat(sprintf("%-24s %-6s %-26s %s\n", "cell", "wave", "all waves", "wave 9 dropped"))
for (k in names(CELLS)) {
  full <- draws_of(k, "m5", "beta_a")
  d9 <- as_draws_df(readRDS(here("models",
          paste0(CELLS[[k]], "_m5_drop9.rds")))$draws("beta_a"))
  for (t in 2:5)
    cat(sprintf("%-24s %-6s %-26s %s\n", if (t == 2) k else "", GRADE_AT[t],
                ci(full[[sprintf("beta_a[%d]", t)]]),
                ci(d9[[sprintf("beta_a[%d]", t)]])))
  cat(sprintf("%-24s %-6s %-26s %s\n", "", "G12",
              ci(full[["beta_a[6]"]]), "not estimated"))
}

cat("\n=== lagged check: response to a change at Grade 12, and the handover ===\n")
cat(sprintf("%-24s %-12s %-26s %s\n", "cell", "quantity",
            "concurrent", "previous wave"))
for (k in names(CELLS)) {
  lg <- function(par) as_draws_df(readRDS(here("models",
          paste0(CELLS[[k]], "_m5_lagged.rds")))$draws(par))
  for (b in c("b_pos", "b_neg")) {
    a <- draws_of(k, "m5", c("beta_a", b))
    l <- lg(c("beta_a", b))
    cat(sprintf("%-24s %-12s %-26s %s\n",
                if (b == "b_pos") k else "",
                if (b == "b_pos") "improving" else "slipping",
                ci(a[["beta_a[6]"]] + a[[paste0(b, "[1]")]]),
                ci(l[["beta_a[6]"]] + l[[paste0(b, "[1]")]])))
  }
  gap <- function(d) (d[["beta_a[6]"]] - d[["beta_sc[6]"]]) -
                     (d[["beta_a[2]"]] - d[["beta_sc[2]"]])
  cat(sprintf("%-24s %-12s %-26s %s\n", "", "handover",
              ci(gap(draws_of(k, "m5", c("beta_a", "beta_sc")))),
              ci(gap(lg(c("beta_a", "beta_sc"))))))
}

