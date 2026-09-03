
source(here::here("script", "00_packages.R"))

dta <- here("data", "SC3_D_13-1-0", "Stata14")
suf <- function(x) file.path(dta, paste0("SC3_", x, "_D_13-1-0.dta"))

ASPIRATION_WAVES <- c(4:9, 11, 13)
SIGNAL_WAVES_A <- 4:9                                  # marks co-observed with an aspiration

SIGNAL_WAVES_SC <- c(5, 9)
INITIAL_WAVE <- 4

WITH_COV <- nzchar(Sys.getenv("AUM_COV"))
SUFFIX <- if (WITH_COV) "_cov" else ""

# Variable names by domain, so German comes for free.
vars_for <- function(domain) {
  stopifnot(domain %in% c("math", "read"))
  if (domain == "math") {
    list(aspiration_col = "math",
         mark_annual = "t724102", mark_half = "t724112", mark_points = "t724602",
         self_concept = "t66001a_g1", interest = "t66201a_g1")
  } else {
    list(aspiration_col = "read",
         mark_annual = "t724101", mark_half = "t724111", mark_points = "t724601",
         self_concept = "t66000a_g1", interest = "t66208a_g1")
  }
}

anchor_const <- function(x, wave, anchor = INITIAL_WAVE) {
  if (all(is.na(x[wave == anchor]))) {
    available <- sort(unique(wave[!is.na(x)]))
    stopifnot("variable is empty in every wave" = length(available) > 0)
    anchor <- available[1]
  }
  ref <- x[wave == anchor]
  s <- sd(ref, na.rm = TRUE)
  # A fallback anchor with one observation gives an NA divisor and turns the
  # whole variable into NA, which the blanking downstream then hides as zeros.
  stopifnot("anchor wave has no usable dispersion" = is.finite(s) && s > 0)
  c(mean = mean(ref, na.rm = TRUE), sd = s, wave = anchor)
}

anchor_scale <- function(x, wave, anchor = INITIAL_WAVE) {
  k <- anchor_const(x, wave, anchor)
  (x - k[["mean"]]) / k[["sd"]]
}

build_domain <- function(domain,
                         aspiration = c("realistic", "idealistic"),
                         coding = c("kldb", "isco")) {
  aspiration <- match.arg(aspiration)
  coding <- match.arg(coding)
  v <- vars_for(domain)
  # The other domain's mark, which model 3 admits into this domain's transition.
  vx <- vars_for(if (domain == "math") "read" else "math")

  suffix <- if (coding == "kldb") "_g2" else "_g4"
  asp_var <- paste0(if (aspiration == "realistic") "t31160a" else "t31060a", suffix)

  onet <- readRDS(here("data", "processed", paste0("onet_", coding, ".rds")))
  code_col <- if (coding == "kldb") "kldb" else "isco08"

  target <- read_dta(suf("pTarget"), col_select = all_of(c(
    "ID_t", "wave", "tx20100", asp_var, "t700031",
    v$mark_annual, v$mark_half, v$mark_points, v$self_concept, v$interest,
    vx$mark_annual, vx$mark_half, vx$mark_points
  ))) %>%
    filter(tx20100 == 1)

  # --- outcome: occupation code to domain importance --------------------------
  target$isco <- neps_na(target[[asp_var]])
  target$b_obs <- onet[[v$aspiration_col]][match(target$isco, onet[[code_col]])]
  # Measurement error on that score, in the same raw units, from the spread of
  # O*NET occupations inside the ISCO code the answer was assigned to.
  target$b_sd <- onet[[paste0(v$aspiration_col, "_sd")]][match(target$isco, onet[[code_col]])]

  mark_signal <- function(vv) -pmin(
    neps_na(target[[vv$mark_annual]]),
    neps_na(target[[vv$mark_half]]),
    points_to_note(neps_na(target[[vv$mark_points]])),
    na.rm = TRUE
  )

  target$interest_val <- neps_na(target[[v$interest]])
  target$s_a <- mark_signal(v)
  target$s_x <- mark_signal(vx)
  target$s_sc <- neps_na(target[[v$self_concept]])

  b_scale <- sd(target$b_obs[target$wave == INITIAL_WAVE], na.rm = TRUE)
  target$v_obs <- (target$b_sd / b_scale)^2

  # Captured before the transform overwrites the columns it was computed from.
  anchoring <- sapply(c("b_obs", "s_a", "s_x", "s_sc"),
                      function(v) anchor_const(target[[v]], target$wave))

  target <- target %>%
    mutate(across(c(b_obs, s_a, s_x, s_sc), ~ anchor_scale(.x, wave)))

  # --- covariates -------------------------------------------------------------
  parent <- read_dta(suf("pParent"),
                     col_select = c(ID_t, p731904_g14, p731954_g14, p400500_g1))
  hisei <- parent %>%
    mutate(across(c(p731904_g14, p731954_g14), neps_na),
           hisei = pmax(p731904_g14, p731954_g14, na.rm = TRUE),
           migrant = as.integer(neps_na(p400500_g1) > 0)) %>%
    group_by(ID_t) %>%
    summarise(hisei = suppressWarnings(max(hisei, na.rm = TRUE)),
              migrant = suppressWarnings(max(migrant, na.rm = TRUE)),
              .groups = "drop") %>%
    mutate(across(c(hisei, migrant), ~ if_else(is.finite(.x), .x, NA_real_)))

  track <- read_dta(suf("CohortProfile"),
                    col_select = c(ID_t, wave, t723080_g1)) %>%
    filter(wave <= INITIAL_WAVE) %>%
    mutate(v = neps_na(t723080_g1)) %>%
    filter(!is.na(v)) %>%
    arrange(ID_t, wave) %>%
    group_by(ID_t) %>%
    summarise(gym = as.numeric(first(v) %in% c(9, 16)), .groups = "drop")

  # first() on an all-NA vector errors or returns NA depending on type.
  first_obs <- function(x) if (all(is.na(x))) NA_real_ else x[!is.na(x)][1]

  modal_obs <- function(x) {
    x <- x[!is.na(x)]
    if (!length(x)) NA_real_ else as.numeric(names(sort(table(x), decreasing = TRUE))[1])
  }

  person <- target %>%
    group_by(ID_t) %>%
    summarise(
      female = as.integer(modal_obs(neps_na(t700031)) == 2),
      # Initial-state predictors, all measured before the first aspiration wave.
      g5_mark = first_obs(s_a[wave == 1]),
      w1_sc = first_obs(s_sc[wave == 1]),
      w2_interest = first_obs(interest_val[wave == 2]),
      .groups = "drop"
    ) %>%
    left_join(hisei, by = "ID_t") %>%
    left_join(track, by = "ID_t")

  # --- keep students with at least one usable aspiration -----------------------
  keep <- target %>%
    filter(wave %in% ASPIRATION_WAVES, !is.na(b_obs)) %>%
    distinct(ID_t) %>%
    pull(ID_t)

  target <- filter(target, ID_t %in% keep)
  person <- filter(person, ID_t %in% keep) %>% arrange(ID_t)

  # --- reshape to person x wave matrices --------------------------------------
  to_matrix <- function(value_col, waves) {
    wide <- target %>%
      filter(wave %in% waves) %>%
      select(ID_t, wave, value = all_of(value_col)) %>%
      group_by(ID_t, wave) %>%
      summarise(value = first_obs(value), .groups = "drop") %>%
      complete(ID_t = person$ID_t, wave = waves) %>%
      arrange(ID_t, wave) %>%
      pivot_wider(names_from = wave, values_from = value)
    as.matrix(wide[match(person$ID_t, wide$ID_t), -1])
  }

  fill_mean <- function(x) {
    x <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
    replace_na(x, 0)
  }
  n_missing_cov <- person %>%
    summarise(across(c(g5_mark, w1_sc, w2_interest, hisei, female, migrant, gym),
                     ~ mean(is.na(.x))))

  fill_prop <- function(x) replace_na(x, mean(x, na.rm = TRUE))

  person <- person %>%
    mutate(across(c(g5_mark, w1_sc, w2_interest, hisei), fill_mean),
           across(c(female, migrant, gym), fill_prop))

  b <- to_matrix("b_obs", ASPIRATION_WAVES)
  a <- to_matrix("s_a", ASPIRATION_WAVES)
  x <- to_matrix("s_x", ASPIRATION_WAVES)
  sc <- to_matrix("s_sc", ASPIRATION_WAVES)
  # Not `v`: that name holds the variable list this function was built from.
  vo <- to_matrix("v_obs", ASPIRATION_WAVES)

  fielded_a <- as.integer(ASPIRATION_WAVES %in% SIGNAL_WAVES_A)
  fielded_sc <- as.integer(ASPIRATION_WAVES %in% SIGNAL_WAVES_SC)
  # The two domains' marks come from the same report card, so they are fielded
  # together and observed together for all but a handful of students.
  fielded_x <- fielded_a

  obs_a <- (!is.na(a)) * rep(fielded_a, each = nrow(a))
  prev <- cbind(NA_real_, a[, -ncol(a), drop = FALSE])
  prev_obs <- cbind(0L, obs_a[, -ncol(obs_a), drop = FALSE])
  obs_news <- obs_a * prev_obs
  news_a <- (a - prev) * obs_news
  news_a[is.na(news_a)] <- 0
  obs_x <- (!is.na(x)) * rep(fielded_x, each = nrow(x))
  obs_sc <- (!is.na(sc)) * rep(fielded_sc, each = nrow(sc))
  observed <- (!is.na(b)) * 1L

  blank <- function(m, keep) { m[keep == 0] <- 0; m }

  Z <- cbind(g5_mark = person$g5_mark, w1_sc = person$w1_sc,
             w2_interest = person$w2_interest, hisei = person$hisei)
  if (WITH_COV)
    Z <- cbind(Z, female = person$female, migrant = person$migrant,
               gym = person$gym)

  list(
    n_person = nrow(b),
    n_wave = length(ASPIRATION_WAVES),
    n_cov = ncol(Z),
    observed = observed,
    y = blank(b, observed),
    v_obs = blank(vo, observed),
    fielded_a = fielded_a,
    fielded_sc = fielded_sc,
    fielded_x = fielded_x,
    obs_a = obs_a,
    obs_sc = obs_sc,
    obs_x = obs_x,
    obs_news = obs_news,
    news_a = news_a,
    s_a = blank(a, obs_a),
    s_sc = blank(sc, obs_sc),
    s_x = blank(x, obs_x),
    z = Z,
    n_mod = 3L,
    w = cbind(female = person$female, hisei = person$hisei,
              migrant = person$migrant),
    ID_t = person$ID_t,
    waves = ASPIRATION_WAVES,
    missing_cov_share = n_missing_cov,
    # Standard deviation of the outcome at the anchor wave, in O*NET points.
    # One standardised unit equals this many points.
    onet_sd = b_scale,
    anchoring = anchoring
  )
}

# --- build and check ----------------------------------------------------------

check <- function(dat, label) {
  stopifnot(
    "no students retained" = dat$n_person > 3000,
    "wave count does not match the Stan data block" =
      dat$n_wave == length(ASPIRATION_WAVES),
    "non-finite values reached the Stan inputs" =
      all(is.finite(dat$y)) && all(is.finite(dat$s_a)) && all(is.finite(dat$s_x)) &&
      all(is.finite(dat$s_sc)) && all(is.finite(dat$z)) && all(is.finite(dat$w)),
    "an unobserved cell carries a value" =
      all(dat$y[dat$observed == 0] == 0) && all(dat$s_a[dat$obs_a == 0] == 0) &&
      all(dat$s_sc[dat$obs_sc == 0] == 0) && all(dat$s_x[dat$obs_x == 0] == 0),
    # The cross-domain signal must not be a copy of the within-domain one.
    "a change score exists where a mark does not" =
      all(dat$obs_news <= dat$obs_a) &&
      all(dat$news_a[dat$obs_news == 0] == 0),
    "the change score is implausibly sparse" =
      sum(dat$obs_news) > 0.5 * sum(dat$obs_a),
    "the two domains' marks are identical" =
      cor(dat$s_a[dat$obs_a == 1 & dat$obs_x == 1],
          dat$s_x[dat$obs_a == 1 & dat$obs_x == 1]) < 0.95,
    "measurement variance missing or negative" =
      all(is.finite(dat$v_obs)) && all(dat$v_obs >= 0),
    "a signal is marked observed at a wave where it was never fielded" =
      all(colSums(dat$obs_a)[dat$fielded_a == 0] == 0) &&
      all(colSums(dat$obs_sc)[dat$fielded_sc == 0] == 0),
    "aspiration series is implausibly sparse" = mean(dat$observed) > 0.2,
    # Anchoring a variable at a wave where it was never measured returns all NA.
    # Self-concept is measured at waves 1, 5 and 9 only, so this check exists.
    "a signal is empty at every wave" = sum(dat$obs_sc) > 0 && sum(dat$obs_a) > 0
  )
  cat(sprintf("%-22s %5d students  aspirations %.0f%%  marks %.0f%%  self-concept %.0f%%\n",
              label, dat$n_person, 100 * mean(dat$observed),
              100 * sum(dat$obs_a) / (dat$n_person * sum(dat$fielded_a)),
              100 * sum(dat$obs_sc) / (dat$n_person * sum(dat$fielded_sc))))
}

# The design is two domains by two aspiration types. Building all four here
# keeps the fitting script free of any provider-specific logic.
for (domain in c("math", "read")) {
  for (aspiration in c("realistic", "idealistic")) {
    label <- paste("aum", domain, aspiration, sep = "_")
    dat <- build_domain(domain, aspiration, coding = "kldb")
    check(dat, label)
    saveRDS(dat, here("data", "processed", paste0(label, SUFFIX, ".rds")))
  }
}

# -----------------------------------------------------------------------------
# RESULT (2026-08-03; re-run 2026-08-23, all four _cov bundles byte-identical)
#
#  cell                  students  aspirations  marks  self-concept
#  aum_math_realistic       6,864       48%      67%       59%
#  aum_math_idealistic      7,119       52%      66%       58%
#  aum_read_realistic       6,864       48%      67%       59%
#  aum_read_idealistic      7,119       52%      66%       58%
#
# Percentages are of the person-wave cells where the quantity could appear.
#
# Reading: the expected-career cells are smaller than the hoped-for ones because
# fewer students name a career they expect than one they would like. Coverage of
# the aspiration is under half, which is why the model integrates the latent
# state rather than requiring complete series.
#
# Note: 46.9% of consecutive aspiration pairs carry a change of exactly zero in
# the outcome and 36.5% name the identical occupation again, so the outcome
# moves rarely and in steps.
# -----------------------------------------------------------------------------

