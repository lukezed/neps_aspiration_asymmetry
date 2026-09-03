
source(here::here("script", "00_packages.R"))

dta <- here("data", "SC3_D_13-1-0", "Stata14")
suf <- function(x) file.path(dta, paste0("SC3_", x, "_D_13-1-0.dta"))

ALL_WAVES <- c(1:9, 11, 13)
GRADE_OF <- c("5", "6", "7", "8", "9", "9", "10", "11", "12", "", "")
YEAR_OF <- c("2010", "2011", "2012", "2013", "2014",
             "2015", "2016", "2016", "2017", "2019", "2021")

cells <- lapply(c(math_realistic  = "aum_math_realistic",
                  read_realistic  = "aum_read_realistic",
                  math_idealistic = "aum_math_idealistic",
                  read_idealistic = "aum_read_idealistic"),
                function(f) readRDS(here("data", "processed", paste0(f, ".rds"))))

SAMPLE <- cells$math_realistic$ID_t                   # the 6,864
N <- length(SAMPLE)

raw <- read_dta(suf("pTarget"), col_select = all_of(c(
  "ID_t", "wave", "tx20100",
  "t724102", "t724112", "t724602",       # maths mark: annual, half-year, points
  "t724101", "t724111", "t724601",       # German mark
  "t66001a_g1", "t66000a_g1",             # self-concept: maths, verbal
  "t66201a_g1", "t66208a_g1"              # interest: maths, German
))) %>%
  filter(tx20100 == 1, ID_t %in% SAMPLE) %>%
  mutate(
    mark_math = -pmin(neps_na(t724102), neps_na(t724112),
                      points_to_note(neps_na(t724602)), na.rm = TRUE),
    mark_read = -pmin(neps_na(t724101), neps_na(t724111),
                      points_to_note(neps_na(t724601)), na.rm = TRUE),
    sc_math = neps_na(t66001a_g1),
    sc_read = neps_na(t66000a_g1),
    int_math = neps_na(t66201a_g1),
    int_read = neps_na(t66208a_g1)
  )

raw_n <- function(col) sapply(ALL_WAVES, function(w)
  sum(!is.na(raw[[col]][raw$wave == w])))

# --- the design table ---------------------------------------------------------

# A built matrix covers the aspiration waves only; place its columns into the
# full wave grid and leave the rest at zero.
on_grid <- function(counts, waves) {
  out <- integer(length(ALL_WAVES))
  out[match(waves, ALL_WAVES)] <- counts
  out
}

ideal <- cells$math_idealistic
in_sample <- match(SAMPLE, ideal$ID_t)
ideal_obs <- colSums(ideal$observed[in_sample[!is.na(in_sample)], , drop = FALSE])

design <- rbind(
  `Realistic aspiration` = on_grid(colSums(cells$math_realistic$observed),
                              cells$math_realistic$waves),
  `Idealistic aspiration` = on_grid(ideal_obs, ideal$waves),
  `Mathematics mark` = raw_n("mark_math"),
  `German mark` = raw_n("mark_read"),
  `Mathematics self-concept` = raw_n("sc_math"),
  `Verbal self-concept` = raw_n("sc_read"),
  `Mathematics interest` = raw_n("int_math"),
  `German interest` = raw_n("int_read")
)

late <- design[c("Mathematics mark", "Mathematics self-concept"),
               match(c(11, 13), ALL_WAVES)]
stopifnot("a late-wave signal is larger than the note claims" =
            all(late["Mathematics mark", ] < 150) &&
            all(late["Mathematics self-concept", ] < 140))

fmt <- function(x) ifelse(x == 0, "", format(x, big.mark = ",", trim = TRUE))
csv_row <- function(name, vals) paste0('"', name, '",',
                                       paste0('"', vals, '"', collapse = ","))

cat("\n=== Table 1, design. Paste into paper.qmd's tbl-design chunk ===\n")
cat('"",', paste0('"W', ALL_WAVES, '"', collapse = ","), "\n", sep = "")
cat(csv_row("School grade", GRADE_OF), "\n")
cat(csv_row("Year", YEAR_OF), "\n")
for (r in rownames(design)) cat(csv_row(r, fmt(design[r, ])), "\n")
cat("\n  late waves, marks W11/W13:", late["Mathematics mark", ],
    " self-concept:", late["Mathematics self-concept", ], "\n")

# --- the descriptives and correlations table ----------------------------------

G <- c(`Grade 9` = 5, `Grade 12` = 9)

ANCHOR_OF <- c(y = "b_obs", s_a = "s_a", s_sc = "s_sc", s_x = "s_x")

DISPLAY_SHIFT <- c(s_a = 7, s_x = 7)

column <- function(dat, mat, obs, w) {
  j <- match(w, dat$waves)
  i <- match(SAMPLE, dat$ID_t)
  v <- dat[[mat]][i, j]
  v[dat[[obs]][i, j] == 0] <- NA
  k <- dat$anchoring[, ANCHOR_OF[[mat]]]
  v * k[["sd"]] + k[["mean"]] + (if (mat %in% names(DISPLAY_SHIFT)) DISPLAY_SHIFT[[mat]] else 0)
}

for (g in names(G)) {
  w <- G[[g]]
  vars <- list(
    "1. Realistic aspiration, mathematics" =
      column(cells$math_realistic, "y", "observed", w),
    "2. Realistic aspiration, reading" =
      column(cells$read_realistic, "y", "observed", w),
    "3. Idealistic aspiration, mathematics" =
      column(cells$math_idealistic, "y", "observed", w),
    "4. Idealistic aspiration, reading" =
      column(cells$read_idealistic, "y", "observed", w),
    "5. Mathematics achievement" = column(cells$math_realistic, "s_a", "obs_a", w),
    "6. German achievement" = column(cells$read_realistic, "s_a", "obs_a", w),
    "7. Mathematics self-concept" =
      column(cells$math_realistic, "s_sc", "obs_sc", w),
    "8. Verbal self-concept" = column(cells$read_realistic, "s_sc", "obs_sc", w)
  )
  M <- do.call(cbind, vars)
  r <- cor(M, use = "pairwise.complete.obs")

  cat("\n=== Table 2, ", g, ". Paste into paper.qmd's tbl-cor chunk ===\n",
      sep = "")
  for (k in seq_along(vars)) {
    x <- vars[[k]]
    # Leading zeros are dropped from correlations, which is the convention the
    # manuscript follows; means and SDs keep theirs.
    below <- if (k == 1) character(0)
             else sub("^(-?)0\\.", "\\1.", sprintf("%.2f", r[k, seq_len(k - 1)]))
    cells_out <- c(format(sum(!is.na(x)), big.mark = ","),
                   sprintf("%.2f", mean(x, na.rm = TRUE)),
                   sprintf("%.2f", sd(x, na.rm = TRUE)),
                   below, rep("", length(vars) - k))
    cat(csv_row(names(vars)[k], cells_out), "\n")
  }
}

check_pattern <- function(w) {
  m <- cbind(exp_m = column(cells$math_realistic, "y", "observed", w),
             exp_r = column(cells$read_realistic, "y", "observed", w),
             mark_m = column(cells$math_realistic, "s_a", "obs_a", w),
             mark_r = column(cells$read_realistic, "s_a", "obs_a", w),
             sc_m = column(cells$math_realistic, "s_sc", "obs_sc", w),
             sc_r = column(cells$read_realistic, "s_sc", "obs_sc", w))
  r <- cor(m, use = "pairwise.complete.obs")
  stopifnot(
    "the realistic aspiration no longer tracks its own domain's mark" =
      r["exp_m", "mark_m"] > 3 * abs(r["exp_m", "mark_r"]),
    "the two self-concepts are no longer uncorrelated" =
      abs(r["sc_m", "sc_r"]) < 0.2,
    "the two marks are no longer correlated" = r["mark_m", "mark_r"] > 0.3,
    "self-concept no longer tracks its own mark more strongly" =
      r["sc_m", "mark_m"] > r["exp_m", "mark_m"]
  )
}
invisible(lapply(G, check_pattern))
cat("\nthe three patterns the prose reads off Table 2 hold at both grades\n")

# --- the numbers quoted in the Sample and Measures prose ---------------------

asp <- read_dta(suf("pTarget"), col_select = all_of(c(
  "ID_t", "wave", "tx20100",
  "t31160a_g2", "t31160a_g4", "t31160a_g14",   # realistic: KldB, ISCO-08, ISEI-08
  "t31060a_g2", "t31060a_g4", "t31060a_g14"    # idealistic: the same three
))) %>% filter(tx20100 == 1)

W49 <- 4:9
kldb <- readRDS(here("data", "processed", "onet_kldb.rds"))

cat("\n=== Numbers quoted in the Sample and Measures prose ===\n")

gap_share <- function(dat) {
  obs <- dat$observed
  apply(obs, 1, function(r) {
    j <- which(r == 1)
    length(j) > 0 && any(diff(j) > 1)
  })
}
for (nm in c("math_realistic", "math_idealistic")) {
  g <- gap_share(cells[[nm]])
  cat(sprintf("  %-18s n %5d   interior gap %5d   unbroken %5d\n",
              nm, length(g), sum(g), sum(!g)))
}

cat(sprintf("  analysis sample, realistic aspiration     %5d\n", N))
cat(sprintf("  analysis sample, idealistic aspiration    %5d\n",
            cells$math_idealistic$n_person))

# Whether the two questions are redundant.
both <- asp %>%
  filter(wave %in% W49) %>%
  transmute(r4 = neps_na(t31160a_g4), i4 = neps_na(t31060a_g4),
            r_isei = neps_na(t31160a_g14), i_isei = neps_na(t31060a_g14)) %>%
  filter(!is.na(r4), !is.na(i4))
cat(sprintf("\n  person-waves W4-9 with both coded         %5d\n", nrow(both)))
cat(sprintf("    same four-digit ISCO-08                 %5.1f%%\n",
            100 * mean(both$r4 == both$i4)))
cat(sprintf("    same major group                        %5.1f%%\n",
            100 * mean(both$r4 %/% 1000 == both$i4 %/% 1000)))
cat(sprintf("    mean ISEI-08, realistic %.1f  idealistic %.1f\n",
            mean(both$r_isei, na.rm = TRUE), mean(both$i_isei, na.rm = TRUE)))

# Coverage and the shape of non-response, on the realistic-aspiration question.
w49 <- asp %>% filter(wave %in% W49)

isei_cor <- function(code_var, isei_var, col) {
  code <- neps_na(w49[[code_var]])
  isei <- neps_na(w49[[isei_var]])
  val <- kldb[[col]][match(code, kldb$kldb)]
  ok <- !is.na(val) & !is.na(isei)
  c(n = sum(ok), r = cor(val[ok], isei[ok]))
}
cat("\n  O*NET score of the named occupation against its own ISEI-08, W4-9:\n")
for (a in list(c("t31160a_g2", "t31160a_g14", "realistic"),
               c("t31060a_g2", "t31060a_g14", "idealistic")))
  for (col in c("math", "read")) {
    x <- isei_cor(a[1], a[2], col)
    cat(sprintf("    %-11s %-5s n %6d   r %+.2f\n", a[3], col, x[1], x[2]))
  }
raw_code <- haven::zap_labels(w49$t31160a_g2)
cat(sprintf("\n  person-waves W4-9, realistic aspiration:\n"))
cat(sprintf("    a valid KldB code                       %5d\n",
            sum(!is.na(raw_code) & raw_code > 0)))
cat(sprintf("    of those, matched to an O*NET occupation%5d\n",
            sum(neps_na(w49$t31160a_g2) %in% kldb$kldb)))
cat(sprintf("    missing by design (-54)                 %5d\n",
            sum(!is.na(raw_code) & raw_code == -54)))
cat(sprintf("    given but not codeable (-55)            %5d\n",
            sum(!is.na(raw_code) & raw_code == -55)))
cat(sprintf("    other missing codes                     %5d\n",
            sum(!is.na(raw_code) & raw_code < 0 & !raw_code %in% c(-54, -55)) +
              sum(is.na(raw_code))))

cat("\n  non-response over W4-9, by question and coding:\n")
cat(sprintf("    %-26s %6s %6s %6s %6s %6s\n",
            "", "valid", "design", "uncod.", "dk", "other"))
for (v in c("t31160a_g2", "t31160a_g4", "t31060a_g2", "t31060a_g4")) {
  z <- haven::zap_labels(w49[[v]])
  lab <- paste(if (startsWith(v, "t31160")) "realistic" else "idealistic",
               if (endsWith(v, "g2")) "KldB" else "ISCO-08")
  cat(sprintf("    %-26s %6d %6d %6d %6d %6d\n", lab,
              sum(!is.na(z) & z > 0), sum(!is.na(z) & z == -54),
              sum(!is.na(z) & z == -55), sum(!is.na(z) & z == -98),
              sum(!is.na(z) & z < 0 & !z %in% c(-54, -55, -98)) + sum(is.na(z))))
}

# Whether an uncodeable answer signals anything about the student. If it did,
# treating it as an unobserved aspiration would bias the outcome.
unc <- asp %>%
  filter(wave %in% W49) %>%
  transmute(ID_t, wave, code = haven::zap_labels(t31160a_g2)) %>%
  filter(!is.na(code)) %>%
  left_join(transmute(raw, ID_t, wave, mark = mark_math), by = c("ID_t", "wave"))
by_wave <- table(unc$wave[unc$code == -55])
cat(sprintf("\n    uncodeable realistic-aspiration answers, W4 %d to W9 %d\n",
            by_wave[["4"]], by_wave[["9"]]))
# The mark is stored negated so that higher is better; flip it back to report a
# German mark, where 1 is best.
cat(sprintf("    mean mark among them %.2f, among coded answers %.2f\n",
            -mean(unc$mark[unc$code == -55], na.rm = TRUE),
            -mean(unc$mark[unc$code > 0], na.rm = TRUE)))
dk <- unc %>% filter(code == -98, !is.na(mark))
cat(sprintf("    'don't know' with an observed mark: %d, of which %d at wave 9\n",
            nrow(dk), sum(dk$wave == 9)))

cross <- readxl::read_excel(here("data", "raw", "KldB_onet_final_matching.xlsx")) %>%
  transmute(kldb = KldB_code, occupation = Onet_label)
kldb_lookup <- readRDS(here("data", "processed", "onet_kldb.rds"))

occ <- asp %>%
  filter(wave %in% ALL_WAVES, !is.na(neps_na(t31160a_g2))) %>%
  transmute(code = neps_na(t31160a_g2)) %>%
  left_join(cross, by = c("code" = "kldb")) %>%
  left_join(select(kldb_lookup, kldb, math, read), by = c("code" = "kldb")) %>%
  filter(!is.na(occupation), !is.na(math)) %>%
  count(occupation, math, read, name = "n")

MIN_REPORTS <- 20

listed <- occ %>%
  filter(n >= MIN_REPORTS) %>%
  arrange(desc(math), desc(n)) %>%
  # O*NET titles carry a disambiguating tail that adds nothing at this width.
  mutate(occupation = sub(",\\s*(Except|All Other).*$", "", occupation))

covered <- sum(listed$n) / sum(occ$n)

# Emitted in two column-blocks side by side. One block of 156 rows runs to
# three pages of mostly white space; two blocks of 78 halve it.
half <- ceiling(nrow(listed) / 2)
pad <- function(x, n) c(x, rep("", n - length(x)))
lhs <- listed[seq_len(half), ]
rhs <- listed[seq(half + 1, nrow(listed)), ]

cat("\n=== Appendix table: occupations named at least", MIN_REPORTS, "times ===\n")
cat('"Occupation","Maths","Reading","$n$","Occupation ","Maths ","Reading ","$n$ "\n')
for (i in seq_len(half)) {
  fields <- c(lhs$occupation[i], sprintf("%.0f", lhs$math[i]),
              sprintf("%.0f", lhs$read[i]), format(lhs$n[i], big.mark = ","),
              pad(rhs$occupation[i], 1),
              pad(if (i <= nrow(rhs)) sprintf("%.0f", rhs$math[i]) else "", 1),
              pad(if (i <= nrow(rhs)) sprintf("%.0f", rhs$read[i]) else "", 1),
              pad(if (i <= nrow(rhs)) format(rhs$n[i], big.mark = ",") else "", 1))
  cat(paste0('"', fields, '"', collapse = ","), "\n", sep = "")
}
cat(sprintf("\n  %d of %d occupations listed, covering %.1f%% of %s reports\n",
            nrow(listed), nrow(occ), 100 * covered,
            format(sum(occ$n), big.mark = ",")))
cat(sprintf("  omitted: %d occupations, %.1f%% of reports\n",
            nrow(occ) - nrow(listed), 100 * (1 - covered)))
cat(sprintf("  the two contents correlate r = %.2f across named occupations\n",
            cor(occ$math, occ$read)))

# Three anchors for the Measures sentence, which carries the calibration in the
# main text so that no float is needed there.
anchors <- c("Police and Sheriff's Patrol Officers", "Secondary School Teachers",
             "Mathematicians")
cat("\n  anchors for the Measures sentence:\n")
for (a in anchors) {
  row <- occ[startsWith(occ$occupation, a), ][1, ]
  cat(sprintf("    %-38s maths %3.0f  reading %3.0f\n", a, row$math, row$read))
}

stopifnot(
  "the listed occupations no longer span the scale" =
    diff(range(listed$math)) > 90,
  "the threshold no longer covers most reports" = covered > 0.9
)

g <- gap_share(cells$math_realistic)
quoted <- c(sample = N, gap = sum(g), unbroken = sum(!g),
            mark_uncodeable = round(-mean(unc$mark[unc$code == -55], na.rm = TRUE), 2),
            mark_coded = round(-mean(unc$mark[unc$code > 0], na.rm = TRUE), 2))
in_paper <- c(sample = 6864, gap = 3774, unbroken = 3090,
              mark_uncodeable = 2.77, mark_coded = 2.80)
drift <- names(quoted)[abs(quoted - in_paper) > 1e-8]
if (length(drift)) {
  print(rbind(computed = quoted, in_paper = in_paper)[, drift, drop = FALSE])
  stop("paper.qmd quotes a number this script no longer produces: ",
       paste(drift, collapse = ", "))
}
cat("\n  the five numbers quoted in Sample and Measures match paper.qmd\n")

# -----------------------------------------------------------------------------
# RESULT (2026-08-03; section 6 added and script re-run 2026-08-22, its
# assertion against paper.qmd passing)
#
# Two stale references stopped this script part-way from 2026-08-22 back to
# whenever the design table's mark and self-concept rows were split by domain.
# late["Report-card mark", ] and raw$mark were both renamed at that point and
# two call sites were not. Neither number reaches the manuscript, so nothing
# printed below had drifted; the script simply could not be run to the end.
#
# 1. DESIGN TABLE, realistic-aspiration analysis sample, N = 6,864.
#
#                          W1     W2     W3     W4     W5     W6     W7     W8     W9    W11    W13
#   Expected career         -      -      -   4,151  3,511  3,763  3,495  2,544  3,777  2,670  2,337
#   Hoped-for career        -      -      -   4,954  4,218  4,298  3,835  2,824  4,066  2,678  2,211
#   Report-card mark      4,278  4,320  5,890 5,898  5,247  5,171  4,649  3,464  3,082    139     15
#   Academic self-concept 4,327    -      -      -   5,293    -      -      -   2,814    133     15
#   Domain interest         -    4,279    -      -      -   5,078    -      -      -       -      -
#
# The wave 11 and 13 counts are reported rather than left blank, so that a blank
# cell means the item was not fielded and nothing else. 139 marks and 133
# self-concept reports cannot identify a wave-specific coefficient, which is why
# neither signal enters the transition after wave 9. The comment in 02_build.R
# gives 136 for wave-11 self-concept because it counts students in either
# aspiration cell; this table runs on the realistic-aspiration sample alone.
#
# 2. DESCRIPTIVES AND CORRELATIONS. Pairwise complete, on the same sample.
#
#   Grade 9 (wave 5)         n       M     SD      1     2     3     4     5
#     1. Expected career, maths   3,511   0.01  0.99
#     2. Expected career, reading 3,511   0.04  0.97   .29
#     3. Mathematics mark         5,247  -0.08  1.03   .17   .20
#     4. Language mark            5,251   0.01  1.00   .04   .27   .42
#     5. Mathematics self-concept 5,293  -0.49  1.05   .22   .06   .64   .07
#     6. Verbal self-concept      5,301  -0.04  0.91   .02   .22   .11   .54  -.05
#
#   Grade 12 (wave 9)
#     1. Expected career, maths   3,777   0.10  0.94
#     2. Expected career, reading 3,777   0.16  0.87   .29
#     3. Mathematics mark         3,082   0.09  1.13   .17   .09
#     4. Language mark            3,087   0.22  1.07   .01   .21   .44
#     5. Mathematics self-concept 2,814  -0.52  1.15   .22  -.01   .71   .08
#     6. Verbal self-concept      2,807  -0.10  1.07  -.10   .20   .05   .61  -.15
#
# These reproduce the values already in paper.qmd exactly, so Table 2 needed no
# correction.
#
# Reading: a realistic aspiration correlates with its own domain's mark and not with
# the other's, .17 against .04 at Grade 9 and .17 against .01 at Grade 12. The
# two self-concepts are uncorrelated, -.05 and -.15, while the two marks
# correlate .42 and .44, which is the internal/external frame of reference
# pattern. Self-concept tracks its own mark more closely, .64 rising to .71,
# than any realistic aspiration tracks anything.
#
# Nothing in the correlation matrix moves between the two grades: the maths mark
# and the realistic aspiration sit at .17 at both. What the model reports is
# therefore not visible in a table of concurrent associations, which is the
# argument for fitting it.
#
# 3. THE OCCUPATIONS BEHIND THE SCALE.
#
#   distinct occupations named                313   over 26,248 reports
#   listed in the appendix                    156   named 20 or more times,
#                                                 set in two column-blocks
#   share of reports they cover             95.3%
#   omitted                                   157   4.7% of reports between them
#   the two contents correlate           r = .46   across named occupations
#
#   Coverage, to show why the tail is dropped rather than listed:
#     50% of reports are covered by  20 occupations
#     80%                            71
#     90%                           111
#     95%                           153
#     100 occupations were named fewer than ten times and carry 1.6% together.
#
#   Anchors quoted in Measures: police officer 31, secondary school teacher 47,
#   mathematician 100.
#
#   The count is person-waves, not students: 1,467 reports of police officer come
#   from 760 students, so a career named at several waves is counted each time.
#   The table's note says so, because "n" over an occupation invites the other
#   reading.
#
# Reading: this is what a coefficient in O*NET points refers to. A police
# officer and a secondary school teacher are sixteen points apart, and the
# reported response to an improving mark is about one. Stating three anchors in
# the Measures sentence carries that calibration without a float; the full list
# is a lookup and belongs in an appendix.
#
# Note: the complete crosswalk is not reproduced. It is the 1,332 row table
# Graf et al. published with their article, used unchanged, so reprinting it
# would duplicate something already public.
#
# Note: this replaced a figure that placed all 313 occupations by their two
# scores. It failed. The named occupations crowd into one corner of the panel
# and no labelling scheme separated them without the labels becoming the figure.
#
# -----------------------------------------------------------------------------

