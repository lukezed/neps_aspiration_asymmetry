
source(here::here("script", "00_packages.R"))

raw_dir <- here("data", "raw")
onet_dir <- file.path(raw_dir, "db_28_3_text")
esco_file <- file.path(raw_dir, "ESCO_to_ONET-SOC.xlsx")
kldb_file <- file.path(raw_dir, "KldB_onet_final_matching.xlsx")

stopifnot(
  "O*NET Skills file missing; see doc/SOURCES.md for the download" =
    file.exists(file.path(onet_dir, "Skills.txt")),
  "ESCO crosswalk missing; see doc/SOURCES.md" = file.exists(esco_file),
  "KldB crosswalk missing; see doc/SOURCES.md" = file.exists(kldb_file)
)

importance <- read_tsv(file.path(onet_dir, "Skills.txt"), show_col_types = FALSE) %>%
  filter(`Element Name` %in% c("Mathematics", "Reading Comprehension"),
         `Scale ID` == "IM") %>%
  transmute(
    onetsoc = `O*NET-SOC Code`,
    domain = if_else(`Element Name` == "Mathematics", "math", "read"),
    score = (`Data Value` - 1) / 4 * 100
  )

# --- main route: KldB, one code to one occupation ----------------------------

onet_kldb <- readxl::read_excel(kldb_file) %>%
  transmute(kldb = KldB_code, onetsoc = Onet_code) %>%
  inner_join(importance, by = "onetsoc", relationship = "many-to-many") %>%
  pivot_wider(id_cols = kldb, names_from = domain, values_from = score) %>%
  # No averaging happens on this route, so there is no coarseness error. The
  # columns exist so both routes present the same interface downstream.
  mutate(math_sd = 0, read_sd = 0)

# --- robustness route: ISCO-08 via O*NET's ESCO crosswalk --------------------

crosswalk_isco <- readxl::read_excel(esco_file, skip = 3) %>%
  set_names(c("esco_code", "esco_title", "onetsoc", "onet_title")) %>%
  filter(!is.na(esco_code), !is.na(onetsoc)) %>%
  transmute(isco08 = as.numeric(str_extract(esco_code, "^[0-9]{4}")), onetsoc) %>%
  filter(!is.na(isco08)) %>%
  distinct()

onet_isco <- importance %>%
  inner_join(crosswalk_isco, by = "onetsoc", relationship = "many-to-many") %>%
  group_by(isco08, domain) %>%
  # Order matters: summarise evaluates sequentially, so computing the mean into
  # `score` first would leave the spread to be taken over a single number.
  summarise(spread = if (dplyr::n() > 1) stats::sd(score) else 0,
            score = mean(score), .groups = "drop") %>%
  pivot_wider(names_from = domain, values_from = c(score, spread)) %>%
  transmute(isco08, math = score_math, read = score_read,
            math_sd = spread_math, read_sd = spread_read)

dir.create(here("data", "processed"), showWarnings = FALSE, recursive = TRUE)
saveRDS(onet_kldb, here("data", "processed", "onet_kldb.rds"))
saveRDS(onet_isco, here("data", "processed", "onet_isco.rds"))

# --- self-check --------------------------------------------------------------

stopifnot(
  "KldB lookup too small" = nrow(onet_kldb) > 1300,
  "ISCO lookup too small" = nrow(onet_isco) > 350,
  "scores outside 0-100" =
    all(dplyr::between(c(onet_kldb$math, onet_kldb$read,
                         onet_isco$math, onet_isco$read), 0, 100)),
  "aggregate ISCO codes are absent by construction" =
    all(onet_isco$isco08 %% 100 != 0),
  # The codes a SOC-2010 crosswalk would silently drop.
  "medical practitioners missing from the ISCO route" =
    all(c(2211, 2212) %in% onet_isco$isco08),
  "ICT professionals missing from the ISCO route" =
    any(onet_isco$isco08 %/% 100 == 25),
  # The one-to-one route must not acquire an averaging error.
  "KldB route should carry no coarseness error" = all(onet_kldb$math_sd == 0)
)

sig <- var(onet_isco$math)
err <- mean(onet_isco$math_sd^2, na.rm = TRUE)

cat("KldB route:", nrow(onet_kldb), "codes, maths SD",
    round(sd(onet_kldb$math), 1),
    "  one code to one occupation, so no coarseness error\n")
cat("ISCO route:", nrow(onet_isco), "codes, maths SD",
    round(sd(onet_isco$math), 1),
    "  reliability", round(sig / (sig + err), 2),
    "(mean within-code SD", paste0(round(sqrt(err), 1), ")"), "\n")

# -----------------------------------------------------------------------------
# RESULT (2026-08-02; re-run 2026-08-23, output byte-identical)
#
# KldB route:  1,327 codes, mathematics SD 15.2   <- USED
#              one code maps to one occupation, so no averaging error
# ISCO route:    425 codes, mathematics SD 10.7
#              reliability 0.50, mean within-code SD 10.6
#
# An ISCO unit group covers a median of nine O*NET occupations, so averaging
# leaves about as much error variance as signal. KldB is one-to-one.
#
# Note: 858 aspiration reports carry an ISCO code ending in 00, an aggregate
# level no four-digit crosswalk can match, and 692 of those are group 2100,
# science and engineering. The same 858 person-waves carry no KldB code either,
# so both routes lose the same students; these are answers too vague to code
# precisely in either taxonomy.
# -----------------------------------------------------------------------------

