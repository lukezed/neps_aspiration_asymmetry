
source(here::here("script","00_packages.R"))
dta <- here("data","SC3_D_13-1-0","Stata14")
suf <- function(x) file.path(dta, paste0("SC3_", x, "_D_13-1-0.dta"))
SAMPLE <- readRDS(here("data","processed","aum_math_realistic_cov.rds"))$ID_t
items <- c(paste0("t66001",letters[1:3]), paste0("t66000",letters[1:3]))
d <- read_dta(suf("pTarget"), col_select=all_of(c("ID_t","wave","tx20100",items))) |>
  filter(tx20100==1, ID_t %in% SAMPLE, wave %in% c(1,5,9))
om <- function(x){ x<-x[complete.cases(x),]; fa<-factanal(x,factors=1)
  l<-as.numeric(fa$loadings); u<-as.numeric(fa$uniquenesses); sum(l)^2/(sum(l)^2+sum(u)) }
al <- function(x){ x<-x[complete.cases(x),]; k<-ncol(x)
  k/(k-1)*(1-sum(apply(x,2,var))/var(rowSums(x))) }
for (dom in c("mathematics","language")) {
  cols <- if (dom=="mathematics") items[1:3] else items[4:6]
  for (w in c(1,5,9)) {
    x <- d |> filter(wave==w) |> select(all_of(cols)) |> mutate(across(everything(), neps_na))
    if (dom=="language") x[[1]] <- 5 - x[[1]]     # first language item is reverse-keyed
    cat(sprintf("  %-12s wave %-2d (Grade %-2s)  n %5d   omega %.2f   alpha %.2f\n",
        dom, w, c("1"="5","5"="9","9"="12")[as.character(w)], sum(complete.cases(x)), om(x), al(x)))
  }
}

# -----------------------------------------------------------------------------
# RESULT (2026-08-23)
#
#   mathematics  Grade 5   n 4,327   omega .87   alpha .86
#   mathematics  Grade 9   n 5,293   omega .89   alpha .89
#   mathematics  Grade 12  n 2,814   omega .90   alpha .89
#   language     Grade 5   n 4,320   omega .75   alpha .73
#   language     Grade 9   n 5,301   omega .83   alpha .82
#   language     Grade 12  n 2,807   omega .87   alpha .87
#
# Longitudinal invariance is not tested here. Graf et al. (2026) establish
# metric, scalar and strict invariance for this scale on this cohort.
# -----------------------------------------------------------------------------

