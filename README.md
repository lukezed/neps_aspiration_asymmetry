# Asymmetric Revision of Career Aspirations

Analysis code for the manuscript *Asymmetric Revision in Adolescent Career
Aspirations: Longitudinal Evidence From School Marks and Self-Concept*
(under review).

The study fits a Bayesian state-space cross-lagged panel model to eight waves of
the German National Educational Panel Study (NEPS), Starting Cohort 3. Occupations
named by students are scored for their mathematics and reading demands, and the
model estimates how far an aspiration moves after a change in a report-card mark,
separately for improvements and declines, alongside domain-specific academic
self-concept.

## What is here

```
stan/
  aum.stan                Model. One program, eight specifications selected by
                          model_id: the reported model is 5 (directions free),
                          7 ties them, 4 drops the change term.
script/
  00_packages.R           Library loads and shared options. Sourced by the rest.
  01_onet_map.R           Builds the occupation-to-O*NET lookup from the KldB
                          and ISCO crosswalks.
  02_build.R              Reads the NEPS files and writes the eight analysis
                          bundles. Every later script starts from these.
  03_fit.R                Fits the Stan model. One call per cell per model_id.
  04_compare.R            Model comparison and every derived contrast the
                          Results report.
  05_wave9_check.R        Sensitivity: refit omitting the final wave.
  06_timing_check.R       Sensitivity: achievement entered lagged.
  07_ar_only.R            Records an autoregression-only fit that does not
                          converge. Kept because the manuscript says so.
  08_tables.R             Builds the design and descriptive tables.
  09_r2.R                 Bayesian R-squared on the one-step-ahead prediction.
  10_level_check.R        Sensitivity: response allowed to depend on the level
                          the mark changed from.
  11_reliability.R        Omega and alpha for the self-concept scales.
  90_prototype.R          Simulation checks. Parameter recovery, and a null
                          simulation with no asymmetry in the generating model.
figure_script/
  transitions.R           Figure 1, observed aspiration transitions.
  coefficients.R          Figures 2 and 3, asymmetry and signal coefficients.
data/                     Empty here. NEPS data are not redistributable.
models/                   Empty here. Fitted objects are written by 03_fit.R.
figure/                   Empty here. Written by figure_script/.
```

Not in this repository: the manuscript source, the working notes, and anything
derived from the data. `data/`, `models/` and `figure/` are present as empty
directories so the paths the scripts expect exist.

## Running it

The data are available from the Leibniz Institute for Educational Trajectories
under a data use agreement, free of charge for scientific use:

> NEPS Network. National Educational Panel Study, Starting Cohort 3, Scientific
> Use File. <https://www.neps-data.de> — doi:10.5157/NEPS:SC3:13.1.0

Request Starting Cohort 3, download the SUF, and put the files in `data/`.
Then, from the repository root:

```r
source("script/00_packages.R")
source("script/01_onet_map.R")
source("script/02_build.R")
source("script/03_fit.R")
source("script/04_compare.R")
```

The remaining numbered scripts are independent of each other and can be run in
any order once `02_build.R` and `03_fit.R` have completed. Fitting needs
cmdstanr and a working CmdStan installation; the four reported cells take a few
hours between them on a laptop.

Each script that produces a number the manuscript quotes ends in a `RESULT`
comment block recording what the run printed and when. Those blocks are a record,
not a source: re-run the script rather than quoting them.

## Licence and attribution

O*NET data are used under CC BY 4.0. O*NET is a trademark of USDOL/ETA. The NEPS
data are governed by the data use agreement above and are not included here.
