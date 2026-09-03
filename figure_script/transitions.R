
source(here::here("script", "00_packages.R"))
library(patchwork)

fig_dir <- here("figure")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

dta <- here("data", "SC3_D_13-1-0", "Stata14")
suf <- function(x) file.path(dta, paste0("SC3_", x, "_D_13-1-0.dta"))

onet <- readRDS(here("data", "processed", "onet_kldb.rds"))
dat <- readRDS(here("data", "processed", "aum_math_realistic.rds"))

# --- shared style -------------------------------------------------------------

BLUE <- "#6497B1"      # stayed
RED  <- "#CC2936"      # moved

theme_paper <- function(base_size = 10) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
      axis.line = element_blank(),
      axis.ticks = element_line(colour = "black", linewidth = 0.3),
      legend.position = "inside",
      legend.justification = c(1, 0),
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      legend.key = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 8, colour = "black",
                                 margin = margin(l = 2)),
      legend.margin = margin(0, 0, 0, 0),
      legend.spacing.y = unit(5, "pt"),
      legend.key.spacing.y = unit(0, "pt"),
      legend.key.height = unit(0.72, "lines"),
      plot.margin = margin(4, 4, 4, 4)
    )
}

# Every figure goes out as PDF for the manuscript, PNG for quick viewing and
# TIFF at 300 dpi because that is what the journals ask for on acceptance.
save_figure <- function(plot, name, width, height) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), plot,
         width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(file.path(fig_dir, paste0(name, ".png")), plot,
         width = width, height = height, units = "in", dpi = 300,
         device = ragg::agg_png)
  ggsave(file.path(fig_dir, paste0(name, ".tiff")), plot,
         width = width, height = height, units = "in", dpi = 300,
         device = ragg::agg_tiff)
}

# --- the transitions ----------------------------------------------------------

codes_of <- function(var) {
  read_dta(suf("pTarget"), col_select = all_of(
    c("ID_t", "wave", "tx20100", var))) %>%
    filter(tx20100 == 1, ID_t %in% dat$ID_t, wave %in% dat$waves) %>%
    transmute(ID_t, wave, code = neps_na(.data[[var]]))
}
CODES <- list(realistic = codes_of("t31160a_g2"),
              idealistic = codes_of("t31060a_g2"))

transitions <- function(domain, aspiration) {
  CODES[[aspiration]] %>%
    mutate(score = onet[[domain]][match(code, onet$kldb)]) %>%
    filter(!is.na(score)) %>%
    arrange(ID_t, wave) %>%
    group_by(ID_t) %>%
    mutate(prev_code = lag(code), prev_score = lag(score)) %>%
    ungroup() %>%
    filter(!is.na(prev_score))
}

transition_plot <- function(domain, aspiration, label) {
  p <- transitions(domain, aspiration)
  cell <- p %>%
    count(prev_score, score, name = "n") %>%
    mutate(moved = factor(if_else(prev_score == score, "Unchanged", "Moved"),
                          levels = c("Unchanged", "Moved")))

  gg <- ggplot(cell, aes(prev_score, score)) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.4,
                linetype = "dotted", colour = "grey55") +
    geom_point(aes(size = n, fill = moved), shape = 21,
               colour = "white", stroke = 0.25, alpha = 0.9) +
    scale_size_area(max_size = 7, breaks = c(10, 100, 500, 1500),
                    limits = c(0, 1600), labels = scales::comma) +
    scale_fill_manual(values = c(Unchanged = BLUE, Moved = RED)) +
    scale_x_continuous(breaks = seq(0, 100, 25)) +
    scale_y_continuous(breaks = seq(0, 100, 25)) +
    coord_cartesian(xlim = c(-6, 106), ylim = c(-6, 106), expand = FALSE) +
    labs(x = NULL, y = NULL, title = label) +
    guides(fill = guide_legend(override.aes = list(size = 3.4), order = 1),
           size = guide_legend(override.aes = list(fill = NA,
                                                   colour = "grey35",
                                                   stroke = 0.5),
                               order = 2)) +
    theme_paper() +
    theme(plot.title = element_text(size = 8, hjust = 0.5,
                                    margin = margin(b = 2)))

  list(plot = gg, n = nrow(p),
       same_score = mean(p$score == p$prev_score),
       same_job = mean(p$code == p$prev_code),
       step = median(diff(sort(unique(onet[[domain]])))),
       values = length(unique(onet[[domain]])))
}

m  <- transition_plot("math", "realistic",  "Mathematics")
r  <- transition_plot("read", "realistic",  "Reading")
mi <- transition_plot("math", "idealistic", "Mathematics")
ri <- transition_plot("read", "idealistic", "Reading")

stopifnot(
  "the reported unchanged share has moved" = abs(m$same_score - 0.469) < 0.005,
  "the reported same-occupation share has moved" = abs(m$same_job - 0.365) < 0.005
)

face <- function(g, keep_y, keep_x) {
  g <- g + theme(legend.position = "right",
                 legend.justification = "center",
                 legend.direction = "vertical",
                 legend.box = "vertical",
                 legend.text = element_text(size = 7.5, colour = "black",
                                            margin = margin(l = 2)))
  g <- g + theme(axis.title.x = if (keep_x)
                   element_text(size = 8.5) else element_blank(),
                 axis.title.y = if (keep_y)
                   element_text(size = 8.5) else element_blank())
  if (!keep_y) g <- g + theme(axis.text.y = element_blank(),
                              axis.ticks.y = element_blank())
  if (!keep_x) g <- g + theme(axis.text.x = element_blank(),
                              axis.ticks.x = element_blank())
  g + labs(x = "Score at report t", y = "Score at report t + 1")
}
combined <- ((face(mi$plot, TRUE, FALSE) | face(ri$plot, FALSE, FALSE)) /
             (face(m$plot,  TRUE, TRUE)  | face(r$plot,  FALSE, TRUE))) +
  patchwork::plot_layout(guides = "collect")

save_figure(combined, "fig1_transitions", 7.0, 6.3)
save_figure(m$plot, "fig1_transitions_math", 5.2, 4.4)
save_figure(r$plot, "fig2_transitions_read", 5.2, 4.4)

for (x in list(list("mathematics", m), list("reading", r))) {
  s <- x[[2]]
  cat(sprintf("%-12s %s pairs, %.1f%% same score, %.1f%% same job, %d values, step %.1f\n",
              x[[1]], format(s$n, big.mark = ","), 100 * s$same_score,
              100 * s$same_job, s$values, s$step))
}

# -----------------------------------------------------------------------------
# RESULT (2026-08-04)
#
# Figures 1 and 2, what the outcome does between one report and the next.
# Consecutive observed reports on the expected-career sample, N = 6,864.
#
#                  pairs   same score   same occupation   attainable   step
#   mathematics   19,384       46.9%             36.5%           33    3.0
#   reading       19,384       47.9%             36.5%           27    3.0
#
# Reading: the mass sits on the diagonal in both domains. Most students name a
# career with the same domain content as the one they named before, and those
# who move mostly move one step. One median step is three points, and that is
# the smallest move the outcome can express: a response smaller than that is
# invisible to the measure whatever the process underneath does.
#
# The ten points between 46.9% and 36.5% are students who named a different
# occupation carrying the identical score. With 33 attainable mathematics values
# over 1,327 occupation codes, collisions are common, so the outcome is coarser
# than the occupation coding it is built from. Reading is coarser still, at 27
# values, and its unchanged share is correspondingly higher.
#
# Note: an earlier version of this script drew a shaded band one step either
# side of the diagonal. It was removed. The points on the diagonal already carry
# the colour that marks them, and a band around them reads as an interval
# estimate, which it is not.
#
# Note: the occupations behind the scale were tried as a third figure, each
# occupation placed by its two scores. It failed. The named occupations crowd
# into a band of moderate mathematics and high reading content, so 313 points
# overlap in a quarter of the panel and no labelling scheme separates them
# without the labels becoming the figure. The information a reader needs from it
# is which occupations sit at which score, and a table answers that better; it
# is in 08_tables.R. The one quantity worth keeping from the attempt is that
# mathematics and reading content correlate r = .46 across the named
# occupations, which is a sentence rather than a panel.
# -----------------------------------------------------------------------------

