
source(here::here("script", "00_packages.R"))
library(patchwork)

SUFFIX <- Sys.getenv("AUM_SUFFIX", "_cov")

CELLS <- sapply(c(math_realistic  = "aum_math_realistic",
                  read_realistic  = "aum_read_realistic",
                  math_idealistic = "aum_math_idealistic",
                  read_idealistic = "aum_read_idealistic"),
                function(x) paste0(x, SUFFIX))

CELLNAME <- c(math_realistic  = "Mathematics, realistic",
              read_realistic  = "Reading, realistic",
              math_idealistic = "Mathematics, idealistic",
              read_idealistic = "Reading, idealistic")

# The figure order is the design's, idealistic above realistic, matching the
# transition figure so the two read the same way round.
PANEL_ORDER <- c("Mathematics, idealistic", "Reading, idealistic",
                 "Mathematics, realistic", "Reading, realistic")

fig_dir <- here("figure")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

draws_of <- function(cell, pars)
  as_draws_df(readRDS(here("models", paste0(CELLS[[cell]], "_m5.rds")))$draws(pars))

DRAWS <- lapply(names(CELLS), function(k)
  draws_of(k, c("beta_a", "beta_sc", "b_pos", "b_neg")))
names(DRAWS) <- names(CELLS)

RED  <- "#CC2936"
BLUE <- "#6497B1"
PURPLE <- "#482677"
GOLD   <- "#E6A817"

theme_paper <- function(base_size = 10) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
      axis.line = element_blank(),
      axis.ticks = element_line(colour = "black", linewidth = 0.3),
      axis.text = element_text(size = 8, colour = "black"),
      axis.title = element_text(size = 8.5),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 8, colour = "black",
                                 margin = margin(l = 2)),
      legend.margin = margin(0, 0, 0, 0),
      legend.key.height = unit(1.0, "lines"),
      legend.key.spacing.y = unit(3, "pt"),
      strip.background = element_blank(),
      strip.text = element_text(size = 8, colour = "black", hjust = 0.5,
                                margin = margin(b = 2)),
      plot.margin = margin(4, 4, 4, 4)
    )
}

LW    <- 0.45   # series line
EBW   <- 0.4    # interval line
CAP   <- 0.10   # interval cap, as a fraction of the x spacing
PT    <- 1.3    # marker
RIB   <- 0.12   # band alpha, low enough that the red does not read as pink

save_figure <- function(plot, name, width, height) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), plot,
         width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(file.path(fig_dir, paste0(name, ".png")), plot,
         width = width, height = height, units = "in", dpi = 300,
         device = ragg::agg_png)
  ggsave(file.path(fig_dir, paste0(name, ".tiff")), plot,
         width = width, height = height, units = "in", dpi = 300,
         device = ragg::agg_tiff, compression = "lzw")
}

qs <- function(x) tibble(m = mean(x), lo = quantile(x, .025), hi = quantile(x, .975))

# --- Figure 2: the two signals, wave by wave ----------------------------------

WAVE <- 2:6                                   # index into beta_a, waves 5 to 9
XLAB <- c("Gr 9\n(W5)", "Gr 9\n(W6)", "Gr 10", "Gr 11", "Gr 12")

signal_rows <- function(k) {
  d <- DRAWS[[k]]
  a <- bind_rows(lapply(WAVE, function(t)
    qs(d[[sprintf("beta_a[%d]", t)]]) |> mutate(x = t)))
  s <- bind_rows(lapply(c(2, 6), function(t)
    qs(d[[sprintf("beta_sc[%d]", t)]]) |> mutate(x = t)))
  bind_rows(mutate(a, signal = "Achievement"),
            mutate(s, signal = "Self-concept")) |>
    mutate(signal = factor(signal, levels = c("Achievement", "Self-concept")),
           cell = factor(CELLNAME[[k]], levels = PANEL_ORDER))
}

sig <- bind_rows(lapply(names(CELLS), signal_rows))

NUDGE <- setNames(c(-0.07, 0.07), levels(sig$signal))
sig <- mutate(sig, xp = x + NUDGE[as.character(signal)])

fig2 <- ggplot(sig, aes(xp, m, colour = signal, fill = signal)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey55",
             linewidth = 0.4) +
  geom_ribbon(data = filter(sig, signal == "Achievement"),
              aes(ymin = lo, ymax = hi), alpha = RIB, colour = NA) +
  geom_line(data = filter(sig, signal == "Achievement"), linewidth = LW) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = CAP, linewidth = EBW) +
  geom_point(size = PT) +
  scale_colour_manual(values = c(Achievement = PURPLE,
                                 `Self-concept` = GOLD)) +
  scale_fill_manual(values = c(Achievement = PURPLE,
                               `Self-concept` = GOLD)) +
  scale_x_continuous(breaks = WAVE, labels = XLAB) +
  facet_wrap(vars(cell), nrow = 2) +
  labs(x = NULL, y = "Coefficient on the signal") +
  theme_paper() +
  theme(legend.position = "right",
        legend.direction = "vertical",
        axis.text.x = element_text(size = 7.5, lineheight = 0.9),
        panel.spacing = unit(0.7, "lines"))

save_figure(fig2, "fig3_signals", 6.5, 4.6)

# --- Figure 3: the response to a change, by direction --------------------------

change_rows <- function(k) {
  d <- DRAWS[[k]]
  bind_rows(lapply(WAVE, function(t) bind_rows(
    tibble(v =  d[[sprintf("beta_a[%d]", t)]] + d[["b_pos[1]"]],
           dir = "Improved", x = t),
    # times minus one: the movement, not the coefficient
    tibble(v = -(d[[sprintf("beta_a[%d]", t)]] + d[["b_neg[1]"]]),
           dir = "Slipped", x = t)))) |>
    mutate(cell = factor(CELLNAME[[k]], levels = PANEL_ORDER))
}

chg <- bind_rows(lapply(names(CELLS), change_rows)) |>
  mutate(dir = factor(dir, levels = c("Improved", "Slipped"))) |>
  group_by(cell, dir, x) |>
  summarise(qs(v), .groups = "drop") |>
  mutate(xp = x + if_else(dir == "Improved", -0.07, 0.07))

fig3 <- ggplot(chg, aes(xp, m, colour = dir, fill = dir)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey55",
             linewidth = 0.4) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = RIB, colour = NA) +
  geom_line(linewidth = LW) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = CAP, linewidth = EBW) +
  geom_point(size = PT) +
  scale_colour_manual(values = c(Improved = RED, Slipped = BLUE)) +
  scale_fill_manual(values = c(Improved = RED, Slipped = BLUE)) +
  scale_x_continuous(breaks = WAVE, labels = XLAB) +
  facet_wrap(vars(cell), nrow = 2) +
  labs(x = NULL, y = "Standardised aspiration movement") +
  theme_paper() +
  theme(legend.position = "right",
        legend.direction = "vertical",
        axis.text.x = element_text(size = 7.5, lineheight = 0.9),
        panel.spacing = unit(0.7, "lines"))

save_figure(fig3, "fig2_asymmetry", 6.5, 4.6)

# What the figures assert, printed so a caption can be checked against it.
for (k in names(CELLS)) {
  d <- DRAWS[[k]]
  gap <- function(t) mean(d[[sprintf("beta_a[%d]", t)]] -
                          d[[sprintf("beta_sc[%d]", t)]] > 0)
  cat(sprintf("%-24s Pr(achievement > self-concept): G9 %.2f  G12 %.2f | ",
              CELLNAME[[k]], gap(2), gap(6)))
  cat(sprintf("Pr(aspiration moves up | improved) %.2f  Pr(moves down | slipped) %.2f\n",
              mean(d[["beta_a[6]"]] + d[["b_pos[1]"]] > 0),
              mean(-(d[["beta_a[6]"]] + d[["b_neg[1]"]]) < 0)))
}

