#!/usr/bin/env Rscript
# 00_plot_settings.R
# ==================================
# Author: Joe Wan
# Shared plot settings for all figures in the manuscript.

# Output directory
out.dir <- file.path(dirname(normalizePath(sys.frame(1)$ofile)), "outputs", "")

# 1. Visual options for toy models (01_toy_models/)
# =================================================
invader.color <- "#64BAFF"
resident.color <- "#FFA154"
draft.theme <- theme_cowplot(font_size = 11) +
    theme(text = element_text(size = 11),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 11),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          strip.text = element_text(size = 11)) +
    theme(strip.background = element_blank())
plot.linewidth <- 0.75
line.types <- c("11", "longdash", "solid")

# Plot units and dimensions
mm.to.in <- 1/25.4
full.width <- 173*mm.to.in
two.column.width <- 110*mm.to.in
one.column.width <- 82*mm.to.in


# 2. Visual options for case studies (02_case_studies/)
# =====================================================
light.green <- "#99C55B"
light.green.faded <- "#D9F2A9"
dark.green <- "#177A55"
dark.green.faded <- "#A3C9B9"
case.study.theme <- theme_cowplot() +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          legend.position = "inside",
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12),
          legend.key.width = unit(2, "cm"))
case.study.linewidth <- 1
case.study.point.size <- 2.5

# Plot export size (in cm)
case.study.width <- 13
case.study.height <- 12
case.study.units <- "cm"
case.study.dpi <- 800
