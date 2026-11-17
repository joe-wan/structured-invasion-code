#!/usr/bin/env Rscript
# 01_toy_models/01_simple_model.R
# =================================
# Author: Joe Wan
# Implements the two-stage toy model (main text eq. 1a--1c) and makes figures:
# Fig. 1b (initial population dynamics) and Fig. 2 (dynamics under invasion)


# 0. Define model and parameters 
# ==============================
# (only this section will be loaded when script is sourced)
simple.stage.model <- function(t, N, params) {
    # Unpack state variables
    N11 <- N[1]
    N12 <- N[2]
    N2  <- N[3]

    # Calculate competition-adjusted self-transition rates
    r1.adj <- params$r1 - params$alpha11 * N12 - params$alpha12 * N2
    r2.adj <- params$r2 - params$alpha21 * N12 - params$alpha22 * N2

    # Calculate derivatives
    dN11 <- r1.adj * N11 + params$f1 * N12
    dN12 <- params$g1 * N11 - params$m1 * N12
    dN2  <- r2.adj * N2

    # Return list of derivatives
    list(c(dN11, dN12, dN2))
}

simple.stage.params <- list(
    r1 = -1,
    f1 = 3, g1 = 1, m1 = 0.5,
    r2 = 1,
    alpha11 = 5, alpha12 = 3,
    alpha21 = 0.5, alpha22 = 1
)


######## BEGIN MAIN CODE BLOCK ########
# Run only if script is directly executed
if (sys.nframe() == 0L) {
#######################################


# General data manipulation
library(dplyr)
library(tidyr)
# ODE solving
library(deSolve)
# Data visualization
library(ggplot2)
library(ggh4x)
library(cowplot)


# 1. Settings and parameters
# ==========================
# Analysis parameters
initial.population <- 1e-3        # Initial total population
adult.proportions <- c(0, 0.5, 1) # Initial adult proportions
adult.proportion.alphas <- setNames(
    seq(0.3, 1, length.out=length(adult.proportions)), adult.proportions)
time.max <- 100                   # Simulation time
time.step <- time.max/2e3         # Time resolution

# Source shared code for visual settings
source("../00_plot_settings.R")
if (!dir.exists(out.dir)) {
    dir.create(out.dir, recursive = T)
}
segment.length <- 0.6     # Length of line segment for initial growth
population.unit <- 10^-3  # Unit of abundance for y-axis
breakpoint <- 2.5         # Breakpoint in x-axis
plot.time.init <- 1       # Maximum time for x-axis in Fig. 1b
plot.time.max <- 80       # Maximum time for x-axis in Fig. 2
# Higher resolution for interpolation to get closer to x-axis crossing
interpolation.times <- seq(0, breakpoint, by=1e-3)


# 2. Simulate population dynamics
# ===============================
# Create solutions for different adult proportions
all.solutions <- data.frame(initial.adult.proportion = adult.proportions) %>%
    rowwise() %>%
    mutate(sol=list(as.data.frame(
        ode(
            y = c(N11 = initial.population*(1-initial.adult.proportion), 
                N12 = initial.population*initial.adult.proportion, 
                # Resident (sp. 2) has carrying capacity of 1
                N2 = 1),
            times = seq(0, time.max, by = time.step),
            func = simple.stage.model,
            parms = simple.stage.params))))

# Reformat into data frame for plotting
solution.df <- all.solutions %>%
    unnest(sol) %>%
    pivot_longer(cols = -c(time, initial.adult.proportion), names_to = "variable", values_to = "abundance")

# Calculate adult proportion and total abundance
stage.df <- solution.df %>%
    select(initial.adult.proportion, time, variable, abundance) %>%
    pivot_wider(names_from = variable, values_from = abundance) %>%
    mutate(N1=N11 + N12) %>%
    mutate(adult.proportion=N12 / N1)


# 3. Plot initial population dynamics (Fig. 1b)
# =============================================
# Helper to calculate effective per-capita growth
total.growth <- function(N11, N12, N2, per.cap=F) {
    dN <- simple.stage.model(0, c(N11, N12, N2), simple.stage.params)[[1]]
    result <- (dN[1] + dN[2])
    if (per.cap) {
        N1 <- N11 + N12
        result <- result / N1
    }
    return(result)
}

# Calculate initial per-capita growth rates
initial.growth.df <- data.frame(initial.adult.proportion = adult.proportions) %>%
    rowwise() %>%
    mutate(gr1=total.growth(
        initial.population * (1 - initial.adult.proportion),
        initial.population * initial.adult.proportion,
        1, T)) %>%
    # Add start and end of a line segment
    mutate(
        theta=atan(gr1),
        y0=initial.population,
        x1=segment.length*cos(theta),
        y1=initial.population + segment.length*sin(theta)*population.unit)

# Make plot
initial.plot <- stage.df %>%
    select(initial.adult.proportion, time, N11, N12, N1) %>%
    pivot_longer(cols = -c(time, initial.adult.proportion), names_to = "variable", values_to = "abundance") %>%
    filter(time <= plot.time.init) %>%
ggplot() +
    geom_line(linewidth = plot.linewidth, aes(x = time, y = abundance/population.unit, linetype=variable), color=invader.color) +
    facet_wrap(~paste0("adult prop.\n= ", initial.adult.proportion),
        # Repeat y axis
        scales="free_y") +
    geom_segment(linewidth = plot.linewidth, aes(x=0, y=y0/population.unit, xend=x1, yend=y1/population.unit),
        data=initial.growth.df,
        arrow=arrow(length=unit(2, "mm")),
        alpha=0.5) +
    coord_cartesian(
        xlim=c(0, plot.time.init), 
        ylim=c(0, 2e-3/population.unit), expand=F) +
    scale_x_continuous(breaks=seq(0, plot.time.init, by=1)) +
    scale_y_continuous(breaks=seq(0, 2, by=1)) +
    scale_linetype_manual("stage",
        breaks=c("N11", "N12", "N1"), 
        labels=c(
            expression("juvenile," ~  N[1]^"(1)"), 
            "adult," ~ N[1]^"(2)", 
            "total"), 
        values=c("11", "longdash", "solid")) +
    ylab(expression("abundance " / 10 ^ -3)) +
    draft.theme
if (interactive()) print(initial.plot)
ggsave(paste0(out.dir, "01_b_initial_population_dynamics.pdf"), initial.plot, width=two.column.width, height=two.column.width/2.5)


# 4. Plot dynamics under invasion (Fig. 2)
# ========================================
# Fig. 2a: Adult proportion
# -------------------------
# Filter the needed data
proportion.df <- stage.df %>%
    select(initial.adult.proportion, time, adult.proportion)
# Data for placing slashes at breakpoint
proportion.breaks.df <- proportion.df %>%
    filter(time==breakpoint) %>%
    select(initial.adult.proportion, time, adult.proportion) %>%
    # Add x-axis breaks
    bind_rows(mutate(., adult.proportion=0))

prop.plot <- proportion.df %>%
    mutate(facet=(time>breakpoint)) %>%
    filter(time<=time.max) %>%
ggplot() +
    geom_line(linewidth = plot.linewidth, aes(x=time, y=adult.proportion, alpha=as.factor(initial.adult.proportion)), color=invader.color) +
    coord_cartesian(ylim=c(0, 1), expand=F, clip="off") +
    facet_wrap(~facet, scales="free_x") +
    # Mark the end of the first facet with a slash
    geom_text(aes(x=breakpoint, y=adult.proportion),
        label="/",  
        data=proportion.breaks.df) +
    geom_blank(aes(x=x, y=y),
        data=data.frame(x=0, y=0, facet=FALSE)) +
    force_panelsizes(cols = c(2, 3)) +
    scale_x_continuous(
        breaks=c(seq(0, breakpoint, by=1),
            seq(time.max, breakpoint, by=-20)),
        labels=scales::number_format(accuracy=1)) +
    scale_y_continuous(breaks = c(0, 0.5, 1),
        labels = c("0.0", "0.5", "1.0")) +
    scale_alpha_manual("initial adult proportion", values = adult.proportion.alphas,
        guide = guide_legend(nrow = 1, title.position = "top", title.hjust = 0.5)) +
    ylab("adult proportion") +
    draft.theme +
    theme(
        legend.position = "bottom",
        legend.justification = "center",
        legend.box.just = "center",
        # No x label
        axis.title.x = element_blank(),
        # No facet labels
        strip.background = element_blank(),
        strip.text = element_blank())
if (interactive()) print(prop.plot)
ggsave(paste0(out.dir, "02_a_adult_proportion.pdf"), prop.plot, 
    width=two.column.width, height=two.column.width/2)

# Fig. 2b: Per-capita growth
# --------------------------
# Calculate per-capita growth rates
growth.df <- stage.df %>%
    rowwise() %>%
    mutate(dN1=total.growth(N11, N12, N2), gr1=dN1/N1)
# Split into pre- and post-breakpoint for plotting
growth.df.pre <- growth.df %>% 
    filter(time <= breakpoint)
growth.df.post <- growth.df %>% 
    filter(time > breakpoint)

# Interpolate pre-breakpoint data to get closer to x-axis crossing
growth.df.pre <- nest(growth.df.pre, data = -initial.adult.proportion) %>%
    rowwise() %>%
    mutate(interp=list(data.frame(
        time=interpolation.times,
        gr1=approx(data$time, data$gr1, xout=interpolation.times)$y
    ))) %>%
    select(-data) %>%
    unnest(interp)

# Data for placing slashes at breakpoint
growth.breaks.df <- growth.df %>%
    filter(time==breakpoint, initial.adult.proportion==0) %>%
    select(time, gr1) %>%
    bind_rows(mutate(., gr1=0))

# Make plot
growth.plot <- bind_rows(growth.df.pre, growth.df.post) %>%
    mutate(facet=(time>breakpoint)) %>%
    filter(time<=plot.time.max, gr1>=0, gr1<=1) %>%
ggplot() +
    geom_line(linewidth = plot.linewidth, aes(x=time, y=gr1, alpha=as.factor(initial.adult.proportion)), color=invader.color) +
    coord_cartesian(ylim=c(0, 1), expand=F, clip="off") +
    scale_y_continuous(breaks = c(0, 0.5, 1),
        labels = c("0.0", "0.5", "1.0")) +
    facet_wrap(~facet, scales="free_x") +
    # Mark the end of the first facet with a slash
    geom_text(aes(x=breakpoint, y=gr1),
        label="/",  
        data=growth.breaks.df) +
    geom_blank(aes(x=x, y=y),
        data=data.frame(x=0, y=0, facet=FALSE)) +
    force_panelsizes(cols = c(2, 3)) +
    # Labeller: no decimals
    scale_x_continuous(
        breaks=c(seq(0, breakpoint, by=1),
            seq(plot.time.max, breakpoint, by=-20)),
        labels=scales::number_format(accuracy=1)) +
    scale_alpha_manual("initial adult proportion", values = adult.proportion.alphas,
        guide = guide_legend(nrow = 1, title.position = "top", title.hjust = 0.5)) +
    ylab("per-cap. growth") +
    draft.theme +
    theme(
        legend.position = "bottom",
        legend.justification = "center",
        legend.box.just = "center",
        # No x label
        axis.title.x = element_blank(),
        # No facet labels
        strip.background = element_blank(),
        strip.text = element_blank())
if (interactive()) growth.plot
ggsave(paste0(out.dir, "02_b_per_capita_growth.pdf"), growth.plot, 
    width=two.column.width, height=two.column.width/2)

# Fig. 2c: Invasion abundance
# ---------------------------
# Get abundance data
abundance.df <- stage.df %>%
    select(initial.adult.proportion, time, N1, N2) %>%
    pivot_longer(cols = -c(time, initial.adult.proportion), names_to = "species", values_to = "abundance")
# Split into pre- and post-breakpoint for plotting
abundance.df.pre <- abundance.df %>% 
    filter(time <= breakpoint)
abundance.df.post <- abundance.df %>% 
    filter(time > breakpoint)

# Data for placing slashes at breakpoint
abundance.breaks.df <- abundance.df.pre %>%
    filter(time==breakpoint) %>%
    select(initial.adult.proportion, species, time, abundance) %>%
    bind_rows(mutate(., abundance=0))

# Make plot
inv.plot <- bind_rows(abundance.df.pre, abundance.df.post) %>%
    filter(time<=plot.time.max) %>%
    mutate(facet=(time>breakpoint)) %>%
ggplot() +
    geom_line(linewidth = plot.linewidth, aes(x = time, y = abundance, alpha=as.factor(initial.adult.proportion), color=species)) +
    facet_wrap(~facet, scales="free_x") +
    force_panelsizes(cols = c(2, 3)) +
    coord_cartesian(ylim=c(0,1), expand=F, clip="off") +
    scale_y_continuous(breaks = c(0, 0.5, 1),
        labels = c("0.0", "0.5", "1.0")) +
    # Mark the end of the first facet with a slash
    geom_text(aes(x=breakpoint, y=abundance),
        label="/",  
        data=abundance.breaks.df) +
    scale_x_continuous(
        breaks=c(seq(0, breakpoint, by=1),
            seq(plot.time.max, breakpoint, by=-20)),
        labels=scales::number_format(accuracy=1)) +
    # Don't show alpha scale
    scale_alpha_manual("initial\nadult prop.", values = adult.proportion.alphas, guide="none") +
    scale_color_manual(values=c(invader.color, resident.color), labels=c("species 1", "species 2")) +
    draft.theme + theme(
        legend.position = "bottom",
        legend.justification = "center",
        legend.box.just = "center",
        # No facet labels
        strip.background = element_blank(),
        strip.text = element_blank())
if (interactive()) print(inv.plot)
ggsave(paste0(out.dir, "02_c_invasion_abundance.pdf"), inv.plot, 
    width=two.column.width, height=two.column.width/2)

######### END MAIN CODE BLOCK #########
# Run only if script is directly executed
}
#######################################
