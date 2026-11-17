#!/usr/bin/env Rscript
# 01_toy_models/02_cyclic_model.R
# =================================
# Author: Joe Wan
# Implements the model with facilitation-mediated cyclic resident dynamics
# (main text eq. 12) and makes figures: Fig. 4b (cyclic resident dynamics),
# Fig. 4c (initial population dynamics), and Fig. 5 (dynamics under invasion)


# 0. Define model and parameters 
# ==============================
# (only this section will be loaded when script is sourced)
cyclic.model <- function(t, N, params) {
    # Unpack state variables
    N11 <- N[1]
    N12 <- N[2]
    N21 <- N[3]
    N22 <- N[4]

    # Calculate competition-adjusted self-transition rates
    r1.adj <- params$r1 - params$alpha11 * N12 - params$alpha12 * N22
    r2.adj <- params$r2 - params$alpha21 * N12 - params$alpha22 * N22

    # Calculate additional facilitative effects (sp. 2 only)
    r2.adj <- r2.adj + params$beta22 * N21
    g2.adj <- params$g2 + params$gamma22 * N21

    # Calculate derivatives
    dN11 <- r1.adj * N11 + params$f1 * N12
    dN12 <- params$g1 * N11 - params$m1 * N12
    dN21 <- r2.adj * N21 + params$f2 * N22
    dN22 <- g2.adj * N21 - params$m2 * N22

    # Return list of derivatives
    list(c(dN11, dN12, dN21, dN22))
}

cyclic.params <- list(
    r1=-1, f1=3, g1=1, m1=0.5,
    r2= 1, f2=3, g2=1, m2=5,
    alpha11=5, alpha12=0.5,
    alpha21=1, alpha22=5,
    beta22=5, gamma22=3
)


######## BEGIN MAIN CODE BLOCK ########
# Run only if script is directly executed
if (sys.nframe() == 0L) {
#######################################


# ODE solving
library(deSolve)
# General data manipulation
library(dplyr)
library(tidyr)
library(stringr)
# Data visualization
library(ggplot2)
library(cowplot)
library(ggh4x)
library(scales)


# 1. Settings and parameters
# ==========================
# Analysis parameters
initial.population <- 1e-2         # Initial total population
adult.proportions <- c(0, 0.5, 1)  # Initial adult proportions
time.resident <- 20                # Time for initial burn-in simulation
time.max <- 40                     # Simulation time
time.step <- time.max/4e3          # Time resolution
# Source shared code for visual settings
source("../00_plot_settings.R")
if (!dir.exists(out.dir)) {
    dir.create(out.dir, recursive=T)
}
# Invasion analysis options
invasion.population <- 1e-3              # Initial invader population
invasion.phases <- c(0, 0.25, 0.5, 0.75) # Phases of resident cycle
invasion.adult.proportion <- 0.4         # Initial adult proportion of invader
invasion.periods <- 40                   # Number of periods to simulate
# Graph options
breakpoint.periods <- 3.25                   # Periods to draw slash mark at
panel.props <- c(1, 1)                       # Ratio of facets before/after break
phase.labels <- c("0", "T/4", "T/2", "3T/4") # Labels for invasion phases
phase.labels <- factor(phase.labels, levels=phase.labels)
phase.alphas <- setNames(seq(0.3, 1, length.out=length(phase.labels)), phase.labels)
min.growth.rate <- -0.7                      # Minimum per-capita growth rate
plot.time.init <- 1                         # Initial dynamics window in periods
segment.length <- 0.3                      # Growth arrow length in plot units
population.unit <- 1e-3                    # Abundance unit for Fig. 4c


# 2. Simulate and plot resident dynamics
# ======================================
resident.solution <- ode(y=c(N11=0, N12=0, N21=initial.population, N22=0), 
    times=seq(0, time.resident, by=time.step), 
    func=cyclic.model, parms=cyclic.params)
resident.solution.df <- as.data.frame(resident.solution) %>%
    mutate(N2=N21 + N22) %>%
    pivot_longer(cols=-time, names_to="variable", values_to="abundance") %>%
    filter(str_starts(variable, "N2"))

# Find local minima and calculate period
min.indices <- which(diff(sign(diff(resident.solution[, "N22"]))) > 0) + 1
min.times <- resident.solution[min.indices, "time"]
period <- min.times[length(min.times)] - min.times[length(min.times) - 1]

# Fig. 4b: Resident dynamics
# --------------------------
resident.plot <- resident.solution.df %>%
    mutate(time=time-1.5) %>%
    filter(str_starts(variable, "N2"), time>=0, time/period < 6) %>%
    ggplot() +
    geom_path(aes(x=time/period, y=abundance, linetype=variable), color=resident.color) +
    scale_linetype_manual("stage",
        breaks=c("N21", "N22", "N2"),
        values=line.types,
        labels=c(
            expression("juvenile," ~ N[2]^"(1)"),
            expression("adult," ~ N[2]^"(2)"),
            "total")) +
    xlab("time / T") +
    coord_cartesian(xlim=c(0, 6), ylim=c(0, 20), expand=FALSE, clip=F) +
    theme_cowplot()
if (interactive()) print(resident.plot)
ggsave(paste0(out.dir, "04_b_cyclic_resident.pdf"), resident.plot, 
    width=two.column.width, height=two.column.width/1.5)

# 3. Invasion simulations
# =======================
# Begin by simulating resident dynamics
start.times <- min.times[length(min.times) - 1] + invasion.phases * period
names(start.times) <- phase.labels
resident.times <- seq(0, max(start.times) + period, by=time.step)
resident.solution.extended <- ode(y=c(N11=0, N12=0, N21=initial.population, N22=0), 
    times=resident.times, func=cyclic.model, parms=cyclic.params)
resident.solution.extended.df <- as.data.frame(resident.solution.extended)

# Function to interpolate resident state for start time
interpolate.resident.state <- function(t) {
    return(sapply(c("N11", "N12", "N21", "N22"), 
        function(v) approx(resident.solution.extended.df$time, 
            resident.solution.extended.df[[v]], xout=t)$y))
}

# Perturbation for adding invader population at start time
perturbation <- c(
    N11=invasion.population*(1-invasion.adult.proportion), 
    N12=invasion.population*invasion.adult.proportion, 
    N21=0, N22=0)

# Function to simulate one invasion
simulate.invasion <- function(t0, label) {
    resident.state <- interpolate.resident.state(t0)
    y0 <- as.numeric(resident.state + perturbation)

    # Integrate forward from t0
    sol <- ode(y=y0, times=seq(t0, t0 + time.max, by=time.step), 
        func=cyclic.model, parms=cyclic.params)

    # Convert to data frame
    result <- as.data.frame(sol)
    colnames(result) <- c("time", "N11", "N12", "N21", "N22")
    result$phase <- label
    result$t0 <- t0
    result$time.relative <- result$time - t0
    return(result)
}

invasions.df <- bind_rows(lapply(seq_along(start.times), function(i) {
    simulate.invasion(start.times[i], names(start.times)[i])
}))
invasions.df.tall <- invasions.df %>% 
    mutate(N1=N11 + N12, N2=N21 + N22) %>%
    pivot_longer(cols=c(N1, N2), names_to="species", values_to="abundance")


# 4. Plot initial invasion dynamics (Fig. 4c)
# ===========================================
# Stage abundance data frame
initial.lines.df <- invasions.df %>%
    filter(time.relative <= plot.time.init * period) %>%
    mutate(N1=N11 + N12, phase=factor(phase, levels=phase.labels)) %>%
    select(phase, time.relative, N11, N12, N1) %>%
    pivot_longer(cols=c(N11, N12, N1),
        names_to="variable", values_to="abundance")

# Calculate arrows for initial growth rate
initial.growth.df <- invasions.df %>%
    group_by(phase) %>%
    arrange(time, .by_group=TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(
        phase=factor(phase, levels=phase.labels),
        N1=N11 + N12,
        r1.adj=cyclic.params$r1 - cyclic.params$alpha11 * N12
            - cyclic.params$alpha12 * N22,
        dN11=r1.adj * N11 + cyclic.params$f1 * N12,
        dN12=cyclic.params$g1 * N11 - cyclic.params$m1 * N12,
        theta=atan((dN11 + dN12) * period / population.unit),
        y0=N1,
        x1=segment.length * cos(theta),
        y1=N1 + segment.length * sin(theta) * population.unit)

# Make four-panel plot
initial.plot <- ggplot(initial.lines.df) +
    geom_line(aes(x=time.relative / period, y=abundance / population.unit,
        linetype=variable), color=invader.color, linewidth=plot.linewidth) +
    facet_wrap(~phase, scales="free_y", nrow=1,
        labeller=as_labeller(c(
            "0"="t[0] == 0", "T/4"="t[0] == T/4",
            "T/2"="t[0] == T/2", "3T/4"="t[0] == 3*T/4"),
            default=label_parsed)) +
    geom_segment(aes(x=0, y=y0 / population.unit,
        xend=x1, yend=y1 / population.unit),
        data=initial.growth.df, linewidth=plot.linewidth,
        arrow=arrow(length=grid::unit(2, "mm")), alpha=0.5) +
    coord_cartesian(xlim=c(0, plot.time.init), ylim=c(0, 1.5), expand=FALSE) +
    scale_x_continuous(breaks=seq(0, plot.time.init, by=0.5)) +
    scale_linetype_manual("stage",
        breaks=c("N11", "N12", "N1"),
        values=setNames(line.types, c("N11", "N12", "N1")),
        labels=c(
            expression("juvenile," ~ N[1]^"(1)"),
            expression("adult," ~ N[1]^"(2)"),
            "total")) +
    xlab(expression("time (periods) since introduction," ~ (t - t[0])/T)) +
    ylab(expression("abundance " / 10^-3)) +
    draft.theme
if (interactive()) print(initial.plot)
ggsave(paste0(out.dir, "04_c_initial_population_dynamics.pdf"), initial.plot,
    width=full.width, height=full.width/2.5)


# 5. Plot invasion dynamics (Fig. 5)
# ==================================
# Fig. 5a: Adult proportions
# --------------------------
proportion.df <- invasions.df %>%
    mutate(N1=N11 + N12,
        adult.proportion=ifelse(N1 > 0, N12 / N1, NA_real_)) %>%
    select(time, phase, adult.proportion)

# Data frames for extra visual elements
proportion.breakpoints.df <- proportion.df %>%
    group_by(phase) %>%
    summarize(adult.proportion=approx(time, adult.proportion, xout=start.times[1] + breakpoint.periods*period, ties="ordered")$y,
        .groups="drop") %>%
    bind_rows(mutate(., adult.proportion=invasion.adult.proportion))
proportion.points.df <- invasions.df %>%
    group_by(phase) %>%
    slice(1) %>%
    mutate(facet=F)

# Plot the adult proportions
proportion.plot <-
    proportion.df %>%
    filter(time <= start.times[1] + invasion.periods*period) %>%
    mutate(facet=(time > start.times[1] + breakpoint.periods*period)) %>%
    ggplot() +
    geom_line(aes(x=(time - start.times[1]) / period, y=adult.proportion, alpha=phase), color=invader.color) +
    facet_wrap(~facet, scales="free_x") +
    force_panelsizes(cols=panel.props) +
    coord_cartesian(ylim=c(0.4, 0.7 + 1e-4), expand=FALSE, clip="off") +
    geom_text(aes(x=breakpoint.periods, y=adult.proportion), data=proportion.breakpoints.df, label="/") +
    geom_point(aes(x=(time - start.times[1]) / period, y=invasion.adult.proportion, alpha=phase),
        shape=16, size=2, color=invader.color,
        data=proportion.points.df) +
    scale_x_continuous(breaks=c(
            seq(0, breakpoint.periods, by=1),
            seq(invasion.periods, breakpoint.periods, by=-10)),
        labels=number_format(accuracy=1)) +
    # Force axis to reach invasion.periods*period
    geom_blank(aes(x=invasion.periods, y=invasion.adult.proportion), data=data.frame(facet=TRUE)) +
    scale_alpha_manual(name="intro. time", values=phase.alphas,
        limits=phase.labels, breaks=phase.labels, labels=phase.labels) +
    theme_cowplot() +
    theme(strip.background=element_blank(),
        strip.text=element_blank(), axis.title.x=element_blank()) +
    ylab("adult proportion")
if (interactive()) print(proportion.plot)
ggsave(paste0(out.dir, "05_a_cyclic_adult_proportion.pdf"), proportion.plot,
    width=two.column.width, height=two.column.width/1.5)

# Fig. 5b.1: Per-capita growth
# ----------------------------
growth.df <- invasions.df %>%
    mutate(
        # competition-adjusted rate for species 1
        r1.adj=(cyclic.params$r1
            - cyclic.params$alpha11 * N12
            - cyclic.params$alpha12 * N22),
        dN11=r1.adj * N11 + cyclic.params$f1 * N12,
        dN12=cyclic.params$g1 * N11 - cyclic.params$m1 * N12,
        N1=N11 + N12,
        growth.rate=ifelse(N1 > 0, (dN11 + dN12) / N1, NA_real_)
    ) %>%
    select(time, phase, growth.rate)

# Centered one-period average using a rounded number of sampled points.
growth.first <- growth.df %>%
    filter(phase == phase.labels[1]) %>%
    arrange(time)
average.n <- max(1L, round(period / median(diff(growth.first$time))))
growth.average <- growth.first %>%
    mutate(growth.rate=as.numeric(stats::filter(growth.rate, rep(1 / average.n, average.n),
        sides=2
    ))) %>%
    filter(!is.na(growth.rate), time >= start.times[1] + period, time <= start.times[1] + invasion.periods*period) %>%
    select(time, growth.rate)
growth.average <- growth.average %>% mutate(facet=time > start.times[1] + breakpoint.periods*period)

# Slash mark at start.times[1] + breakpoint.periods*period, drawn at x=breakpoint.periods (relative axis)
growth.breakpoints.df <- growth.df %>%
    group_by(phase) %>%
    summarize(
        growth.rate=approx(time, growth.rate, xout=start.times[1] + breakpoint.periods*period, ties="ordered")$y,
        .groups="drop"
    ) %>%
    bind_rows(mutate(., growth.rate=min.growth.rate))

# Plot the per-capita growth rates
growth.plot <-
    growth.df %>%
    filter(time <= start.times[1] + invasion.periods*period) %>%
    mutate(facet=(time > start.times[1] + breakpoint.periods*period)) %>%
    filter(!is.na(growth.rate)) %>%
    ggplot() +
    geom_line(aes(x=(time - start.times[1]) / period, y=growth.rate, alpha=phase), color=invader.color) +
    facet_wrap(~facet, scales="free_x") +
    force_panelsizes(cols=panel.props) +
    coord_cartesian(ylim=c(min.growth.rate, 1), expand=FALSE, clip="off") +
    geom_text(aes(x=breakpoint.periods, y=growth.rate), data=growth.breakpoints.df, label="/") +
    geom_point(aes(x=(time - start.times[1]) / period, y=growth.rate, alpha=phase),
        shape=16, size=2, color=invader.color,
        show.legend=FALSE,
        data=mutate(slice(group_by(growth.df, phase), 1), facet=F)
    ) +
    scale_x_continuous(
        breaks=c(
            seq(0, breakpoint.periods, by=1),
            seq(invasion.periods, breakpoint.periods, by=-10)
        ),
        labels=number_format(accuracy=1)
    ) +
    scale_alpha_manual(
        name="intro. time", values=phase.alphas,
        # Force order
        limits=phase.labels,
        breaks=phase.labels,
        labels=phase.labels
    ) +
    geom_hline(yintercept=0, linetype="dashed", color="gray") +
    # Force axis to reach invasion.periods*period
    geom_blank(aes(x=invasion.periods, y=0), data=data.frame(facet=TRUE)) +
    theme_cowplot() +
    theme(
        strip.background=element_blank(),
        strip.text=element_blank(),
        axis.title.x=element_blank()
    ) +
    ylab("per-cap. growth") +
    geom_line(aes(x=(time - start.times[1]) / period, y=growth.rate),
        data=growth.average, inherit.aes=FALSE,
        color="black", linetype="dashed", linewidth=0.75
    )
if (interactive()) print(growth.plot)
ggsave(paste0(out.dir, "05_b_1_cyclic_per_capita_growth.pdf"), growth.plot,
    width=two.column.width, height=two.column.width/1.5)

# Fig. 5b.2: Invasion abundances
# -------------------------------
# Get invader abundances at introduction for the start-time points.
invasion.points.df <- invasions.df.tall %>% 
    filter(time <= start.times[1] + invasion.periods*period, species == "N1") %>%
    group_by(phase) %>%
    slice(1) %>%
    mutate(facet=F)
# Values at the breakpoint.periods*period for the little "/" mark
invasion.breakpoints.df <- invasions.df.tall %>%
    group_by(phase, species) %>%
    summarize(
        abundance=approx(time, abundance, xout=start.times[1] + breakpoint.periods*period, ties="ordered")$y,
        .groups="drop"
    ) %>%
    bind_rows(mutate(., abundance=0)) # draw a "/" from value to 0 as in Sample 1

# Plot the invasion dynamics
abundance.plot <- invasions.df.tall %>%
    filter(time <= start.times[1] + invasion.periods*period) %>%
    mutate(facet=(time > start.times[1] + breakpoint.periods*period)) %>%
ggplot() +
    geom_line(aes(x=(time-start.times[1])/period, 
        y=abundance, color=species, group=interaction(species, phase))) +
    facet_wrap(~facet, scales="free_x") +
    force_panelsizes(cols=panel.props) +
    coord_cartesian(ylim=c(0, 12), expand=FALSE, clip="off") +
    # Add slashes at breakpoint.periods*period
    geom_text(aes(x=breakpoint.periods, y=abundance), 
        data=invasion.breakpoints.df, label="/") +
    geom_point(aes(x=(time-start.times[1])/period, y=abundance, alpha=phase), 
        shape=16, size=2, color=invader.color, data=invasion.points.df) +
    scale_x_continuous(breaks=c(seq(0, breakpoint.periods, by=1), 
        seq(invasion.periods, breakpoint.periods, by=-10)),
        labels=number_format(accuracy=1)) +
    # Force axis to go to tmax
    geom_blank(aes(x=invasion.periods*period/period, y=0),
        data=data.frame(facet=T)) +
    scale_alpha_manual(values=phase.alphas, limits=phase.labels, guide="none") +
    scale_color_manual(values=c("N1"=invader.color, "N2"=resident.color),
        labels=c("species 1", "species 2")) +
    theme_cowplot() +
    theme(strip.background=element_blank(), strip.text=element_blank(), 
        axis.title.x=element_blank()) +
    ylab("abundance")
if (interactive()) print(abundance.plot)
ggsave(paste0(out.dir, "05_b_2_cyclic_invasion_abundance.pdf"), abundance.plot, 
    width=two.column.width, height=two.column.width/1.5)


######### END MAIN CODE BLOCK #########
# Run only if script is directly executed
}
#######################################
