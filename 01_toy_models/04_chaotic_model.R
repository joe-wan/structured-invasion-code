#!/usr/bin/env Rscript
# 01_toy_models/04_chaotic_model.R
# ================================
# Author: Joe Wan
# Implements the appendix model coupling the simple stage-structured invader to
# a four-species chaotic resident community, and makes figures:
# Appendix S8, Fig. S8.1 to Fig. S8.6


# 0. Define model and parameters 
# ==============================
# (only this section will be loaded when script is sourced)

# Implements a two-stage invader coupled to a four-species resident community
chaotic.model <- function(t, N, chaotic.params) {
    # Unpack state variables
    N11 <- N[1]
    N12 <- N[2]
    N2 <- N[3]
    N3 <- N[4]
    N4 <- N[5]
    N5 <- N[6]

    # Unpack parameters
    r <- chaotic.params$r
    A <- chaotic.params$A
    f1 <- chaotic.params$f1
    g1 <- chaotic.params$g1
    m1 <- chaotic.params$m1
    r1 <- r[1]
    a11 <- A[1, 1]
    a12 <- A[1, 2]
    a13 <- A[1, 3]
    a14 <- A[1, 4]
    a15 <- A[1, 5]

    # LV species
    pop_vec <- c(N11+N12, N2, N3, N4, N5)
    V <- r - A %*% pop_vec

    # Calculate derivatives
    dN11 <- (r1 - a11*N12 - a12*N2 - a13*N3 - a14*N4 - a15*N5) * N11 + f1 * N12
    dN12 <- g1 * N11 - m1 * N12
    dN2 <- N2 * V[2]
    dN3 <- N3 * V[3]
    dN4 <- N4 * V[4]
    dN5 <- N5 * V[5]

    # Return list of derivatives
    return(list(c(dN11, dN12, dN2, dN3, dN4, dN5)))
}

chaotic.params <- list(
    r=c(-1, 5, 5, 5, 5),
    A=matrix(c(5, 5, 1, 3, 3,
                0.6, 5, 10, 6.5, 0,
                0.7, 0, 5, 5, 5,
                0.8, 10, 0, 5, 0,
                0.9, 5, 0, 0, 5),
                nrow=5, byrow=T),
    f1=3, g1=1, m1=0.5)


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
resident.time <- 1000           # Time to simulate before invader introduction
invader.time <- 1500            # Time to simulate after invader introduction
time.step <- resident.time/1e4  # Time resolution
invader.population <- 1e-12     # Initial invader abundance
# Source shared code for visual settings
source("../00_plot_settings.R")
if (!dir.exists(out.dir)) {
    dir.create(out.dir, recursive=T)
}
# Graph settings
smooth.time <- 50                 # Time window for averaging population structure
smooth.ignore <- 5                # Time to exclude from averaging
invasion.time.limits <- c(0, 150) # Time limits for invasion plots
# Labels and scales
invasion.time.breaks <- seq(0, 150, by=30)
relative.abundance.label <- expression("abundance"~"/"~N[0])
stage.linetypes <- c(N11="11", N12="longdash", N1="solid")
stage.labels <- c("juvenile", "adult", "total")
species.colors <- c("N1"=invader.color, "N2"=resident.color,
    N3="#CD9EF7", N4="#9BDB4D", N5="#D4D4D4")
invader.gradient <- c("#EEEFEF", "#E3EFF1", "#D8EEF3", "#CDECF5", "#C1E9F7", "#B5E6F9", "#A9E2FA", "#9DDDFC", "#91D7FE", "#84D0FF", "#76C7FF", "#69BEFF", "#5BB2FC", "#4CA2F5", "#3C90EB", "#2C7DDF", "#1D6AD2", "#1158C4", "#0949B7", "#043DAA", "#01329E", "#002992", "#002285", "#001B79", "#00156C")


# 2. Plot resident (Fig. S8.1)
# ============================
# Initial condition near the fixed point of species 2–5
N0 <- c(0, 0, round(solve(chaotic.params$A[-1, -1], chaotic.params$r[-1]), 2))
names(N0) <- c("N11", "N12", "N2", "N3", "N4", "N5")

# Simulate resident until chaos
sol <- ode(y=N0, times=seq(0, resident.time, by=time.step), 
        func=chaotic.model, parms=chaotic.params)
resident.df <- sol %>%
    as.data.frame() %>%
    pivot_longer(cols=c(N11:N5), names_to="variable", values_to="abundance")

# Fig. S8.1a: Resident dynamics
# -----------------------------
chaotic.resident.plot <- resident.df %>%
        filter(!str_starts(variable, "N1")) %>%
ggplot(aes(x=time, y=abundance, color=variable)) +
    geom_line(linewidth=0.75) +
    coord_cartesian(xlim=c(0, 150),
        ylim=c(0, 0.9),
        expand=F) +
    scale_color_manual(values=species.colors, breaks=names(species.colors)) +
    scale_y_continuous(breaks=c(0, 0.3, 0.6, 0.9)) +
    theme_cowplot()
if (interactive()) print(chaotic.resident.plot)
ggsave(paste0(out.dir, "S8_1_a_chaotic_resident.pdf"), chaotic.resident.plot, width=two.column.width, height=two.column.width/2)

# Fig. S8.1b: Chaotic attractor
# -----------------------------
chaotic.resident.phase.plot <- resident.df %>%
        pivot_wider(names_from=variable, values_from=abundance) %>%
        filter(time > 25) %>%
ggplot(aes(x=N4, y=N3+N5)) +
    geom_path(aes(color=time), alpha=0.3) +
    scale_x_continuous(breaks=c(0, 0.25, 0.5)) +
    scale_y_continuous(breaks=c(0, 0.3, 0.6, 0.9)) +
    coord_cartesian(xlim=c(0, NA), ylim=c(0, 1.2), expand=FALSE) +
    scale_color_gradientn(colors=c("#666666", "#676767")) +
    theme_cowplot() +
    # Remove the legend for the phase plot, since it is not informative
    theme(legend.position="none")
if (interactive()) print(chaotic.resident.phase.plot)
ggsave(paste0(out.dir, "S8_1_b_chaotic_attractor.pdf"), chaotic.resident.phase.plot, width=two.column.width/2, height=two.column.width/2)

# 3. Plot invasion trajectory (Figs. S8.2 to S8.4)
# ================================================
# Simulate after invasion
invasion.df <- ode(y=as.numeric(sol[nrow(sol),-1]) + c(invader.population, 0, 0, 0, 0 ,0), 
                        times=seq(0, invader.time, by=time.step), 
                        func=chaotic.model, parms=chaotic.params)
colnames(invasion.df) <- colnames(sol)
invason.df.long <- invasion.df %>%
    as.data.frame() %>%
    mutate(N1=N11 + N12) %>%
    pivot_longer(cols=c(N11:N5, N1), names_to="variable", values_to="abundance")

# Fig. S8.2: Invasion abundance trajectory
# ----------------------------------------
chaos.invasion.plot <- invason.df.long %>%
    filter(str_length(variable) == 2) %>%
    mutate(variable=forcats::fct_relevel(variable, "N1", after=Inf)) %>%
    ggplot(aes(x=time, y=abundance, color=variable, group=variable)) +
    geom_line(aes(
        alpha=ifelse(str_starts(variable, "N1"), 1, 0.5),
        linewidth=ifelse(str_starts(variable, "N1"), 0.75, 0.2))) +
    expand_limits(y=0) +
    coord_cartesian(xlim=c(0, invader.time), expand=F) +
    scale_alpha_identity() +
    scale_linewidth_identity() +
    scale_color_manual(values=species.colors, breaks=names(species.colors)) +
    theme_cowplot()
if (interactive()) print(chaos.invasion.plot)
ggsave(paste0(out.dir, "S8_2_chaotic_invasion_abundance.pdf"), chaos.invasion.plot, width=two.column.width, height=two.column.width/2)

# Fig. S8.3a: Population structure during invasion
# ------------------------------------------------
# Calculate population structure
population.structure.df <- as.data.frame(invasion.df) %>%
    filter(time > smooth.ignore, time < invader.time) %>%
    mutate(N1=N11 + N12, adult.prop=N12/N1)

# Add centered average to smooth out chaotic fluctuations
window.length <- round(smooth.time/time.step/2)*2
population.structure.df <- population.structure.df %>%
    mutate(adult.proportion.average=
        (lead(cumsum(adult.prop), window.length/2) - lag(cumsum(adult.prop), window.length/2))
        / window.length)
# Smooth the average with a spline
adult.prop.spline <- population.structure.df %>%
    filter(!is.na(adult.proportion.average)) %>%
    (function (.) return(smooth.spline(.$time, .$adult.proportion.average, all.knots=F, nknots=20)))
# Hold the smoothed average at its boundary values outside the full window
population.structure.df$adult.prop.smooth <- predict(adult.prop.spline,
    pmax(min(adult.prop.spline$x), pmin(max(adult.prop.spline$x), population.structure.df$time)))$y

adult.prop.limits <- c(0.59, 0.69)
adult.prop.breaks <- c(0.60, 0.64, 0.68)

chaos.proportion.plot <- population.structure.df %>%
    filter(row_number() %% 2 == 0) %>%
    ggplot(aes(x=time, y=adult.prop)) +
    scale_y_continuous(breaks=adult.prop.breaks) +
    # geom_point(size=0.3, alpha=0.1) +
    geom_line(aes(color=time), linewidth=0.2) +
    geom_line(aes(y=adult.prop.smooth), linetype="dashed", na.rm=TRUE) +
    scale_color_gradientn(colors=invader.gradient, limits=c(0, invader.time)) +
    coord_cartesian(xlim=c(0, invader.time), ylim=adult.prop.limits, expand=F) +
    theme_cowplot()
if (interactive()) print(chaos.proportion.plot)
ggsave(paste0(out.dir, "S8_3_a_chaotic_adult_proportion.pdf"), chaos.proportion.plot, width=two.column.width, height=two.column.width/2)

# Fig. S8.3b: Population structure attractors
# -------------------------------------------
chaos.proportion.phase.plot <- as.data.frame(invasion.df) %>%
    filter(time > 25) %>%
    mutate(N1=N11 + N12, adult.prop=N12/N1) %>%
    ggplot(aes(x=N4, y=adult.prop)) +
    geom_path(aes(color=time), alpha=0.3) +
    scale_color_gradientn(colors=invader.gradient, limits=c(0, invader.time)) +
    scale_x_continuous(breaks=c(0, 0.25, 0.5)) +
    scale_y_continuous(breaks=adult.prop.breaks) +
    coord_cartesian(xlim=c(0, 0.55), ylim=adult.prop.limits, expand=F) +
    theme_cowplot()
if (interactive()) print(chaos.proportion.phase.plot)
ggsave(paste0(out.dir, "S8_3_b_chaotic_population_structure_attractor.pdf"), chaos.proportion.phase.plot, width=one.column.width, height=one.column.width/1.5)

# Fig. S8.4: Per-capita growth rate trajectory
# --------------------------------------------
# Calculate instantaneous growth directly from the ODE derivatives
calculate.per.capita.growth <- function(i, trajectory) {
    state <- as.numeric(trajectory[i, -1])
    derivatives <- chaotic.model(trajectory$time[i], state, chaotic.params)[[1]]
    return(sum(derivatives[1:2]) / sum(state[1:2]))
}
growth.df <- as.data.frame(invasion.df)
growth.df$per.capita.growth <- vapply(seq_len(nrow(growth.df)),
    calculate.per.capita.growth, numeric(1), trajectory=growth.df)

# Exclude the initial stage transient and add a centered average
growth.df <- growth.df %>%
    filter(time > smooth.ignore, time < invader.time) %>%
    mutate(growth.average=
        (lead(cumsum(per.capita.growth), window.length/2) - lag(cumsum(per.capita.growth), window.length/2))
        / window.length)
growth.spline <- growth.df %>%
    filter(!is.na(growth.average)) %>%
    (function (.) return(smooth.spline(.$time, .$growth.average, all.knots=F, nknots=20)))
# Hold the smoothed average at its boundary values outside the full window
growth.df$growth.smooth <- predict(growth.spline,
    pmax(min(growth.spline$x), pmin(max(growth.spline$x), growth.df$time)))$y

chaos.per.capita.growth.plot <- growth.df %>%
    ggplot(aes(x=time)) +
    geom_hline(yintercept=0, color="grey70", linewidth=0.2) +
    geom_line(aes(y=per.capita.growth), color=species.colors[["N1"]], linewidth=0.2) +
    geom_line(aes(y=growth.smooth), linetype="dashed", linewidth=0.5) +
    coord_cartesian(xlim=c(0, invader.time), expand=F) +
    labs(x=expression(time~t), y="per-capita growth rate") +
    theme_cowplot()
if (interactive()) print(chaos.per.capita.growth.plot)
ggsave(paste0(out.dir, "S8_4_chaotic_per_capita_growth.pdf"), chaos.per.capita.growth.plot, width=one.column.width, height=one.column.width/2)


# 4. Plot algorithm visualizations (Figs. S8.5 and S8.6)
# ======================================================
# Algorithm S3.2: Rescaling algorithm for Lyapunov exponent
# ---------------------------------------------------------
# Note: follows pseudocode in Appendix S3, but adds options to save trajectory or growth rates
rescaling.IGR <- function(end.time, time.step, initial.abundance,
    burn.in,
    return.growth.rates=F, return.trajectory=F,
    output.step=time.step, max.step=NULL) {
    state <- c(initial.abundance, 0, as.numeric(sol[nrow(sol), -(1:3)]))
    growth.rates <- numeric(floor(end.time/time.step))
    growth.times <- numeric(length(growth.rates))
    count <- 0
    cycle <- 0
    time <- 0

    output.times <- c(0, time.step)
    if (return.trajectory) {
        state.names <- colnames(sol)[-1]
        output.times <- seq(0, time.step, by=output.step)
        output.times <- c(output.times[output.times < time.step], time.step)
        initial <- as.data.frame(as.list(c(time=0, setNames(state, state.names))))
        initial$N1 <- initial.abundance
        initial$cycle <- cycle
        trajectory.df <- initial
    }

    while (time + time.step <= end.time) {
        solution <- ode(y=state, times=time + output.times,
            func=chaotic.model, parms=chaotic.params, hmax=max.step)
        state <- as.numeric(solution[nrow(solution), -1])
        total.abundance <- sum(state[1:2])
        growth.rate <- log(total.abundance/initial.abundance)/time.step
        if (time > burn.in) {
            count <- count + 1
            growth.rates[count] <- growth.rate
            growth.times[count] <- time + time.step
        }

        if (return.trajectory) {
            trajectory <- as.data.frame(solution[-1, , drop=F])
            names(trajectory) <- c("time", state.names)
            trajectory$N1 <- trajectory$N11 + trajectory$N12
            trajectory$cycle <- cycle
            trajectory.df <- bind_rows(trajectory.df, trajectory)
        }

        state[1:2] <- (initial.abundance/total.abundance) * state[1:2]
        time <- time + time.step
        cycle <- cycle + 1

        if (return.trajectory) {
            reset <- initial
            reset$time <- time
            reset[1, state.names] <- state
            reset$N1 <- sum(state[1:2])
            reset$cycle <- cycle
            trajectory.df <- bind_rows(trajectory.df, reset)
        }
    }

    growth.rates <- growth.rates[seq_len(count)]
    IGR <- mean(growth.rates)
    if (!return.growth.rates && !return.trajectory) return(IGR)

    result <- list(IGR=IGR)
    if (return.growth.rates) {
        result$history <- data.frame(time=growth.times[seq_len(count)], growth.rate=growth.rates)
    }
    if (return.trajectory) result$trajectory <- trajectory.df
    return(result)
}

# Fig. S8.5a: Rescaling algorithm
# --------------------------------
rescaling.demonstration <- rescaling.IGR(
    end.time=150, time.step=30, initial.abundance=invader.population, burn.in=5,
    return.trajectory=T, output.step=0.5, max.step=0.5)
rescaling.trajectory.df <- rescaling.demonstration$trajectory %>%
    pivot_longer(cols=c(N11, N12, N1), names_to="variable", values_to="abundance")

chaos.algorithm.plot <- rescaling.trajectory.df %>%
    ggplot(aes(x=time, y=abundance/invader.population, linetype=variable,
        group=interaction(variable, cycle))) +
    geom_path(color=species.colors[["N1"]], linewidth=0.75) +
    geom_hline(yintercept=1, linetype="dotted") +
    scale_linetype_manual(NULL, breaks=names(stage.linetypes),
        values=stage.linetypes, labels=stage.labels) +
    coord_cartesian(xlim=invasion.time.limits, expand=F) +
    scale_x_continuous(breaks=invasion.time.breaks) +
    labs(x=expression(time~t), y=relative.abundance.label) +
    theme_cowplot()
if (interactive()) print(chaos.algorithm.plot)
ggsave(paste0(out.dir, "S8_5_a_chaotic_rescaling_algorithm.pdf"), chaos.algorithm.plot, width=two.column.width, height=two.column.width/2)

# Fig. S8.5b: Exponential growth comparison
# -----------------------------------------
rescaling.results <- rescaling.IGR(
    end.time=invader.time, time.step=time.step, initial.abundance=invader.population,
    burn.in=5,
    return.growth.rates=T)
rescaling.df <- rescaling.results$history
rescaling.result <- rescaling.results$IGR

exponential.fit.df <- as.data.frame(invasion.df) %>%
    filter(time == 50) %>%
    mutate(N1=N11 + N12) %>%
    select(time, N1, N11, N12) %>%
    pivot_longer(cols=c(N1, N11, N12), names_to="variable", values_to="abundance") %>%
    transmute(variable, initial.abundance=abundance/exp(rescaling.result*time)) %>%
    crossing(time=seq(0, 300, by=time.step)) %>%
    transmute(time, variable, exponential=initial.abundance*exp(rescaling.result*time))
exponential.comparison.df <- as.data.frame(invasion.df) %>%
    filter(time <= 300) %>%
    mutate(N1=N11 + N12) %>%
    select(time, N1, N11, N12) %>%
    pivot_longer(cols=c(N1, N11, N12), names_to="variable", values_to="abundance") %>%
    full_join(exponential.fit.df, by=c("time", "variable")) %>%
    pivot_longer(cols=c(exponential, abundance), names_to="type", values_to="value")

chaos.exponential.plot <- exponential.comparison.df %>%
    filter(time <= invasion.time.limits[2]) %>%
    ggplot(aes(x=time, y=value/invader.population,
        linetype=if_else(type == "exponential", "N1", variable), alpha=type,
        group=interaction(variable, type))) +
    geom_line(color=species.colors[["N1"]], linewidth=0.75) +
    scale_linetype_manual(NULL, breaks=names(stage.linetypes),
        values=stage.linetypes, labels=stage.labels) +
    scale_alpha_manual(values=c(abundance=1, exponential=0.5)) +
    coord_cartesian(xlim=invasion.time.limits, ylim=c(0, 40), expand=F) +
    scale_x_continuous(breaks=invasion.time.breaks) +
    labs(x=expression(time~t), y=relative.abundance.label) +
    guides(alpha="none") +
    theme_cowplot()
if (interactive()) print(chaos.exponential.plot)
ggsave(paste0(out.dir, "S8_5_b_chaotic_exponential_comparison.pdf"), chaos.exponential.plot, width=two.column.width, height=two.column.width/2)

# Fig. S8.6: Algorithm comparison
# --------------------------------
# Interpolate the resident-only trajectory with a spline
resident.solution <- ode(y=as.numeric(sol[nrow(sol), -1]),
    times=seq(0, invader.time, by=time.step), func=chaotic.model, parms=chaotic.params)
resident.splines <- lapply(4:7, function(i) {
    return(splinefun(resident.solution[, 1], resident.solution[, i]))
})
# Define an ODE simulating invader using interpolated resident-only trajectory
linearized.invader.model <- function(t, N, params) {
    residents <- c(0, vapply(resident.splines, function(f) return(f(t)), numeric(1)))
    juvenile.growth <- params$r[1] - sum(params$A[1, ] * residents)
    return(list(c(juvenile.growth*N[1] + params$f1*N[2], params$g1*N[1] - params$m1*N[2])))
}

# Under periodic approximation, calculate Floquet exponent at each time point
comparison.times <- seq(0, invader.time, by=1)
juvenile.simulation <- ode(y=c(1, 0), times=comparison.times,
        func=linearized.invader.model, parms=chaotic.params)
adult.simulation <- ode(y=c(0, 1), times=comparison.times,
        func=linearized.invader.model, parms=chaotic.params)
# Function to get exponent at a given time index
calculate.floquet.exponent <- function(i) {
    # Per usual numerical recipe, paste together columns to get M
    monodromy <- matrix(c(juvenile.simulation[i, 2], adult.simulation[i, 2],
        juvenile.simulation[i, 3], adult.simulation[i, 3]), nrow=2)
    rho <- Re(eigen(monodromy)$values[1])
    lambda <- log(rho)/comparison.times[i]
    return(lambda)
}
# Store Floquet results
floquet.df <- data.frame(time=comparison.times,
    IGR=vapply(seq_along(comparison.times), calculate.floquet.exponent, numeric(1)))

# Reformat the rescaling results to match
scaling.df <- rescaling.df %>%
    transmute(time, IGR=cumsum(growth.rate)/row_number())
# Assemble plotting data frame
performance.colors <- c(Floquet="#0D52BF", rescaling=invader.color)
performance.df <- bind_rows(
    mutate(floquet.df, method="Floquet"),
    mutate(scaling.df, method="rescaling")) %>%
    filter(is.finite(IGR), (time >= 0 & time <= 110) | (time >= 1440 & time <= 1500)) %>%
    mutate(panel=factor(if_else(time <= 110, "early", "late"), levels=c("early", "late")))
performance.bounds <- data.frame(
    panel=factor(c("early", "early", "late", "late"), levels=c("early", "late")),
    time=c(0, 110, 1440, 1500), IGR=0)

# Plot the performance comparison
chaos.performance.plot <- performance.df %>%
    ggplot(aes(x=time, y=IGR, color=method)) +
    geom_blank(data=performance.bounds, aes(x=time, y=IGR), inherit.aes=F) +
    geom_hline(yintercept=rescaling.result, linetype="dotted", linewidth=0.4) +
    geom_line(linewidth=0.6) +
    facet_grid(cols=vars(panel), scales="free_x", space="free_x") +
    scale_color_manual(NULL, values=performance.colors, breaks=names(performance.colors)) +
    scale_x_continuous(breaks=c(0, 50, 100, 1450, 1500)) +
    scale_y_continuous(breaks=seq(0, 0.06, by=0.02)) +
    coord_cartesian(ylim=c(0, 0.07), expand=F) +
    labs(x=expression("simulation time,"~t), y=expression("estimated"~IGR[1])) +
    theme_cowplot() +
    theme(strip.background=element_blank(), strip.text=element_blank(),
        panel.spacing.x=grid::unit(8, "mm"),
        legend.position="inside", legend.position.inside=c(0.98, 0.05),
        legend.justification=c(1, 0), legend.direction="vertical",
        legend.background=element_rect(fill="white", color=NA),
        legend.margin=margin(2, 2, 2, 2), legend.key.width=grid::unit(7, "mm"))
if (interactive()) print(chaos.performance.plot)
ggsave(paste0(out.dir, "S8_6_chaotic_algorithm_comparison.pdf"), chaos.performance.plot, width=two.column.width, height=two.column.width/2)

######### END MAIN CODE BLOCK #########
# Run only if script is directly executed
}
#######################################
