#!/usr/bin/env Rscript
# 01_toy_models/03_floquet_theory.R
# =================================
# Author: Joe Wan
# Implements Floquet analysis of the facilitation-mediated cyclic resident model
# (main text eq. 12) and makes Fig. 6 (Floquet decomposition)

# Matrix logarithm
library(expm)
# General data manipulation
library(dplyr)
library(tidyr)
# ODE solving
library(deSolve)
# Data visualization
library(ggplot2)
library(ggh4x)
library(cowplot)
# Import the cyclic model code
source("02_cyclic_model.R")


# 1. Settings and parameters
# ==========================
# Analysis parameters
period.rough <- 0.78
optim.tolerance <- 0.5
spline.resolution <- 1e5
initial.stages <- c(N11=1, N12=0)
# Floquet decomposition plot
invasion.population <- 1e-3
plot.periods <- 4
simulation.periods <- 10
points.per.period <- 100
# Floquet M and B comparison plot
trajectory.periods <- 4.5
# Source shared code for visual settings
source("../00_plot_settings.R")
if (!dir.exists(out.dir)) {
    dir.create(out.dir, recursive = T)
}


# 2. High-quality period and minimum of resident cycle
# ====================================================
# Optimize period and initial state on limit cycle
# ------------------------------------------------
# Rough estimate of resident state on limit cycle
N2.rough <- ode(y=c(N11=0, N12=0, N21=1e-3, N22=0), times=c(0, 10), 
    func=cyclic.model, parms=cyclic.params)[2, 4:5]
# Function (squared distance to initial point) to be minimized
period.objective <- function(N21, N22, period) {
    initial.state <- c(N11=0, N12=0, N21=N21, N22=N22)
    solution <- ode(y=initial.state, times=c(0, period), func=cyclic.model, parms=cyclic.params)
    N <- solution[2, 4:5]
    return(sum((N - N2.rough)^2))
}
# Optimize the period and initial state to find a periodic solution
par.start <- c(N22=unname(N2.rough[2]), period=period.rough)
par.cycle <- optim(par=par.start, fn=function(par) period.objective(N2.rough[1], par[1], par[2]), method="L-BFGS-B", 
    lower=par.start + c(-optim.tolerance, -optim.tolerance), 
    upper=par.start + c(optim.tolerance, optim.tolerance))$par
N2.cycle <- c(N21=unname(N2.rough[1]), N22=unname(par.cycle[1]))
period <- unname(par.cycle['period'])

# Optimize to find adult minimum (N22) on the cycle
# -------------------------------------------------
# Function (adult abundance at time t) to be minimized
minimum.objective <- function(time) {
    initial.state <- c(N11=0, N12=0, N2.cycle)
    solution <- ode(y=initial.state, 
        times=c(0, time), func=cyclic.model, parms=cyclic.params)
    return(unname(solution[2, "N22"]))
}
# Find rough minimum
times <- seq(0, period, length.out=100)
min.time.rough <- times[which.min(sapply(times, minimum.objective))]
# Optimize to find the minimum
min.time <- optim(par=min.time.rough, fn=minimum.objective, method="L-BFGS-B", lower=0, upper=period)$par
N2.min <- ode(y=c(N11=0, N12=0, N2.cycle), times=c(0, min.time), func=cyclic.model, parms=cyclic.params)[2, 4:5]

# # Make spline of the resident dynamics
# # ------------------------------------
# spline.times <- seq(0, period*1.1, by=period/spline.resolution)
# sol.fine <- ode(y=c(N11=0, N12=0, N2.min), times=spline.times, func=cyclic.model, parms=cyclic.params)
# spline.N21 <- splinefun(sol.fine[, "time"], sol.fine[, "N21"], method="natural")
# spline.N22 <- splinefun(sol.fine[, "time"], sol.fine[, "N22"], method="natural")
# N2.spline <- function(time) {
#     return(c(N21=spline.N21(time %% period), N22=spline.N22(time %% period)))
# }

# 3. Define ODE system for Floquet analysis
# =========================================
# Resident (N21, N22) variables advance as if the invasion is absent using the
# original function, while the invader (N11, N12) advances according to Jacobian
# C at the current resident state.
floquet.ode <- function(t, state, parameters) {
    # Unpack state:
    y <- matrix(state[1:2], nrow = 2)
    N21 <- state[3]
    N22 <- state[4]

    # Implement derivatives from original function
    dNs <- cyclic.model(t, c(0, 0, N21, N22), parameters)[[1]]
    dN21 <- dNs[3]
    dN22 <- dNs[4]

    # Get transition matrix
    C <- matrix(c(
        parameters$r1 - parameters$alpha12 * N22,
        parameters$f1,
        parameters$g1,
        -parameters$m1
    ), nrow = 2, byrow = TRUE)

    # Matrix derivatives
    dy <- C %*% y

    # Unpack matrix derivatives into vector
    list(c(as.vector(dy), dN21, dN22))
}

# 4. General helpers for Floquet analysis
# =======================================
# Helper to integrate a matrix initial condition for resident (part of Algorithm S3.7)
# ------------------------------------------------------------------------------------
integrate.matrix.ODE <- function(initial.matrix, time, resident=NULL) {
    n <- ncol(initial.matrix)
    N2 <- if (is.null(resident)) N2.min else resident
    result <- matrix(NA, nrow=n, ncol=n)
    for (i in 1:n) {
        y0 <- c(initial.matrix[, i])
        sol <- ode(y=c(y0, N2), times=c(0, time), 
            func=floquet.ode, parms=cyclic.params)
        result[, i] <- sol[2, 2:3]
    }
    return(result)
}

# Algorithm S3.5: Helper to find monodromy matrix for arbitrary initial time
# --------------------------------------------------------------------------
# Note: this code is simpler than the pseudocode in S3.5 because the concatenation
# and simulation are implemented in the integrate.matrix.ODE helper
monodromy.matrix <- function(t0=0) {
    # Advance the resident to t0
    N2 <- N2.min
    if (t0 > 0) {
        sol <- ode(y=c(0, 0, N2), times=c(0, t0), func=cyclic.model, parms=cyclic.params)
        N2 <- sol[2, 4:5]
    }
    # Integrate for one period
    identity.matrix <- diag(1, 2, 2)
    M <- integrate.matrix.ODE(identity.matrix, period, resident=N2)
    return(M)
}

# Algorithm S3.6: Floquet calculation of invasion growth
# ------------------------------------------------------
# Note: different t0 values give consistent IGRs (within numerical error)
floquet.IGR <- function(t0=0) {
    M <- monodromy.matrix(t0)
    spectrum <- eigen(M)$values
    rho <- max(Re(spectrum))
    return(log(rho) / period)
}

# Algorithm S3.7: Finding the Floquet exponent matrix B
# -----------------------------------------------------
floquet.B <- function(t0=0) {
    # For illustration, we follow pseudocode by using inverse to compute M
    identity.matrix <- diag(1, 2, 2)
    Phi0 <- integrate.matrix.ODE(identity.matrix, t0)
    PhiT <- integrate.matrix.ODE(identity.matrix, t0 + period)
    M <- PhiT %*% solve(Phi0)
    # # Equivalent to above, directly reusing the monodromy matrix helper:
    # M <- monodromy.matrix(t0)
    B <- logm(M) / period
    return(B)
}

# Algorithm S3.8: Periodic Floquet transformation
# -----------------------------------------------
floquet.transformation <- function(times, t0=0) {
    identity.matrix <- diag(1, 2, 2)
    B <- floquet.B(t0)

    # Advance the resident to t0
    N2 <- N2.min
    if (t0 > 0) {
        sol <- ode(y=c(0, 0, N2), times=c(0, t0), func=cyclic.model, parms=cyclic.params)
        N2 <- sol[2, 4:5]
    }
    
    # Back-calculate transformation for each specified time
    result.list <- vector("list", length(times))
    for (time.index in seq_along(times)) {
        t <- times[time.index]
        Phit <- integrate.matrix.ODE(identity.matrix, t, resident=N2)
        Yt <- expm(B * (t - t0))
        Pt <- Phit %*% solve(Yt)
        result.list[[time.index]] <- Pt
    }
    # Unpack result if only one time was requested
    if (length(times) == 1) {
        return(result.list[[1]])
    }
    return(result.list)
}


# 5. Make Floquet decomposition plot (Fig. 6a)
# ============================================
# Get the Floquet exponent matrix
B <- floquet.B()
invasion.growth <- floquet.IGR()

# Eigenvalue analysis of B
spectrum <- eigen(B)
leading.index <- which.max(Re(spectrum$values))
# Calculate leading eigenvector (\vec{w})
leading.eigenvector <- spectrum$vectors[, leading.index]
leading.eigenvector <- leading.eigenvector / sum(leading.eigenvector)
# # We can confirm that associated eigenvalue is consistent with helper's IGR
# Re(spectrum$values[leading.index]) - invasion.growth
# # = 0

# Calculate the periodic transformation, P(t)
cycle.times <- seq(0, period, length.out=points.per.period + 1)
transformation.matrices <- floquet.transformation(seq(0, period, 
    length.out=points.per.period + 1))
# P(t) * w gives the periodic stage structure variation
periodic.fluctuation <- vapply(transformation.matrices, 
    function(P) as.vector(P %*% leading.eigenvector), 
    numeric(2))

# Simulate the linearized invasion, starting with juveniles only.
invasion.times <- seq(0, simulation.periods * period,
    length.out=simulation.periods * points.per.period + 1)
invasion.solution <- ode(y=c(invasion.population * initial.stages, N2.min),
    times=invasion.times, func=floquet.ode, parms=cyclic.params)
invasion.df <- as.data.frame(invasion.solution)

# Match the asymptotic solution at the final simulated juvenile abundance.
cycle.indices <- (seq_along(invasion.times) - 1L) %% points.per.period + 1L
periodic.solution <- periodic.fluctuation[, cycle.indices]
exponential.solution <- exp(invasion.growth * invasion.times)
final.index <- length(invasion.times)
combined.scaling <- invasion.df$N11[final.index] /
    (exponential.solution[final.index] * periodic.solution[1, final.index])
# Normalize the periodic stages so their mean total abundance is one.
periodic.scaling <- 1 / mean(colSums(periodic.solution))

# Combine all components into a single data frame for plotting
decomposition.wide.df <- invasion.df %>%
    transmute(time,
        N11.obs=N11, N12.obs=N12,
        N11.by=periodic.scaling * periodic.solution[1, ],
        N12.by=periodic.scaling * periodic.solution[2, ],
        N11.exp=exponential.solution * combined.scaling / periodic.scaling,
        N12.exp=0)
decomposition.df <- decomposition.wide.df %>%
    pivot_longer(cols=-time, names_to=c("variable", "type"),
        names_sep="\\.", values_to="abundance")
decomposition.df <- bind_rows(decomposition.df,
    decomposition.df %>%
        group_by(time, type) %>%
        summarize(variable="N1", abundance=sum(abundance), .groups="drop")) %>%
    filter(type != "exp" | variable == "N1")

# Invisible points to keep the y-axis scales nice
component.types <- c("obs", "exp", "by")
component.points.df <- tribble(
    ~type, ~abundance,
    "obs", 3,
    "exp", 0,
    "exp", 3,
    "by", 0.4,
    "by", 0.6 + 1e-3) %>%
    mutate(type=factor(type, levels=component.types), variable="N1")
component.labels <- c(
    obs="full dynamics",
    exp="exponential growth",
    by="periodic variation")

# Fig. 6a: Floquet decomposition
# -------------------------------
floquet.plot.full <- decomposition.df %>%
    filter(type %in% component.types, time <= plot.periods * period) %>%
    filter(type != "by" | variable != "N1") %>%
    mutate(type=factor(type, levels=component.types),
        abundance=if_else(type == "by", abundance, abundance / invasion.population)) %>%
    ggplot() +
    geom_line(aes(x=time / period, y=abundance, linetype=variable),
        color=invader.color, linewidth=plot.linewidth) +
    facet_wrap(~type, scales="free_y", ncol=3,
        labeller=as_labeller(component.labels)) +
    scale_linetype_manual("stage",
        breaks=c("N11", "N12", "N1"),
        values=setNames(line.types, c("N11", "N12", "N1")),
        labels=c(
            expression("juvenile," ~ N[1]^"(1)"),
            expression("adult," ~ N[1]^"(2)"),
            "total")) +
    geom_point(aes(x=0, y=abundance), size=0, data=component.points.df) +
    coord_cartesian(expand=FALSE, clip="off") +
    xlab("time / T") +
    ylab("component") +
    theme_cowplot() +
    theme(strip.background=element_blank(), strip.clip="off")
if (interactive()) print(floquet.plot.full)
ggsave(paste0(out.dir, "06_a_floquet_decomposition.pdf"), floquet.plot.full,
    width=full.width, height=full.width/3)


# 6. Plot M and B matrices from Floquet theory (Fig. 6b)
# ======================================================
# Helpers to project using M and B matrices
project.M <- function(periods) {
    as.vector((M %^% periods) %*% initial.stages)
}
project.B <- function(periods) {
    as.vector(expm(B * periods * period) %*% initial.stages)
}

# Get monodromy matrix
M <- monodromy.matrix()
# Calculate the Floquet exponent matrix from the monodromy matrix
B <- logm(M) / period
# (equivalent to using the helper floquet.B() above)

# Calculate discrete-time trajectory generated by monodromy matrix
M.times <- seq(0, floor(trajectory.periods))
M.stages <- vapply(M.times, project.M, numeric(2))
M.trajectory.df <- data.frame(time=M.times,
    N11=M.stages[1, ], N12=M.stages[2, ]) %>%
    pivot_longer(cols=c(N11, N12), names_to="variable", values_to="abundance")

# Calculate continuous-time trajectory generated by Floquet exponent matrix
B.times <- seq(0, trajectory.periods, by=1 / points.per.period)
B.stages <- vapply(B.times, project.B, numeric(2))
B.trajectory.df <- data.frame(time=B.times,
    N11=B.stages[1, ], N12=B.stages[2, ]) %>%
    pivot_longer(cols=c(N11, N12), names_to="variable", values_to="abundance")

# Format original invasion trajectory by reusing simulation result from above
trajectory.panels <- c("discrete", "continuous")
trajectory.labels <- c(
    discrete="atop('discrete-time reduction', '(monodromy matrix ' * bold(M) * ')')",
    continuous="atop('continuous-time reduction', '(exponent matrix ' * bold(B) * ')')")
trajectory.lines.df <- bind_rows(
    invasion.df %>%
        transmute(time=time / period,
            N11=N11 / invasion.population, N12=N12 / invasion.population) %>%
        pivot_longer(cols=c(N11, N12), names_to="variable", values_to="abundance") %>%
        mutate(panel=trajectory.panels[1]),
    B.trajectory.df %>% mutate(panel=trajectory.panels[2])) %>%
    filter(time <= trajectory.periods) %>%
    mutate(panel=factor(panel, levels=trajectory.panels))
trajectory.points.df <- bind_rows(
    M.trajectory.df %>% mutate(panel=trajectory.panels[1]),
    M.trajectory.df %>% mutate(panel=trajectory.panels[2])) %>%
    mutate(panel=factor(panel, levels=trajectory.panels))

# Plot discrete- and continuous-time trajectories
trajectory.plot <- ggplot(trajectory.lines.df, aes(x=time, y=abundance)) +
    geom_line(aes(linetype=variable), color=invader.color, linewidth=plot.linewidth) +
    # Draw filled adults first, then white juveniles over them.
    geom_point(aes(shape=variable),
        data=trajectory.points.df %>% arrange(variable == "N11"),
        size=2, color=invader.color, fill="white") +
    facet_wrap(~panel, nrow=1, axes="all_y", axis.labels="all_y",
        labeller=as_labeller(trajectory.labels, default=label_parsed)) +
    scale_linetype_manual("stage",
        breaks=c("N11", "N12"),
        values=setNames(line.types[1:2], c("N11", "N12")),
        labels=c(
            expression("juvenile," ~ N[1]^"(1)"),
            expression("adult," ~ N[1]^"(2)"))) +
    scale_shape_manual("stage",
        breaks=c("N11", "N12"),
        values=c(N11=21, N12=16),
        labels=c(
            expression("juvenile," ~ N[1]^"(1)"),
            expression("adult," ~ N[1]^"(2)"))) +
    scale_y_continuous(breaks=c(0, 0.5, 1)) +
    coord_cartesian(xlim=c(0, trajectory.periods), ylim=c(0, NA),
        expand=FALSE, clip="off") +
    xlab("time / T") +
    ylab("abundance") +
    theme_cowplot() +
    theme(strip.background=element_blank(), strip.clip="off")
if (interactive()) print(trajectory.plot)
ggsave(paste0(out.dir, "06_b_floquet_M_B_comparison.pdf"), trajectory.plot,
    width=full.width, height=full.width/2)
