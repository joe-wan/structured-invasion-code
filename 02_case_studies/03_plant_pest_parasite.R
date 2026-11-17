#!/usr/bin/env Rscript
# 02_case_studies/03_plant_pest_parasite.R
# ========================================
# Authors: Guan-Yu Chen, Hsiang-Chih Lo
# Implements the plant-pest-parasite simulation (Appendix S6) and makes
# Figure 3f, showing invasion of a pest species. Also performs an extra
# analysis for parasite invasion (Appendix S6, Figure S6.1).

library(tidyverse)
library(deSolve)
library(ggplot2)
library(cowplot)

# Source the helper function for IGR calculation
source("00_simulation_igr.R")

# Source shared code for visual settings
source("../00_plot_settings.R")
if (!dir.exists(out.dir)) {
    dir.create(out.dir, recursive=T)
}


# 1. Define the ODE models and their parameters
# =============================================
# The full plant-pest-parasite model with two pest species
tritrophic.model <- function(times, state, parms) {
    with(as.list(c(state, parms)), {
        dH_dt = h1*o1*P1H + h2*o2*P2H - (c1*b1*P1 + c2*b2*P2)*H - d*H
        dP1H_dt = b1*P1*H - (o1 + m1)*P1H
        dP2H_dt = b2*P2*H - (o2 + m2)*P2H
        dP1_dt = e1*a1*P1*S - b1*P1*H + e1H*a1*P1H*S - m1*P1
        dP2_dt = e2*a2*P2*S - b2*P2*H + e2H*a2*P2H*S - m2*P2
        dS_dt = r*S*(1 - S/K) - (a1*P1 + a2*P2 + a1*P1H + a2*P2H)*S
        return(list(c(dH_dt, dP1H_dt, dP2H_dt, dP1_dt, dP2_dt, dS_dt)))
    })
}

# Above model, with only one pest species (P1)
tritrophic.model.P1 <- function(times, state, parms) {
    with(as.list(c(state, parms)), {
        dH_dt = h1*o1*P1H - c1*b1*P1*H - d*H
        dP1H_dt = b1*P1*H - (o1 + m1)*P1H
        dP1_dt = e1*a1*P1*S - b1*P1*H + e1H*a1*P1H*S - m1*P1
        dS_dt = r*S*(1 - S/K) - (a1*P1 + a1*P1H)*S 
        return(list(c(dH_dt, dP1H_dt, dP1_dt, dS_dt)))
    })
}

# Template state variable vectors
state.vector <- c(H = 0.1, P1H = 0, P2H = 0, P1 = 0.01, P2 = 0.01, S = 0.1) # Multi-pest invasion
state.vector.P1 <- c(H = 0.1, P1H = 0, P1 = 0.01, S = 0.1) # Parasite invasion
# States classified as invader or resident
state.type <- c("res", "inv", "res", "inv", "res", "res") # Multi-pest invasion
state.type.P1 <- c("inv", "inv", "res", "res") # Parasite invasion

# Default model parameters
base.params <- list(
    r = 1, K = 10,
    a1 = 0.38, a2 = 0.51, e1 = 0.5, e2 = 0.5,
    b1 = 0.2, b2 = 0.42, m1 = 0.06, m2 = 0.05, e1H = 0.5, e2H = 0.5,
    o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.02)

# Values of focal parameter (b1) for continuous analytic and pointwise numeric IGR
b1.values <- seq(0.05, 0.35, by = 0.025)
b1.focal <- seq(0.15, 0.35, by = 0.025)
# Values of a1 for variant analysis
a1.values <- c(0.3, 0.4)

# 2. Calculate analytic IGR
# =========================
# Invasion analysis for pest 1 (P1) invading pest 2 (P2)/parasite (H) system
# --------------------------------------------------------------------------
analytic.pest.IGR <- function(a1, b1) {
    params <- base.params
    params$a1 <- a1
    params$b1 <- b1
    with(params, {
        # Calculate roots to find state variables at invasion
        P2 <- d*(o2 + m2)/(b2*h2*o2 - c2*b2*(o2 + m2))
        A <- -(e2H*a2^2*b2^2*P2*K / (r*(o2 + m2)^2))
        B <- ((e2H*a2*b2 - e2H*a2^2*b2*P2) / (o2 + m2) - 
            e2*a2^2*b2*P2 / (r*(o2 + m2)))*K - b2
        C <- (e2*a2 - e2*a2^2*P2/r)*K - m2
        E <- B^2 - 4*A*C
        if (E < 0 || is.na(E)) {
            H <- 0
            S <- m2/(e2*a2)
            P2 <- r/a2*(1 - S/K)
        } else {
            H1 <- (-B + sqrt(B^2 - 4*A*C)) / (2*A)
            H2 <- (-B - sqrt(B^2 - 4*A*C)) / (2*A)
            H_12 <- c(H1, H2)
            H_12 <- H_12[H_12 > 0 & !is.na(H_12) & is.finite(H_12)]
            if (length(H_12) == 0) {
                # If no positive roots, set H = 0 because the S* from E1 = S*
                H <- 0
                S <- m2/(e2*a2)
                P2 <- r/a2*(1 - S/K)
            } else {
                # Smallest root is the biologically meaningful one
                H <- min(H_12)
                S <- (b2*H + m2)*(m2 + o2) / (e2*a2*(m2 + o2) + e2H*a2*b2*H)
            }
        }

        # Construct transition matrix (C, the invader Jacobian at invasion)
        transition.matrix <- matrix(
            data = c(-(m1 + o1), b1*H, e1H*a1*S, e1*a1*S - b1*H - m1),
            nrow = 2, ncol = 2, byrow = T)

        # Return the leading eigenvalue of the matrix
        return(max(eigen(transition.matrix)$values))
    })
}

# Set up data frame for analytic IGR calculation
pest.IGR.df <- expand.grid(b1 = b1.values, a1 = a1.values)
# Apply the IGR calculation function
pest.IGR.df$analytic.IGR <- mapply(analytic.pest.IGR, 
    pest.IGR.df$a1, pest.IGR.df$b1)

# Invasion analysis for parasite (H) invading pest 1 (P1)
# -------------------------------------------------------
analytic.parasite.IGR <- function(a1, b1) {
    params <- base.params
    params$a1 <- a1
    params$b1 <- b1
    with(params, {
        # Calculate the resident pest population at invasion
        P1 <- r/a1*(1 - m1/(e1*a1*K))

        # Construct transition matrix (C, the invader Jacobian at invasion)
        transition.matrix <- matrix(
            data = c(-(c1*b1*P1 + d), h1*o1, b1*P1, -(m1 + o1)),
            nrow = 2, ncol = 2, byrow = T)

        # Return the leading eigenvalue of the matrix
        return(max(eigen(transition.matrix)$values))
    })
}

# Set up data frame for analytic IGR calculation
parasite.IGR.df <- expand.grid(b1 = b1.focal, a1 = a1.values)
# Apply the IGR calculation function
parasite.IGR.df$analytic.IGR <- mapply(analytic.parasite.IGR, 
    parasite.IGR.df$a1, parasite.IGR.df$b1)


# 3. Numerical IGR calculation using the helper function
# ======================================================
# Initial invader abundance for IGR calculation
initial.invader.abundance <- 10^-10
# Vector of output times for IGR calculation
time <- seq(0, 1250, by = 0.1)
# Initial second-stage proportions
initial.proportions <- c(0.4, 0.5, 0.6)

# Pest 1 (P1) invading pest 2 (P2)/parasite (H) system
# ----------------------------------------------------
# Helper function to apply the simulation.IGR function
simulation.pest.IGR.helper <- function(a1, b1) {
    params <- base.params
    params$a1 <- a1
    params$b1 <- b1
    return(simulation.IGR(
        state.vector = state.vector,
        state.type = state.type,
        model = tritrophic.model,
        params = params,
        initial.invader.abundance = initial.invader.abundance,
        time = time,
        initial.resident.vector = c(H = 1, P2H = 1, P2 = 1, S = 1),
        initial.proportions = initial.proportions))
}
# Apply with mapply
pest.IGR.df$numeric.IGR <- mapply(simulation.pest.IGR.helper,
    pest.IGR.df$a1, pest.IGR.df$b1)

# Parasite (H) invading pest 1 (P1)
# ---------------------------------
# Helper function to apply the simulation.IGR function
simulation.parasite.IGR.helper <- function(a1, b1) {
    params <- base.params
    params$a1 <- a1
    params$b1 <- b1
    return(simulation.IGR(
        state.vector = state.vector.P1,
        state.type = state.type.P1,
        model = tritrophic.model.P1,
        params = params,
        initial.invader.abundance = initial.invader.abundance,
        time = time,
        initial.resident.vector = c(P1 = 1, S = 1),
        initial.proportions = initial.proportions))
}
# Apply with mapply
parasite.IGR.df$numeric.IGR <- mapply(simulation.parasite.IGR.helper,
    parasite.IGR.df$a1, parasite.IGR.df$b1)


# 4. Plot the results
# ===================
# Main pest (P1) analysis (Fig. 3f)
# ---------------------------------
ggplot(pest.IGR.df) +
    geom_line(mapping = aes(x = b1, y = analytic.IGR, color = factor(a1)),
        linewidth = case.study.linewidth) +
    geom_point(mapping = aes(x = b1, y = analytic.IGR, color = factor(a1)),
        size = case.study.point.size) +
    labs(x = expression(beta[1]),
        y = "IGR",
        color = expression(alpha[1])) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0.159, linetype = "dotted") +
    geom_vline(xintercept = 0.262, linetype = "dotted") +
    scale_color_manual(values = c(dark.green, light.green)) +
    case.study.theme +
    theme(legend.position.inside = c(0.72, 0.8))
ggsave(paste0(out.dir, "03_f_plant_pest_parasite.pdf"),
    width = case.study.width, height = case.study.height,
    units = case.study.units, dpi = case.study.dpi)

# Supplementary parasite (H) analysis (Appendix S6, Figure S6.1)
# --------------------------------------------------------------
ggplot(parasite.IGR.df) +
    geom_line(mapping = aes(x = b1, y = analytic.IGR, color = factor(a1)),
        linewidth = case.study.linewidth) +
    geom_point(mapping = aes(x = b1, y = analytic.IGR, color = factor(a1)),
        size = case.study.point.size) +
    labs(x = expression(beta), y = "IGR",
        color = expression(alpha)) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0.207, linetype = "dotted") +
    geom_vline(xintercept = 0.273, linetype = "dotted") +
    scale_color_manual(values = c(dark.green, light.green)) +
    case.study.theme +
    theme(legend.position.inside = c(0.68, 0.2))
ggsave(paste0(out.dir, "S6_1_parasite_invasion.pdf"),
    width = case.study.width, height = case.study.height,
    units = case.study.units, dpi = case.study.dpi)
