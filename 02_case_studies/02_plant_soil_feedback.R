#!/usr/bin/env Rscript
# 02_case_studies/02_plant_soil_feedback.R
# ========================================
# Authors: Shiang-Chi Huang, Jia-Zhen Liu, Chin-Te Tsai, Hsiang-Chih Lo
# Implements the plant-soil feedback simulation (Appendix S5) and
# makes Figure 3d.

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


# 1. Define the ODE model and its parameters
# ==========================================
psf.model <- function(times, state, parms) {
    with(as.list(c(state, parms)), {
        # Combined abundance variables for colonization
        P1 <- sigma_f11*P11 + sigma_f12*P12
        P2 <- sigma_f22*P22 + sigma_f21*P21
        # Calculate derivatives
        dP00_dt <- d1*P01 + d2*P02  - r1*P1*P00 - r2*P2*P00
        dP01_dt <- m1*sigma_m11*P11 + m2*sigma_m21*P21 - 
            (r1*sigma_c11*P1 + r2*sigma_c21*P2)*P01 - d1*P01
        dP02_dt <- m1*sigma_m12*P12 + m2*sigma_m22*P22 - 
            (r1*sigma_c12*P1 + r2*sigma_c22*P2)*P02 - d2*P02
        dP11_dt <- r1*P1*P00 + r1*sigma_c11*P1*P01 + c1*P12 - m1*sigma_m11*P11
        dP12_dt <- r1*sigma_c12*P1*P02 - c1*P12 - m1*sigma_m12*P12
        dP22_dt <- r2*P2*P00 + r2*sigma_c22*P2*P02 + c2*P21 - m2*sigma_m22*P22
        dP21_dt <- r2*sigma_c21*P2*P01 - c2*P21 - m2*sigma_m21*P21
        # Return list of derivatives
        return(list(c(dP00_dt, dP01_dt, dP02_dt, dP11_dt, dP12_dt, dP22_dt, dP21_dt)))
    })
}

# Template state variable vector
state.vector <- c(P00 = 1, P01 = 1, P02 = 1, P11 = 1, P12 = 1, P22 = 1, P21 = 1)
# States classified as invader or resident
state.type <- c("res", "inv", "res", "inv", "inv", "res", "inv")

# Default model parameters
base.params <- list(
  r1 = 0.5, r2 = 0.5, 
  m1 = 0.1, m2 = 0.1, 
  sigma_m11 = 1, sigma_m22 = 1,
  sigma_m12 = 1, sigma_m21 = 1,
  sigma_f11 = 1, sigma_f22 = 1,
  sigma_f12 = 1, sigma_f21 = 1,
  sigma_c11 = 1, sigma_c22 = 0.5,
  sigma_c12 = 1, sigma_c21 = 1,
  d1 = 0.5, d2 = 0.1,
  c1 = 0.5, c2 = 0.5
)

# Values of focal parameter (sigma_c22) for continuous analytic and pointwise numeric IGR
sigma_c22.values <- seq(0, 3, length.out = 1000)
sigma_c22.focal <- seq(0, 3, length.out = 11)
# Values of d2 for variant analyses
d2.values <- c(0.1, 0.9)


# 2. Calculate analytic IGR (Appendix S5)
# =======================================
# Helper function for analytic IGR calculation
analytic.IGR <- function(parms) {
    # Calculate transformed parameters
    parms.setting <- with(as.list(parms), {
        c(beta1 = sigma_m11/sigma_m12,
          beta2 = sigma_m22/sigma_m21,
          gamma1 = sigma_f11/sigma_f12,
          gamma2 = sigma_f22/sigma_f21,
          k1 = c1/(r1*sigma_f12),
          k2 = c2/(r2*sigma_f21),
          delta1 = d1/(r1*sigma_f12),
          delta2 = d2/(r2*sigma_f21),
          mu1 = m1*sigma_m12/(r1*sigma_f12),
          mu2 = m2*sigma_m21/(r2*sigma_f21))
    })

    # Remaining calculations should have access to all parameters
    with(as.list(c(parms, parms.setting)), {
        # Calculate equilibrium of P02
        G <- 1/2*(delta2/(gamma2*sigma_c22*(1 - sigma_c22)) +
            beta2*mu2/(gamma2*sigma_c22) +
            (gamma2 - beta2*mu2)/(gamma2*(1 - sigma_c22)))
        sqrt.term <- sqrt(G^2 - beta2*mu2*(1 - beta2*mu2/gamma2)/
            (gamma2*sigma_c22*(1 - sigma_c22)))
        P02_star <- G + ifelse(sigma_c22 < 1, -1, 1) * sqrt.term
        # If no equilibrium exists, return NA
        if (P02_star < 0 | is.na(P02_star)) return(NA)
        
        # Calculate associated stage structure
        P00_star <- (beta2*mu2 - gamma2*sigma_c22*P02_star)/gamma2
        A <- -gamma1*sigma_c12*P02_star*r1*sigma_f12
        B <- (gamma1*P00_star - mu1*beta1 - sigma_c12*P02_star + k1 + mu1) * 
            r1*sigma_f12
        C <- (P00_star + k1)*r1*sigma_f12
        Q_star <- (-B - sqrt(B^2 - 4*A*C))/(2*A)
        # If no equilibrium exists, return NA
        if (is.na(Q_star)) return(NA)

        # Final invasion growth rate calculation
        b <- (beta1*mu1 - gamma1*P00_star)/(gamma1*sigma_c12*P02_star)
        IGR1 <- r1*sigma_f12/(Q_star+b)*
            (P00_star + k1 + (beta1*mu1 - gamma1*P00_star)/gamma1 -
            (beta1*mu1 - gamma1*P00_star)/(gamma1*sigma_c12*P02_star)*(k1 + mu1))

        return(IGR1)
    })
}

# Helper function for calculating analytic IGR
analytic.IGR.helper <- function(sigma_c22, d2) {
    params <- base.params
    params$sigma_c22 <- sigma_c22
    params$d2 <- d2
    return(analytic.IGR(parms = params))
}

# Calculate analytic IGR for each parameter combination
analytic.IGR.df <- expand.grid(
    sigma_c22 = sigma_c22.values, d2 = d2.values)
analytic.IGR.df$IGR <- mapply(analytic.IGR.helper,
    sigma_c22 = analytic.IGR.df$sigma_c22,
    d2 = analytic.IGR.df$d2)


# 3. Calculate numeric IGR
# ========================
# Initial invader abundance for IGR calculation
initial.invader.abundance <- 10^-6
# Vector of output times for IGR calculation
time <- seq(0, 300, by = 0.1)

# Helper function for numeric IGR
numeric.IGR.helper <- function(sigma_c22, d2) {
    params <- base.params
    params$sigma_c22 <- sigma_c22
    params$d2 <- d2
    return(simulation.IGR(
        state.vector = state.vector,
        state.type = state.type,
        model = psf.model,
        params = params,
        initial.invader.abundance = initial.invader.abundance,
        time = time,
        initial.resident.vector = c(P00 = 1/3, P02 = 1/3, P22 = 1/3),
        initial.proportions = c(0.4, 0.5, 0.6)))
}

# Calculate numeric IGR for each parameter combination
numeric.IGR.df <- expand.grid(
    sigma_c22 = sigma_c22.focal, d2 = d2.values)
numeric.IGR.df$IGR <- mapply(numeric.IGR.helper,
    sigma_c22 = numeric.IGR.df$sigma_c22,
    d2 = numeric.IGR.df$d2)


# 4. Plot the results (Fig. 3d)
# =============================
# Labels for d2 values in the plot
d2.legend <- c("d2 = 0.1" = 0.1, "d2 = 0.9" = 0.9)

# Change d2 to final legend name
numeric.IGR.df <- numeric.IGR.df %>%
    mutate(d2 = factor(d2, levels = d2.legend, labels = names(d2.legend)))
analytic.IGR.df <- analytic.IGR.df %>%
    mutate(d2 = factor(d2, levels = d2.legend, labels = names(d2.legend)))

# Make and export figure
ggplot() +
    geom_line(data = analytic.IGR.df,
        aes(x = sigma_c22, y = IGR, color = d2),
        linewidth = case.study.linewidth) +
    geom_point(data = numeric.IGR.df,
        aes(x = sigma_c22, y = IGR, color = d2),
        size = case.study.point.size) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 1, linetype = "dotted") +
    scale_color_manual(
        values = c("d2 = 0.1" = light.green, "d2 = 0.9" = dark.green),
        labels = c("d2 = 0.1" = 0.1, "d2 = 0.9" = 0.9)) +
    labs(x = expression(sigma[c*","*22]),
        y = "IGR",
        color = expression(d[2])) +
    case.study.theme +
    scale_y_continuous(breaks = c(-0.05, 0.00, 0.05, 0.10, 0.15, 0.20)) +
    theme(legend.position.inside = c(0.6, 0.6))
ggsave(paste0(out.dir, "03_d_plant_soil_feedback.pdf"),
    width = case.study.width, height = case.study.height,
    units = case.study.units, dpi = case.study.dpi)
