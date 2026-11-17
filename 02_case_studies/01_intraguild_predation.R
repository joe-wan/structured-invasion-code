#!/usr/bin/env Rscript
# 02_case_studies/01_intraguild_predation.R
# ===================================
# Author: Hsiang-Chih Lo
# Implements the intraguild predation/cannibalism simulation (Appendix S4) and
# makes Figure 3b.

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
intraguild.predation.model <- function(times, state, parms) {
    with(as.list(c(state, parms)), {
        # Intraguild predation term
        mu2_igp <- mu2*(aF1F2*F1 + aM1F2*M1 + u2*R) / (iF1*aF1F2*F1 + iM1*aM1F2*M1 + u2*R)
        # Calculate derivatives
        dF1_dt <- f*AF - u1*F1*R - aF1F2*F1*F2 - mu1*F1
        dF2_dt <- u1*F1*R - u2*F2*R - mu2_igp*F2
        dAF_dt <- u2*F2*R - muA*AF
        dM1_dt <- f*AM - u1*M1*R - aM1F2*M1*F2 - mu1*M1
        dM2_dt <- u1*M1*R - u2*M2*R - mu2*M2
        dAM_dt <- u2*R*M2 - muA*AM
        dR_dt  <- rR*R*(1 - R/K) - u1*F1*R - u2*F2*R - u1*M1*R - u2*M2*R
        # Return list of derivatives
        return(list(c(dF1_dt, dF2_dt, dAF_dt, dM1_dt, dM2_dt, dAM_dt, dR_dt)))
    })
}

# Template state variable vector
state.vector <- c(F1 = 1, F2 = 1, AF = 1, M1 = 1, M2 = 1, AM = 1, R = 1)
# States classified as invader or resident
state.type <- c("inv", "inv", "inv", "res", "res", "res", "res")

# Default model parameters
base.params <- list(
    K = 5, rR = 0.5, f = 0.5,
    u1 = 0.0001, u2 = 0.015,
    mu1 = 0.0005, mu2 = 0.00025, muA = 0.1,
    aF1F2 = 0.001, aM1F2 = 0.001,
    iF1 = 0.9, iM1 = 0.5
)

# Values of focal parameter (b21) for continuous analytic and pointwise numeric IGR
b21.values <- round(c(50:150)*0.01, 2)
b21.focal <- round(c(5:14)/10+0.05, 2)
# Values for variant analyses
a21N1.values <- c(0.00001, 0.0001)
cannibalism.values <- c(0, 1)

# 2. Calculate analytic IGR
# =========================
# Set up parameter combinations for analysis
IGR.df <- crossing(!!!base.params,
        b21 = b21.values) %>%
    mutate(iM1 = b21) %>%
    select(-aM1F2, -aF1F2) %>%
    crossing(a21N1 = a21N1.values) %>%
    mutate(aM1F2 = a21N1) %>%
    crossing(cannibalism = cannibalism.values) %>%
    mutate(aF1F2 = cannibalism*aM1F2) %>%
    mutate(params = pmap(., ~c(...)))

# Implement formulae (see Appendix S4) under timescale separation
IGR.df <- IGR.df %>%
    mutate(a = u1*u2* (f / muA -1),
           b = -(u1*mu2 + u2*mu1),
           c = -(mu1 * mu2),
           R = (-b + sqrt(b^2 - 4*a*c)) / (2*a)) %>%
    mutate(N22 = rR*(1 - R/K) / ((u1*u2*f*R) / (muA * (u1*R + mu1)) + u2)) %>%
    mutate(N21 = u2*f*R*N22 / (muA*(u1*R + mu1))) %>%
    mutate(N2A = u2*R*N22 / muA) %>%
    mutate(m12 = mu2 * (a21N1 * N21 + u2 * R) / (b21 * a21N1 * N21 + u2 * R)) %>%
    mutate(A = u1 * R,
           B = (u1 * R + mu1) - (u2 * R + m12),
           C = f * u2 * R / muA) %>%
    mutate(Q = (-B + sqrt(B^2 + 4*A*C)) / (2*A)) %>%
    mutate(analytic.IGR = u1 * R * Q - (u2 * R + m12)) %>%
    select(aM1F2, aF1F2, cannibalism, b21, analytic.IGR, params)


# 3. Numerical IGR calculation using the helper function
# ======================================================
# Initial invader abundance for IGR calculation
initial.invader.abundance <- 10^-8
# Vector of output times for IGR calculation
time <- seq(0, 10000000, by = 1000)

# Map helper onto the parameter combination data frame
IGR.df <- IGR.df %>% 
    mutate(numeric.IGR = map2_dbl(params, b21,
        ~ if(.y %in% b21.focal){
            simulation.IGR(
                state.vector = state.vector,
                state.type = state.type,
                model = intraguild.predation.model,
                params = .x,
                initial.invader.abundance = initial.invader.abundance,
                time = time,
                initial.resident.vector = c(M1 = 1, M2 = 1, AM = 1, R = 1),
                initial.proportions = c(0, 0.5, 1))
        } else {
            NA_real_
        })) %>%
  mutate(error = paste((numeric.IGR - analytic.IGR)/analytic.IGR*100, "%"))


# 4. Plot the results (Fig. 3b)
# =============================
IGR.df %>%
    mutate(type = as.factor(-log10(aM1F2) - 2*cannibalism - 1)) %>%
    ggplot(aes(color = as.factor(type))) +
    geom_line(aes(b21, analytic.IGR*10^6, linetype = as.factor(cannibalism)),
        linewidth = case.study.linewidth) +
    scale_linetype_manual(values = c("solid", "dashed"),
        labels = c("IGP only", "IGP + cannibalism")) +
    scale_color_manual(
        values = c("1" = dark.green.faded, "2" = light.green.faded,
                   "3" = dark.green, "4" = light.green),
        breaks = c(3, 4), labels = c("High", "Low")) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 1, linetype = "dotted") +
    geom_point(aes(b21, numeric.IGR*10^6,
        size = as.factor(cannibalism), shape = as.factor(cannibalism)),
        data = ~ filter(.x, b21 %in% b21.focal)) +
    scale_size_manual(values = c(4, 3),
        labels = c("IGP only", "IGP + cannibalism")) +
    scale_shape_manual(values = c(17, 19),
        labels = c("IGP only", "IGP + cannibalism")) +
    labs(title = "",
        color = "Aggression level",
        linetype = "Aggression type",
        size = "Aggression type",
        shape = "Aggression type",
        x = expression(b[21]),
        y = expression(IGR/10^{-6})) +
    case.study.theme +
    theme(legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.position.inside = c(0.55, 0.3))
# Export PDF
ggsave(paste0(out.dir, "03_b_intraguild_predation.pdf"),
    width = case.study.width, height = case.study.height,
    units = case.study.units, dpi = case.study.dpi)
