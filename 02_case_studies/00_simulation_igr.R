#!/usr/bin/env Rscript
# 02_case_studies/00_simulation_igr.R
# ===================================
# Author: Hsiang-Chih Lo
# Implements the direct-simulation IGR function (equivalent to Algorithm S3.1
# in Appendix S3).

# Given a model and invasion analysis, use the direct simulation approach to
# identify convergence of population structure and return calculated IGR.
#
# Arguments:
#   state.vector: Template state variable vector in the format used by ODE.
#     Specific values do not matter, as they will be discarded, but the vector
#     should have the correct names if they are required by the model function.
#   state.type: Character vector: for each stage in state.vector, specify 
#     "inv" for invader stages and "res" for resident.
#   model: ODE model function in deSolve format.
#   params: Model parameters passed to the ODE model.
#   initial.invader.abundance: Initial total invader abundance.
#   time: Output time vector for ODE solver.
#   initial.resident.vector: Initial state for simulating resident equilibrium;
#     provide only the resident stages in the same order as state.vector.
#   initial.proportions: Vector of second-stage proportions to use in different
#     trajectories.
#   maximum.abundance: Maximum total invader abundance per trajectory when
#     identifying convergence of trajectories (default: 10^-3).
#   resident.time: Time to simulate resident equilibrium (default: 100000).

simulation.IGR <- function(
        state.vector, state.type, model, params, initial.invader.abundance, time, 
        initial.resident.vector, initial.proportions, maximum.abundance = 10^-3,
        resident.time = 100000) {
    # 1. Calculate resident's equilibrium state
    # -----------------------------------------
    # Fill in the initial resident states, setting invader states to zero
    initial.state <- state.vector
    initial.state[state.type=="res"] <- initial.resident.vector
    initial.state[state.type=="inv"] <- 0
    # Simulate the resident equilibrium
    resident.equilibrium <- ode(y = initial.state,
            times = seq(0, resident.time, by = 1),
            func = model,
            parms = params) %>%
        # Get final state vector
        tail(1) %>%
        as.data.frame() %>%
        select(-time) %>%
        as.numeric()
    
    # 2. Run parallel invasions with different starting conditions
    # ------------------------------------------------------------
    # Helper function which makes a starting state
    get.starting.state <- function(p2) {
        x <- state.vector * 0
        x[state.type=="inv"][1:2] <- initial.invader.abundance * c(1 - p2, p2)
        x[state.type=="res"] <- resident.equilibrium[state.type=="res"]
        return(x)
    }
    # Simulate all trajectories
    all.solutions <- data.frame(initial.proportions = initial.proportions) %>%
        mutate(initial.state = map(initial.proportions, get.starting.state)) %>%
        mutate(sol = map(initial.state, ~as.data.frame(ode(
            y = .x, times = time, func = model, parms = params))))

    # 3. Calculate per-capita growth rate of the invader population
    # -------------------------------------------------------------
    # Helper function to calculate per-capita growth of the total population
    total.growth <- function(time, state) {
        dN <- model(time, state, params)[[1]]
        total.derivative <- sum(dN[state.type=="inv"])
        total.population <- sum(state[state.type=="inv"])
        return(total.derivative / total.population)
    }
    # Calculate per-capita growth rate at each time step
    growth.rate.df <- all.solutions %>%
        unnest(sol) %>%
        mutate(state = pmap(select(., -initial.proportions, -initial.state, -time), ~ c(...))) %>%
        mutate(gr1 = map2_dbl(time, state, total.growth))

    # 4. Identify the time of convergence of the trajectories
    # -------------------------------------------------------
    # Find the standard deviation between trajectories
    deviation.df <- growth.rate.df %>%
        group_by(time) %>%
        # # Require all invader stages in every trajectory to be below the cutoff
        # filter(all(if_all(
        #     all_of(names(state.vector)[state.type=="inv"]),
        #     ~ .x < maximum.abundance
        # ))) %>%
        # Require total invader abundance to be below the cutoff
        filter(all(rowSums(across(
            all_of(names(state.vector)[state.type=="inv"])
        )) < maximum.abundance)) %>%
        # Calculate standard deviation
        summarise(SD = sd(gr1), .groups = "drop")
    # Find time of minimum
    minimum.time <- deviation.df$time[which.min(deviation.df$SD)]
    # Extract per capita growth rate at this time
    IGR <- growth.rate.df %>%
        select(initial.proportions, time, gr1) %>% 
        group_by(time) %>% 
        summarise(growth_rate = mean(gr1)) %>%
        filter(time == minimum.time) %>%
        pull(growth_rate)
    
    # Return the IGR
    return(IGR)
}
