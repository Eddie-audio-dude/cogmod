graphics.off()
rm(list = ls(all = TRUE))

library(tidyverse)
library(patchwork)
library(WienR)

set.seed(12222)

# This could be derived using MDS, but here I'm just specifying them myself (and in 2 dimensions)
stim_coords <- matrix(
    c(
        -1, -1,
        -1,  1,
         1, -1,
         1,  1
    ),
    nrow = 4, ncol = 2, byrow = TRUE
)

# The following creates every possible combination of 2-item study lists and probe items
study_items <- c()
probe_item <- c()

for (i in 1:(nrow(stim_coords) - 1)) {
    for (j in (i + 1):nrow(stim_coords)) {
        for (k in 1:nrow(stim_coords)) {
            study_items <- rbind(
                study_items,
                c(i, j)
            )
            
            probe_item <- c(
                probe_item,
                k
            )
        }
    }
}

# Function arguments:
# par: this is a named vector of parameter values
# stim_coords: this is a matrix of the coordinates of the stimuli, where each row is a stimulus and each column is a dimension
# study_items: this is a matrix where each row is a trial and each column indicates the items that were studied on that trial
# probe_item: this is a vector giving the index of the probe item on each trial
# response: this is a vector where each value is 2 or 1, depending on whether the participant responsed "yes" (2) or "no" (1) on that trial
# rt: this is a vector of the response times from each trial
# n_sims: this indicates whether or not simulations are meant to be returned instead of the NLL
ebrw_nll <- function(par, stim_coords, study_items, probe_item, response, rt, n_sims = 0) {
    # 1. Compute the mean and SD of the drift for each trial
    
    stim_dists <- as.matrix(dist(stim_coords))
    
    evidence_mean <- rep(0, length(probe_item))
    evidence_sd <- rep(0, length(probe_item))
    
    for (i in 1:length(probe_item)) {
        summed_sim <- sum(exp(-par["specificity"] * stim_dists[probe_item[i], study_items[i,]]))
        p <- summed_sim / (summed_sim + par["criterion"])
        
        evidence_mean[i] <- par["retrieval_rate"] * (2 * p - 1)
        evidence_sd[i] <- 2 * sqrt(par["retrieval_rate"] * p * (1 - p))
    }
    
    if (is.infinite(n_sims)) {
        # 2a. Compute predicted response probabilities and mean RT's for each trial
        result <- c()
        
        evidence_var <- evidence_sd^2
        
        for (i in 1:length(probe_item)) {
            trial_p_yes <- (1 - exp(-2 * par["a"] * evidence_mean[i] * par["w"] / evidence_var[i])) / (1 - exp(-2 * par["a"] * evidence_mean[i] / evidence_var[i]))
            trial_mean_rt_yes <- par["t0"] + (par["a"] / evidence_mean[i]) * (1 / tanh(par["a"] * evidence_mean[i] / evidence_var[i]) - par["w"] / tanh(par["a"] * evidence_mean[i] * par["w"] / evidence_var[i]))
            trial_mean_rt_no <- par["t0"] + (par["a"] / evidence_mean[i]) * (1 / tanh(par["a"] * evidence_mean[i] / evidence_var[i]) - (1 - par["w"]) / tanh(par["a"] * evidence_mean[i] * (1 - par["w"]) / evidence_var[i]))
            
            result <- rbind(
                result,
                tibble(trial = i, p_yes = trial_p_yes, mean_rt_yes = trial_mean_rt_yes, mean_rt_no = trial_mean_rt_no)
            )
        }
        
        # 3a. Return final result
        return(result)
    } else if (n_sims > 0) {
        # 2b. Simulate a response for each trial
        result <- c()
        
        # For each trial...
        for (i in 1:length(probe_item)) {
            # Simulate `n_sims` replications of this trial, using its corresponding parameters
            trial_sim <- sampWiener(
                N = n_sims,
                a = par["a"] / evidence_sd[i],
                w = par["w"],
                v = evidence_mean[i] / evidence_sd[i],
                t0 = par["t0"]
            )
            
            # Append the new set of simulations to our `result`, including an indicator of which trial we are simulating and an index for each simulation
            result <- rbind(
                result,
                tibble(trial = i, sim_index = 1:n_sims, rt = trial_sim$q, response = trial_sim$response)
            )
        }
        
        # 3b. Return final result
        return(result)
    } else {
        # 2c. Calculate the log-likelihood of each observed response/RT on each trial
        result <- try(WienerPDF(
            t = rt,
            response = response,
            a = par["a"] / evidence_sd,
            w = par["w"],
            v = evidence_mean / evidence_sd,
            t0 = par["t0"]
        ))
        
        # 3c. Return final result
        if (class(result) == "try-error") {
            return(Inf)
        } else {
            return(-sum(result$logvalue))
        }
    }
}

# Parameters for simulation
sim_par <- c(
    "retrieval_rate" = 3,     # This is the "nu" parameter
    "a" = 2,                  # Response caution
    "w" = 0.5,                # Response bias
    "t0" = 0.2,               # Residual time
    "specificity" = 1,        # Specificity of memory representations (the "c" parameter)
    "criterion" = 1           # Criterion (the "kappa" parameter)
)

# Simulate trials
sim_trials <- ebrw_nll(
    par = sim_par,
    stim_coords = stim_coords,
    study_items = study_items,
    probe_item = probe_item,
    n_sims = 1
)

# Fit the model
init_par <- c(
    "retrieval_rate" = 3,     # This is the "nu" parameter
    "a" = 2,                  # Response caution
    "w" = 0.5,                # Response bias
    "t0" = 0,                 # Residual time
    "specificity" = 1,        # Specificity of memory representations (the "c" parameter)
    "criterion" = 1           # Criterion (the "kappa" parameter)
)

lower <- c(
    "retrieval_rate" = 0,     # This is the "nu" parameter
    "a" = 0,                  # Response caution
    "w" = 0,                  # Response bias
    "t0" = 0,                 # Residual time
    "specificity" = 0,        # Specificity of memory representations (the "c" parameter)
    "criterion" = 0           # Criterion (the "kappa" parameter)
)

upper <- c(
    "retrieval_rate" = Inf,     # This is the "nu" parameter
    "a" = Inf,                  # Response caution
    "w" = 1,                    # Response bias
    "t0" = min(sim_trials$rt),  # Residual time
    "specificity" = Inf,        # Specificity of memory representations (the "c" parameter)
    "criterion" = Inf           # Criterion (the "kappa" parameter)
)

fit <- nlminb(
    start = init_par,               # Need to provide initial guess of parameter values
    objective = ebrw_nll,           # Tell R the name of the function to optimize
    lower = lower,                  # The lower bounds on each parameter
    upper = upper,                  # The upper bounds on each parameter
    stim_coords = stim_coords,      # The coordinates of each stimulus
    rt = sim_trials$rt,             # The vector of RT's on each trial
    response = sim_trials$response, # The vector of responses on each trial
    study_items = study_items,      # The study items on each trial
    probe_item = probe_item         # The probe item on each trial
)

# Compare fitted to observed
math_pred <- ebrw_nll(
    par = fit$par,
    stim_coords = stim_coords,
    study_items = study_items,
    probe_item = probe_item,
    n_sims = Inf
)

obs_pred <- full_join(sim_trials, math_pred, by = "trial")

obs_pred$probe_type <- NA

for (i in 1:length(probe_item)) {
    # Note the use of "%in%"
    if (probe_item[i] %in% study_items[i,]) {
        obs_pred$probe_type[i] <- "Target"
    } else {
        obs_pred$probe_type[i] <- "Foil"
    }
}

obs_pred$probe_type <- factor(obs_pred$probe_type, levels = c("Target", "Foil"))

resp_plot <- obs_pred %>%
    ggplot(aes(x = 1 * (response == "upper"), y = p_yes, color = factor(response))) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    geom_point() +
    coord_equal() +
    scale_color_manual(values = c("upper" = "#377eb8", "lower" = "#e41a1c")) +
    labs(x = "Observed response", y = "Predicted P(\"Yes\")", color = "Response")

rt_plot <- obs_pred %>%
    mutate(pred_rt = if_else(response == 1, mean_rt_yes, mean_rt_no)) %>%
    ggplot(aes(x = rt, y = pred_rt, color = factor(response))) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    geom_point() +
    coord_equal() +
    scale_color_manual(values = c("upper" = "#377eb8", "lower" = "#e41a1c")) +
    labs(x = "Observed RT", y = "Predicted mean RT", color = "Response")

print(
    resp_plot + rt_plot + plot_layout(nrow = 1, guides = "collect") + plot_annotation(title = "Individual trials")
)

agg_resp_plot <- obs_pred %>%
    group_by(probe_type) %>%
    summarize(mean_obs = mean(response == "upper"), mean_pred = mean(p_yes)) %>%
    ggplot(aes(x = mean_obs, y = mean_pred, color = probe_type)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    geom_point() +
    coord_equal() +
    scale_color_manual(values = c("Target" = "#377eb8", "Foil" = "#e41a1c")) +
    labs(x = "Observed response proportion", y = "Predicted P(\"Yes\")", color = NULL)

agg_rt_plot <- obs_pred %>%
    mutate(pred_rt = if_else(response == 1, mean_rt_yes, mean_rt_no)) %>%
    group_by(probe_type) %>%
    summarize(mean_obs = mean(rt), mean_pred = mean(pred_rt)) %>%
    ggplot(aes(x = mean_obs, y = mean_pred, color = probe_type)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    geom_point() +
    coord_equal() +
    scale_color_manual(values = c("Target" = "#377eb8", "Foil" = "#e41a1c")) +
    labs(x = "Observed mean RT", y = "Predicted mean RT", color = NULL)

print(
    agg_resp_plot + agg_rt_plot + plot_layout(nrow = 1, guides = "collect") + plot_annotation(title = "Aggregated over target/foil trials")
)
