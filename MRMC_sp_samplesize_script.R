.libPaths("/users/cpennington/R/rocker/4.3.0")
print(.libPaths())
source("MRMC_sp_samplesize_functions_clean.R")

library(dplyr)
library(MRMCaov)

############################
### Observer Variability ###
############################

# From Table 1 in Obuchowski 2000
var_table <- list(
  S = c(0.005, 0.01),
  M = c(0.025, 0.05),
  L = c(0.05, 0.10)
)

#########################
### Define Scenarios ###
#########################

# All combinations with ratio = 1:1 only (Table 2 in Obuchowski 2000)
build_OR_scenarios <- function() {
  expand.grid(
    readers = c(4, 6, 10),
    observer_var = c("S", "M", "L"),
    accuracy_level = c(0.75, 0.90),
    delta = c(0.05, 0.10, 0.15),
    stringsAsFactors = FALSE
  ) %>%
    filter(!(readers == 4 & observer_var == "L"))  # Infeasible per paper
}

####################################
### Convert OR -> RMH parameters ###
####################################

OR_scenario_to_RMH <- function(readers, observer_var, accuracy_level, delta) {
  AUC1 <- accuracy_level
  AUC2 <- accuracy_level - delta
  n0 <- 1  # control
  n1 <- 1  # case
  
  inter_var <- var_table[[observer_var]][2]
  intra_var <- var_table[[observer_var]][1]
  
  rmh <- OR_to_RMH(
    AUC1 = AUC1,
    AUC2 = AUC2,
    var_R = inter_var,
    var_TR = intra_var,
    corr1 = 0.47,  # average across Rockette studies
    corr2 = 0,
    corr3 = 0,
    n0 = n0,
    n1 = n1
  )
  
  return(list(delta = delta, sigma_r = rmh$sigma_r, sigma_tr = rmh$sigma_tr))
}

################
### Main Sim ###
################

args <- as.integer(commandArgs(trailing = TRUE))
print(paste("Running scenario index:", args))

# Build full grid
scenarios <- build_OR_scenarios()
selected <- scenarios[args, ]

# Convert to RMH
rmh <- OR_scenario_to_RMH(
  readers = selected$readers,
  observer_var = selected$observer_var,
  accuracy_level = selected$accuracy_level,
  delta = selected$delta
)

# Estimate sample size via simulation
ss_result <- find_min_sample_size_uniroot(
  readers_per_block = selected$readers,
  blocks = 1,
  sigma_r = rmh$sigma_r,
  sigma_tr = rmh$sigma_tr,
  n_sim = 200,
  target_power = 0.80
)

# Save results
sim_results <- data.frame(
  readers = selected$readers,
  observer_var = selected$observer_var,
  accuracy_level = selected$accuracy_level,
  delta = selected$delta,
  sigma_r = rmh$sigma_r,
  sigma_tr = rmh$sigma_tr,
  simulated_ss = ss_result[[1]]
)

write.csv(sim_results, file = paste0("scenario", args, "_sp_samplesize_results.csv"), row.names = FALSE)
cat("Scenario", args, "completed.\n")
