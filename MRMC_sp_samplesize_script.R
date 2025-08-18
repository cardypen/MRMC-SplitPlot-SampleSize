# .libPaths("/users/cpennington/R/rocker/4.3.0")
# print(.libPaths())
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
  AUC1 <- accuracy_level - delta
  AUC2 <- accuracy_level
  n0 <- 100  # control
  n1 <- 100  # case
  
  inter_var <- var_table[[observer_var]][2]
  intra_var <- var_table[[observer_var]][1]
  
  rmh <- OR_to_RMH(
    AUC1 = AUC1,
    AUC2 = AUC2,
    var_R = inter_var,
    var_TR = intra_var,
    corr1 = 0.47,  # average across Rockette studies
    corr2 = 0.35, 
    corr3 = 0.35,
    # corr2 = 0,
    # corr3 = 0,
    n0 = n0,
    n1 = n1,
    b_method = "specified",
    b_input=1
  )
  
  return(list(delta = delta, sigma_r = rmh$var_R, sigma_tr = rmh$var_TR, 
              sigma_C = rmh$var_C, sigma_TC = rmh$var_TC, sigma_RC = rmh$var_RC, 
              sigma_trc = rmh$var_error, delta1 = rmh$delta1, delta2 = rmh$delta2,
              AUC1 = AUC1, AUC2 = AUC2))
}

################
### Main Sim ###
################

#args <- as.integer(commandArgs(trailing = TRUE))
args <- 3
print(paste("Running scenario index:", args))

# Build full grid
scenarios <- build_OR_scenarios()
selected <- scenarios[args, ]

print(selected)

# Convert to RMH
rmh <- OR_scenario_to_RMH(
  readers = selected$readers,
  observer_var = selected$observer_var,
  accuracy_level = selected$accuracy_level,
  delta = selected$delta
)

# check parameter conversion
RMH_to_OR(
  n0 = 100, 
  n1 = 100, 
  b = 1,
  delta1 = rmh$delta1,
  delta2 = rmh$delta2,
  var_R = rmh$sigma_r,
  var_TR = rmh$sigma_tr,
  var_C = rmh$sigma_C,
  var_TC = rmh$sigma_TC,
  var_RC = rmh$sigma_RC,
  var_error = rmh$sigma_trc
)

# Estimate sample size via simulation
ss_result <- find_min_sample_size_uniroot(
  readers_per_block = selected$readers,
  blocks = 1,
  sigma_r = rmh$sigma_r,
  sigma_tr = rmh$sigma_tr,
  delta = rmh$delta, 
  rangeb = var_table[[selected$observer_var]][2], 
  rangew = var_table[[selected$observer_var]][1],
  theta = selected$accuracy_level,
  mu_nondisease = 0,  # not sure where to get these (but I think they should correspond to certain accuracy levels (theta))
  mu_disease = rmh$delta1, # not sure where to get these (but I think they should correspond to certain accuracy levels (theta))
  tau = rmh$delta2 - rmh$delta1,
  sigma_c = rmh$sigma_C, # var_C
  sigma_rc = rmh$sigma_RC, # var_RC
  sigma_tc = rmh$sigma_TC,  # var_TC
  sigma_trc = rmh$sigma_trc, # var_error
  n_sim = 50,
  target_power = 0.80
)

# ss_result_optim <- find_min_sample_size_optim(
#   readers_per_block = selected$readers,
#   blocks = 1,
#   sigma_r = rmh$sigma_r,
#   sigma_tr = rmh$sigma_tr,
#   delta = rmh$delta, 
#   rangeb = var_table[[selected$observer_var]][2], 
#   rangew = var_table[[selected$observer_var]][1],
#   theta = selected$accuracy_level,
#   mu_nondisease = 0,  # not sure where to get these (but I think they should correspond to certain accuracy levels (theta))
#   mu_disease = rmh$delta1, # not sure where to get these (but I think they should correspond to certain accuracy levels (theta))
#   tau = rmh$delta2 - rmh$delta1,
#   sigma_c = rmh$sigma_C, # var_C
#   sigma_rc = rmh$sigma_RC, # var_RC
#   sigma_tc = rmh$sigma_TC,  # var_TC
#   sigma_trc = rmh$sigma_trc, # var_error
#   n_sim = 50,
#   target_power = 0.80
# )

### feed RMH parameters back into RMH_to_OR to make sure they match the original OR parameters 
######### looks right
### check average AUCs for two treatments in the simulation (using iMRMC output) 
######## now prints average AUCs for each treatment on each run
######## AUCs look correct
### then if those match, check power calculation (or see if variance is too large)

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
print(paste("Scenario", args, "completed"))
