## validate moments by comparing to iMRMC output from large simulated dataset

source("functions_moments.R")

# Parameters from scenario 1
target_power <- 0.80
args <- 1
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


# Inputs:
N0 <- 100  # you used n0=100 in OR_to_RMH
N1 <- 100  # you used n1=100 in OR_to_RMH
NR <- selected$readers

# Deltas from OR_to_RMH (these are the Δ_A and Δ_B in Eq. 8/9)
DeltaA <- rmh$delta1
DeltaB <- rmh$delta2

# Build sigma sums
sigmas <- build_sigma_sums(rmh)

# Expected AUCs (you also have AUC1/AUC2, but recompute here from deltas & sigmas for consistency)
AUCA <- expected_auc(DeltaA, sigmas$sigma_Omega, sigmas$sigma_A)
AUCB <- expected_auc(DeltaB, sigmas$sigma_Omega, sigmas$sigma_B)

moments_df <- compute_moments_df(DeltaA, DeltaB, sigmas, AUCA, AUCB)
print(moments_df)


####### compute moments using iMRMC on large simulated dataset ########

# simulate data and put it in iMRMC format

scen1_data <- sim_one_splitplot_cardy(mu_nondisease = 0, 
                                             mu_disease = rmh$delta1, 
                                             hypothesis = "alt",
                                             tau = rmh$delta2-rmh$delta1,
                                             sigma_r = rmh$sigma_r, sigma_c = rmh$sigma_C,
                                             sigma_rc = rmh$sigma_RC, sigma_tr = rmh$sigma_tr,
                                             sigma_tc = rmh$sigma_TC, sigma_trc = rmh$sigma_trc, 
                                             block_ss = 1000,
                                             readers_per_block = 20,
                                             blocks = 1)

scen1_data_for_imrmc<-createIMRMCdf(scen1_data,
                                       keyColumns = list(readerID = "reader", 
                                                         caseID = "case", 
                                                         modalityID = "modality",
                                                         score = "score", 
                                                         truth = "truth_binary"),
                                       truePositiveFactor = 1)

scen1_imrmc_results<-doIMRMC(scen1_data_for_imrmc)

scen1_imrmc_results$VarDecomp$BDG
print(moments_df)















