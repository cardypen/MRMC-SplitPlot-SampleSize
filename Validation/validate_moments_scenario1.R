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
                                             block_ss = 2000,
                                             readers_per_block = 30,
                                             blocks = 1)

scen1_data_for_imrmc<-createIMRMCdf(scen1_data,
                                       keyColumns = list(readerID = "reader", 
                                                         caseID = "case", 
                                                         modalityID = "modality",
                                                         score = "score", 
                                                         truth = "truth_binary"),
                                       truePositiveFactor = 1)

scen1_imrmc_results<-doIMRMC(scen1_data_for_imrmc)
# should save this analysis as it's a bit slow

scen1_imrmc_results$varDecomp$BDG$Ustat$comp
print(moments_df)



### check with table 2 from Chen, Gong, Gallas 2018 paper

deltaA_true <- AUCB - AUCA

v_fc  <- var_deltaA(MA, MB, MAB, N0, N1, NR, design="FC")
pow_fc <- power_two_sided(deltaA_true, v_fc, alpha = 0.05)
cat(sprintf("FC design: Var(ΔA)=%.6g, SD=%.6g, Power=%.3f\n", v_fc, sqrt(v_fc), pow_fc))


# ΔA is the true effect you’re targeting power for (e.g., AUC_B - AUC_A)
# Chen, Gong, Gallas scenarios from table 1
tab1_row1 <- list(delta=0.05, sigma_r=0.011, sigma_tr=0.11, sigma_C=0.03, 
                  sigma_TC = 0.3, sigma_RC=0.2, sigma_trc=0.2, AUC1 = 0.65, AUC2 = 0.70)

# Inputs:
N0 <- 80  # you used n0=100 in OR_to_RMH
N1 <- 60  # you used n1=100 in OR_to_RMH
NR <- 16

# Deltas from OR_to_RMH (these are the Δ_A and Δ_B in Eq. 8/9)
# DeltaA <- rmh$delta1
# DeltaB <- rmh$delta2
# need to figure out how to get delta1 and delta2 when we have AUC1 and AUC2 (OR params) and RMH variances

# Build sigma sums
sigmas <- build_sigma_sums(tab1_row1)

# Expected AUCs (you also have AUC1/AUC2, but recompute here from deltas & sigmas for consistency)
# AUCA <- expected_auc(DeltaA, sigmas$sigma_Omega, sigmas$sigma_A)
# AUCB <- expected_auc(DeltaB, sigmas$sigma_Omega, sigmas$sigma_B)

moments_df <- compute_moments_df(DeltaA, DeltaB, sigmas, AUCA, AUCB)
print(moments_df)













