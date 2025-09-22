## validate moments by comparing to iMRMC output from large simulated dataset

source("functions_moments.R")


library(purrr)
library(gt)
library(tidyr)

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
tab1_row1 <- list(delta=0.05, sigma_r=0.011, sigma_tr=0.011, sigma_C=0.3, 
                  sigma_TC = 0.3, sigma_RC=0.2, sigma_trc=0.2, AUC1 = 0.65, AUC2 = 0.70)

tab1_row2 <- list(delta=0.05, sigma_r=0.030, sigma_tr=0.030, sigma_C=0.3, 
                  sigma_TC = 0.3, sigma_RC=0.2, sigma_trc=0.2, AUC1 = 0.80, AUC2 = 0.85)


# Inputs:
N0 <- 80
N1 <- 60  
NR <- 16


# Build sigma sums
sigmas <- build_sigma_sums(tab1_row1)

# epected AUC is calculated like this: pnorm( Delta / sqrt(sigma_Omega + sigma_mod) )
# back into delta
delta1 <- qnorm(tab1_row1$AUC1) * sqrt(sigmas$sigma_Omega + sigmas$sigma_A)
delta2 <- qnorm(tab1_row1$AUC2) * sqrt(sigmas$sigma_Omega + sigmas$sigma_B)

# Expected AUCs (you also have AUC1/AUC2, but recompute here from deltas & sigmas for consistency)
 AUCA <- expected_auc(delta1, sigmas$sigma_Omega, sigmas$sigma_A)
 AUCB <- expected_auc(delta2, sigmas$sigma_Omega, sigmas$sigma_B)

moments_df <- compute_moments_df(delta1, delta2, sigmas, AUCA, AUCB)
print(moments_df)


### check with table 2 from Chen, Gong, Gallas 2018 paper

deltaA_true <- AUCB - AUCA

# Extract the rows as numeric vectors
MA  <- as.numeric(moments_df[moments_df$modality == "A", paste0("M", 1:8)])
MB  <- as.numeric(moments_df[moments_df$modality == "B", paste0("M", 1:8)])
MAB <- as.numeric(moments_df[moments_df$modality == "Cross", paste0("M", 1:8)])

v_fc  <- var_deltaA(MA, MB, MAB, N0, N1, NR, design="FC")
pow_fc <- power_two_sided(deltaA_true, v_fc, alpha = 0.05)
cat(sprintf("FC design: Var(ΔA)=%.6g, SD=%.6g, Power=%.3f\n", v_fc, sqrt(v_fc), pow_fc))






# ---- full Table 1 scenarios (same as before) ----
scenarios <- list(
  # HH
  list(structure="HH", AUC1=0.65, AUC2=0.70, sigma_r=0.011, sigma_tr=0.011,
       sigma_C=0.3, sigma_TC=0.3, sigma_RC=0.2, sigma_trc=0.2),
  list(structure="HH", AUC1=0.80, AUC2=0.85, sigma_r=0.030, sigma_tr=0.030,
       sigma_C=0.3, sigma_TC=0.3, sigma_RC=0.2, sigma_trc=0.2),
  list(structure="HH", AUC1=0.90, AUC2=0.95, sigma_r=0.056, sigma_tr=0.056,
       sigma_C=0.3, sigma_TC=0.3, sigma_RC=0.2, sigma_trc=0.2),
  # HL
  list(structure="HL", AUC1=0.65, AUC2=0.70, sigma_r=0.0055, sigma_tr=0.0055,
       sigma_C=0.3, sigma_TC=0.3, sigma_RC=0.2, sigma_trc=0.2),
  list(structure="HL", AUC1=0.80, AUC2=0.85, sigma_r=0.0055, sigma_tr=0.0055,
       sigma_C=0.3, sigma_TC=0.3, sigma_RC=0.2, sigma_trc=0.2),
  list(structure="HL", AUC1=0.90, AUC2=0.95, sigma_r=0.0055, sigma_tr=0.0055,
       sigma_C=0.3, sigma_TC=0.3, sigma_RC=0.2, sigma_trc=0.2),
  # LH
  list(structure="LH", AUC1=0.65, AUC2=0.70, sigma_r=0.011, sigma_tr=0.011,
       sigma_C=0.1, sigma_TC=0.1, sigma_RC=0.2, sigma_trc=0.6),
  list(structure="LH", AUC1=0.80, AUC2=0.85, sigma_r=0.030, sigma_tr=0.030,
       sigma_C=0.1, sigma_TC=0.1, sigma_RC=0.2, sigma_trc=0.6),
  list(structure="LH", AUC1=0.90, AUC2=0.95, sigma_r=0.056, sigma_tr=0.056,
       sigma_C=0.1, sigma_TC=0.1, sigma_RC=0.2, sigma_trc=0.6),
  # LL
  list(structure="LL", AUC1=0.65, AUC2=0.70, sigma_r=0.0055, sigma_tr=0.0055,
       sigma_C=0.1, sigma_TC=0.1, sigma_RC=0.2, sigma_trc=0.6),
  list(structure="LL", AUC1=0.80, AUC2=0.85, sigma_r=0.0055, sigma_tr=0.0055,
       sigma_C=0.1, sigma_TC=0.1, sigma_RC=0.2, sigma_trc=0.6),
  list(structure="LL", AUC1=0.90, AUC2=0.95, sigma_r=0.0055, sigma_tr=0.0055,
       sigma_C=0.1, sigma_TC=0.1, sigma_RC=0.2, sigma_trc=0.6)
)



### convert original RMH parameters to OR parameters to match how it will be used in practice
tab1_OR_params <- list()
for (i in seq_along(scenarios)) {
  scen <- scenarios[[i]]
  OR_params <- convert_to_OR(AUC1 = scen$AUC1, AUC2 = scen$AUC2, 
                             var_R = scen$sigma_r,
                             var_TR = scen$sigma_tr, var_C = scen$sigma_C,
                             var_TC = scen$sigma_TC, var_RC = scen$sigma_RC,
                             var_error = scen$sigma_trc, n0 = 80, n1 = 60)
  tab1_OR_params[[i]] <- OR_params
}


# scen1_OR_params<-convert_to_OR(AUC1 = scenarios[[1]]$AUC1, AUC2 = scenarios[[1]]$AUC2, 
#               var_R = scenarios[[1]]$sigma_r,
#               var_TR = scenarios[[1]]$sigma_tr, var_C = scenarios[[1]]$sigma_C,
#               var_TC = scenarios[[1]]$sigma_TC, var_RC = scenarios[[1]]$sigma_RC,
#               var_error = scenarios[[1]]$sigma_trc, n0 = 80, n1 = 60)

converted_RMH_params<- map_df(tab1_OR_params, OR_to_RMH)


# ---- run all scenarios ----
results <- map_df(scenarios, evaluate_power_all, N0=80, N1=60, NR=16)

converted_RMH_params2<-converted_RMH_params %>% select(n0, n1, delta1, delta2, 
                                                       var_R, var_TR, var_C,
                                                       var_TC, var_RC, var_error)


df_list <- split(converted_RMH_params2, seq(nrow(converted_RMH_params2)))
df_list <- lapply(df_list, as.list)

results2 <- map_df(df_list, evaluate_power_new, NR=16)


# Prepare data: separate variance and power columns for each design
results_wide <- results %>%
  mutate(var_x1e3 = 1000 * var_delta,
         power_pct = 100 * power) %>%
  select(structure, AUC1, AUC2, design, var_x1e3, power_pct) %>%
  pivot_wider(
    names_from = design,
    values_from = c(var_x1e3, power_pct)
  )

# Build gt table with top-level spanners
results_gt <- results_wide %>%
  gt(rowname_col = "structure") %>%
  tab_spanner(
    label = "Variance (×1000)",
    columns = starts_with("var_x1e3")
  ) %>%
  tab_spanner(
    label = "Power (%)",
    columns = starts_with("power_pct")
  ) %>%
  cols_label(
    var_x1e3_FC = "FC", var_x1e3_PSP2 = "PSP2",
    var_x1e3_PSP4 = "PSP4", var_x1e3_PSP8 = "PSP8",
    power_pct_FC = "FC", power_pct_PSP2 = "PSP2",
    power_pct_PSP4 = "PSP4", power_pct_PSP8 = "PSP8"
  ) %>%
  fmt_number(
    columns = starts_with("var_x1e3"),
    decimals = 2
  ) %>%
  fmt_number(
    columns = starts_with("power_pct"),
    decimals = 1
  )

results_gt

print(results_table)




