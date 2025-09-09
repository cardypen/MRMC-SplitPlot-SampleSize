# calculate moments

library(iMRMC)
library(dplyr)
library(MRMCaov)

# Specify OR parameters

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


# Use OR to RMH

####################################
### Convert OR -> RMH parameters ###
####################################

OR_scenario_to_RMH <- function(readers, observer_var, accuracy_level, delta) {
  AUC1 <- accuracy_level - (delta/2)
  AUC2 <- accuracy_level + (delta/2)  # make AUC1 and AUC2 symmetric around the accuracy level
  n0 <- 100  # control
  n1 <- 100  # case
  
  inter_var <- var_table[[observer_var]][2]
  intra_var <- var_table[[observer_var]][1]
  
  rmh <- OR_to_RMH(
    AUC1 = AUC1,
    AUC2 = AUC2,
    var_R = inter_var*4, #### multiply by 4 to convert from range to var
    var_TR = intra_var*4, #### multiply by 4 to convert from range to var
    corr1 = 0.47,  # average across Rockette studies
    corr2 = 0.3, 
    corr3 = 0.3,
    # corr2 = 0,
    # corr3 = 0,
    n0 = n0,
    n1 = n1,
    b_method = "specified",
    b_input=1
  )
  
  print(rmh)
  
  return(list(delta = delta, sigma_r = rmh$var_R, sigma_tr = rmh$var_TR, 
              sigma_C = rmh$var_C, sigma_TC = rmh$var_TC, sigma_RC = rmh$var_RC, 
              sigma_trc = rmh$var_error, delta1 = rmh$delta1, delta2 = rmh$delta2,
              AUC1 = AUC1, AUC2 = AUC2))
}


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


# Calculate moments

#############################
### Moments & Variances  ####
#############################

# Helper: normal PDF/CDF wrappers (base R pnorm/dnorm are fine)

# Given RMH params from OR_to_RMH (equal across truth & modalities),
# map them to the equal-variance special case of the generalized R&M model.
# In this equal-component case:
#   σ_R0=σ_R1=var_R;  σ_C0=σ_C1=var_C;  σ_RC0=σ_RC1=var_RC
#   σ_AR0=σ_AR1=var_TR; σ_AC0=σ_AC1=var_TC; σ_ARC0=σ_ARC1=var_error
# Then:
#   σΩ      = 2*(var_R + var_C + var_RC)
#   σA = σB = 2*(var_TR + var_TC + var_error)

build_sigma_sums <- function(rmh) {
  var_R     <- rmh$sigma_r
  var_TR    <- rmh$sigma_tr
  var_C     <- rmh$sigma_C
  var_TC    <- rmh$sigma_TC
  var_RC    <- rmh$sigma_RC
  var_error <- rmh$sigma_trc
  
  # assuming equal variance for both truth states and both modalities
  sigma_Omega <- 2*(var_R + var_C + var_RC)
  sigma_mod   <- 2*(var_TR + var_TC + var_error)  
  
  list(
    sigma_Omega = sigma_Omega,
    sigma_A     = sigma_mod,
    sigma_B     = sigma_mod,
    # For Table 3 rows (equal-components simplifications):
    sigmaA_l = list(
      `2` = (var_TC + var_error),
      `3` = (var_TC + var_error),
      `4` = 2*(var_TC + var_error),
      `5` = 2*(var_TR + var_error),
      `6` = 2*(var_TR + var_TC + var_error),
      `7` = 2*(var_TR + var_TC + var_error)
    ),
    sigmaOmega_l = list(
      # For covariance rows in Table 3; these also appear inside single-modality integrals (Eq. 12)
      `2` = (var_C + var_RC),
      `3` = (var_C + var_RC),
      `4` = 2*(var_C + var_RC),
      `5` = 2*(var_R + var_RC),
      `6` = 2*(var_R + var_C + var_RC),
      `7` = 2*(var_R + var_C + var_RC)
    )
  )
}

# Coefficients c for fully-crossed design (Table 2)
coeff_vector <- function(N0, N1, NR) {
  c1 <- 1/(N0*N1*NR)
  c2 <- (N0-1)/(N0*N1*NR)
  c3 <- (N1-1)/(N0*N1*NR)
  c4 <- (N0-1)*(N1-1)/(N0*N1*NR)
  c5 <- (NR-1)/(N0*N1*NR)
  c6 <- (N0-1)*(NR-1)/(N0*N1*NR)
  c7 <- (N1-1)*(NR-1)/(N0*N1*NR)
  c8 <- (N0-1)*(N1-1)*(NR-1)/(N0*N1*NR) - 1  # note the "− 1" for the eighth entry
  c(c1,c2,c3,c4,c5,c6,c7,c8)
}

# Expected AUC from Δ and totals (Eq. 9)
expected_auc <- function(Delta, sigma_Omega, sigma_mod) {
  pnorm( Delta / sqrt(sigma_Omega + sigma_mod) )
}

# Single-modality moment M_l for l=2..7 (Eq. 12)
moment_single_l <- function(l, Delta, sigmas) {
  sOm   <- sigmas$sigma_Omega
  sMod  <- sigmas$sigma_A   # same form if called for B in this equal-variance case
  sOm_l <- sigmas$sigmaOmega_l[[as.character(l)]]
  sMod_l<- sigmas$sigmaA_l[[as.character(l)]]
  
  s_common <- sOm + sMod - sOm_l - sMod_l
  denom   <- sqrt(sOm_l + sMod_l)
  
  integrand <- function(x) {
    p <- pnorm( (Delta + x*sqrt(s_common)) / denom )
    p*p * dnorm(x)
  }
  integrate(integrand, lower = -Inf, upper = Inf, subdivisions = 400L, rel.tol = 1e-7)$value
}

# Cross-modality moment MAB_l for l=1..7 (Eq. 15)
moment_cross_l <- function(l, DeltaA, DeltaB, sigmas) {
  sOm   <- sigmas$sigma_Omega
  sOm_l <- sigmas$sigmaOmega_l[[as.character(l)]]
  # In equal-variance case σA=σB
  sA <- sigmas$sigma_A
  sB <- sigmas$sigma_B
  
  s_common <- sOm - sOm_l
  denomA   <- sqrt(sA + sOm_l)
  denomB   <- sqrt(sB + sOm_l)
  
  integrand <- function(x) {
    pA <- pnorm( (DeltaA + x*sqrt(s_common)) / denomA )
    pB <- pnorm( (DeltaB + x*sqrt(s_common)) / denomB )
    pA * pB * dnorm(x)
  }
  integrate(integrand, lower = -Inf, upper = Inf, subdivisions = 400L, rel.tol = 1e-7)$value
}

# Build all 8 moments for a single modality (A or B)
moments_single <- function(Delta, sigmas, AUC) {
  # M1 = AUC ; M8 = AUC^2 ; M2..M7 via integral
  M <- numeric(8)
  M[1] <- AUC
  for (l in 2:7) M[l] <- moment_single_l(l, Delta, sigmas)
  M[8] <- AUC*AUC
  M
}

# Build all 8 cross-modality moments MAB
moments_cross <- function(DeltaA, DeltaB, sigmas, AUCA, AUCB) {
  M <- numeric(8)
  # l = 1..7 via integral with σΩ(l) (Table 3)
  for (l in 1:7) M[l] <- moment_cross_l(l, DeltaA, DeltaB, sigmas)
  # Special case: l=8 => AUCA * AUCB
  M[8] <- AUCA * AUCB
  M
}

# ---- Driver given your 'rmh' and design sizes ----
# Inputs you already have/know:
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

# Moments
MA   <- moments_single(DeltaA, sigmas, AUCA)
MB   <- moments_single(DeltaB, sigmas, AUCB)
MAB  <- moments_cross(DeltaA, DeltaB, sigmas, AUCA, AUCB)

# Coefficients
cvec <- coeff_vector(N0, N1, NR)

# Variances & covariance (Eqs. 11, 13, 10)
var_A   <- sum(cvec * MA)
var_B   <- sum(cvec * MB)
cov_AB  <- sum(cvec * MAB)
var_diff <- var_A + var_B - 2*cov_AB

cat(sprintf("\nAUC_A = %.6f,  sd(AUC_A) = %.6f\n", AUCA, sqrt(var_A)))
cat(sprintf("AUC_B = %.6f,  sd(AUC_B) = %.6f\n", AUCB, sqrt(var_B)))
cat(sprintf("cov(AUC_A, AUC_B) = %.6f\n", cov_AB))
cat(sprintf("sd(AUC_A - AUC_B) = %.6f\n\n", sqrt(var_diff)))


















































