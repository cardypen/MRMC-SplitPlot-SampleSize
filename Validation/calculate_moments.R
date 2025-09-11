# calculate moments
source("functions_moments.R")


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

# simgle modality moments
MA   <- moments_single(DeltaA, sigmas, AUCA)
MB   <- moments_single(DeltaB, sigmas, AUCB)

# cross moments
MAB  <- moments_cross(DeltaA, DeltaB, sigmas, AUCA, AUCB)

# coefficients
cvec <- coeff_vector(N0, N1, NR)

# variances & covariance (Eqs. 11, 13, 10)
var_A   <- sum(cvec * MA)
var_B   <- sum(cvec * MB)
cov_AB  <- sum(cvec * MAB)
var_diff <- var_A + var_B - 2*cov_AB

cat(sprintf("\nAUC_A = %.6f,  sd(AUC_A) = %.6f\n", AUCA, sqrt(var_A)))
cat(sprintf("AUC_B = %.6f,  sd(AUC_B) = %.6f\n", AUCB, sqrt(var_B)))
cat(sprintf("cov(AUC_A, AUC_B) = %.6f\n", cov_AB))
cat(sprintf("sd(AUC_A - AUC_B) = %.6f\n\n", sqrt(var_diff)))


















































