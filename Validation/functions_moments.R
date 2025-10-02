### functions


library(iMRMC)
library(dplyr)
library(MRMCaov)
library(tidyverse)
library(gt)

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
    var_R = (inter_var/4)^2, #### divide by 4 to convert from range to var
    var_TR = (intra_var/4)^2, #### divide by 4 to convert from range to var
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
      `6` = var_TC + 2*(var_TR + var_error),
      `7` = var_TC + 2*(var_TR + var_error)
    ),
    sigmaOmega_l = list(
      # For covariance rows in Table 3; these also appear inside single-modality integrals (Eq. 12)
      `1` = 0,
      `2` = (var_C + var_RC),
      `3` = (var_C + var_RC),
      `4` = 2*(var_C + var_RC),
      `5` = 2*(var_R + var_RC),
      `6` = var_C + 2*(var_R + var_RC),
      `7` = var_C + 2*(var_R + var_RC)
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



# ####### original
# # Single-modality moment M_l for l=2..7 (Eq. 12)
# moment_single_l <- function(l, Delta, sigmas) {
#   sOm   <- sigmas$sigma_Omega
#   sMod  <- sigmas$sigma_A   # same form if called for B in this equal-variance case
#   sOm_l <- sigmas$sigmaOmega_l[[as.character(l)]]
#   sMod_l<- sigmas$sigmaA_l[[as.character(l)]]
# 
#   s_common <- sOm + sMod - sOm_l - sMod_l
#   denom   <- sqrt(sOm_l + sMod_l)
# 
#   integrand <- function(x) {
#     p <- pnorm( (Delta + x*sqrt(s_common)) / denom )
#     p*p * dnorm(x)
#   }
#   #integrate(integrand, lower = -Inf, upper = Inf, subdivisions = 400L, rel.tol = 1e-7)$value
#   integrate(integrand, lower = -Inf, upper = Inf, subdivisions = 2000L, rel.tol = 1e-12, abs.tol = 1e-12)$value
# }
# ### try "In our software, we simply sample the 1-D integral at 256 points on the interval (−10;10) and use the midpoint rule (rectangle rule)."
# ### Appendix A Gallas, Hillis under eq 24
# 
# 
# 
# 
# # Cross-modality moment MAB_l for l=1..7 (Eq. 15)
# moment_cross_l <- function(l, DeltaA, DeltaB, sigmas) {
#   sOm   <- sigmas$sigma_Omega
#   sOm_l <- sigmas$sigmaOmega_l[[as.character(l)]]
#   # In equal-variance case σA=σB
#   sA <- sigmas$sigma_A
#   sB <- sigmas$sigma_B
# 
#   s_common <- sOm - sOm_l
#   denomA   <- sqrt(sA + sOm_l)
#   denomB   <- sqrt(sB + sOm_l)
# 
#   integrand <- function(x) {
#     pA <- pnorm( (DeltaA + x*sqrt(s_common)) / denomA )
#     pB <- pnorm( (DeltaB + x*sqrt(s_common)) / denomB )
#     pA * pB * dnorm(x)
#   }
#   #integrate(integrand, lower = -Inf, upper = Inf, subdivisions = 400L, rel.tol = 1e-7)$value
#   integrate(integrand, lower = -Inf, upper = Inf, subdivisions = 2000L, rel.tol = 1e-12, abs.tol = 1e-12)$value
# }








# ### try "In our software, we simply sample the 1-D integral at 256 points on the interval (−10;10) and use the midpoint rule (rectangle rule)."
# ### Appendix A Gallas, Hillis under eq 24

midpoint_integral <- function(f, a = -10, b = 10, n = 256) {
  h <- (b - a) / n
  midpoints <- seq(a + h/2, b - h/2, length.out = n)
  h * sum(f(midpoints))
}


moment_single_l <- function(l, Delta, sigmas, n_points = 256) {
  sOm   <- sigmas$sigma_Omega
  sMod  <- sigmas$sigma_A
  sOm_l <- sigmas$sigmaOmega_l[[as.character(l)]]
  sMod_l<- sigmas$sigmaA_l[[as.character(l)]]

  s_common <- sOm + sMod - sOm_l - sMod_l
  denom   <- sqrt(sOm_l + sMod_l)

  f <- function(x) {
    p <- pnorm((Delta + x * sqrt(s_common)) / denom)
    p * p * dnorm(x)
  }

  midpoint_integral(f, a = -10, b = 10, n = n_points)
}
#
#
moment_cross_l <- function(l, DeltaA, DeltaB, sigmas, n_points = 256) {
  sOm   <- sigmas$sigma_Omega
  sOm_l <- sigmas$sigmaOmega_l[[as.character(l)]]
  sA <- sigmas$sigma_A
  sB <- sigmas$sigma_B

  s_common <- sOm - sOm_l
  denomA   <- sqrt(sA + sOm_l)
  denomB   <- sqrt(sB + sOm_l)

  f <- function(x) {
    pA <- pnorm((DeltaA + x * sqrt(s_common)) / denomA)
    pB <- pnorm((DeltaB + x * sqrt(s_common)) / denomB)
    pA * pB * dnorm(x)
  }

  midpoint_integral(f, a = -10, b = 10, n = n_points)
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



#########################################
### Return all 8 moments as dataframe ###
#########################################

compute_moments_df <- function(DeltaA, DeltaB, sigmas, AUCA, AUCB) {
  # build vectors
  MA  <- moments_single(DeltaA, sigmas, AUCA)
  MB  <- moments_single(DeltaB, sigmas, AUCB)
  MAB <- moments_cross(DeltaA, DeltaB, sigmas, AUCA, AUCB)
  
  # assemble tidy dataframe
  df <- data.frame(
    modality = c("A", "B", "Cross"),
    M1  = c(MA[1],  MB[1],  MAB[1]),
    M2  = c(MA[2],  MB[2],  MAB[2]),
    M3  = c(MA[3],  MB[3],  MAB[3]),
    M4  = c(MA[4],  MB[4],  MAB[4]),
    M5  = c(MA[5],  MB[5],  MAB[5]),
    M6  = c(MA[6],  MB[6],  MAB[6]),
    M7  = c(MA[7],  MB[7],  MAB[7]),
    M8  = c(MA[8],  MB[8],  MAB[8])
  )
  return(df)
}



#############################################################
### Moments -> Var(ΔA) -> Power (Chen, Gong, Gallas 2018) ###
#############################################################

# Build c1..c4 for given N0, N1 (used in VR^Δ and VC^Δ)
.c_coeffs <- function(N0, N1) {
  c1 <- 1/(N0*N1)
  c2 <- (N0 - 1)/(N0*N1)
  c3 <- (N1 - 1)/(N0*N1)
  c4 <- (N0 - 1)*(N1 - 1)/(N0*N1)
  list(c1=c1, c2=c2, c3=c3, c4=c4)
}

# Given MA, MB, MAB (each length-8), construct Δ-moments (length-8)
delta_moments <- function(MA, MB, MAB) {
  # Mk^Δ = Mk^(1) + Mk^(2) - 2 Mk^(1×2)
  MA + MB - 2*MAB
}

# Decompose Var(ΔA) into VR^Δ and VC^Δ using Mk^Δ and c1..c4
VR_VC_delta <- function(MD, N0, N1) {
  stopifnot(length(MD) == 8)
  cc <- .c_coeffs(N0, N1)
  with(cc, {
    VRd <- c1*(MD[1]-MD[5]) + c2*(MD[2]-MD[6]) + c3*(MD[3]-MD[7]) + c4*(MD[4]-MD[8])
    VCd <- c1*MD[5] + c2*MD[6] + c3*MD[7] - (1 - c4)*MD[8]
    list(VR_delta = VRd, VC_delta = VCd)
  })
}

# Var(ΔA) for different designs
# - FC: Var = (1/NR)*VR^Δ + VC^Δ
# - PSP with equal reader-group sizes (G groups): Var = (1/NR)*VR^Δ + (1/G)*VC^Δ
# - PSP with unequal group sizes: Var = (1/NR)*VR^Δ + (sum NR_g^2 / NR^2) * VC^Δ
var_deltaA <- function(MA, MB, MAB, N0, N1, NR,
                       design = c("FC","PSP_equal","PSP_unequal"),
                       G = NULL, NR_groups = NULL) {
  design <- match.arg(design)
  MD <- delta_moments(MA, MB, MAB)
  parts <- VR_VC_delta(MD, N0, N1)
  VRd <- parts$VR_delta
  VCd <- parts$VC_delta
  
  if (design == "FC") {
    return( (1/NR)*VRd + VCd )
  }
  if (design == "PSP_equal") {
    if (is.null(G)) stop("Provide G (number of reader groups) for PSP_equal.")
    return( (1/NR)*VRd + (1/G)*VCd )
  }
  # PSP_unequal
  if (is.null(NR_groups)) stop("Provide NR_groups (vector of reader counts per group) for PSP_unequal.")
  w <- sum(NR_groups^2) / (sum(NR_groups)^2)
  (1/NR)*VRd + w*VCd
}

# Power calculators
# Two-sided superiority: α (default 0.05)
power_two_sided <- function(deltaA, var_deltaA, alpha = 0.05) {
  zcrit <- round(qnorm(1 - alpha/2), 2)
  sdD <- sqrt(var_deltaA)
  pnorm(deltaA/sdD - zcrit) + pnorm(-deltaA/sdD - zcrit)
}

# One-sided non-inferiority:
# H0: A2 - A1 <= -δ  vs  H1: A2 - A1 > -δ
# Using Chen+2018 convention: Power = Φ( (ΔA + δ)/sd - z_{α/2} )
# (You can change to z_{α} if you prefer one-sided critical value.)
power_noninferiority <- function(deltaA, var_deltaA, delta_margin, alpha = 0.05, use_two_sided_z = TRUE) {
  zcrit <- if (use_two_sided_z) qnorm(1 - alpha/2) else qnorm(1 - alpha)
  sdD <- sqrt(var_deltaA)
  pnorm( (deltaA + delta_margin)/sdD - zcrit )
}








# ---- helper to evaluate one scenario across designs ----
evaluate_power_all <- function(params, N0, N1, NR, alpha = 0.05) {
  sigmas <- build_sigma_sums(params)
  
  delta1 <- qnorm(params$AUC1) * sqrt(sigmas$sigma_Omega + sigmas$sigma_A)
  delta2 <- qnorm(params$AUC2) * sqrt(sigmas$sigma_Omega + sigmas$sigma_B)
  
  # AUCA <- expected_auc(delta1, sigmas$sigma_Omega, sigmas$sigma_A)
  # AUCB <- expected_auc(delta2, sigmas$sigma_Omega, sigmas$sigma_B)
  # deltaA <- AUCB - AUCA
  AUCA <- params$AUC1
  AUCB <- params$AUC2
  deltaA <- AUCB - AUCA
  
  moments_df <- compute_moments_df(delta1, delta2, sigmas, AUCA, AUCB)
  
  MA  <- as.numeric(moments_df[moments_df$modality == "A", paste0("M", 1:8)])
  MB  <- as.numeric(moments_df[moments_df$modality == "B", paste0("M", 1:8)])
  MAB <- as.numeric(moments_df[moments_df$modality == "Cross", paste0("M", 1:8)])
  
  # full variance decomposition for FC
  v_fc <- var_deltaA(MA, MB, MAB, N0, N1, NR, design="FC")
  
  # Δ-moments and VR/VC decomposition
  MD <- delta_moments(MA, MB, MAB)
  parts <- VR_VC_delta(MD, N0, N1)
  vR <- parts$VR_delta
  vC <- parts$VC_delta
  
  # design-specific variances
  variances <- c(
    FC   = (1/NR) * vR + vC,         # Eq. (11)
    PSP2 = (1/NR) * vR + (1/2) * vC, # Eq. (12) with G=2
    PSP4 = (1/NR) * vR + (1/4) * vC, # G=4
    PSP8 = (1/NR) * vR + (1/8) * vC  # G=8
  )
  
  # compute power for each design
  powers <- map_dbl(variances, ~ power_two_sided(deltaA, .x, alpha))
  
  tibble(
    structure = params$structure,
    AUC1 = params$AUC1,
    AUC2 = params$AUC2,
    delta_true = deltaA,
    design = names(variances),
    var_delta = variances,
    sd_delta = sqrt(variances),
    power = powers
  )
}





###################################################################
### check power by plugging in variance of delta a from table 2 ###
###################################################################
manual_vars <- data.frame(
  structure = c("HH","HH","HH","HL","HL","HL",
                "LH","LH","LH","LL","LL","LL"),
  AUC1      = c(0.65,0.80,0.90,
                0.65,0.80,0.90,
                0.65,0.80,0.90,
                0.65,0.80,0.90),
  AUC2      = c(0.70,0.85,0.95,
                0.70,0.85,0.95,
                0.70,0.85,0.95,
                0.70,0.85,0.95),
  var_FC    = c(0.00142,0.00096,0.00042,  # Table 2 values here
                0.00134,0.00077,0.00028,
                0.00071,0.00055,0.00027, 
                0.00063,0.00035,0.00012),
  var_PSP2  = c(0.00083,0.00063,0.00030,  # Table 2 values here
                0.00075,0.00043,0.00015,
                0.00052,0.00045,0.00023, 
                0.00044,0.00025,0.000083),
  var_PSP4  = c(0.00053,0.00046,0.00024,  # Table 2 values here
                0.00045,0.00026,0.00009,
                0.00043,0.00040,0.00022, 
                0.00034,0.00020,0.000067),
  var_PSP8  = c(0.00038,0.00038,0.00021,  # Table 2 values here
                0.00030,0.00017,0.00006,
                0.00038,0.00038,0.00021, 
                0.00030,0.00017,0.000059)
)


manual_long <- manual_vars %>%
  tidyr::pivot_longer(cols = starts_with("var_"),
                      names_to = "design",
                      values_to = "var_delta") %>%
  dplyr::mutate(design = sub("var_", "", design),
                delta_true = AUC2 - AUC1,
                sd_delta = sqrt(var_delta),
                power = power_two_sided(delta_true, var_delta, alpha = 0.05))


manual_table2 <- manual_long %>%
  mutate(var_x1e3 = 1000 * var_delta,
         power_pct = 100 * power_two_sided(AUC2 - AUC1, var_delta)) %>%
  select(structure, AUC1, AUC2, design, var_x1e3, power_pct) %>%
  pivot_wider(names_from = design, 
              values_from = c(var_x1e3, power_pct))

power_from_manual_var<-manual_table2 %>%
  gt(rowname_col = "structure") %>%
  tab_spanner(label = "Variance (×1000)",
              columns = starts_with("var_x1e3")) %>%
  tab_spanner(label = "Power (%)",
              columns = starts_with("power_pct")) %>%
  cols_label(
    var_x1e3_FC   = "FC", var_x1e3_PSP2 = "PSP2",
    var_x1e3_PSP4 = "PSP4", var_x1e3_PSP8 = "PSP8",
    power_pct_FC   = "FC", power_pct_PSP2 = "PSP2",
    power_pct_PSP4 = "PSP4", power_pct_PSP8 = "PSP8"
  ) %>%
  fmt_number(columns = starts_with("var_x1e3"), decimals = 3) %>%
  fmt_number(columns = starts_with("power_pct"), decimals = 3)



power_two_sided_check <- function(deltaA, var_deltaA, alpha = 0.05) {
  zcrit <- qnorm(1 - alpha/2)
  sdD <- sqrt(var_deltaA)
  pnorm(deltaA/sdD - zcrit) + pnorm(-deltaA/sdD - zcrit)
}

evaluate_power_check <- function(params, N0, N1, NR, alpha = 0.05) {
  sigmas <- build_sigma_sums(params)
  
  delta1 <- qnorm(params$AUC1) * sqrt(sigmas$sigma_Omega + sigmas$sigma_A)
  delta2 <- qnorm(params$AUC2) * sqrt(sigmas$sigma_Omega + sigmas$sigma_B)
  
  # AUCA <- expected_auc(delta1, sigmas$sigma_Omega, sigmas$sigma_A)
  # AUCB <- expected_auc(delta2, sigmas$sigma_Omega, sigmas$sigma_B)
  # deltaA <- AUCB - AUCA
  AUCA <- params$AUC1
  AUCB <- params$AUC2
  deltaA <- AUCB - AUCA
  
  moments_df <- compute_moments_df(delta1, delta2, sigmas, AUCA, AUCB)
  
  MA  <- as.numeric(moments_df[moments_df$modality == "A", paste0("M", 1:8)])
  MB  <- as.numeric(moments_df[moments_df$modality == "B", paste0("M", 1:8)])
  MAB <- as.numeric(moments_df[moments_df$modality == "Cross", paste0("M", 1:8)])
  
  # full variance decomposition for FC
  v_fc <- var_deltaA(MA, MB, MAB, N0, N1, NR, design="FC")
  
  # Δ-moments and VR/VC decomposition
  MD <- delta_moments(MA, MB, MAB)
  parts <- VR_VC_delta(MD, N0, N1)
  vR <- parts$VR_delta
  vC <- parts$VC_delta
  
  # design-specific variances
  # variances <- c(
  #   FC   = (1/NR) * vR + vC,         # Eq. (11)
  #   PSP2 = (1/NR) * vR + (1/2) * vC, # Eq. (12) with G=2
  #   PSP4 = (1/NR) * vR + (1/4) * vC, # G=4
  #   PSP8 = (1/NR) * vR + (1/8) * vC  # G=8
  # )
  
  variances_check<- c(FC = 1.42, PSP2 = )
  
  # compute power for each design
  powers <- map_dbl(variances, ~ power_two_sided(deltaA, .x, alpha))
  
  tibble(
    structure = params$structure,
    AUC1 = params$AUC1,
    AUC2 = params$AUC2,
    delta_true = deltaA,
    design = names(variances),
    var_delta = variances,
    sd_delta = sqrt(variances),
    power = powers
  )
}











### power function that takes in OR parameters

evaluate_power_new <- function(params, NR, alpha = 0.05) {
  N0<-params$n0
  N1<-params$n1
  
  delta1 <- params$delta1
  delta2 <- params$delta2
  
  rename_lookup <- c(
    var_R = "sigma_r",
    var_TR = "sigma_tr",
    var_C = "sigma_C",
    var_TC = "sigma_TC",
    var_RC = "sigma_RC",
    var_error = "sigma_trc"
  )
  
  params_for_sigma_sums <- params %>%
    set_names(rename_lookup[names(params)])
  
  sigmas <- build_sigma_sums(params_for_sigma_sums)
  
  # delta1 <- qnorm(params$AUC1) * sqrt(sigmas$sigma_Omega + sigmas$sigma_A)
  # delta2 <- qnorm(params$AUC2) * sqrt(sigmas$sigma_Omega + sigmas$sigma_B)
  
  AUCA <- expected_auc(delta1, sigmas$sigma_Omega, sigmas$sigma_A)
  AUCB <- expected_auc(delta2, sigmas$sigma_Omega, sigmas$sigma_B)
  deltaA <- AUCB - AUCA
  # AUCA <- params$AUC1
  # AUCB <- params$AUC2
  # deltaA <- AUCB - AUCA
  
  moments_df <- compute_moments_df(params$delta1, params$delta2, sigmas, AUCA, AUCB)
  
  MA  <- as.numeric(moments_df[moments_df$modality == "A", paste0("M", 1:8)])
  MB  <- as.numeric(moments_df[moments_df$modality == "B", paste0("M", 1:8)])
  MAB <- as.numeric(moments_df[moments_df$modality == "Cross", paste0("M", 1:8)])
  
  # full variance decomposition for FC
  v_fc <- var_deltaA(MA, MB, MAB, N0, N1, NR, design="FC")
  
  # Δ-moments and VR/VC decomposition
  MD <- delta_moments(MA, MB, MAB)
  parts <- VR_VC_delta(MD, N0, N1)
  vR <- parts$VR_delta
  vC <- parts$VC_delta
  
  # design-specific variances
  variances <- c(
    FC   = (1/NR) * vR + vC,         # Eq. (11)
    PSP2 = (1/NR) * vR + (1/2) * vC, # Eq. (12) with G=2
    PSP4 = (1/NR) * vR + (1/4) * vC, # G=4
    PSP8 = (1/NR) * vR + (1/8) * vC  # G=8
  )
  
  # compute power for each design
  powers <- map_dbl(variances, ~ power_two_sided(deltaA, .x, alpha))
  
  tibble(
    structure = params$structure,
    AUC1 = AUCA,
    AUC2 = AUCB,
    delta_true = deltaA,
    design = names(variances),
    var_delta = variances,
    sd_delta = sqrt(variances),
    power = powers
  )
}



###### convert to OR ######

convert_to_OR <- function(AUC1, AUC2, var_R, var_TR, var_C, var_TC, var_RC, var_error, n0=100, n1=100) {
  sigma_Omega <- 2*(var_R + var_C + var_RC)
  sigma_mod   <- 2*(var_TR + var_TC + var_error)
  
  delta1 <- qnorm(AUC1) * sqrt(sigma_Omega + sigma_mod)
  delta2 <- qnorm(AUC2) * sqrt(sigma_Omega + sigma_mod)
  
  or_params<- RMH_to_OR(
    n0 = n0,
    n1 = n1,
    b = 1, # this is true under the diseased/nondiseased equal variance assumption
    delta1 = delta1,
    delta2 = delta2,
    var_R = var_R,
    var_TR = var_TR,
    var_C = var_C,
    var_TC = var_TC,
    var_RC = var_RC,
    var_error = var_error)
  
  or_params
}




#############################################
### functions for sample size calculation ###
#############################################


# ---- Build S_k and constant term for analytic Var(ΔA) ----
# weight_w: design weight on V_C^Δ
#   FC:        w = 1
#   PSP_equal: w = 1/G
#   PSP_unequal: w = sum(NR_g^2)/(NR^2)  (pass this as weight_w)
build_var_coeffs <- function(MA, MB, MAB, NR, weight_w) {
  # Δ-moments:
  MD <- delta_moments(MA, MB, MAB)  # length-8
  
  # Coeffs for VR^Δ and VC^Δ per c1..c4:
  # VR^Δ = c1*(MD1-MD5) + c2*(MD2-MD6) + c3*(MD3-MD7) + c4*(MD4-MD8)
  T_R <- c(MD[1]-MD[5], MD[2]-MD[6], MD[3]-MD[7], MD[4]-MD[8])
  
  # VC^Δ = c1*MD5 + c2*MD6 + c3*MD7 - (1 - c4)*MD8
  #      = -MD8 + (c1*MD5 + c2*MD6 + c3*MD7 + c4*MD8)
  T_C <- c(MD[5], MD[6], MD[7], MD[8])  # the part multiplied by c's
  const_term <- weight_w * (-MD[8])     # the -MD8 part
  
  # Combine into S_k = (1/NR)*T_Rk + weight_w*T_Ck, for k=1..4
  S <- (1/NR) * T_R + weight_w * T_C   # length-4: S1..S4
  
  # Express Var = const_term + S1*c1 + S2*c2 + S3*c3 + S4*c4
  # with c1..c4 rewritten as:
  # c1 = 1/(N0 N1)
  # c2 = 1/N1 - 1/(N0 N1)
  # c3 = 1/N0 - 1/(N0 N1)
  # c4 = 1 - 1/N0 - 1/N1 + 1/(N0 N1)
  
  # Collect coefficients:
  const      <- const_term + S[4]
  coeff_1N0  <- S[3] - S[4]
  coeff_1N1  <- S[2] - S[4]
  coeff_1N0N1<- S[1] - S[2] - S[3] + S[4]
  
  list(const = const,
       A = coeff_1N0,
       B = coeff_1N1,
       C = coeff_1N0N1)
}

# ---- Required variance threshold from (alpha, power, delta_A) ----
required_variance <- function(deltaA, alpha = 0.05, power = 0.80) {
  zcrit  <- qnorm(1 - alpha/2)
  zbeta  <- qnorm(power)
  (deltaA / (zcrit + zbeta))^2
}

# ---- Closed-form N0 (then N1=r*N0) given ratio r ----
solve_cases_closed_form <- function(const, A, B, C, r, Vreq) {
  # Var(N0) = const + (A + B/r)/N0 + C/(r*N0^2) <= Vreq
  # Let x = 1/N0:   q x^2 + p x + (const - Vreq) <= 0
  p <- A + B / r
  q <- C / r
  a <- q
  b <- p
  c <- const - Vreq
  disc <- b*b - 4*a*c
  
  if (disc < 0) return(NULL) # infeasible at given r and parameters
  
  # For a>0 (typically), inequality holds between the two roots
  x1 <- (-b - sqrt(disc)) / (2*a)
  x2 <- (-b + sqrt(disc)) / (2*a)
  xr <- max(x1, x2)          # rightmost root (largest x)
  
  if (!is.finite(xr) || xr <= 0) return(NULL)
  
  N0 <- ceiling(1 / xr)
  if (!is.finite(N0) || N0 < 2) N0 <- 2L
  N0
}

# ---- Public API: analytic sample size for cases ----
# params: list with AUC1, AUC2 and sigma_* (your OR_to_RMH -> build_sigma_sums path)
# NR: number of readers; design: "FC", "PSP_equal", or "PSP_unequal"
# r: case ratio N1/N0; for PSP_equal set G; for PSP_unequal pass weight_w directly.
analytic_case_ss <- function(params, NR, design = c("FC","PSP_equal","PSP_unequal"),
                             r = 1, G = NULL, weight_w = NULL,
                             alpha = 0.05, target_power = 0.80) {
  design <- match.arg(design)
  # Build sigma sums
  sigmas <- build_sigma_sums(params)
  
  # Work with AUC targets directly (your code already does this consistently)
  AUCA   <- params$AUC1
  AUCB   <- params$AUC2
  deltaA <- AUCB - AUCA
  
  # Moments
  delta1 <- qnorm(AUCA) * sqrt(sigmas$sigma_Omega + sigmas$sigma_A)
  delta2 <- qnorm(AUCB) * sqrt(sigmas$sigma_Omega + sigmas$sigma_B)
  mdf    <- compute_moments_df(delta1, delta2, sigmas, AUCA, AUCB)
  
  MA  <- as.numeric(mdf[mdf$modality == "A",     paste0("M",1:8)])
  MB  <- as.numeric(mdf[mdf$modality == "B",     paste0("M",1:8)])
  MAB <- as.numeric(mdf[mdf$modality == "Cross", paste0("M",1:8)])
  
  # Design weight w for V_C^Δ
  w <- switch(design,
              FC = 1.0,
              PSP_equal = {
                if (is.null(G)) stop("Provide G for PSP_equal.")
                1.0 / G
              },
              PSP_unequal = {
                if (is.null(weight_w)) stop("Provide weight_w for PSP_unequal (sum NR_g^2 / NR^2).")
                weight_w
              })
  
  # Build coefficients for Var(ΔA) = const + A/N0 + B/N1 + C/(N0*N1)
  K <- build_var_coeffs(MA, MB, MAB, NR = NR, weight_w = w)
  
  # Required variance for target power
  Vreq <- required_variance(deltaA, alpha = alpha, power = target_power)
  
  # Closed-form solution in N0 with N1 = r*N0
  N0 <- solve_cases_closed_form(const = K$const, A = K$A, B = K$B, C = K$C, r = r, Vreq = Vreq)
  if (is.null(N0)) {
    return(list(feasible = FALSE,
                message = "No finite solution for the given (NR, r, AUCs, variances) and power. Try increasing NR, relaxing power, or changing r."))
  }
  N1 <- ceiling(r * N0)
  
  # Report achieved power at the rounded integers (sanity check)
  # Build c-coeffs at these N0,N1 and compute analytic var & power
  cc <- .c_coeffs(N0, N1) # your helper
  # Reconstruct Var via VR/VC using your existing pathway for verification
  parts <- VR_VC_delta(delta_moments(MA, MB, MAB), N0, N1)
  vR <- parts$VR_delta
  vC <- parts$VC_delta
  v_design <- switch(design,
                     FC = (1/NR) * vR + 1.0   * vC,
                     PSP_equal = (1/NR) * vR + (1.0/G) * vC,
                     PSP_unequal = (1/NR) * vR + weight_w * vC)
  pow <- power_two_sided(deltaA, v_design, alpha = alpha)
  
  list(feasible = TRUE,
       design = design,
       NR = NR,
       r = r,
       AUC1 = AUCA, AUC2 = AUCB, deltaA = deltaA,
       N0 = N0, N1 = N1,
       var_delta = v_design,
       achieved_power = pow,
       target_power = target_power,
       alpha = alpha)
}


uniroot_case_ss <- function(params, NR, 
                             design = c("FC","PSP_equal","PSP_unequal"),
                             r = 1, G = NULL, NR_groups = NULL,
                             alpha = 0.05, target_power = 0.80,
                             search_min = 10, search_max = 5000) {
  design <- match.arg(design)
  
  # Extract deltas and AUCs
  delta1 <- qnorm(params$AUC1) * sqrt(2*(params$sigma_r + params$sigma_C + params$sigma_RC) +
                                        2*(params$sigma_tr + params$sigma_TC + params$sigma_trc))
  delta2 <- qnorm(params$AUC2) * sqrt(2*(params$sigma_r + params$sigma_C + params$sigma_RC) +
                                        2*(params$sigma_tr + params$sigma_TC + params$sigma_trc))
  
  AUCA <- params$AUC1
  AUCB <- params$AUC2
  deltaA <- AUCB - AUCA
  
  # Precompute sigma and moments (independent of N0/N1)
  sigmas <- build_sigma_sums(params)
  moments_df <- compute_moments_df(delta1, delta2, sigmas, AUCA, AUCB)
  MA  <- as.numeric(moments_df[moments_df$modality == "A", paste0("M", 1:8)])
  MB  <- as.numeric(moments_df[moments_df$modality == "B", paste0("M", 1:8)])
  MAB <- as.numeric(moments_df[moments_df$modality == "Cross", paste0("M", 1:8)])
  
  # Power function in terms of N0
  f_power <- function(N0) {
    N1 <- ceiling(r * N0)
    varD <- var_deltaA(MA, MB, MAB, N0, N1, NR, design = design, G = G, NR_groups = NR_groups)
    power_two_sided(deltaA, varD, alpha = alpha)
  }
  
  # Function for root finding: f(N0) - target_power = 0
  f_root <- function(N0) f_power(N0) - target_power
  
  # Check feasibility
  if (f_root(search_max) < 0) {
    return(list(feasible = FALSE, design = design, NR = NR, r = r,
                AUC1 = AUCA, AUC2 = AUCB, deltaA = deltaA,
                N0 = NA, N1 = NA, var_delta = NA,
                achieved_power = f_power(search_max),
                target_power = target_power, alpha = alpha))
  }
  
  # Solve for N0 using uniroot
  sol <- uniroot(f_root, lower = search_min, upper = search_max)
  
  N0_opt <- ceiling(sol$root)
  N1_opt <- ceiling(r * N0_opt)
  
  varD_opt <- var_deltaA(MA, MB, MAB, N0_opt, N1_opt, NR, design = design, G = G, NR_groups = NR_groups)
  power_opt <- power_two_sided(deltaA, varD_opt, alpha = alpha)
  
  list(feasible = TRUE,
       design = design,
       NR = NR,
       r = r,
       AUC1 = AUCA,
       AUC2 = AUCB,
       deltaA = deltaA,
       N0 = N0_opt,
       N1 = N1_opt,
       var_delta = varD_opt,
       achieved_power = power_opt,
       target_power = target_power,
       alpha = alpha)
}









