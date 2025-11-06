### functions


library(iMRMC)
library(dplyr)
library(MRMCaov)
library(tidyverse)
library(gt)
library(fpow)

# Specify OR parameters

library(dplyr)

# --- Table 1 (Obuchowski 2000): rangeb and rangew ---
# Each entry: c(rangeb, rangew)
var_table <- list(
  S = c(rangew = 0.005, rangeb = 0.01),
  M = c(rangew = 0.025, rangeb = 0.05),
  L = c(rangew = 0.05,  rangeb = 0.10)
)

# --- Build all OR-style scenarios (Table 2) ---
build_OR_scenarios <- function() {
  
  base_grid <- expand.grid(
    readers = c(4, 6, 10),
    observer_var = c("S", "M", "L"),
    accuracy_level = c(0.75, 0.90),
    delta = c(0.05, 0.10, 0.15),
    ratio = c(1, 2, 4),   # non-diseased : diseased ratio
    stringsAsFactors = FALSE
  ) %>%
    arrange(ratio, readers, observer_var, accuracy_level, delta) %>%
    mutate(r1 = 0.47,
           r2 = 0.30,
           r3 = 0.30,
           rb = 0.80)
  
  # Convert var_table list to data frame for joining
  var_df <- tibble(
    observer_var = names(var_table),
    rangeb = sapply(var_table, function(x) x["rangeb"]),
    rangew = sapply(var_table, function(x) x["rangew"])
  )
  
  left_join(base_grid, var_df, by = "observer_var")
}



####################################
### Convert OR -> RMH parameters ###
####################################

# wrapper for OR_to_RMH from MRMCaov

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
# (Gallas & Hillis 2014)

build_sigma_sums <- function(rmh) {
  var_R     <- rmh$sigma_r
  var_TR    <- rmh$sigma_tr
  var_C     <- rmh$sigma_C
  var_TC    <- rmh$sigma_TC
  var_RC    <- rmh$sigma_RC
  var_error <- rmh$sigma_trc
  
  # assuming equal variance for both truth states and both modalities
  sigma_Omega <- 2*(var_R + var_C + var_RC) #(Gallas & Hillis 2014)(eq. 3)
  sigma_mod   <- 2*(var_TR + var_TC + var_error)  #(Gallas & Hillis 2014)(eq 4 & 5)
  
  list(
    sigma_Omega = sigma_Omega,
    sigma_A     = sigma_mod,
    sigma_B     = sigma_mod,
    # For Table 3 rows (equal-components simplifications): (Gallas & Hillis 2014)
    sigmaA_l = list(
      `2` = (var_TC + var_error),
      `3` = (var_TC + var_error),
      `4` = 2*(var_TC + var_error),
      `5` = 2*(var_TR + var_error),
      `6` = var_TC + 2*(var_TR + var_error),
      `7` = var_TC + 2*(var_TR + var_error)
    ),
    sigmaOmega_l = list(
      # appear inside single-modality integrals (Eq. 12) (Gallas & Hillis 2014)
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

#Coefficients c for fully-crossed design (Table 2) (Gallas & Hillis 2014)
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


# Expected AUC from Δ and totals (Eq. 9) (Gallas & Hillis 2014)
expected_auc <- function(Delta, sigma_Omega, sigma_mod) {
  pnorm( Delta / sqrt(sigma_Omega + sigma_mod) )
}



# ### try "In our software, we simply sample the 1-D integral at 256 points on the interval (−10;10) and use the midpoint rule (rectangle rule)."
# ### Appendix A (Gallas & Hillis 2014) under eq 24

midpoint_integral <- function(f, a = -10, b = 10, n = 256) {
  h <- (b - a) / n
  midpoints <- seq(a + h/2, b - h/2, length.out = n)
  h * sum(f(midpoints))
}


moment_single_l <- function(l, Delta_AB, sigmas, n_points = 256) {
  sOm   <- sigmas$sigma_Omega
  sMod  <- sigmas$sigma_A
  sOm_l <- sigmas$sigmaOmega_l[[as.character(l)]]
  sMod_l<- sigmas$sigmaA_l[[as.character(l)]]

  s_common <- sOm + sMod - sOm_l - sMod_l
  denom   <- sqrt(sOm_l + sMod_l)
  
  #

  f <- function(x) {
    p <- pnorm((Delta_AB + x * sqrt(s_common)) / denom)
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


###### try putting this^ in the moments_cross and moments_single functions directly?








# Build all 8 moments for a single modality (A or B) (Gallas & Hillis 2014)
moments_single <- function(Delta_AB, sigmas, AUC) {
  # M1 = AUC ; M8 = AUC^2 ; M2..M7 via integral
  M <- numeric(8)
  M[1] <- AUC
  for (l in 2:7) {
    sOm   <- sigmas$sigma_Omega
    sMod  <- sigmas$sigma_A
    sOm_l <- sigmas$sigmaOmega_l[[as.character(l)]]
    sMod_l<- sigmas$sigmaA_l[[as.character(l)]]
    
    s_common <- sOm + sMod - sOm_l - sMod_l
    denom   <- sqrt(sOm_l + sMod_l)
    
    #
    
    f <- function(x) {
      p <- pnorm((Delta_AB + x * sqrt(s_common)) / denom)
      p * p * dnorm(x)
    }
    
    M[l] <- midpoint_integral(f, a = -10, b = 10, n = 256)
  }
    
  #M[l] <- moment_single_l(l, Delta_AB2, sigmas)
  M[8] <- AUC*AUC
  M
}

# Build all 8 cross-modality moments MAB (Gallas & Hillis 2014)
moments_cross <- function(DeltaA, DeltaB, sigmas, AUCA, AUCB) {
  M <- numeric(8)
  # l = 1..7 via integral with σΩ(l) (Table 3)
  for (l in 1:7) {
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
    
    M[l] <- midpoint_integral(f, a = -10, b = 10, n = 256)
  }
    
    
  # M[l] <- moment_cross_l(l, DeltaA, DeltaB, sigmas)
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

# Build c1..c4 for given N0, N1 (used in VR^Δ and VC^Δ) (from Chen, Gong, Gallas 2018)
.c_coeffs <- function(N0, N1) {
  c1 <- 1/(N0*N1)
  c2 <- (N0 - 1)/(N0*N1)
  c3 <- (N1 - 1)/(N0*N1)
  c4 <- (N0 - 1)*(N1 - 1)/(N0*N1)
  list(c1=c1, c2=c2, c3=c3, c4=c4)
}

# Given MA, MB, MAB (each length-8), construct Δ-moments (length-8) (Chen, Gong, Gallas 2018) (Section 2.2)
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
  
  #design-specific variances
  variances <- c(
    FC   = (1/NR) * vR + vC,         # Eq. (11)
    PSP2 = (1/NR) * vR + (1/2) * vC, # Eq. (12) with G=2
    PSP4 = (1/NR) * vR + (1/4) * vC, # G=4
    PSP8 = (1/NR) * vR + (1/8) * vC  # G=8
  )
  
  #variances_check<- c(FC = 1.42, PSP2 = )
  
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
# inputs: list of parameters including n0, n1, delta1, delta2, var_R, var_TR, var_C, var_TC, var_RC, var_error
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


# takes in RMH params and gives sample size via uniroot
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






# takes MRMCsamplesize style inputs and gives sample size via uniroot
mrmc_ss <- function(readers, 
                    accuracy_level, 
                    delta, 
                    ratio, 
                    rangeb, 
                    rangew, 
                    r1, r2, r3, rb, 
                            design = c("FC","PSP_equal","PSP_unequal"),
                            G = NULL, NR_groups = NULL,
                            alpha = 0.05, target_power = 0.80,
                            search_min = 10, search_max = 5000) {
  design <- match.arg(design)
  
  
  #### needs to take in MRMC sample size style params (readers, accuracy_level, delta, ratio, rangeb, rangew, r1, r2, r3, rb)
  
  ### then these are converted to OR parameters internally
  
  sb <- rangeb/4
  sw <- rangew/4
  
  lambda <- ncparamF(alpha, 1-target_power, nu1 = 1, nu2 = 1*(readers-1)) #ndf (nu1) is always number of treatments - 1.
  K<-1
  
  varTR <- sb^2*(1-rb)
  varR <- sb^2
  varW <- sw^2
  num <- ((readers*delta^2)/ (2*lambda)) -  (varTR + sw^2/K)
  den <- (1-r1) + (readers-1)*(r2-r3)

  sigma.square.c = num/den

  varE <-  sigma.square.c + sw^2
  Cov1 <- r1*varE
  Cov2 <- r2*varE
  Cov3 <- r3*varE
  
  N_total <- 200
  n1 <- round(N_total / (1 + ratio))     # diseased
  n0 <- N_total - n1 # non-diseased
  
  
  vars_to_print<-data.frame(varR=varR, varTR=varTR, varE=varE, Cov1=Cov1, Cov2=Cov2, Cov3=Cov3)
  
  #print(vars_to_print)
  
  ### then those are converted to RMH parameters, also internally
  
  RMH_params <- OR_to_RMH(
    AUC1 = accuracy_level,
    AUC2 = accuracy_level + delta,
    var_R = varR,
    var_TR = varTR,
    corr1 = r1,
    corr2 = r2,
    corr3 = r3,
    n0 = n0,
    n1 = n1,
    var_error = varE,
    b_method = "specified",
    b_input = 1 # this is true under the diseased/nondiseased equal variance assumption
    # b_method = "mean_to_sigma", mean_sig_input = 2
  )
  

  
  params <- RMH_params %>%
    mutate(AUC1 = accuracy_level,
           AUC2 = accuracy_level + delta) %>%
    rename(sigma_r = var_R,
           sigma_tr = var_TR,
           sigma_C = var_C,
           sigma_TC = var_TC,
           sigma_RC = var_RC,
           sigma_trc = var_error)
  
  ### then the sample size logic below (from moments), can be used
  
  # Extract deltas and AUCs
  delta1 <- qnorm(params$AUC1) * sqrt(2*(params$sigma_r + params$sigma_C + params$sigma_RC) +
                                        2*(params$sigma_tr + params$sigma_TC + params$sigma_trc))
  delta2 <- qnorm(params$AUC2) * sqrt(2*(params$sigma_r + params$sigma_C + params$sigma_RC) +
                                        2*(params$sigma_tr + params$sigma_TC + params$sigma_trc))
  
  # # check delta1 and delta2 from both calculations are equal
  # print(paste0("delta1: ", round(delta1,3), ", delta2: ", round(delta2,3)))
  # print(paste0("delta1 from or_to_rmh: ", round(params$delta1,3), ", delta2 from or_to_rmh: ", round(params$delta2,3)))
  
  #print(params)
  
  AUCA <- params$AUC1
  AUCB <- params$AUC2
  deltaA <- AUCB - AUCA
  
  # Precompute sigma and moments (independent of N0/N1)
  sigmas <- build_sigma_sums(params)
  moments_df <- compute_moments_df(delta1, delta2, sigmas, AUCA, AUCB)
  MA  <- as.numeric(moments_df[moments_df$modality == "A", paste0("M", 1:8)])
  MB  <- as.numeric(moments_df[moments_df$modality == "B", paste0("M", 1:8)])
  MAB <- as.numeric(moments_df[moments_df$modality == "Cross", paste0("M", 1:8)])
  
  coeffs_v1 <- coeff_vector(N0 = n0, N1 = n1, NR = readers)
  coeffs_v1_to_print <- as.data.frame(t(coeffs_v1))        
  colnames(coeffs_v1_to_print) <- paste0("c", 1:8) 
  
  
  # coeffs from Chen, Gong, Gallas 2018 (not including reader)
  coeffs_v2 <- .c_coeffs(N0 = n0, N1 = n1)
  coeffs_v2_to_print <- as.data.frame(coeffs_v2)
  
  # Power function in terms of N0
  f_power <- function(N1) {
    N0 <- ceiling(ratio * N1)
    varD <- var_deltaA(MA, MB, MAB, N0, N1, readers, design = design, G = G, NR_groups = NR_groups)
    power_two_sided(deltaA, varD, alpha = alpha)
  }
  
  # Function for root finding: f(N1) - target_power = 0
  f_root <- function(N1) f_power(N1) - target_power
  
  # Check feasibility
  if (f_root(search_max) < 0) {
    return(list(feasible = FALSE, design = design, readers = readers, ratio = ratio,
                AUC1 = AUCA, AUC2 = AUCB, deltaA = deltaA,
                N0 = NA, N1 = NA, var_delta = NA,
                achieved_power = f_power(search_max),
                target_power = target_power, alpha = alpha))
  }
  
  # Solve for N0 using uniroot
  sol <- uniroot(f_root, lower = search_min, upper = search_max)
  
  N1_opt <- ceiling(sol$root)
  N0_opt <- ceiling(ratio * N1_opt)
  
  varD_opt <- var_deltaA(MA, MB, MAB, N0_opt, N1_opt, readers, design = design, G = G, NR_groups = NR_groups)
  power_opt <- power_two_sided(deltaA, varD_opt, alpha = alpha)
  
  list(feasible = TRUE,
       design = design,
       readers = readers,
       ratio = ratio,
       AUC1 = AUCA,
       AUC2 = AUCB,
       deltaA = deltaA,
       N0_pilot = n0,
       N1_pilot = n1,
       N0 = N0_opt,
       N1 = N1_opt,
       N_total = N0_opt + N1_opt,
       var_deltaA = varD_opt,
       achieved_power = power_opt,
       target_power = target_power,
       alpha = alpha,
       OR_variances = vars_to_print,
       moments = moments_df,
       coeffs_8 = coeffs_v1_to_print,
       coeffs_4 = coeffs_v2_to_print)
}





#########################################
### Format results to check with Java ###
#########################################

write_imrmc_summary <- function(x,
                                input_csv_path = "this_is_a_placeholder.csv",
                                out_path = "imrmc_summary_test.omrmc",
                                version = "4.0.3",
                                modality_names = c(A = "TestA", B = "TestB")) {
  stopifnot(is.list(x), is.character(input_csv_path), is.character(out_path))
  # pull core fields (with simple safety)
  nR <- x$readers %||% x[["NReader"]] %||% stop("Missing readers")
  n0 <- x$N0_pilot %||% stop("Missing N0")
  n1 <- x$N1_pilot %||% stop("Missing N1")
  AUC_A <- x$AUC1 %||% stop("Missing AUC1")
  AUC_B <- x$AUC2 %||% stop("Missing AUC2")
  
  # moments: expect a data.frame with columns: modality, M1..M8
  mm <- x$moments
  if (is.null(mm) || !all(c("modality", paste0("M", 1:8)) %in% names(mm))) {
    stop("x$moments must have columns: modality, M1..M8")
  }
  # enforce ordering A, B, Cross to match the desired three lines
  want_order <- c("A", "B", "Cross")
  if (!all(want_order %in% mm$modality)) {
    stop("moments$modality must include 'A', 'B', and 'Cross'")
  }
  mm <- mm[match(want_order, mm$modality), paste0("M", 1:8), drop = FALSE]
  
  # format rules to mimic your example:
  # - AUCs with 2 decimals
  # - M1..M7 with 7 decimals
  # - M8 with 4 decimals
  fmt_auc <- function(z) sprintf("%.2f", as.numeric(z))
  fmt_m <- function(v) {
    v <- as.numeric(v)
    c(sprintf("%.7f", v[1]),
      sprintf("%.7f", v[2]),
      sprintf("%.7f", v[3]),
      sprintf("%.7f", v[4]),
      sprintf("%.7f", v[5]),
      sprintf("%.7f", v[6]),
      sprintf("%.7f", v[7]),
      sprintf("%.4f", v[8]))
  }
  A_line_vals <- fmt_m(as.numeric(mm[1, ]))
  B_line_vals <- fmt_m(as.numeric(mm[2, ]))
  X_line_vals <- fmt_m(as.numeric(mm[3, ]))
  
  # modality names (shown in the small header block)
  modA_name <- modality_names[["A"]] %||% "TestA"
  modB_name <- modality_names[["B"]] %||% "TestB"
  
  # build lines exactly as shown
  lines <- c(
    sprintf("MRMC summary statistics from iMRMC Version %s", version),
    "Summary statistics based on input file named:",
    input_csv_path,
    "",
    "BEGIN SUMMARY",
    sprintf("NReader=  %d", as.integer(nR)),
    sprintf("Nnormal=  %d", as.integer(n0)),
    sprintf("NDisease= %d", as.integer(n1)),
    "",
    sprintf("Modality A = %s", modA_name),
    sprintf("Modality B = %s", modB_name),
    "",
    "Reader-Averaged AUCs",
    sprintf("AUC_A = %s", fmt_auc(AUC_A)),
    sprintf("AUC_B = %s", fmt_auc(AUC_B)),
    "",
    "",
    "**********************BDG Moments***************************",
    "         Moments,         M1,         M2,         M3,         M4,         M5,         M6,         M7,         M8",
    sprintf("Modality1(AUC_A), %s, %s, %s, %s, %s, %s, %s, %s,",
            A_line_vals[1], A_line_vals[2], A_line_vals[3], A_line_vals[4],
            A_line_vals[5], A_line_vals[6], A_line_vals[7], A_line_vals[8]),
    sprintf("Modality2(AUC_B), %s, %s, %s, %s, %s, %s, %s, %s,",
            B_line_vals[1], B_line_vals[2], B_line_vals[3], B_line_vals[4],
            B_line_vals[5], B_line_vals[6], B_line_vals[7], B_line_vals[8]),
    sprintf("    comp product, %s, %s, %s, %s, %s, %s, %s, %s,",
            X_line_vals[1], X_line_vals[2], X_line_vals[3], X_line_vals[4],
            X_line_vals[5], X_line_vals[6], X_line_vals[7], X_line_vals[8]),
    "",
    "END SUMMARY"
  )
  
  # write
  con <- file(out_path, open = "w", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(lines, con = con, sep = "\n")
  
  invisible(out_path)
}

# a tiny helper for null-coalescing (so this file is self-contained)
`%||%` <- function(a, b) if (!is.null(a)) a else b



