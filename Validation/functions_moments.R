### functions


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




#####################
### simulate data ###
#####################

#simulates one split-plot MRMC data set
# default variance components from table 4 (from those originally chosen by Roe and Metz)
sim_one_splitplot_cardy <- function( mu_nondisease = 0, 
                                     mu_disease = 1.53, 
                                     hypothesis = "alt",
                                     tau = 0.25,
                                     sigma_r = 0.011, sigma_c = 0.1,
                                     sigma_rc = 0.2, sigma_tr = 0.03,
                                     sigma_tc = 0.1, sigma_trc = 0.2, 
                                     block_ss = 60,
                                     readers_per_block = 3,
                                     blocks = 2,
                                     modalities = 2) {
  
  total_ss<-block_ss*blocks
  readers<-readers_per_block*blocks
  
  # set tau based null vs. alternative hypothesis selection
  if (hypothesis == "alt") {
    tau_nondisease_m1<- 0
    tau_nondisease_m2<- 0
    tau_disease_m1<- 0
    tau_disease_m2<- tau
  } else {
    tau_nondisease_m1<- 0
    tau_disease_m1<- 0
    tau_nondisease_m2<- 0
    tau_disease_m2<- 0
  }
  
  
  #simulate each random effects term in the model (equation A2)
  #random-reader effect
  R<-rnorm(readers, 0, sqrt(sigma_r))
  
  #random-case effect
  C<- rnorm(total_ss, 0, sqrt(sigma_c))
  
  #random modality-reader effect
  TR<- rnorm(readers*modalities, 0, sqrt(sigma_tr))
  
  #random modality-case effect
  TC<- rnorm(total_ss*modalities, 0, sqrt(sigma_tc))
  
  #random reader-case effect
  RC<- rnorm(readers*total_ss, 0, sqrt(sigma_rc))
  
  #random effect due to pure error (also includes 3 way interaction per this paper https://pmc.ncbi.nlm.nih.gov/articles/PMC9497942/)
  TRC<- rnorm(readers*total_ss*modalities, 0, sqrt(sigma_trc))
  
  
  # # even diseased/non-diseased split
  # disease_cases <- total_ss / 2  
  # nondisease_cases <- total_ss / 2 
  # 
  # # cases per block (rounded up to ensure every case is assigned a truth state)
  # #block_ss <- ceiling(total_ss / blocks)
  # 
  # base_disease_per_block <- floor(disease_cases / blocks)
  # base_nondisease_per_block <- floor(nondisease_cases / blocks)
  # 
  # extra_disease_cases <- disease_cases %% blocks
  # extra_nondisease_cases <- nondisease_cases %% blocks
  # 
  # disease_distribution <- rep(base_disease_per_block, blocks)
  # nondisease_distribution <- rep(base_nondisease_per_block, blocks)
  # 
  # disease_distribution[sample(1:blocks, extra_disease_cases)] <- disease_distribution[sample(1:blocks, extra_disease_cases)] + 1
  # nondisease_distribution[sample(1:blocks, extra_nondisease_cases)] <- nondisease_distribution[sample(1:blocks, extra_nondisease_cases)] + 1
  # 
  # # truth assignments by block
  # truth_blocks <- unlist(lapply(1:blocks, function(b) {
  #   c(rep("disease", modalities * disease_distribution[b]),
  #     rep("nondisease", modalities * nondisease_distribution[b]))
  # }))
  # 
  # 
  # truth <- rep(truth_blocks, times = readers)
  
  # target number of disease/nondisease cases overall
  disease_total <- total_ss / 2
  nondisease_total <- total_ss / 2
  
  # start with floor(block_ss / 2) disease per block
  disease_per_block <- rep(floor(block_ss / 2), blocks)
  nondisease_per_block <- rep(block_ss, blocks) - disease_per_block
  
  # adjust to ensure total matches disease_total
  disease_adjustment <- disease_total - sum(disease_per_block)
  
  if (disease_adjustment != 0) {
    adjust_blocks <- sample(1:blocks, abs(disease_adjustment))
    disease_per_block[adjust_blocks] <- disease_per_block[adjust_blocks] + sign(disease_adjustment)
    nondisease_per_block[adjust_blocks] <- nondisease_per_block[adjust_blocks] - sign(disease_adjustment)
  }
  
  # create truth values for each block, expanding by modalities
  truth_blocks <- unlist(lapply(1:blocks, function(b) {
    c(rep("disease", modalities * disease_per_block[b]),
      rep("nondisease", modalities * nondisease_per_block[b]))
  }))
  
  # expand to full reader-modality-case structure
  truth <- rep(truth_blocks, times = readers)
  
  
  reader <- rep(1:readers, each = total_ss * modalities)
  modality <- rep(c(1, 2), times = total_ss * readers)
  case <- rep(rep(1:total_ss, each = modalities), times = readers)
  
  reader_block <- rep(rep(1:blocks, each = readers_per_block), each = block_ss * modalities)
  case_block <- rep(rep(1:blocks, each = block_ss * modalities), times = readers)
  
  
  
  #remove data to make split plot
  # print(paste("length reader:",length(reader)))
  # print(paste("length case:",length(case)))
  # print(paste("length truth:",length(truth)))
  # print(paste("length modality:",length(modality)))
  # print(paste("length case_block:",length(case_block)))
  # print(paste("length reader_block:",length(reader_block)))
  
  id_df <- data.frame(reader, case, truth, modality, case_block, reader_block)
  
  id_df <- id_df %>%
    mutate(truth_binary = ifelse(truth == "disease", 1, 0)) %>%
    filter(case_block == reader_block)
  
  id_df <- id_df %>%
    mutate(block = case_block) %>%
    select(-reader_block, -case_block)
  
  
  
  
  #generate scores
  dat<- id_df %>%
    mutate(mu = ifelse(truth == "disease", mu_disease, mu_nondisease),
           tau = case_when(truth == "disease" & modality == 1 ~ tau_disease_m1,
                           truth == "disease" & modality == 2 ~ tau_disease_m2,
                           truth == "nondisease" & modality == 1 ~ tau_nondisease_m1,
                           truth == "nondisease" & modality == 2 ~ tau_nondisease_m2),
           R = R[reader],
           C = C[case],
           TR = TR[reader + (modality-1)*readers],
           TC = TC[case+(modality-1)*(total_ss)],
           RC = RC[case + (reader-1)*total_ss],
           TRC = TRC[case + (modality-1)*total_ss + (reader-1)*total_ss*modalities],
           score = mu + tau + R + C + TR + TC + RC + TRC)
  
  
  dat_clean<-dat %>% 
    select(reader, case, modality, truth_binary, score, block) %>%
    mutate(reader = as.factor(paste0("reader", reader)),
           case = as.factor(paste0("case", case)),
           modality = as.factor(ifelse(modality == 1, "TestA", "TestB")),
           block = as.factor(paste0("block", block)))
  
  
  return(dat_clean)
}

########################################
### Return all 8 moments as dataframe ###
########################################

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
  zcrit <- qnorm(1 - alpha/2)
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












###########################
### Simulation Function ###
###########################


#simulates one split-plot MRMC data set
# default variance components from table 4 (from those originally chosen by Roe and Metz)
sim_one_splitplot_cardy <- function( mu_nondisease = 0, 
                                     mu_disease = 1.53, 
                                     hypothesis = "alt",
                                     tau = 0.25,
                                     sigma_r = 0.011, sigma_c = 0.1,
                                     sigma_rc = 0.2, sigma_tr = 0.03,
                                     sigma_tc = 0.1, sigma_trc = 0.2, 
                                     block_ss = 60,
                                     readers_per_block = 3,
                                     blocks = 2,
                                     modalities = 2) {
  
  total_ss<-block_ss*blocks
  readers<-readers_per_block*blocks
  
  # set tau based null vs. alternative hypothesis selection
  if (hypothesis == "alt") {
    tau_nondisease_m1<- 0
    tau_nondisease_m2<- 0
    tau_disease_m1<- 0
    tau_disease_m2<- tau
  } else {
    tau_nondisease_m1<- 0
    tau_disease_m1<- 0
    tau_nondisease_m2<- 0
    tau_disease_m2<- 0
  }
  
  
  #simulate each random effects term in the model (equation A2)
  #random-reader effect
  R<-rnorm(readers, 0, sqrt(sigma_r))
  
  #random-case effect
  C<- rnorm(total_ss, 0, sqrt(sigma_c))
  
  #random modality-reader effect
  TR<- rnorm(readers*modalities, 0, sqrt(sigma_tr))
  
  #random modality-case effect
  TC<- rnorm(total_ss*modalities, 0, sqrt(sigma_tc))
  
  #random reader-case effect
  RC<- rnorm(readers*total_ss, 0, sqrt(sigma_rc))
  
  #random effect due to pure error (also includes 3 way interaction per this paper https://pmc.ncbi.nlm.nih.gov/articles/PMC9497942/)
  TRC<- rnorm(readers*total_ss*modalities, 0, sqrt(sigma_trc))
  
  
  # # even diseased/non-diseased split
  # disease_cases <- total_ss / 2  
  # nondisease_cases <- total_ss / 2 
  # 
  # # cases per block (rounded up to ensure every case is assigned a truth state)
  # #block_ss <- ceiling(total_ss / blocks)
  # 
  # base_disease_per_block <- floor(disease_cases / blocks)
  # base_nondisease_per_block <- floor(nondisease_cases / blocks)
  # 
  # extra_disease_cases <- disease_cases %% blocks
  # extra_nondisease_cases <- nondisease_cases %% blocks
  # 
  # disease_distribution <- rep(base_disease_per_block, blocks)
  # nondisease_distribution <- rep(base_nondisease_per_block, blocks)
  # 
  # disease_distribution[sample(1:blocks, extra_disease_cases)] <- disease_distribution[sample(1:blocks, extra_disease_cases)] + 1
  # nondisease_distribution[sample(1:blocks, extra_nondisease_cases)] <- nondisease_distribution[sample(1:blocks, extra_nondisease_cases)] + 1
  # 
  # # truth assignments by block
  # truth_blocks <- unlist(lapply(1:blocks, function(b) {
  #   c(rep("disease", modalities * disease_distribution[b]),
  #     rep("nondisease", modalities * nondisease_distribution[b]))
  # }))
  # 
  # 
  # truth <- rep(truth_blocks, times = readers)
  
  # target number of disease/nondisease cases overall
  disease_total <- total_ss / 2
  nondisease_total <- total_ss / 2
  
  # start with floor(block_ss / 2) disease per block
  disease_per_block <- rep(floor(block_ss / 2), blocks)
  nondisease_per_block <- rep(block_ss, blocks) - disease_per_block
  
  # adjust to ensure total matches disease_total
  disease_adjustment <- disease_total - sum(disease_per_block)
  
  if (disease_adjustment != 0) {
    adjust_blocks <- sample(1:blocks, abs(disease_adjustment))
    disease_per_block[adjust_blocks] <- disease_per_block[adjust_blocks] + sign(disease_adjustment)
    nondisease_per_block[adjust_blocks] <- nondisease_per_block[adjust_blocks] - sign(disease_adjustment)
  }
  
  # create truth values for each block, expanding by modalities
  truth_blocks <- unlist(lapply(1:blocks, function(b) {
    c(rep("disease", modalities * disease_per_block[b]),
      rep("nondisease", modalities * nondisease_per_block[b]))
  }))
  
  # expand to full reader-modality-case structure
  truth <- rep(truth_blocks, times = readers)
  
  
  reader <- rep(1:readers, each = total_ss * modalities)
  modality <- rep(c(1, 2), times = total_ss * readers)
  case <- rep(rep(1:total_ss, each = modalities), times = readers)
  
  reader_block <- rep(rep(1:blocks, each = readers_per_block), each = block_ss * modalities)
  case_block <- rep(rep(1:blocks, each = block_ss * modalities), times = readers)
  
  
  
  #remove data to make split plot
  # print(paste("length reader:",length(reader)))
  # print(paste("length case:",length(case)))
  # print(paste("length truth:",length(truth)))
  # print(paste("length modality:",length(modality)))
  # print(paste("length case_block:",length(case_block)))
  # print(paste("length reader_block:",length(reader_block)))
  
  id_df <- data.frame(reader, case, truth, modality, case_block, reader_block)
  
  id_df <- id_df %>%
    mutate(truth_binary = ifelse(truth == "disease", 1, 0)) %>%
    filter(case_block == reader_block)
  
  id_df <- id_df %>%
    mutate(block = case_block) %>%
    select(-reader_block, -case_block)
  
  
  
  
  #generate scores
  dat<- id_df %>%
    mutate(mu = ifelse(truth == "disease", mu_disease, mu_nondisease),
           tau = case_when(truth == "disease" & modality == 1 ~ tau_disease_m1,
                           truth == "disease" & modality == 2 ~ tau_disease_m2,
                           truth == "nondisease" & modality == 1 ~ tau_nondisease_m1,
                           truth == "nondisease" & modality == 2 ~ tau_nondisease_m2),
           R = R[reader],
           C = C[case],
           TR = TR[reader + (modality-1)*readers],
           TC = TC[case+(modality-1)*(total_ss)],
           RC = RC[case + (reader-1)*total_ss],
           TRC = TRC[case + (modality-1)*total_ss + (reader-1)*total_ss*modalities],
           score = mu + tau + R + C + TR + TC + RC + TRC)
  
  
  dat_clean<-dat %>% 
    select(reader, case, modality, truth_binary, score, block) %>%
    mutate(reader = as.factor(paste0("reader", reader)),
           case = as.factor(paste0("case", case)),
           modality = as.factor(ifelse(modality == 1, "TestA", "TestB")),
           block = as.factor(paste0("block", block)))
  
  
  return(dat_clean)
}



#simulate multiple split-plot MRMC studies
# n is the desired number of simulated datasets (set to 2000 in the paper)
sim_splitplot_cardy<- function(n=10,
                               mu_nondisease = 0, 
                               mu_disease = 1.5, 
                               hypothesis = "alt",
                               tau = 0.25,
                               sigma_r = 0.011, sigma_c = 0.1,
                               sigma_rc = 0.2, sigma_tr = 0.03,
                               sigma_tc = 0.1, sigma_trc = 0.2, 
                               #readers = 6,
                               #total_ss = 120,
                               block_ss = 60,
                               readers_per_block = 3,
                               blocks = 2,
                               modalities = 2) {
  sim_data<-list()
  for (i in 1:n) {
    sim_data[[i]]<- sim_one_splitplot_cardy(mu_nondisease = mu_nondisease, 
                                            mu_disease = mu_disease, 
                                            hypothesis = hypothesis,
                                            tau = tau,
                                            sigma_r = sigma_r, sigma_c = sigma_c,
                                            sigma_rc = sigma_rc, sigma_tr = sigma_tr,
                                            sigma_tc = sigma_tc, sigma_trc = sigma_trc, 
                                            block_ss = block_ss,
                                            readers_per_block = readers_per_block,
                                            blocks = blocks,
                                            modalities = modalities)
    #print(paste0("Simulated dataset ", i, " of ", n))
  }
  
  return(sim_data)
}

