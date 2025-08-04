####### Split-plot MRMC Sample Size Simulation Functions #######

# These functions can be used to simulate split-plot MRMC datasets with different numbers of readers and blocks, as well as different variance parameters.


library(dplyr)
library(MRMCaov)
library(iMRMC)
library(MRMCsamplesize)
library(fpow)
library(purrr)
library(parallel)


#split-plot tools

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
  R<-list(rnorm(readers, 0, sqrt(sigma_r)), 
          rnorm(readers, 0, sqrt(sigma_r)))
  
  #random-case effect
  C<- list(rnorm(total_ss, 0, sqrt(sigma_c)),
           rnorm(total_ss, 0, sqrt(sigma_c)))
  
  #random modality-reader effect
  TR<- list(rnorm(readers*modalities, 0, sqrt(sigma_tr)),
            rnorm(readers*modalities, 0, sqrt(sigma_tr)))
  
  #random modality-case effect
  TC<- list(rnorm(total_ss*modalities, 0, sqrt(sigma_tc)),
            rnorm(total_ss*modalities, 0, sqrt(sigma_tc)))
  
  #random reader-case effect
  RC<- list(rnorm(readers*total_ss, 0, sqrt(sigma_rc)),
            rnorm(readers*total_ss, 0, sqrt(sigma_rc)))
  
  #random effect due to pure error (also includes 3 way interaction per this paper https://pmc.ncbi.nlm.nih.gov/articles/PMC9497942/)
  TRC<- list(rnorm(readers*total_ss*modalities, 0, sqrt(sigma_trc)),
             rnorm(readers*total_ss*modalities, 0, sqrt(sigma_trc)))
  
  
  # even diseased/non-diseased split
  disease_cases <- total_ss / 2  
  nondisease_cases <- total_ss / 2 
  
  # cases per block (rounded up to ensure every case is assigned a truth state)
  block_ss <- ceiling(total_ss / blocks)
  
  base_disease_per_block <- floor(disease_cases / blocks)
  base_nondisease_per_block <- floor(nondisease_cases / blocks)
  
  extra_disease_cases <- disease_cases %% blocks
  extra_nondisease_cases <- nondisease_cases %% blocks
  
  disease_distribution <- rep(base_disease_per_block, blocks)
  nondisease_distribution <- rep(base_nondisease_per_block, blocks)
  
  disease_distribution[sample(1:blocks, extra_disease_cases)] <- disease_distribution[sample(1:blocks, extra_disease_cases)] + 1
  nondisease_distribution[sample(1:blocks, extra_nondisease_cases)] <- nondisease_distribution[sample(1:blocks, extra_nondisease_cases)] + 1
  
  # truth assignments by block
  truth_blocks <- unlist(lapply(1:blocks, function(b) {
    c(rep("disease", modalities * disease_distribution[b]),
      rep("nondisease", modalities * nondisease_distribution[b]))
  }))
  

  truth <- rep(truth_blocks, times = readers)
  reader <- rep(1:readers, each = total_ss * modalities)
  modality <- rep(c(1, 2), times = total_ss * readers)
  case <- rep(rep(1:total_ss, each = modalities), times = readers)
  
  reader_block <- rep(rep(1:blocks, each = readers_per_block), each = block_ss * modalities)
  case_block <- rep(rep(1:blocks, each = block_ss * modalities), times = readers)
  
  
  
  #remove data to make split plot
  
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
           R = ifelse(truth == "disease", R[[1]][reader], R[[2]][reader]),
           C = ifelse(truth == "disease", C[[1]][case], C[[2]][case]),
           TR =  ifelse(truth == "disease", TR[[1]][reader + (modality-1)*readers], 
                        TR[[2]][reader + (modality-1)*readers]),
           TC = ifelse(truth == "disease", TC[[1]][case+(modality-1)*(total_ss)], 
                       TC[[2]][case+(modality-1)*(total_ss)]),
           RC = ifelse(truth == "disease", RC[[1]][case + (reader-1)*total_ss], 
                       RC[[2]][case + (reader-1)*total_ss]),
           TRC = ifelse(truth == "disease", TRC[[1]][case + (modality-1)*total_ss + (reader-1)*total_ss*modalities], 
                        TRC[[2]][case + (modality-1)*total_ss + (reader-1)*total_ss*modalities]),
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
    print(paste0("Simulated dataset ", i, " of ", n))
  }
  
  return(sim_data)
}



################################
#### verification functions ####
################################

pull_ustat_params <- function(fac_list) {
  suppressMessages({
    ures <- lapply(fac_list, function(sim_data) {
      # Prepare truth rows
      truth_rows <- sim_data %>%
        select(case, truth_binary) %>%
        distinct() %>%
        mutate(reader = "truth", modality = "truth", score = truth_binary) %>%
        select(-truth_binary)
      
      # Prepare data for iMRMC
      df_pivoted <- bind_rows(sim_data %>%
                                mutate(reader = as.character(reader),
                                       modality = as.character(modality)) %>%
                                select(-truth_binary, -block),
                              truth_rows) %>%
        rename(readerID = reader, modalityID = modality, caseID = case)
      
      # Run iMRMC
      imrmc_analysis <- doIMRMC(df_pivoted)
      
      # Per-reader AUCs (from index offset)
      readers <- length(unique(df_pivoted$readerID)) - 1
      auc_reader <- imrmc_analysis$perReader[(readers*2 + 1):(readers*2 + readers), ]
      
      # U-stat AUCs
      ustat <- imrmc_analysis$Ustat[3, ]
      theta_control <- ustat["AUCA"][[1]]
      theta_treat   <- ustat["AUCB"][[1]]
      delta         <- theta_treat - theta_control
      
      # Per-reader mean AUC (per modality)
      reader_means <- rowMeans(auc_reader[, c("AUCA", "AUCB")])
      
      # MRMCsampleSize-aligned theta: average across readers
      theta <- mean(reader_means)
      
      # Between-reader range of means
      rangeb <- diff(range(reader_means))
      
      # Inter-modality correlation across readers
      r1 <- cor(auc_reader$AUCA, auc_reader$AUCB)
      
      # Null rejection from U-stat
      rejectNull <- ifelse(ustat["pValueBDG"][[1]] < 0.05, 1, 0)
      
      # Return result
      data.frame(
        theta_control = theta_control,
        theta_treat   = theta_treat,
        theta         = theta,
        delta         = delta,
        rangeb        = rangeb,
        r1            = r1,
        rejectNull    = rejectNull
      )
    })
    
    bind_rows(ures)
  })
}


# pull_ustat_params<- function(fac_list) {
#   
#   suppressMessages({
#     
#     # consider Ustat results first because it is fastest and works for unbalanced designs
#     ures<- lapply(fac_list, function(sim_data) {
#       #use package iMRMC
#       #must format data based on iMRMC conventions
#       # Create the "truth" rows
#       truth_rows <- sim_data %>%
#         select(case, truth_binary) %>%
#         distinct() %>%
#         mutate(
#           reader = "truth",
#           modality = "truth",
#           score = truth_binary,
#         ) %>%
#         select(-truth_binary)
#       
#       
#       # Bind the original data with the "truth" rows
#       df_pivoted <- bind_rows(sim_data %>% 
#                                 mutate(modality = as.character(modality)) %>%
#                                 select(-truth_binary), truth_rows) %>%
#         mutate(readerID = reader,
#                modalityID = modality,
#                caseID = case) %>% select(-reader, -modality, -case)
# 
#       
#       
#       imrmc_analysis<-doIMRMC(df_pivoted)
#       
#       readers<-length(unique(df_pivoted$readerID))-1
#       
#       calc_r1<-imrmc_analysis$perReader[(readers*2 + 1):(readers*2 + readers),]
# 
#       
#       calc_r1_rangeb<- calc_r1 %>%
#         rowwise() %>%
#         mutate(
#           mean_AUC = (AUCA + AUCB) / 2,
#           cov_AUC = abs(((AUCA - mean_AUC) * (AUCB - mean_AUC)) / 2),
#           correlation = cov_AUC / (sqrt(varAUCA) * sqrt(varAUCB))
#         )
#       
#       
#       r1<- calc_r1_rangeb %>%
#         select(-N0, -N1, -modalityA, -modalityB, -AUCAminusAUCB, -varAUCAminusAUCB) %>%
#         as.data.frame() %>% summarise(r1 = mean(correlation)) %>% pull(r1)
#       
#       
#       theta_control<-imrmc_analysis$Ustat$AUCA[3]
#       theta_treat<-imrmc_analysis$Ustat$AUCB[3]
#       theta <- mean(c(imrmc_analysis$Ustat$AUCA[3], imrmc_analysis$Ustat$AUCB[3]))
#       
#       #rangeb<-max(calc_r1_rangeb$mean_AUC) - min(calc_r1_rangeb$mean_AUC)
#       rangeb <- range(calc_r1$AUCA + calc_r1$AUCB) / 2
#       
#       
#       
#       ustat<-imrmc_analysis$Ustat[3,]
#       
# 
#       
#       ustat_ci<-c(ustat["botCIBDG"][[1]], ustat["topCIBDG"][[1]])
#       
#       ustat_results<-data.frame(test_statistic = ustat["AUCAminusAUCB"], 
#                                 p_value = ustat["pValueBDG"], ci_lower=ustat_ci[1], 
#                                 ci_upper=ustat_ci[2], theta = theta, 
#                                 theta_control = theta_control, theta_treat = theta_treat,
#                                 rangeb = rangeb, r1 = r1)
#       
#       
#       return(ustat_results)
#     } )
#     
#     # Convert list to dataframe
#     df <- list_rbind(ures)
#     
# 
#     
#     # power
#     # proportion of rejections of the null in data simulated under the null hypothesis
#     # combine results for easier comparison
#     ustat_res<- df %>%
#       mutate(test_statistic = AUCAminusAUCB,
#              p_value = pValueBDG,
#              rejectNull = ifelse(p_value < 0.05, 1, 0)) %>%
#       select(-AUCAminusAUCB, -pValueBDG)
#     
#   })
#   
#   return(ustat_res)
# }





#function to get RMH (ROC-level like simulation inputs) parameters back from OR MRMCaov results
test_ORparams<- function(mrmc_res, cases = 120) {
  
  ### annoyingly, this won't work in a function (because of eval(substitute()) I think)
  # Error in eval(substitute(new_mrmc(response, test, reader, case, data,  : object 'sim_data' not found
  #so need to run mrmc() outside of function and pass it in
  # mrmc_res<-mrmc(response = empirical_auc(sim_data$truth, sim_data$score), 
  #                test = modality, 
  #                reader = reader,
  #                case = case,
  #                data = sim_data,
  #                cov = DeLong())
  
  mrmc_sum<-summary(mrmc_res)
  v_cov_comps<-mrmc_sum$vcov_comps
  aucs<-mrmc_sum$test_means$Estimate
  
  pars<-data.frame(AUC1 = aucs[1], AUC2 = aucs[2], 
                   var_R = abs(v_cov_comps$Estimate[1]), 
                   var_TR = abs(v_cov_comps$Estimate[2]), 
                   corr1 = abs(v_cov_comps$Correlation[4]), 
                   corr2 = abs(v_cov_comps$Correlation[5]),
                   corr3 = abs(v_cov_comps$Correlation[6]),
                   var_error = v_cov_comps$Estimate[3],
                   n0 = cases/2,
                   n1 = cases/2)
  
  
  OR_pars<-OR_to_RMH(params = pars, b_method = "specified", b_input = 1)
  
  return(OR_pars)
  
}




#function to get summary from simulation parameters
#AUC level parameters
summarise_ss_results<- function(ss_results,
                                readers,
                                total_ss) {
  res<-ss_results %>% 
    summarise(sim_ss = total_ss,
              readers = readers,
              mean_theta = mean(theta),
              mean_theta_control = mean(theta_control),
              mean_theta_treat = mean(theta_treat),
              mean_theta_diff = mean(theta_treat - theta_control),
              mean_rangeb = mean(rangeb),
              mean_r1 = mean(r1),
              power = sum(rejectNull)/n())
  return(res)
}




#using the mean of the OR parameters, calculate sample size from modified MRMCsamplesize calculation
#don't need rangew
#uses AUC-level parameters
# check_mrmc_sample_size<-function(ss_pars_mean,
#                                  J=20,
#                                  power = 0.8) {
#   R<-1
#   K<-1
# 
#   theta<-ss_pars_mean$theta
#   delta<-ss_pars_mean$delta
#   varTR<-ss_pars_mean$varTR
#   varE<-ss_pars_mean$varE
#   r1<-ss_pars_mean$r1
#   r2<-ss_pars_mean$r2
#   r3<-ss_pars_mean$r3
# 
#   # test12 <- summary(mrmc_res)
#   # test12$test_means %>% select(Estimate,StdErr) %>% mutate(var = StdErr^2,.keep='unused')
#   # test12$test_diffs %>% select(Comparison,Estimate,StdErr,`p-value`) %>% mutate(var = StdErr^2,.keep='unused')
# 
# 
#   #theta = mean(test12$test_means$Estimate)
#   A <-  stats::qnorm(theta)*1.414
#   var.theta <- ((0.0099 * exp(-A^2/2)) * ( (5*A^2 +8) + (A^2 +8)/R))
# 
#   # delta <- test12$test_diffs$Estimate %>% abs
#   #
#   # varTR <- test12$vcov_comps$Estimate[2]
#   # varE <- test12$vcov_comps$Estimate[3] #variance of error
#   # r1 = test12$vcov_comps$Correlation[4]
#   # r2 = test12$vcov_comps$Correlation[5]
#   # r3 = test12$vcov_comps$Correlation[6]
# 
# 
#   # lambda <- fpow::ncparamF(0.05, 1-power, nu1 = 1, nu2 = 1*(J-1)) #ndf (nu1) is always number of treatments - 1.
#   # #num <- ((J*delta^2)/ (2*lambda)) -  (varTR + sw^2/K)
#   # den <- (1-r1) + (J-1)*(r2-r3)
#   #
#   # num1 <- ((K*J*delta^2)/ (2*lambda)) - K*varTR - varE
#   # den1 <- K*(1-r1) + K*(J-1)*(r2-r3) -1
#   #
#   # sigma.square.c = num1/(den1)
# 
#   lambda <- fpow::ncparamF(0.05, 1-power, nu1 = 1, nu2 = 1*(J-1)) #ndf (nu1) is always number of treatments - 1.
#   #num <- ((J*delta^2)/ (2*lambda)) -  (varTR + sw^2/K)
#   den <- (1-r1) + (J-1)*(r2-r3)
# 
#   num1 <- ((J*delta^2)/ (2*lambda)) - varTR - varE
# 
#   sigma.square.c = num1/(den-1)
# 
#   #nUnits_i
#   ss<-ceiling(var.theta/sigma.square.c)
# 
#   #find rangew to compare...
#   sw2<-varE-sigma.square.c
#   rangew<-4*sqrt(abs(sw2))
# 
#   #print(c(num1, sigma.square.c, var.theta, ss))
# 
#   # if(sigma.square.c <=0){
#   #   stop("Number of readers (J) are not enough to power the MRMC study for the
#   #        given assumptions. Increase J and re-estimate the sample size OR be
#   #        less conservative in the current assumptions")}
# 
#   ss_res<-data.frame(lambda = lambda,
#                      numerator = num1,
#                      sigma.square.c = sigma.square.c,
#                      var.theta = var.theta,
#                      rangew = rangew,
#                      n_units = ss,
#                      total_ss = ss*2)
# 
#   return(ss_res)
# 
# }



############### Analysis Functions ####################

ustat_analysis_forsim<- function(sim_data) {
  suppressMessages({
    
    ###########################
    ###### U statistic ########
    ###########################
    
    #use package iMRMC
    #must format data based on iMRMC conventions
    # Create the "truth" rows
    truth_rows <- sim_data %>%
      select(case, truth_binary) %>%
      distinct() %>%
      mutate(
        reader = "truth",
        modality = "truth",
        score = truth_binary,
      ) %>%
      select(-truth_binary)
    
    # Bind the original data with the "truth" rows
    df_pivoted <- bind_rows(sim_data %>% 
                              mutate(reader = as.character(reader),
                                     modality = as.character(modality)) %>%
                              select(-truth_binary, -block), truth_rows) %>%
      rename(readerID = reader, 
             modalityID = modality, 
             caseID = case, 
             score = score)
    
    
    imrmc_analysis<-doIMRMC(df_pivoted)
    ustat<-imrmc_analysis$Ustat[3,]
    
    ustat_ci<-c(ustat["botCIBDG"][[1]], ustat["topCIBDG"][[1]])
    
    ustat_results<-data.frame(test_statistic = ustat["AUCAminusAUCB"], 
                              p_value = ustat["pValueBDG"], ci_lower=ustat_ci[1], 
                              ci_upper=ustat_ci[2])
    
    
  })
  
  return(ustat_results)
}


### function to run ustat analysis on multiple data sets
get_ustat_results<- function(list_data, total_readers = 6, modalities = 2, hypothesis = "alt") {
  results<- mapply(function(sim_dat) ustat_analysis_forsim(sim_dat), 
                   list_data, SIMPLIFY = FALSE)
  result_analysis <- bind_rows(lapply(results, 
                                      function(dat) dat %>%
                                        summarise(reject_null = sum(typeIerror),
                                                  coverage_count = sum(coverage), .groups = "drop")))
  final_res<- result_analysis %>%
    summarise(reject_null = sum(reject_null), coverage_total = sum(coverage_count))
  
  return(list(sim_metrics = final_res, results = results))
}



### function to find the minimum sample size using uniroot given readers per block, number of blocks, and other parameters
find_min_sample_size_uniroot <- function(readers_per_block = 3, 
                                         blocks = 2,
                                         sigma_r = 0.011, sigma_tr = 0.03,
                                         delta = 0.05, rangeb = 0.1, rangew = 0.05,
                                         theta = 0.75,
                                         mu_nondisease = 0, 
                                         mu_disease = 1.53, 
                                         hypothesis = "alt",
                                         tau = 0.25, # RMH params
                                         sigma_c = 0.1, # var_C
                                         sigma_rc = 0.2, # var_RC
                                         sigma_tc = 0.1,  # var_TC
                                         sigma_trc = 0.2, # var_error
                                         n_sim = 200, target_power = 0.80) {
  
  total_readers<- readers_per_block * blocks
  
  
  starting_point <- tryCatch({
    sampleSize_MRMC(endpoint = "AUC",
                    J = total_readers,
                    delta = delta,
                    rangeb = rangeb,
                    rangew = rangew,
                    theta = theta,
                    r1 = 0.47,
                    power = target_power)$ORSampleSizeResults$nTotal
  }, error = function(e) {
    1000  # conservative fallback if sampleSize_MRMC fails
    #should use largest possible sample size from sampleSize_MRMC here
  })
  
  
  
  ending_point <- tryCatch({
    sampleSize_MRMC(endpoint = "AUC",
                    J = readers_per_block,
                    delta = delta,
                    rangeb = rangeb,
                    rangew = rangew,
                    theta = theta,
                    r1 = 0.47,
                    power = target_power)$ORSampleSizeResults$nTotal
  }, error = function(e) {
    starting_point+500
  })
  
  if (ending_point == starting_point) {
    # If the ending point is the same as the starting point, increase it by 100
    #starting_point <- max(starting_point - 100,10)
    starting_point <- starting_point - 10
    ending_point <- ending_point + 100
  }
  
  print(paste("Starting point for sample size search:", starting_point))
  print(paste("Ending point for sample size search:", ending_point))
  
  
  power_diff <- function(block_ss) {
    block_ss <- round(block_ss)  # ensure it's an integer
    
    cat("Trying block sample size:", block_ss, "\n")
    
    data <- sim_splitplot_cardy(n = n_sim,
                                readers_per_block = readers_per_block,
                                block_ss = block_ss,
                                sigma_r = sigma_r,
                                sigma_tr = sigma_tr,
                                mu_nondisease = mu_nondisease, 
                                mu_disease = mu_disease,
                                tau = tau,
                                sigma_c = sigma_c,
                                sigma_rc = sigma_rc,
                                sigma_tc = sigma_tc, sigma_trc = sigma_trc,
                                hypothesis = "alt",
                                blocks = blocks)
    
    # Extract OR parameters from sim
    #ustat_df <- pull_ustat_params(data)
    
    # summary_stats <- ustat_df %>%
    #   summarise(
    #     mean_theta = mean(theta),
    #     mean_delta = mean(delta),
    #     mean_rangeb = mean(rangeb),
    #     mean_r1 = mean(r1),
    #     power = mean(rejectNull)
    #   )
    
    
    # Parallel processing
    parms_list1 <- mclapply(data, ustat_analysis_forsim, mc.cores = detectCores() - 1)
    pvals <- sapply(parms_list1, function(x) x$pValueBDG)
    power <- mean(pvals < 0.05)
    
    
    cat("Power at block_ss =", block_ss, "is", power, "\n")
    
    return(power - target_power)
  }
  
  
  result <- tryCatch({
    uniroot(power_diff, lower = starting_point, upper = ending_point, tol = 0.01, trace = 1)$root
  }, error = function(e) {
    warning("Could not find root within given bounds.")
    return(NA)
  })
  
  return(ceiling(result))
}




### function to find the minimum sample size using optim given readers per block, number of blocks, and other parameters
find_min_sample_size_optim <- function(readers_per_block = 3, 
                                       blocks = 2,
                                       sigma_r = 0.011, sigma_tr = 0.03,
                                       n_sim = 200, target_power = 0.80) {
  
  total_readers<- readers_per_block * blocks
  
  
  starting_point <- tryCatch({
    sampleSize_MRMC(endpoint = "AUC",
                    J = total_readers,
                    delta = 0.05,
                    rangeb = 0.1,
                    rangew = 0.05,
                    theta = 0.85,
                    r1 = 0.47)$ORSampleSizeResults$nTotal
  }, error = function(e) {
    1000
    #should use largest possible sample size from sampleSize_MRMC here
  })
  
  print(paste("Starting point for sample size search:", starting_point))
  
  
  ending_point <- tryCatch({
    sampleSize_MRMC(endpoint = "AUC",
                    J = readers_per_block,
                    delta = 0.05,
                    rangeb = 0.1,
                    rangew = 0.05,
                    theta = 0.85,
                    r1 = 0.47)$ORSampleSizeResults$nTotal
  }, error = function(e) {
    starting_point+500
  })

  ending_point <- ifelse(ending_point == starting_point, 
         starting_point + 100, 
         ending_point)
  
  print(paste("Ending point for sample size search:", ending_point))
  
  
  power_diff <- function(block_ss) {
    block_ss <- round(block_ss)  # ensure it's an integer
    
    cat("Trying block sample size:", block_ss, "\n")
    
    data <- sim_splitplot_cardy(n = n_sim,
                                readers_per_block = readers_per_block,
                                block_ss = block_ss,
                                sigma_r = sigma_r,
                                sigma_tr = sigma_tr,
                                hypothesis = "alt",
                                blocks = blocks)
    
    # Parallel processing
    parms_list1 <- mclapply(data, ustat_analysis_forsim, mc.cores = detectCores() - 1)
    pvals <- sapply(parms_list1, function(x) x$pValueBDG)
    power <- mean(pvals < 0.05)
    
    
    cat("Power at block_ss =", block_ss, "is", power, "\n")
    
    return(power - target_power)
  }
  
  
  result <- tryCatch({
    optim(fn = power_diff, 
          par = (starting_point + ending_point) / 2, 
          method = "Brent", 
          lower = starting_point, 
          upper = ending_point, 
          control = list(trace = 1, abstol = 0.02, reltol=1))
    
  }, error = function(e) {
    warning("Could not find root within given bounds.")
    return(NA)
  })
  
  return(ceiling(result))
}




