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
    corr2 = 0.35, 
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

################
### Main Sim ###
################

#args <- as.integer(commandArgs(trailing = TRUE))
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


test_data_for_gui <- sim_one_splitplot_cardy(mu_nondisease = 0, 
                                             mu_disease = rmh$delta1, 
                                             hypothesis = "alt",
                                             tau = rmh$delta2-rmh$delta1,
                                             sigma_r = rmh$sigma_r, sigma_c = rmh$sigma_C,
                                             sigma_rc = rmh$sigma_RC, sigma_tr = rmh$sigma_tr,
                                             sigma_tc = rmh$sigma_TC, sigma_trc = rmh$sigma_trc, 
                                             block_ss = 1000,
                                             readers_per_block = 200,
                                             blocks = 1)

test_data_for_gui_imrmc<-createIMRMCdf(test_data_for_gui,
              keyColumns = list(readerID = "reader", 
                                 caseID = "case", 
                                 modalityID = "modality",
                                 score = "score", 
                                 truth = "truth_binary"),
              truePositiveFactor = 1) %>%
  mutate(readerID = case_when(readerID == "truth" ~ "-1",
                              .default = readerID),
         modalityID = case_when(modalityID == "truth" ~ "0",
                                .default = modalityID))


write.csv(test_data_for_gui_imrmc, paste0("test_data_for_gui",args,".csv"), row.names = FALSE)
#must use a text editor to modify csv to get it in the right format for Java
#delete column names, remove quotation marks, add BEGIN DATA:, add header with NR: N0: N1: NM:

sampleSize_MRMC(endpoint = "AUC",
                J = selected$readers,
                theta = 0.825,
                delta = selected$delta,
                rangeb = var_table[[selected$observer_var]][2],
                rangew = var_table[[selected$observer_var]][1],
                r1 = 0.47,
                power = target_power,
                reader_var_estimation_method = 1)$ORSampleSizeResults$nTotal







