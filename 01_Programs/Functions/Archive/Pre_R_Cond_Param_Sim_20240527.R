#####################################################
#####################################################
#####################################################
### PRE-RESPONSE CONDITIONAL PARAMETER SIMULATION ###
#####################################################
#####################################################
#####################################################

#### PURPOSE: This file contains code to run pre-response (t=0,1) simulations
####          and simulate response given driver file input. 
####          Additionally, this file contains code which analyzes the simulated
####          data in order to derive response-conditional error variance and
####          mean parameters


#### DATE CREATED: 23 OCT 2023
#### PROGRAMMER: GABRIEL DURHAM (GJD)
#### EDITS: 06 NOV 2023 (GJD) -  Changed sigma_02/sigma_12 input column names in
####                            Process_Driver_Row_Pre_R_MC read-in
####        09 JAN 2024 (GJD) - Changed back to sigma2_0/sigma2_1
####        24 MAY 2024 (GJD) - Added stored_run variable to driver and changed
####                  Execute_Raw_Driver_File_Pre_R() function to only run on 
####                  settings without stored run
####                  Also changed output to use cond_param_setting_pre_r instead
####                  of sim_id


# Expand Driver File
## Some driver parameters except lists (e.g., sample size),
## expand the file so there's a row for each unique simulation setting
### raw_driver_row = Row of raw driver file
#### Returns expanded driver file with one row for each simulation setting
Expand_Driver_File <- function(raw_driver_file){
  #Columns that can contain lists
  list_columns <- c("sample_sizes")
  static_data <- raw_driver_file[, !(names(raw_driver_file)%in%list_columns)]
  list_data <- raw_driver_file[, c("cond_param_setting_pre_r", list_columns)]
  
  #Need to add any other list variables manually if we update the driver
  temp_expanded_data <- data.frame(cond_param_setting_pre_r = numeric(),
                                   sample_size = numeric())
  
  
  for (row in rownames(raw_driver_file)) {
    #Update if we add any list variables
    temp_expanded_data_row <- 
      expand.grid(sample_size = eval(parse(text=raw_driver_file[row,"sample_sizes"])))
    temp_expanded_data_row$cond_param_setting_pre_r <- 
      as.numeric(raw_driver_file[row,"cond_param_setting_pre_r"][1])
    temp_expanded_data <- rbind(temp_expanded_data, 
                                temp_expanded_data_row)
  }
  
  Output <- merge(static_data, temp_expanded_data, by="cond_param_setting_pre_r", all=TRUE)
  
  return(Output)
}




# Process a Driver Row - 
### driver_row = Row of driver file to process
#### Returns a list with simulation parameters
Process_Driver_Row_Pre_R_MC <- function(driver_row){
  Output <- NULL
  Output[["n_iter"]] <- driver_row$n_iter[1]
  Output[["r_model"]] <- driver_row$r_model[1]
  Output[["p_r"]] <- driver_row$p_r[1]
  Output[["sigma2_0"]] <- driver_row$sigma2_0[1]
  Output[["sigma2_1"]] <- driver_row$sigma2_1[1]
  Output[["rho_0"]] <- driver_row$rho_0[1]
  Output[["rho_1"]] <- driver_row$rho_1[1]
  Output[["rho_01"]] <- driver_row$rho_01[1]
  Output[["phi_01"]] <- driver_row$phi_01[1]
  Output[["n_i"]] <- driver_row$sample_size[1]
  Output[["stored_run"]] <- driver_row$stored_run[1]
  return(Output)
}






# Construct Relevant Parameters for Pre-Response Outcome Simulation
## Constructs parameters for residual simulation and longitudinal trajectory building
### parms = Output of Process_Driver_Row_Pre_R_MC() function
#### Returns list with: Sigma_0, Sigma_1, P_1
Construct_Pre_Resp_Helper_Parms <- function(parms){
  
  P_1 <- ifelse(parms[["rho_0"]]==0, 0,
                parms[["rho_01"]]/parms[["rho_0"]])
  V1_ni <- as.matrix(integer(parms[["n_i"]])+1, ncol=1)
  
  Sigma_0 <- (parms[["sigma2_0"]]-parms[["rho_0"]])*diag(parms[["n_i"]]) + 
    parms[["rho_0"]]*V1_ni%*%t(V1_ni)
  
  Sigma_1 <- (parms[["sigma2_1"]]-parms[["rho_1"]])*diag(parms[["n_i"]]) + 
    parms[["rho_1"]]*V1_ni%*%t(V1_ni) - (P_1^2)*Sigma_0 -
    2*P_1*(parms[["phi_01"]]-(P_1*parms[["sigma2_0"]]))*diag(parms[["n_i"]])
  
  Output <- NULL
  Output[["P_1"]] <- P_1
  Output[["Sigma_0"]] <- Sigma_0
  Output[["Sigma_1"]] <- Sigma_1
  
  return(Output)
}








# Build Pre-Response Error Variance Matrix
### parms = Output of Process_Driver_Row_Pre_R_MC() function
### pre_r_helper_parms = Output of Construct_Pre_Resp_Helper_Parms() function
Build_Pre_Resp_Error_Var_Matrix <- function(parms, pre_r_helper_parms){
  P_1 <- pre_r_helper_parms[["P_1"]]
  
  off_diag <- (parms[["phi_01"]] - P_1*parms[["sigma2_0"]])*diag(parms[["n_i"]])
  
  top_block <- cbind(pre_r_helper_parms[["Sigma_0"]], off_diag)
  bottom_block <- cbind(off_diag, pre_r_helper_parms[["Sigma_1"]])
  
  Output <- rbind(top_block, bottom_block)
  return(Output)
}



# Generate Model Errors
### error_variance = Output of Build_Pre_Resp_Error_Var_Matrix() function
Generate_Model_Errors_Pre_R <- function(error_variance){
  errors <- mvrnorm(n = 1, mu=integer(nrow(error_variance)), Sigma=error_variance)
  eps_0 <- errors[(1:(nrow(error_variance)/2))]
  eps_1 <- errors[(((nrow(error_variance)/2)+1):nrow(error_variance))]
  Output <- NULL
  Output[["e_0"]] <- eps_0
  Output[["e_1"]] <- eps_1
  
  return(Output)
}





# Generate Outcome Deviations Given Errors
## We want to calculate the deviation from the marginal model
### model_errors = Output from Generate_Model_Errors_Pre_R() function
### pre_r_helper_parms = Output of Construct_Pre_Resp_Helper_Parms() function
Generate_Outcome_Devs_Pre_R <- function(model_errors, pre_r_helper_parms){
  eps_0 <- model_errors[["e_0"]]
  eps_1 <- model_errors[["e_1"]]
  P_1 <- pre_r_helper_parms[["P_1"]]
  
  
  Y_0 <- eps_0
  Y_1 <- P_1*eps_0 + eps_1
  
  Output <- NULL
  Output[["Y_0"]] <- Y_0
  Output[["Y_1"]] <- Y_1
  return(Output)
}



################################
### SHELL RESPONSE FUNCTIONS ###
################################

# These functions are scaled-down response-generation functions
# just for conditional parameter MC simulation. I.e., they assume
# outcomes are already centered so they don't waste time with
# centering - still need to scale

# Response Function 1
## Response based on inverse CDF transform using Beta(p/1-p, 1) distribution
### sim_outcomes_devs = Output of Generate_Outcome_Devs_Pre_R() function
### parms = Output of Process_Driver_Row_Pre_R_MC() function - Will use to calculate scaling factor
### scaling_factor = Option of putting in scaling factor manually to avoid repeated calculation
###                  (Will divide by scaling factor)
Response_Function_1_Shell <- function(sim_outcomes_devs, parms, scaling_factor=NULL){
  # Pull sample size and marginal probability of response
  n_i <- parms[["n_i"]]
  marg_p_r <- parms[["p_r"]]
  # Calculate scaling factor if not provided
  if (is.null(scaling_factor)) {
    scaling_factor <- sqrt((1/n_i)*(parms[["sigma2_1"]] + (n_i-1)*parms[["rho_1"]]))
  }
  beta_input <- pnorm(mean(sim_outcomes_devs[["Y_1"]])/scaling_factor)
  
  p_response <- qbeta(beta_input, shape1=(marg_p_r/(1-marg_p_r)), shape2=1, 
                      ncp = 0, log = FALSE)
  
  Output <- rbinom(n=1, size=1, prob=p_response)
  return(Output)
}



# Save Observations
## Saves observations to allow for component calculations down the line
### sim_outcomes_devs = Output of Generate_Outcome_Devs_Pre_R() function
### saved_outcomes = List with outcomes saved from past simulations
###                  If NULL, will create a new one
Save_Observations <- function(sim_outcomes_devs, saved_outcomes=NULL){
  if (is.null(saved_outcomes)) {
    temp_output <- NULL
    for (var in names(sim_outcomes_devs)) {
      temp_output[[var]] <- NULL
      for (i in (1:length(sim_outcomes_devs[[var]]))) {
        temp_output[[var]][[toString(i)]] <- c(sim_outcomes_devs[[var]][i])
      }
    }
  } else{
    temp_output <- saved_outcomes
    for (var in names(sim_outcomes_devs)) {
      for (i in (1:length(sim_outcomes_devs[[var]]))) {
        temp_output[[var]][[toString(i)]] <- c(temp_output[[var]][[toString(i)]],
                                               sim_outcomes_devs[[var]][i])
      }
    }
  }
  Output <- temp_output
  return(Output)
}



# Process Saved Output - Pre Response
## Processes stored output and returns empirical mean/(co)variance values
## of pre-response outcomes
### saved_obs = Output from Save_Observations() function
Collapse_Saved_Obs_Pre_R <- function(saved_obs){
  temp_mean_vals <- NULL
  for (var in names(saved_obs)) {
    temp_mean_vals[[var]] <- c()
    for (i in names(saved_obs[[var]])) {
      temp_mean_vals[[var]] <- c(temp_mean_vals[[var]],
                                 mean(saved_obs[[var]][[i]]))
    }
  }
  temp_varnc_vals <- NULL
  for (var in names(saved_obs)) {
    temp_varnc_vals[[var]] <- c()
    for (i in names(saved_obs[[var]])) {
      temp_varnc_vals[[var]] <- c(temp_varnc_vals[[var]],
                                  var(saved_obs[[var]][[i]]))
    }
  }
  temp_rho_t_vals <- NULL
  for (var in names(saved_obs)) {
    temp_rho_t_vals[[var]] <- c()
    for (i in (1:(length(saved_obs[[var]])-1))) {
      vec_i <- saved_obs[[var]][[toString(i)]]
      for (j in ((i+1):length(saved_obs[[var]]))) {
        vec_j <- saved_obs[[var]][[toString(j)]]
        temp_rho_t_vals[[var]] <- c(temp_rho_t_vals[[var]],
                                    cov(vec_i, vec_j))
      }
    }
  }
  
  temp_rho_c_vals <- NULL
  temp_phi_vals <- NULL
  for (k in (1:(length(saved_obs)-1))) {
    var_1 <- names(saved_obs)[k]
    for (m in ((k+1):length(saved_obs))) {
      var_2 <- names(saved_obs)[m]
      lab <- paste0(var_1," x ", var_2)
      temp_rho_c_vals[[lab]] <- c()
      temp_phi_vals[[lab]] <- c()
      for (i in (1:length(saved_obs[[var_1]]))) {
        vec_i <- saved_obs[[var_1]][[toString(i)]]
        for (j in (1:length(saved_obs[[var_2]]))) {
          vec_j <- saved_obs[[var_2]][[toString(j)]]
          if (j!=i) {
            temp_rho_c_vals[[lab]] <- c(temp_rho_c_vals[[lab]],
                                        cov(vec_i, vec_j)) 
          } else  {
            temp_phi_vals[[lab]] <- c(temp_phi_vals[[lab]],
                                      cov(vec_i, vec_j)) 
          }
        }
      }
    }
    
    
  }
  
  Output <- NULL
  Output[["mean"]] <- NULL
  for (var in names(temp_mean_vals)) {Output[["mean"]][[var]] <- mean(temp_mean_vals[[var]])}
  Output[["Sigma^2"]] <- NULL
  for (var in names(temp_varnc_vals)) {Output[["sigma^2"]][[var]] <- mean(temp_varnc_vals[[var]])}
  Output[["rho"]] <- NULL
  for (var in names(temp_rho_t_vals)) {Output[["rho"]][[var]] <- mean(temp_rho_t_vals[[var]])}
  for (var in names(temp_rho_c_vals)) {Output[["rho"]][[var]] <- mean(temp_rho_c_vals[[var]])}
  Output[["phi"]] <- NULL
  for (var in names(temp_phi_vals)) {Output[["phi"]][[var]] <- mean(temp_phi_vals[[var]])}
  Output[["n"]] <- length(saved_obs[[names(saved_obs)[1]]][[1]])
  return(Output)
}



# Simulate Outcomes Y0 and Y1 and Response
## Simulate a pre-response outcome trajectory and response
### parms = Output of Process_Driver_Row_Pre_R_MC() function 
### pre_r_helper_parms = Output of Construct_Pre_Resp_Helper_Parms() function
### error_variance = Output of Build_Pre_Resp_Error_Var_Matrix() function
### response_function = Function used to simulate response given pre-response
###                     Only argument can be centered outcome trajectories
Simulate_Outcomes_01_and_Response <- function(parms, pre_r_helper_parms, 
                                              error_variance, response_function){
  #Simulate model errors
  model_errors <- Generate_Model_Errors_Pre_R(error_variance)
  #Simulate outcome trajectories (without mean components)
  centered_outcomes <- Generate_Outcome_Devs_Pre_R(model_errors, pre_r_helper_parms)
  #Generate response
  response <- response_function(centered_outcomes)
  
  Output <- NULL
  Output[["Y"]] <- centered_outcomes
  Output[["R"]] <- response
  
  return(Output)
}




# Build Single-Input Response Function
## Processes parameters to tailor response function to exact simulation parameters
### parms = Output of Process_Driver_Row_Pre_R_MC() function
Build_Single_Input_R_Function <- function(parms){
  n_i <- parms[["n_i"]]
  if (parms[["r_model"]]==1) {
    scaling_factor <- sqrt((1/n_i)*(parms[["sigma2_1"]] + (n_i-1)*parms[["rho_1"]]))
    Output <- function(input){return(Response_Function_1_Shell(sim_outcomes_devs=input,
                                                               parms=parms,
                                                               scaling_factor=scaling_factor))}
  }
  return(Output)
}



# Simulate Iteratively and Store Output
## Runs repeated simulations and stores output
### parms = Output of Process_Driver_Row_Pre_R_MC() function 
### pre_r_helper_parms = Output of Construct_Pre_Resp_Helper_Parms() function
### error_variance = Output of Build_Pre_Resp_Error_Var_Matrix() function
### Response_function = Function used to simulate response given pre-response
###                     Only argument can be centered outcome trajectories
Simulate_Multi_Store_Output <- function(parms){
  helper_parms <- Construct_Pre_Resp_Helper_Parms(parms=parms)
  response_function <- Build_Single_Input_R_Function(parms=parms)
  error_variance <- Build_Pre_Resp_Error_Var_Matrix(parms=parms,
                                                    pre_r_helper_parms=helper_parms)
  stored_output_total <- NULL
  stored_output_r <- NULL
  stored_output_nr <- NULL
  for (iter in (1:parms[["n_iter"]])) {
    iter_results <- Simulate_Outcomes_01_and_Response(parms=parms, 
                                                      pre_r_helper_parms=helper_parms, 
                                                      error_variance=error_variance, 
                                                      response_function=response_function)
    
    stored_output_total <- Save_Observations(sim_outcomes_devs=iter_results[["Y"]], 
                                             saved_outcomes=stored_output_total)
    if (iter_results[["R"]]==1) {
      stored_output_r <- Save_Observations(sim_outcomes_devs=iter_results[["Y"]],
                                           saved_outcomes=stored_output_r)
    } else if (iter_results[["R"]]==0) {
      stored_output_nr <- Save_Observations(sim_outcomes_devs=iter_results[["Y"]],
                                            saved_outcomes=stored_output_nr)
    }
  }
  Output <- NULL
  Output[["Total"]] <- stored_output_total
  Output[["Responders"]] <- stored_output_r
  Output[["Nonresponders"]] <- stored_output_nr
  return(Output)
}


# Process Simulation Results - Pre-Response
## Consolidates the results from processing one row of a driver file
### sim_output = Output from Simulate_Multi_Store_Output() function
### sim_id = ID of simulation setting
### n_i = Cluster size
### n_iter = Number of iterations
### saved_results = Output from other driver rows
###                 If NULL, new dataframe created
#### Returns list of dataframes tracking results (Total, Responders, Nonresponders)
Process_Simulation_Results_Pre_R <- function(sim_output, sim_id, n_i, saved_results=NULL){
  sim_results <- NULL
  for (pop in names(sim_output)) {
    sim_results[[pop]] <- data.frame(sim_id = sim_id,
                                     n_i = n_i,
                                     mean_e_0 = sim_output[[pop]][["mean"]][["Y_0"]],
                                     mean_e_1 = sim_output[[pop]][["mean"]][["Y_1"]],
                                     sigma2_0 = sim_output[[pop]][["sigma^2"]][["Y_0"]],
                                     sigma2_1 = sim_output[[pop]][["sigma^2"]][["Y_1"]],
                                     rho_0 = sim_output[[pop]][["rho"]][["Y_0"]],
                                     rho_1 = sim_output[[pop]][["rho"]][["Y_1"]],
                                     rho_01 = sim_output[[pop]][["rho"]][["Y_0 x Y_1"]],
                                     phi_01 = sim_output[[pop]][["phi"]][["Y_0 x Y_1"]],
                                     n_obs = sim_output[[pop]][["n"]])
  }
  if (is.null(saved_results)) {
    Output <- sim_results
  } else {
    Output <- NULL
    for (pop in names(sim_output)) {
      Output[[pop]] <- rbind(saved_results[[pop]], sim_results[[pop]])
    }
  }
  return(Output)
}







# Export Driver Results
## Exports Driver Results to Excel
### driver = Driver file (allows for recording design parameters)
### driver_results = Output of Execute_Driver_File_Pre_R() 
### path = Where to save output
### file_name = What to name file
Export_Driver_Results <- function(driver, driver_results, path, file_name){
  #Check if file name contains .xlsx suffix, if not add it
  suffix <- substr(file_name, nchar(file_name)-4, nchar(file_name))
  if (suffix!=".xlsx") {file_name <- paste0(file_name, ".xlsx")}
  #Check if path has a / at the end, if not then add it
  if (substr(file_name, nchar(file_name), nchar(file_name))!="/") {path <- paste0(path, "/")}
  full_file_location <- paste0(path, file_name)
  
  #Create object to export
  export_obj <- NULL
  for (pop in names(driver_results)) {
    export_obj[[pop]] <- driver_results[[pop]]
  }
  export_obj[["driver"]] <- driver
  
  write.xlsx(export_obj,
             file=full_file_location,
             rowNames=FALSE)
}

# Execute Driver File Rows for Pre Response Simulation
## Executes a driver file rows from start to finish
### driver_rows = Dataframe containing expanded driver file rows to process
Execute_Driver_Rows_Pre_R <- function(driver_rows){
  driver <- driver_rows
  stored_sim_results <- NULL
  for (row in rownames(driver)) {
    driver_row <- driver[row,]
    parms <- Process_Driver_Row_Pre_R_MC(driver_row=driver_row)
    driver_row_observations <- Simulate_Multi_Store_Output(parms=parms)
    driver_row_results <- NULL
    for (pop in names(driver_row_observations)) {
      driver_row_results[[pop]] <- Collapse_Saved_Obs_Pre_R(saved_obs=driver_row_observations[[pop]])
    }
    stored_sim_results <- 
      Process_Simulation_Results_Pre_R(sim_output=driver_row_results, 
                                       cond_param_setting_pre_r=driver_row$cond_param_setting_pre_r, 
                                       n_i=parms[["n_i"]],
                                       saved_results=stored_sim_results)
  }
  Output <- stored_sim_results
  return(Output)
}


# Execute Raw Driver File for Pre Response Simulation
## Executes a driver file from start to finish
### raw_driver = Dataframe containing raw driver file
### do_par = Boolean indicating whether to parallelize
### n_threads = Number of parallel threads to run, used if do_par==TRUE
###             Note: Parallelization done using doRNG package, ensuring
###                   replicability with respect to random seed
Execute_Raw_Driver_File_Pre_R <- function(raw_driver, do_par=FALSE, n_threads=1){
  driver <- Expand_Driver_File(raw_driver)
  # Only run settings without stored data
  driver <- driver[is.na(driver$stored_run),]
  
  #Don't mess around with overhead if only one thread
  if (n_threads==1) {do_par <- FALSE}
  if (do_par) {
    row_partitions <- Partition_Vector(vec=rownames(driver),k=n_threads)
    partitioned_output <- foreach(i=1:n_threads) %dorng% {
      Execute_Driver_Rows_Pre_R(driver_rows=driver[row_partitions[[i]],])
    }
    #Merge output
    temp_output <- NULL
    for (pop in names(partitioned_output[[1]])) {
      temp_output[[pop]] <- partitioned_output[[1]][[pop]]
      for (i in (2:n_threads)) {
        temp_output[[pop]] <- rbind(temp_output[[pop]],
                                    partitioned_output[[i]][[pop]])
      }
    }
    Output <- temp_output
  } else {
    #If not parallelizing, just throw all rows into one execution call
    Output <- Execute_Driver_Rows_Pre_R(driver_rows=driver)
  }
  
  return(Output)
}

