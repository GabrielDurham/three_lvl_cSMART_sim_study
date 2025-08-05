#################################
#################################
#################################
### ALT_SIMULATE_AND_FIT_DATA ###
#################################
#################################
#################################

#### PURPOSE: This file contains code that serves to support simulating and
####          fitting data given a sheet of a driver file


#### DATE CREATED:  10 MAY 2024
#### PROGRAMMER:    GABRIEL DURHAM (GJD)
#### EDITS:         16 MAY 2024 - Added functionality for alternate model fitting
####                              and static model fitting in particular based
####                              on function_v3.5
#### EDITS:         20 MAY 2024 - Added functionality to Execute_Driver_File_Main()
####                              so it only simulates rows with run_simulation==1
####                              (And made it save the driver)
####                08 JUL 2024 - Replaced call of Obtain_All_Sim_Params() to
####                              a Derive_Sim_Parms() call
####                              Also stored pre_r_mc_parms in output
####                18 SEP 2024 - Fixed static alt fit incorporation of covariates
####                10 FEB 2025 - Added line to retry the run if there's a singular matrix
####                              in the fitting (happens sometimes for complicated)
####                              WV models in small samples
####                24 FEB 2025 - Set seed prior to starting each driver row for
####                              reproducibility
####                25 FEB 2025 - Changed to driver row seed
####                26 FEB 2025 - Changed max_attempts to driver row input
####                02 JUN 2025 - Incorporated crit_t, t_max, t_0 as driver parameters



# Create Design Matrices for Marginal Mean Model of Prototypical SMART
## Partition into coefficients of 1, a_1, a_2nr and a_1*a_2nr
### driver_parms = Output of Process_Driver_Row_Main() function
### sim_data = Simulated SMART data - Output of Sim_SMART_Data() function
Create_Partial_Design_Mats_Prototypical <- function(driver_parms, sim_data) {
  N <- nrow(sim_data)
  # Create some time placeholder variables 
  crit_t <- driver_parms[["crit_t"]]
  t_diff <- sim_data$t-crit_t
  t_ind <- ifelse(sim_data$t>crit_t, 1, 0)
  
  gamma_0 <- integer(N)+1
  gamma_1 <- (1-t_ind)*sim_data$t + (t_ind*crit_t)
  #gamma_2 <- ((1-t_ind)*sim_data$t + (t_ind*crit_t))*sim_data$a_1
  gamma_2 <- ((1-t_ind)*sim_data$t + (t_ind*crit_t))
  gamma_3 <- t_ind*t_diff
  #gamma_4 <- t_ind*t_diff*sim_data$a_1
  gamma_4 <- t_ind*t_diff
  #gamma_5 <- t_ind*t_diff*sim_data$a_2
  gamma_5 <- t_ind*t_diff
  #gamma_6 <- t_ind*t_diff*sim_data$a_1*sim_data$a_2
  gamma_6 <- t_ind*t_diff
  covars <- as.matrix(sim_data[,driver_parms[["covars"]]])
  
  Output <- NULL
  Output[["ones"]] <- as.matrix(cbind(gamma_0, gamma_1, gamma_3, covars))
  Output[["a_1"]] <- as.matrix(cbind(gamma_2, gamma_4))
  Output[["a_2nr"]] <- as.matrix(gamma_5, ncol=1)
  colnames(Output[["a_2nr"]]) <- c("gamma_5")
  Output[["a_1*a_2nr"]] <- as.matrix(gamma_6, ncol=1)
  colnames(Output[["a_1*a_2nr"]]) <- c("gamma_6")
  
  return(Output)
}





# Create Partial Design Matrices for Marginal Mean Model SMART
## Partition into coefficients of main decisions
### driver_parms = Output of Process_Driver_Row_Main() function
### sim_data = Simulated SMART data - Output of Sim_SMART_Data() function
Create_Partial_Design_Mats <- function(driver_parms, sim_data) {
  if (driver_parms[["SMART_structure"]]=="prototypical") {
    Output <- Create_Partial_Design_Mats_Prototypical(driver_parms=driver_parms, 
                                                      sim_data=sim_data)
  }
  return(Output)
}





# Create Indicator Matrix for Prototypical SMART
### driver_parms = Output of Process_Driver_Row_Main() function
Create_Ind_Mat <- function(driver_parms) {
  SMART_structure_df <- Define_SMART_Pathways(SMART_structure=driver_parms[["SMART_structure"]])
  if (driver_parms[["SMART_structure"]]=="prototypical") {
    Output <- data.frame(
      "ones"=integer(4)+1,
      "a_1"=SMART_structure_df$a_1,
      "a_2nr"=SMART_structure_df$a_2nr,
      "a_1*a_2nr"=SMART_structure_df$a_1*SMART_structure_df$a_2nr
    )
  }
  rownames(Output) <- SMART_structure_df$dtr
  return(Output)
}



# Create Matrix Indicating Consistency with DTR
### driver_parms = Output of Process_Driver_Row_Main() function
### sim_data = Simulated SMART data - Output of Sim_SMART_Data() function
Create_Consistency_Mat <- function(driver_parms, sim_data) {
  if (driver_parms[["SMART_structure"]]=="prototypical") {
    Output <- Create_Consistency_Mat_Prototypical(sim_data=sim_data)
  }
  return(Output)
}


# Create a Consistency Matrix for a Prototypical SMART
### sim_data = Simulated SMART data - Output of Sim_SMART_Data() function
Create_Consistency_Mat_Prototypical <- function(sim_data) {
  temp_output <- data.frame(d_1_1=numeric(), d_1_m1=numeric(),
                            d_m1_1=numeric(), d_m1_m1=numeric())
  clusters <- unique(sim_data$cluster_id)
  for (id in clusters) {
    cluster_data <- sim_data[sim_data$cluster_id==id,]
    a_1 <- cluster_data$a_1[1]
    a_2 <- ifelse(is.na(cluster_data$a_2[1]), 0, cluster_data$a_2[1])
    r <- cluster_data$r[1]
    
    # Code Consistent DTRs
    c_d_1_1 <- ifelse(((a_1==1)&(r==1))|((a_1==1)&(a_2==1)), 1, 0)
    c_d_1_m1 <- ifelse(((a_1==1)&(r==1))|((a_1==1)&(a_2==-1)), 1, 0)
    c_d_m1_1 <- ifelse(((a_1==-1)&(r==1))|((a_1==-1)&(a_2==1)), 1, 0)
    c_d_m1_m1 <- ifelse(((a_1==-1)&(r==1))|((a_1==-1)&(a_2==-1)), 1, 0)
    
    temp_output <- rbind(temp_output, 
                         data.frame(d_1_1=c_d_1_1,
                                    d_1_m1=c_d_1_m1,
                                    d_m1_1=c_d_m1_1,
                                    d_m1_m1=c_d_m1_m1))
  }
  Output <- temp_output
  return(Output)
}


# Obtain Model Fitting Parameters
### driver_parms = Output of Process_Driver_Row_Main() function
### fit_number = The fit number in the driver row to use
### fit_settings_df = Dataframe of driver sheet of model fit settings
Obtain_Fit_Parms <- function(driver_parms, fit_number, fit_settings_df) {
  fit_setting <- driver_parms[[paste0("fit_setting_", fit_number)]]
  fit_parm_row <- fit_settings_df[fit_settings_df$fit_setting==fit_setting,]
  
  Output <- NULL
  for (parm in colnames(fit_parm_row)) {
    Output[[parm]] <- fit_parm_row[[parm]][1]
  }
  return(Output)
}


# Derive Arguments for Model Fitting Function
### driver_parms = Output of Process_Driver_Row_Main() function
### fit_number = The fit number in the driver row to use
### fit_settings_df = Dataframe of driver sheet of model fit settings
### sim_data = Simulated SMART data - Output of Sim_SMART_Data() function
Pull_Function_Args <- function(driver_parms, fit_number, fit_settings_df,
                               sim_data) {
  fit_parms <- Obtain_Fit_Parms(driver_parms=driver_parms, 
                                fit_number=fit_number, 
                                fit_settings_df=fit_settings_df)
  # Restrict data if we request it
  if ("t_to_model" %in% names(fit_parms)) {
    if (fit_parms[["t_to_model"]]=="saturated") {
      t_to_model <- 
        c(driver_parms[["t_0"]], driver_parms[["crit_t"]], driver_parms[["t_max"]])
      sim_data <- sim_data[sim_data$t %in% t_to_model, ]
    }
  }

  Output <- NULL
  Output[["N"]] <- length(unique(sim_data$cluster_id))
  Output[["M"]] <- as.vector(table(sim_data[sim_data$t==0, "cluster_id"]))
  Output[["max_T"]] <- length(unique(sim_data$t))
  Output[["Ind"]] <- Create_Ind_Mat(driver_parms=driver_parms)
  Output[["D"]] <- nrow(Output[["Ind"]])
  Output[["Y"]] <- sim_data$y
  Output[["X"]] <- Create_Partial_Design_Mats(driver_parms=driver_parms, 
                                              sim_data=sim_data)
  Output[["within"]] <- Create_Consistency_Mat(driver_parms=driver_parms, 
                                               sim_data=sim_data)
  Output[["weight"]] <- 1/rowSums(Output[["within"]])
  Output[["var_homo_across_time"]] <- fit_parms[["var_homo_over_time"]]==1
  Output[["var_homo_across_AI"]] <- fit_parms[["var_homo_over_AI"]]==1
  Output[["diagonal_structure"]] <- fit_parms[["correlation_structure_wi_unit"]]
  Output[["diagonal_homo_across_AI"]] <- fit_parms[["cov_homo_over_AI_wi_unit"]]
  Output[["diagonal_ICC_lower_thresh"]] <- fit_parms[["ICC_lower_thresh_wi_unit"]]
  Output[["off_diagonal_structure"]] <- fit_parms[["correlation_structure_ic"]]
  Output[["off_diagonal_homo_across_AI"]] <- fit_parms[["cov_homo_over_AI_ic"]]
  Output[["off_diagonal_ICC_lower_thresh"]] <- fit_parms[["ICC_lower_thresh_ic"]]
  Output[["max_iter"]] <- fit_parms[["max_iter"]]
  Output[["dof_adjustment"]] <- fit_parms[["dof_adjustment"]]
  Output[["use_t"]] <- fit_parms[["use_t"]]
  return(Output)
  
  
}




# Reformat Summary Parameters
### model_fit_output = Output of solve_SMART_Multilayer() function
# Reformat_Output <- function(model_fit_output) {
#   temp_output <- model_fit_output[["summary_paras"]]
#   for (row in rownames(temp_output)) {
#     broken_string <- strsplit(temp_output[row, "Parameter"], split=" * ")
#     temp_output[row, "Parameter"] <- broken_string[[1]][length(broken_string[[1]])]
#   }
#   Output <- temp_output[order(temp_output$Parameter),]
#   rownames(Output) <- 1:nrow(Output)
#   return(Output)
# }

# Reformat Summary Parameters and Variance Estimator Labels
### model_fit_output = Output of solve_SMART_Multilayer() function
Reformat_Output <- function(model_fit_output) {
  temp_output_paras <- model_fit_output[["summary_paras"]]
  temp_output_var <- model_fit_output[["var_estimator"]]
  for (row in rownames(temp_output_paras)) {
    broken_string <- strsplit(temp_output_paras[row, "Parameter"], split=" * ")
    temp_output_paras[row, "Parameter"] <- broken_string[[1]][length(broken_string[[1]])]
  }
  for (i in (1:nrow(temp_output_var))) {
    rowname <- rownames(temp_output_var)[i]
    colname <- colnames(temp_output_var)[i]
    broken_string_row <- strsplit(rowname, split=" * ")
    rownames(temp_output_var)[i] <-
      broken_string_row[[1]][length(broken_string_row[[1]])]
    broken_string_col <- strsplit(colname, split=" * ")
    colnames(temp_output_var)[i] <-
      broken_string_col[[1]][length(broken_string_col[[1]])]
  }
  params_in_new_order <- order(temp_output_paras$Parameter)
  Output <- list(summary_paras=temp_output_paras[params_in_new_order,],
                 var_estimator=temp_output_var[params_in_new_order, params_in_new_order])
  rownames(Output[["summary_paras"]]) <- 1:nrow(Output[["summary_paras"]])
  return(Output)
}

# Fit a Model Given Data
### driver_parms = Output of Process_Driver_Row_Main() function
### fit_number = The fit number in the driver row to use
### fit_settings_df = Dataframe of driver sheet of model fit settings
### sim_data = Simulated SMART data - Output of Sim_SMART_Data() function
Fit_Model <- function(driver_parms, fit_number, fit_settings_df, sim_data) {
  fit_args <- Pull_Function_Args(driver_parms=driver_parms, 
                                 fit_number=fit_number, 
                                 fit_settings_df=fit_settings_df,
                                 sim_data=sim_data)
  Output <- solve_SMART_Multilayer(N=fit_args[["N"]],
                                   M=fit_args[["M"]],
                                   max_T=fit_args[["max_T"]],
                                   D=fit_args[["D"]],
                                   Ind=fit_args[["Ind"]],
                                   Y=fit_args[["Y"]],
                                   X=fit_args[["X"]],
                                   within=fit_args[["within"]],
                                   weight=fit_args[["weight"]],
                                   var_homo_across_time=fit_args[["var_homo_across_time"]],
                                   var_homo_across_AI=fit_args[["var_homo_across_AI"]],
                                   diagonal_structure=fit_args[["diagonal_structure"]],
                                   diagonal_homo_across_AI=fit_args[["diagonal_homo_across_AI"]],
                                   diagonal_ICC_lower_thresh=fit_args[["diagonal_ICC_lower_thresh"]],
                                   off_diagonal_structure=fit_args[["off_diagonal_structure"]],
                                   off_diagonal_homo_across_AI=fit_args[["off_diagonal_homo_across_AI"]],
                                   off_diagonal_ICC_lower_thresh=fit_args[["off_diagonal_ICC_lower_thresh"]],
                                   max_iter=fit_args[["max_iter"]],
                                   dof_adjustment=fit_args[["dof_adjustment"]],
                                   use_t=fit_args[["use_t"]],
                                   verbose=0)
  reformatted_components <- Reformat_Output(model_fit_output=Output)
  Output[["summary_paras"]] <- reformatted_components[["summary_paras"]]
  Output[["var_estimator"]] <- reformatted_components[["var_estimator"]]
  return(Output)
}


# Run a Single Iteration of a Driver Row
### driver_parms = Output of Process_Driver_Row_Main() function
### dist_data = Output from Build_All_Distribution_Mats() function
### fit_settings_df = Dataframe of driver sheet of model fit settings
### alt_fit_settings_dfs = List of dataframe of driver sheets of alternate model 
###                        fit settings (indexed by alternate fit types)
Run_Single_Iteration_Alt <- function(driver_parms, dist_data, fit_settings_df,
                                     alt_fit_settings_dfs) {
  sim_data <- Sim_SMART_Data_Alt(driver_parms=driver_parms, 
                                 dist_data=dist_data)
  Output <- NULL
  if (driver_parms[["n_fit"]]>0) {
    for (fit_i in (1:driver_parms[["n_fit"]])) {
      out_label <- paste0("fit_", fit_i)
      Output[[out_label]] <- Fit_Model(driver_parms=driver_parms, 
                                       fit_number=fit_i, 
                                       fit_settings_df=fit_settings_df, 
                                       sim_data=sim_data)
      
    }
  }
  
  if (driver_parms[["n_alt_fit"]]>0) {
    for (alt_fit_i in (1:driver_parms[["n_alt_fit"]])) {
      out_label <- paste0("alt_fit_", alt_fit_i)
      if (driver_parms[[paste0("alt_fit_type_", alt_fit_i)]]=="static") {
        Output[[out_label]] <- 
          Fit_Alt_Model_Static(driver_parms=driver_parms, 
                               alt_fit_number=alt_fit_i, 
                               static_alt_fit_settings_df=alt_fit_settings_dfs[["static"]], 
                               sim_data=sim_data)
      }
    }
  }
  
  return(Output)
}


# Run a Multiple Iterations of a Driver Row
### driver_parms = Output of Process_Driver_Row_Main() function
### dist_data = Output from Build_All_Distribution_Mats() function
### fit_settings_df = Dataframe of driver sheet of model fit settings
### alt_fit_settings_dfs = List of dataframe of driver sheets of alternate model 
###                        fit settings (indexed by alternate fit types)
### iter_list = List of iteration numbers (to serve as output labels)
Run_Multiple_Iterations_Alt <- function(driver_parms, dist_data, fit_settings_df, 
                                        alt_fit_settings_dfs, iter_list) {
  Output <- NULL
  max_attempts <- driver_parms[["max_attempts"]]
  for (iter in (1:length(iter_list))) {
    iter_output <- NULL
    # Allow multiple runs before a fatal error arises
    for (attempt in (1:max_attempts)) {
      tryCatch({
        iter_output <- 
          Run_Single_Iteration_Alt(driver_parms=driver_parms, 
                                   dist_data=dist_data, 
                                   fit_settings_df=fit_settings_df,
                                   alt_fit_settings_dfs=alt_fit_settings_dfs)
        break  # Exit the loop if successful
      }, error = function(e) {
        if (attempt == max_attempts) {
          stop("All attempts to run Run_Single_Iteration_Alt() have failed.")
        }
      })
    }
    for (storage_label in names(iter_output)) {
      Output[[storage_label]][[iter]] <- iter_output[[storage_label]]
    }
  }
  return(Output)
}


# Execute a Row of a Driver File
### alt_driver = Output of Read_Driver_File() function
### driver_rowname = Rowname of driver file to execute
### pre_r_mc_parms = Output of Import_Pre_R_MC_Results() function
### do_par = Boolean indicating whether to parallelize
### n_threads = Number of parallel threads to run, used if do_par==TRUE
###             Note: Parallelization done using doRNG package, ensuring
###                   replicability with respect to random seed
Execute_Driver_Row_Alt <- function(alt_driver, driver_rowname,
                                   do_par=FALSE, n_threads=1) {
  driver_parms <- 
    Process_Driver_Row_Main(driver_row=alt_driver[["driver"]][driver_rowname,])
  cp_settings_df <- alt_driver[["cp_settings"]]
  cp_settings <- 
    cp_settings_df[cp_settings_df$alt_cp_setting==driver_parms[["cond_parm_setting"]], ]
  dist_data <- 
    Build_All_Distribution_Mats(driver_parms=driver_parms, 
                                cp_settings=cp_settings)
  Output <- NULL
  Output[["driver_parms"]] <- driver_parms
  Output[["dist_data"]] <- dist_data
  # Set a seed prior to each row being run
  set.seed(driver_parms[["seed"]])
  #Parallelize
  if (n_threads==1) {do_par <- FALSE}
  if (do_par) {
    iter_partitions <- Partition_Vector(vec=(1:driver_parms[["n_iter"]]), k=n_threads)
    partitioned_output <- foreach(i=1:n_threads) %dorng% {
      Run_Multiple_Iterations_Alt(driver_parms=driver_parms, 
                                  dist_data=dist_data, 
                                  fit_settings_df=alt_driver[["fit_settings"]], 
                                  alt_fit_settings_dfs=alt_driver[["alt_fit_settings"]],
                                  iter_list=iter_partitions[[i]])
    }
    #Merge Output
    temp_output <- NULL
    for (name in names(partitioned_output[[1]])) {
      temp_output[[name]] <- partitioned_output[[1]][[name]]
      for (i in (2:n_threads)) {
        temp_output[[name]] <- c(temp_output[[name]],
                                 partitioned_output[[i]][[name]])
      }
    }
    Output <- temp_output
  } else {
    Output <- 
      Run_Multiple_Iterations_Alt(driver_parms=driver_parms, 
                                  dist_data=dist_data, 
                                  fit_settings_df=alt_driver[["fit_settings"]], 
                                  alt_fit_settings_dfs=alt_driver[["alt_fit_settings"]],
                                  iter_list=(1:driver_parms[["n_iter"]]))
  }
  
  return(Output)
}


# Execute a Main Sheet of a Driver File - Alternate
### alt_driver = Output of Read_Driver_File() function
### do_par = Boolean indicating whether to parallelize
### n_threads = Number of parallel threads to run, used if do_par==TRUE
###             Note: Parallelization done using doRNG package, ensuring
###                   replicability with respect to random seed
### path = Path with past output stored
### folder_name = Folder name to store output
Execute_Driver_File_Main_Alt <- function(alt_driver, do_par=FALSE, n_threads=1,
                                         path, folder_name) {
  # Create folder to store output
  dir.create(file.path(path, folder_name))
  for (row in rownames(alt_driver[["driver"]])) {
    if (alt_driver[["driver"]][row, "run_simulation"]==1) {
      row_output <- Execute_Driver_Row_Alt(alt_driver=alt_driver,
                                           driver_rowname=row,
                                           do_par=do_par,
                                           n_threads=n_threads)
      saveRDS(object=row_output, 
              file=paste0(path, folder_name, 
                          paste0("/sim_", alt_driver[["driver"]][row,"sim_label"])))
    }
  }
  saveRDS(object=alt_driver, 
          file=paste0(path, folder_name, "/Driver"))
}




#### Alternate model fitting

# Fit an Alternate Model Given Data
### driver_parms = Output of Process_Driver_Row_Main() function
### alt_fit_number = The alternate fit number in the driver row to use
### static_alt_fit_settings_df = Dataframe of driver sheet of model fit settings
### sim_data = Simulated SMART data - Output of Sim_SMART_Data() function
Fit_Alt_Model_Static <- function(driver_parms, 
                                 alt_fit_number, 
                                 static_alt_fit_settings_df, 
                                 sim_data) {
  
  fit_args <- Pull_Function_Args_Alt_Static(driver_parms=driver_parms, 
                                            alt_fit_number=alt_fit_number, 
                                            static_alt_fit_settings_df=static_alt_fit_settings_df,
                                            sim_data=sim_data)
  
  Output <- solve_SMART(Y=fit_args[["Y"]],
                        X=fit_args[["X"]],
                        cluster_id=fit_args[["id"]],
                        A1=fit_args[["A1"]],
                        R=fit_args[["R"]],
                        A2=fit_args[["A2"]],
                        aimed_comparison=matrix(c(1,1,1,0,0,0,1,0,0,1,1,0,0,1,0,1,0,1,0,0,1,0,1,1),6,4),
                        variance_structure=fit_args[["variance_structure"]],
                        correlation_structure=fit_args[["correlation_structure"]],
                        ICC_lower_thresh=fit_args[["ICC_lower_thresh"]],
                        max_iter=fit_args[["max_iter"]],
                        convergence_thresh=fit_args[["convergence_thresh"]],
                        alpha=fit_args[["alpha"]],
                        estimate_weight=fit_args[["estimate_weight"]],
                        dof_adjustment=fit_args[["dof_adjustment"]],
                        use_t=fit_args[["use_t"]],
                        bias_correction=fit_args[["bias_correction"]],
                        verbose=0)
  #Output[["summary_paras"]] <- Reformat_Output(model_fit_output=Output)
  return(Output)
}


# Derive Arguments for Static Model Fitting Function
### driver_parms = Output of Process_Driver_Row_Main() function
### alt_fit_number = The fit number in the driver row to use
### static_alt_fit_settings_df = Dataframe of driver sheet of model fit settings
### sim_data = Simulated SMART data - Output of Sim_SMART_Data() function
Pull_Function_Args_Alt_Static <- function(driver_parms, alt_fit_number, 
                                          static_alt_fit_settings_df, sim_data) {
  alt_fit_setting <- driver_parms[[paste0("alt_fit_setting_", alt_fit_number)]]
  fit_settings <- 
    static_alt_fit_settings_df[static_alt_fit_settings_df$alt_fit_setting==alt_fit_setting,]
  Output <- NULL
  static_data <- Collapse_Sim_Data(sim_data=sim_data, 
                                   outcome_type=fit_settings$outcome_type)
  
  Output[["Y"]] <- static_data$y
  if (driver_parms[["n_covar"]]>0) {
    Output[["X"]] <- static_data[,driver_parms[["covars"]]]
  } else {Output[["X"]] <- NULL}
  Output[["id"]] <- static_data$cluster_id
  Output[["A1"]] <- static_data$a_1
  Output[["R"]] <- static_data$r
  Output[["A2"]] <- static_data$a_2
  Output[["variance_structure"]] <- 
    ifelse(fit_settings$var_homo_over_AI==1, 0, 1)
  # Can only accept "Independent" or "Exchangeable"
  Output[["correlation_structure"]] <- 
    ifelse(fit_settings$correlation_structure=="Independent", 0,
           ifelse(fit_settings$cov_homo_over_AI==1, 1, 2))
  # Pull in numeric tuning parameters  
  for (var in c("ICC_lower_thresh", "max_iter", "convergence_thresh", "alpha")) {
    Output[[var]] <- fit_settings[[var]]
  }
  # Pull in Boolean tuning parameters
  for (var in c("estimate_weight", "dof_adjustment", "use_t", "bias_correction")) {
    Output[[var]] <- fit_settings[[var]]==1
  }
  
  return(Output)
}


# Collapse Simulated Data for Alternate Fits
### sim_data = Simulated SMART data - Output of Sim_SMART_Data() function
### outcome_type = Outcome construction to analyze
Collapse_Sim_Data <- function(sim_data, outcome_type) {
  if (outcome_type=="EOS") {
    temp_output <- sim_data[sim_data$t==max(sim_data$t),]
  } else if (outcome_type=="Total") {
    temp_outcome <- aggregate(x=list(y=sim_data[,"y"]),
                              by=list(person_id=sim_data$person_id,
                                      cluster_id=sim_data$cluster_id),
                              FUN=sum)
    temp_output <- sim_data[sim_data$t==max(sim_data$t),!(names(sim_data)=="y")]
    temp_output <- merge(temp_output, temp_outcome, by=c("cluster_id", "person_id"))
  }
  
  Output <- temp_output[order(temp_output$cluster_id, temp_output$person_id),
                        !(names(sim_data)=="t")]
  
  return(Output)
}
