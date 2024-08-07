#############################
#############################
#############################
### SIMULATE_AND_FIT_DATA ###
#############################
#############################
#############################

#### PURPOSE: This file contains code that serves to support simulating and
####          fitting data given a sheet of a driver file


#### DATE CREATED: 10 MAY 2024
#### PROGRAMMER: GABRIEL DURHAM (GJD)
#### EDITS: 







# Create Design Matrices for Marginal Mean Model of Prototypical SMART
## Partition into coefficients of 1, a_1, a_2nr and a_1*a_2nr
### driver_parms = Output of Process_Driver_Row_Main() function
### sim_data = Simulated SMART data - Output of Sim_SMART_Data() function
### crit_t = Critical t value (time of second decision). Default is 1
Create_Partial_Design_Mats_Prototypical <- function(driver_parms, sim_data, crit_t=1) {
  N <- nrow(sim_data)
  # Create some time placeholder variables 
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
### crit_t = Critical t value (time of second decision). Default is 1
Create_Partial_Design_Mats <- function(driver_parms, sim_data, crit_t=1) {
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
### crit_t = Critical t value (time of second decision). Default is 1
Pull_Function_Args <- function(driver_parms, fit_number, fit_settings_df,
                               sim_data, crit_t=1) {
  fit_parms <- Obtain_Fit_Parms(driver_parms=driver_parms, 
                                fit_number=fit_number, 
                                fit_settings_df=fit_settings_df)
  
  Output <- NULL
  Output[["N"]] <- length(unique(sim_data$cluster_id))
  Output[["M"]] <- as.vector(table(sim_data[sim_data$t==0, "cluster_id"]))
  Output[["max_T"]] <- 3
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
Reformat_Output <- function(model_fit_output) {
  temp_output <- model_fit_output[["summary_paras"]]
  for (row in rownames(temp_output)) {
    broken_string <- strsplit(temp_output[row, "Parameter"], split=" * ")
    temp_output[row, "Parameter"] <- broken_string[[1]][length(broken_string[[1]])]
  }
  Output <- temp_output[order(temp_output$Parameter),]
  rownames(Output) <- 1:nrow(Output)
  return(Output)
}

# Fit a Model Given Data
### driver_parms = Output of Process_Driver_Row_Main() function
### fit_number = The fit number in the driver row to use
### fit_settings_df = Dataframe of driver sheet of model fit settings
### sim_data = Simulated SMART data - Output of Sim_SMART_Data() function
### crit_t = Critical t value (time of second decision). Default is 1
Fit_Model <- function(driver_parms, fit_number, fit_settings_df, sim_data, crit_t=1) {
  fit_args <- Pull_Function_Args(driver_parms=driver_parms, 
                                 fit_number=fit_number, 
                                 fit_settings_df=fit_settings_df,
                                 sim_data=sim_data, 
                                 crit_t=crit_t)
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
  Output[["summary_paras"]] <- Reformat_Output(model_fit_output=Output)
  return(Output)
}


# Run a Single Iteration of a Driver Row
### driver_parms = Output of Process_Driver_Row_Main() function
### sim_parms = Output of Obtain_All_Sim_Params() function
### fit_settings_df = Dataframe of driver sheet of model fit settings
Run_Single_Iteration <- function(driver_parms, sim_parms, fit_settings_df) {
  sim_data <- Sim_SMART_Data(driver_parms=driver_parms,
                             sim_parms=sim_parms)
  Output <- NULL
  for (fit_i in (1:driver_parms[["n_fit"]])) {
    out_label <- paste0("fit_", fit_i)
    Output[[out_label]] <- Fit_Model(driver_parms=driver_parms, 
                                     fit_number=fit_i, 
                                     fit_settings_df=fit_settings_df, 
                                     sim_data=sim_data, 
                                     crit_t=crit_t)
    
  }
  if (driver_parms[["n_alt_fit"]]>0) {
    print(".")
  }
  return(Output)
}


# Run a Multiple Iterations of a Driver Row
### driver_parms = Output of Process_Driver_Row_Main() function
### sim_parms = Output of Obtain_All_Sim_Params() function
### fit_settings_df = Dataframe of driver sheet of model fit settings
### iter_list = List of iteration numbers (to serve as output labels)
Run_Multiple_Iterations <- function(driver_parms, sim_parms, fit_settings_df, iter_list) {
  Output <- NULL
  for (iter in (1:length(iter_list))) {
    iter_output <- Run_Single_Iteration(driver_parms=driver_parms, 
                                        sim_parms=sim_parms, 
                                        fit_settings_df=fit_settings_df)
    for (storage_label in names(iter_output)) {
      Output[[storage_label]][[iter]] <- iter_output[[storage_label]]
    }
  }
  return(Output)
}



# Execute a Row of a Driver File
### full_driver = Output of Read_Driver_File() function
### driver_rowname = Rowname of driver file to execute
### pre_r_mc_parms = Output of Import_Pre_R_MC_Results() function
### do_par = Boolean indicating whether to parallelize
### n_threads = Number of parallel threads to run, used if do_par==TRUE
###             Note: Parallelization done using doRNG package, ensuring
###                   replicability with respect to random seed
Execute_Driver_Row <- function(full_driver, driver_rowname, pre_r_mc_parms,
                               do_par=FALSE, n_threads=1) {
  driver_parms <- Process_Driver_Row_Main(driver_row=full_driver[["driver"]][driver_rowname,])
  sim_parms <- Obtain_All_Sim_Params(driver_parms=driver_parms,
                                     cp_settings_df=full_driver[["cp_settings"]],
                                     pre_r_marg_parms=full_driver[["var_parms_pre_r"]],
                                     pre_r_cond_var_parms=pre_r_mc_parms,
                                     post_r_var_parms=full_driver[["var_parms_post_r"]])
  Output <- NULL
  Output[["driver_parms"]] <- driver_parms
  Output[["sim_parms"]] <- sim_parms
  
  #Parallelize
  if (n_threads==1) {do_par <- FALSE}
  if (do_par) {
    iter_partitions <- Partition_Vector(vec=(1:driver_parms[["n_iter"]]), k=n_threads)
    partitioned_output <- foreach(i=1:n_threads) %dorng% {
      Run_Multiple_Iterations(driver_parms=driver_parms, 
                              sim_parms=sim_parms, 
                              fit_settings_df=full_driver[["fit_settings"]], 
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
    Output <- Run_Multiple_Iterations(driver_parms=driver_parms, 
                                      sim_parms=sim_parms, 
                                      fit_settings_df=full_driver[["fit_settings"]], 
                                      iter_list=(1:driver_parms[["n_iter"]]))
  }
  
  return(Output)
}


# Execute a Main Sheet of a Driver File
### full_driver = Output of Read_Driver_File() function
### pre_r_mc_parms = Output of Import_Pre_R_MC_Results() function
### do_par = Boolean indicating whether to parallelize
### n_threads = Number of parallel threads to run, used if do_par==TRUE
###             Note: Parallelization done using doRNG package, ensuring
###                   replicability with respect to random seed
### path = Path with past output stored
### folder_name = Folder name to store output
Execute_Driver_File_Main <- function(full_driver, pre_r_mc_parms, 
                                     do_par=FALSE, n_threads=1,
                                     path, folder_name) {
  # Create folder to store output
  dir.create(file.path(path, folder_name))
  for (row in rownames(full_driver[["driver"]])) {
    row_output <- Execute_Driver_Row(full_driver=full_driver,
                                     driver_rowname=row,
                                     pre_r_mc_parms=pre_r_mc_parms,
                                     do_par=do_par,
                                     n_threads=n_threads)
    saveRDS(object=row_output, 
            file=paste0(path, folder_name, 
                        paste0("/sim_", full_driver[["driver"]][row,"sim_label"])))
  }
}
