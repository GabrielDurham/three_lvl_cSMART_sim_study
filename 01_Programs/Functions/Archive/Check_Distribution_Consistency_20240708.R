######################################
######################################
######################################
### CHECK DISTRIBUTION CONSISTENCY ###
######################################
######################################
######################################

#### PURPOSE: This file contains code to run check the empirical distribution
####          of a given simulation setting
####          and other key simulation parameters for a given simulation setting.
####          The Obtain_Univ_Sim_Params() function will return a list of list of ....
####          Obtain_Univ_Sim_Params()[["n_x"]] - List of parameters for sample size x
####          ...[[d_a1_a2]] - List of parameters for dtr (a1, a2)
####          ...[["var_params"]] = Variance parameters
####              pre_r: Pre-response
####                marg: Marginal variance parameters
####                cond: Conditional variance parameters
####              post_r: Post-response
####                cond: Conditional variance parameters (need mean parameters to derive marg)
####                      (R and NR)
####          ...[["coefs"]] = Simulation coefficients
####              Sigma_01: Pre-response Error variance matrix (t=0,1)
####              Sigma_2: Post-response Error variance matrix (t=2) 
####                       (R and NR)
####              P_vec: Vectors of P coefficients
####                       (R and NR)


#### DATE CREATED:  04 JUN 2024
#### PROGRAMMER:    GABRIEL DURHAM (GJD)
#### EDITS:         


# Process a Driver Row for Distribution Check
### driver = Main sheet of a driver file
### sim_label = sim_label of row to emulate
### cluster_size = Desired cluster size
### n_clusters = Number of clusters to simulate
Process_Driver_Row_for_Check <- function(driver, sim_label, cluster_size, 
                                         n_clusters) {
  driver_row <- driver[driver$sim_label==sim_label, ]
  driver_row$cluster_sizes <- c(cluster_size)
  driver_row$p_cluster_size <- c(1)
  driver_row$n_covar <- 0
  driver_row$n_clusters <- n_clusters
  Output <- Process_Driver_Row_Main(driver_row=driver_row)
  return(Output)
}


# Make Simulated Data Wide (One Row Per Cluster)
### sim_data = Sim_SMART_Data() output
### driver_parms = Output of Process_Driver_Row_Main() function
Make_Sim_Data_Wide <- function(sim_data, driver_parms) {
  temp_output <- data.frame()
  for (cluster in unique(sim_data$cluster_id)) {
    cluster_data_long <- sim_data[sim_data$cluster_id==cluster, ]
    cluster_data_wide <- data.frame(cluster_id=cluster,
                                    a_1=cluster_data_long$a_1[1],
                                    r=cluster_data_long$r[1],
                                    a_2=cluster_data_long$a_2[1])
    for (person in unique(cluster_data_long$person_id)) {
      for (t in c(0,1,2)) {
        var_label <- paste0("y_p_", person, "_t_", t)
        cluster_data_wide[[var_label]] <- 
          cluster_data_long[cluster_data_long$person_id==person&cluster_data_long$t==t, "y"]
      }
    }
    if (driver_parms[["SMART_structure"]]=="prototypical") {
      cluster_data_wide$path <- 
        ifelse(cluster_data_wide$a_1==1&cluster_data_wide$r==1, 1,
               ifelse(cluster_data_wide$a_1==1&cluster_data_wide$a_2==1, 2, 
                      ifelse(cluster_data_wide$a_1==1&cluster_data_wide$a_2==-1, 3,
                             ifelse(cluster_data_wide$a_1==-1&cluster_data_wide$r==1, 4,
                                    ifelse(cluster_data_wide$a_1==-1&cluster_data_wide$a_2==1, 5,
                                           ifelse(cluster_data_wide$a_1==-1&cluster_data_wide$a_2==-1, 6, NA))))))
      cluster_data_wide$c_d_1_1 <- ifelse(cluster_data_wide$path %in% c(1,2), 1, 0)
      cluster_data_wide$c_d_1_m1 <- ifelse(cluster_data_wide$path %in% c(1,3), 1, 0)
      cluster_data_wide$c_d_m1_1 <- ifelse(cluster_data_wide$path %in% c(4,5), 1, 0)
      cluster_data_wide$c_d_m1_m1 <- ifelse(cluster_data_wide$path %in% c(4,6), 1, 0)
    }
    
    temp_output <- rbind(temp_output, cluster_data_wide)
  }
  Output <- temp_output
  return(Output)
}

# Simulate SMART Data and Convert to Wide Format
### driver_parms = Output of Process_Driver_Row_Main() function
### sim_parms = Output of Obtain_All_Sim_Params() function
Sim_SMART_Data_Wide <- function(driver_parms, sim_parms) {
  sim_data_long <- Sim_SMART_Data(driver_parms=driver_parms, 
                                  sim_parms=sim_parms)
  Output <- Make_Sim_Data_Wide(sim_data=sim_data_long, driver_parms=driver_parms)
  return(Output)
}


# Grab the Target Variance for a Single Setting
## E.g., one conditional variance or marginal variance
### pre_r_var_comps = Portion of Obtain_All_Sim_Params() output containing 
###                   applicable pre-response variance components
### post_r_var_comps = Portion of Obtain_All_Sim_Params() output containing 
###                    applicable post-response variance components
### cluster_size = Individuals in a cluster
Grab_Target_Var_Single_Setting <- function(pre_r_var_comps, post_r_var_comps,
                                           cluster_size) {
  # Hard code times
  times <- c(0,1,2)
  pre_r_times <- c(0,1)
  post_r_times <- c(2)
  
  #Grab variance components for diagonal and off-diagonal blocks
  d_block <- matrix(NA, nrow=3, ncol=3)
  od_block <- matrix(NA, nrow=3, ncol=3)
  for (t_1 in times) {
    if (t_1 %in% pre_r_times) {
      d_block[t_1+1, t_1+1] <- pre_r_var_comps[[paste0("sigma2_", t_1)]]
      od_block[t_1+1, t_1+1] <- pre_r_var_comps[[paste0("rho_", t_1)]]
    } else if (t_1 %in% post_r_times) {
      d_block[t_1+1, t_1+1] <- post_r_var_comps[[paste0("sigma2_", t_1)]]
      od_block[t_1+1, t_1+1] <- post_r_var_comps[[paste0("rho_", t_1)]]
    }
    if (t_1<max(times)) {for (t_2 in times[times>t_1]) {
      if (t_2 %in% pre_r_times) {
        d_block[t_1+1, t_2+1] <- pre_r_var_comps[[paste0("phi_", t_1, t_2)]]
        d_block[t_2+1, t_1+1] <- pre_r_var_comps[[paste0("phi_", t_1, t_2)]]
        
        od_block[t_1+1, t_2+1] <- pre_r_var_comps[[paste0("rho_", t_1, t_2)]]
        od_block[t_2+1, t_1+1] <- pre_r_var_comps[[paste0("rho_", t_1, t_2)]]
      } else {
        d_block[t_1+1, t_2+1] <- post_r_var_comps[[paste0("phi_", t_1, t_2)]]
        d_block[t_2+1, t_1+1] <- post_r_var_comps[[paste0("phi_", t_1, t_2)]]
        
        od_block[t_1+1, t_2+1] <- post_r_var_comps[[paste0("rho_", t_1, t_2)]]
        od_block[t_2+1, t_1+1] <- post_r_var_comps[[paste0("rho_", t_1, t_2)]]
      }
    }}
  }
  # Expand blocks into a full variance matrix
  var_mat <- d_block
  for (i in (2:cluster_size)) {var_mat <- rbind(var_mat, od_block)}
  for (i in (2:cluster_size)) {
    new_cols <- od_block
    for (ii in (2:cluster_size)) {
      if (i==ii) {new_cols <- rbind(new_cols, d_block)
      } else {new_cols <- rbind(new_cols, od_block)}
    }
    var_mat <- cbind(var_mat, new_cols)
  }
  Output <- list(mat_full=var_mat,
                 blocks=list(diag=d_block, off_diag=od_block))
  return(Output)
}

# Grab Target Conditional and Marginal Distributions
### sim_parms = Output of Obtain_All_Sim_Params() function
### cluster_size = Number of units in cluster
Grab_Target_Dists <- function(sim_parms, cluster_size) {
  # Hard code times
  times <- c(0,1,2)
  sim_parms_n_i <- sim_parms[[paste0("n_", cluster_size)]]
  temp_output <- NULL
  for (dtr in names(sim_parms_n_i)) {
    sim_parms_dtr_n_i <- sim_parms_n_i[[dtr]]
    
    # Grab means
    marg_means <- c()
    cond_means <- list(R=c(), NR=c())
    
    for (t in times) {
      marg_means <- 
        c(marg_means, 
          sim_parms_dtr_n_i[["target_means"]][["marg"]][[paste0("mean_", t)]])
      cond_means[["R"]] <- 
        c(cond_means[["R"]], 
          sim_parms_dtr_n_i[["target_means"]][["cond"]][["R"]][[paste0("delta_", t)]])
      cond_means[["NR"]] <- 
        c(cond_means[["NR"]], 
          sim_parms_dtr_n_i[["target_means"]][["cond"]][["NR"]][[paste0("delta_", t)]])
    }
    # Grab variances
    pre_r_target_var_parms <- sim_parms_dtr_n_i[["target_var"]][["pre_r"]]
    post_r_target_var_parms <- sim_parms_dtr_n_i[["target_var"]][["post_r"]]
    
    marg_var <- 
      Grab_Target_Var_Single_Setting(pre_r_var_comps=pre_r_target_var_parms[["marg"]], 
                                     post_r_var_comps=post_r_target_var_parms[["marg"]],
                                     cluster_size=cluster_size)
    cond_var <- list(
      R=Grab_Target_Var_Single_Setting(pre_r_var_comps=pre_r_target_var_parms[["cond"]][["R"]], 
                                       post_r_var_comps=post_r_target_var_parms[["cond"]][["R"]],
                                       cluster_size=cluster_size),
      NR=Grab_Target_Var_Single_Setting(pre_r_var_comps=pre_r_target_var_parms[["cond"]][["NR"]], 
                                        post_r_var_comps=post_r_target_var_parms[["cond"]][["NR"]],
                                        cluster_size=cluster_size)
    )
    temp_output[[dtr]] <- list(
      cond = list(R=list(means=cond_means[["R"]],
                         var=cond_var[["R"]]),
                  NR=list(means=cond_means[["NR"]],
                          var=cond_var[["NR"]])),
      marg = list(means=marg_means,
                  var=marg_var)
    )
  }
  Output <- temp_output
  return(Output)
}

# Calculate Empirical Distributions (Cond/Marg) for all DTRs in a Prototypical SMART
### driver_parms = Output of Process_Driver_Row_Main() function
### wide_data = Output of Make_Sim_Data_Wide() function
### sim_parms = Output of Obtain_All_Sim_Params() function
### cluster_size = Number of units in cluster
Grab_Emp_Dists_Proto <- function(wide_data, driver_parms, sim_parms, 
                                 cluster_size) {
  # Hard-code times
  times <- c(0,1,2)
  # Calculate weights
  p_dtrs <- driver_parms[["p_dtr"]]
  p_1 <- p_dtrs[1] + p_dtrs[2]
  p_2_1 <- p_dtrs[1]/(p_dtrs[1] + p_dtrs[2])
  p_2_m1 <- p_dtrs[3]/(p_dtrs[3] + p_dtrs[4])
  wgt <- list("1"=(1/p_1)*1,
              "2"=(1/p_1)*(1/p_2_1),
              "3"=(1/p_1)*(1/(1-p_2_1)),
              "4"=(1/(1-p_1))*1,
              "5"=(1/(1-p_1))*(1/p_2_m1),
              "6"=(1/(1-p_1))*(1/(1-p_2_m1)))
  wide_data$wgt <- NA
  for (path in (1:6)) {
    wgt_val <- wgt[[toString(path)]]
    wide_data$wgt <- ifelse(wide_data$path==path, wgt_val, wide_data$wgt)
  }
  temp_output <- NULL
  for (dtr in c("d_1_1", "d_1_m1", "d_m1_1", "d_m1_m1")) {
    # Calculate empirical distributions
    consistent_clusters <- wide_data[wide_data[[paste0("c_", dtr)]]==1, ]
    emp_dists <- list(
      marg=Summarize_Emp_Data(in_wide_data=consistent_clusters,
                              cluster_size=cluster_size,
                              times=times,
                              weight_col="wgt"),
      cond=list(
        R=Summarize_Emp_Data(in_wide_data=consistent_clusters[consistent_clusters$r==1, ],
                             cluster_size=cluster_size,
                             times=times,
                             weight_col="wgt"),
        NR=Summarize_Emp_Data(in_wide_data=consistent_clusters[consistent_clusters$r==0, ],
                              cluster_size=cluster_size,
                              times=times,
                              weight_col="wgt")
      )
    )
    # Re-organize mean vectors
    emp_dists[["marg"]][["mean"]][["parms"]] <- 
      c(emp_dists[["marg"]][["mean"]][["parms"]][["y_0"]],
        emp_dists[["marg"]][["mean"]][["parms"]][["y_1"]],
        emp_dists[["marg"]][["mean"]][["parms"]][["y_2"]])
    for (r in c("R", "NR")) {
      emp_dists[["cond"]][[r]][["mean"]][["parms"]] <- 
        c(emp_dists[["cond"]][[r]][["mean"]][["parms"]][["y_0"]],
          emp_dists[["cond"]][[r]][["mean"]][["parms"]][["y_1"]],
          emp_dists[["cond"]][[r]][["mean"]][["parms"]][["y_2"]])
    }
    temp_output[[dtr]] <- emp_dists
  }
  
  
  
  Output <- temp_output
  return(Output)
}


# Calculate Empirical Distributions (Cond/Marg) for all DTR for Simulated Data
### driver_parms = Output of Process_Driver_Row_Main() function
### wide_data = Output of Make_Sim_Data_Wide() function
### sim_parms = Output of Obtain_All_Sim_Params() function
### cluster_size = Number of units in cluster
Grab_Emp_Dists <- function(driver_parms, wide_data, sim_parms, cluster_size) {
  if (driver_parms[["SMART_structure"]]=="prototypical") {
    Output <- Grab_Emp_Dists_Proto(wide_data=wide_data, 
                                   driver_parms=driver_parms, 
                                   sim_parms=sim_parms, 
                                   cluster_size=cluster_size)
  }
  return(Output)
}




# Check Empirical/Target Distribution Consistency
### driver_parms = Output of Process_Driver_Row_Main() function
### wide_data = Output of Make_Sim_Data_Wide() function
### sim_parms = Output of Obtain_All_Sim_Params() function
### cluster_size = Number of units in cluster
Check_Distribution_Consistency <- function(wide_data, driver_parms, sim_parms, 
                                           cluster_size) {
  # Hard-code times
  times <- c(0,1,2)
  # Pull empirical and target distributions
  emp_dists_all <- Grab_Emp_Dists(driver_parms=driver_parms, 
                                  wide_data=wide_data, 
                                  sim_parms=sim_parms, 
                                  cluster_size=cluster_size)
  target_dists_all <- Grab_Target_Dists(sim_parms=sim_parms, 
                                        cluster_size=cluster_size)
  
  # Calculate differences and restructure output
  temp_output <- NULL
  for (dtr in names(emp_dists_all)) {
    
    emp_dists <- emp_dists_all[[dtr]]
    target_dists <- target_dists_all[[dtr]]
    
    # Calculate differences
    # Marginal
    mean_diffs_marg <- emp_dists[["marg"]][["mean"]][["parms"]] - 
      target_dists[["marg"]][["means"]]
    #Variance
    var_diffs_marg <- list(
      mat_full=emp_dists[["marg"]][["var"]][["mat_full"]] - 
        target_dists[["marg"]][["var"]][["mat_full"]],
      blocks=list(
        diag=emp_dists[["marg"]][["var"]][["mat_means"]][["diag"]] - 
          target_dists[["marg"]][["var"]][["blocks"]][["diag"]],
        off_diag=emp_dists[["marg"]][["var"]][["mat_means"]][["off_diag"]] - 
          target_dists[["marg"]][["var"]][["blocks"]][["off_diag"]]
      ))
    # Conditional
    mean_diffs_cond <- NULL
    var_diffs_cond <- NULL
    for (r in c("R", "NR")) {
      mean_diffs_cond[[r]] <- emp_dists[["cond"]][[r]][["mean"]][["parms"]] - 
        target_dists[["cond"]][[r]][["means"]]
      var_diffs_cond[[r]] <- list(
        mat_full=emp_dists[["cond"]][[r]][["var"]][["mat_full"]] - 
          target_dists[["cond"]][[r]][["var"]][["mat_full"]],
        blocks=list(
          diag=emp_dists[["cond"]][[r]][["var"]][["mat_means"]][["diag"]] - 
            target_dists[["cond"]][[r]][["var"]][["blocks"]][["diag"]],
          off_diag=emp_dists[["cond"]][[r]][["var"]][["mat_means"]][["off_diag"]] - 
            target_dists[["cond"]][[r]][["var"]][["blocks"]][["off_diag"]]
        ))
    }
    n_dtr <- nrow(wide_data[wide_data[[paste0("c_", dtr)]]==1, ])
    n_r <- nrow(wide_data[(wide_data[[paste0("c_", dtr)]]==1&wide_data$r==1), ])
    n_nr <- nrow(wide_data[(wide_data[[paste0("c_", dtr)]]==1&wide_data$r==0), ])
    temp_output[[dtr]] <- NULL
    temp_output[[dtr]][["dists"]] <- list(emp=emp_dists, target=target_dists)
    temp_output[[dtr]][["diffs"]] <- list(marg=list(mean=mean_diffs_marg, 
                                                    var=var_diffs_marg), 
                                          cond=list(R=list(mean=mean_diffs_cond[["R"]],
                                                           var=var_diffs_cond[["R"]]),
                                                    NR=list(mean=mean_diffs_cond[["NR"]],
                                                            var=var_diffs_cond[["NR"]])))
    temp_output[[dtr]][["n"]] <- list(marg=n_dtr, cond=list(R=n_r, NR=n_nr))
  }
  
  Output <- temp_output
  return(Output)
  
}





### consistency_check_dtr = Single DTR element of 
###                         Check_Distribution_Consistency() output
### setting = Type of setting ("marg", "cond_r", "cond_nr")
### digits = Number of digits to display
Format_Consistency_Single_Setting <- function(consistency_check_dtr, setting,
                                              digits=4) {
  if (setting=="marg") {
    emp_dists <- consistency_check_dtr[["dists"]][["emp"]][["marg"]]
    tgt_dists <- consistency_check_dtr[["dists"]][["target"]][["marg"]]
    diffs <- consistency_check_dtr[["diffs"]][["marg"]]
    n <- consistency_check_dtr[["n"]][["marg"]]
  } else if (setting=="cond_r") {
    emp_dists <- consistency_check_dtr[["dists"]][["emp"]][["cond"]][["R"]]
    tgt_dists <- consistency_check_dtr[["dists"]][["target"]][["cond"]][["R"]]
    diffs <- consistency_check_dtr[["diffs"]][["cond"]][["R"]]
    n <- consistency_check_dtr[["n"]][["cond"]][["R"]]
  } else if (setting=="cond_nr") {
    emp_dists <- consistency_check_dtr[["dists"]][["emp"]][["cond"]][["NR"]]
    tgt_dists <- consistency_check_dtr[["dists"]][["target"]][["cond"]][["NR"]]
    diffs <- consistency_check_dtr[["diffs"]][["cond"]][["NR"]]
    n <- consistency_check_dtr[["n"]][["cond"]][["NR"]]
  }
  
  mean_emp <- format(round(emp_dists[["mean"]][["parms"]], digits), nsmall=digits)
  mean_tgt <- format(round(tgt_dists[["means"]], digits), nsmall=digits)
  mean_diff <- format(round(diffs[["mean"]], digits), nsmall=digits)
  
  
  
  var_diag_emp <- 
    format(round(emp_dists[["var"]][["mat_means"]][["diag"]], digits), nsmall=digits)
  var_diag_tgt <- 
    format(round(tgt_dists[["var"]][["blocks"]][["diag"]], digits), nsmall=digits)
  var_diag_diff <- 
    format(round(diffs[["var"]][["blocks"]][["diag"]], digits), nsmall=digits)
  
  var_off_diag_emp <- 
    format(round(emp_dists[["var"]][["mat_means"]][["off_diag"]], digits), nsmall=digits)
  var_off_diag_tgt <- 
    format(round(tgt_dists[["var"]][["blocks"]][["off_diag"]], digits), nsmall=digits)
  var_off_diag_diff <- 
    format(round(diffs[["var"]][["blocks"]][["off_diag"]], digits), nsmall=digits)
  
  Output <- list(
    mean=list(emp=mean_emp, 
              target=mean_tgt, 
              diff=mean_diff),
    var=list(diag=list(emp=var_diag_emp, 
                       target=var_diag_tgt, 
                       diff=var_diag_diff),
             off_diag=list(emp=var_off_diag_emp, 
                           target=var_off_diag_tgt, 
                           diff=var_off_diag_diff)),
    n=n
  )
  return(Output)
}





# Display Consistency Comparisons for a Single DTR
### consistency_check = Output of Check_Distribution_Consistency() function
### dtr = DTR to analyze
### title = Title of printed output
### verbose = Whether or not to print output
Display_Topline_Consistency <- function(consistency_check, dtr, title, verbose=TRUE) {
  # Hard-code times
  times <- c(0,1,2)
  T <- length(times)
  
  # Format output
  consistency_check_dtr <- consistency_check[[dtr]]
  fmt_output <- list(
    marg=Format_Consistency_Single_Setting(consistency_check_dtr=consistency_check_dtr, 
                                           setting="marg"),
    cond_r=Format_Consistency_Single_Setting(consistency_check_dtr=consistency_check_dtr, 
                                             setting="cond_r"),
    cond_nr=Format_Consistency_Single_Setting(consistency_check_dtr=consistency_check_dtr, 
                                              setting="cond_nr")
  )
  if (verbose) {
    divider <- "\n"
    for (i in (1:(8+nchar(title)))) {divider <- paste0(divider, "*")}
    
    cat(divider)
    cat()
    cat(paste0("\n*** ", title, " ***"))
    cat(divider)
    # Print out output for individual settings
    for (name in names(fmt_output)) {
      output_data <- fmt_output[[name]]
      beginner_text <- 
        ifelse(name=="marg", "\nMarginal Distribution: ",
               ifelse(name=="cond_r", "\nConditional Distribution (R): ", 
                      ifelse(name=="cond_nr", "\nConditional Distribution (NR): ", NA)))
      cat(beginner_text)
      cat(paste0("\n*Empirical results based off of ", output_data[["n"]], 
                 " simulated clusters."))
      cat("\nMean: ")
      cat("\nEmpirical: ")
      cat('\n', output_data[["mean"]][["emp"]])
      cat("\nTarget: ")
      cat('\n', output_data[["mean"]][["target"]])
      cat("\nDifference: ")
      cat('\n', output_data[["mean"]][["diff"]])
      cat('\n')
      cat("\nWithin-Person Variance: ")
      cat("\nEmpirical: ")
      for (t in (1:T)) {cat('\n', output_data[["var"]][["diag"]][["emp"]][t,])}
      cat("\nTarget: ")
      for (t in (1:T)) {cat('\n', output_data[["var"]][["diag"]][["target"]][t,])}
      cat("\nDifference: ")
      for (t in (1:T)) {cat('\n', output_data[["var"]][["diag"]][["diff"]][t,])}
      cat('\n')
      cat("\nBetween-Person Variance: ")
      cat("\nEmpirical: ")
      for (t in (1:T)) {cat('\n', output_data[["var"]][["off_diag"]][["emp"]][t,])}
      cat("\nTarget: ")
      for (t in (1:T)) {cat('\n', output_data[["var"]][["off_diag"]][["target"]][t,])}
      cat("\nDifference: ")
      for (t in (1:T)) {cat('\n', output_data[["var"]][["off_diag"]][["diff"]][t,])}
      cat('\n\n')
    }
  }
  Output <- list(marg=fmt_output[["marg"]],
                 cond=list(R=fmt_output[["cond_r"]], NR=fmt_output[["cond_nr"]]))
  return(Output)
}