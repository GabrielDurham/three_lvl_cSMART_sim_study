##########################################
##########################################
##########################################
### GENERATE ALL SIMULATION PARAMETERS ###
##########################################
##########################################
##########################################

#### PURPOSE: This file contains code to collect all data-generative simulation 
####          parameters for a row of the main driver file.
####          The Obtain_All_Sim_Params() function will return a list of list of ....
####          Obtain_All_Sim_Params()[["n_x"]] - List of parameters for sample size x
####          ...[[d_a1_a2]] - List of parameters for dtr (a1, a2)
####          ...[["var_params"]] = Variance parameters
####              pre_r: Pre-response
####                marg: Marginal variance parameters
####                cond: Conditional variance parameters
####              post_r: Post-response
####                marg: Marginal variance parameters
####                cond: Conditional variance parameters (need mean parameters to derive marg)
####                      (R and NR)
####          ...[["coefs"]] = Simulation coefficients
####              Sigma_01: Pre-response Error variance matrix (t=0,1)
####              Sigma_2: Post-response Error variance matrix (t=2) 
####                       (R and NR)
####              P_vec: Vectors of P coefficients
####                       (R and NR)
####          ...[["mean_parms"]] = List of mean parameters for simulation
####          ...[["response_function"]] = Response function


#### DATE CREATED:  08 JAN 2024
#### PROGRAMMER:    GABRIEL DURHAM (GJD)
#### EDITS:         27 MAY 2024 (GJD): Added adjustment to mean parameter derivation
####                  via Derive_Cond_Mean_Adj() function which enables specifying
####                  the path-conditional t=2 means
####                04 JUN 2024 (GJD): Added target mean to output and changed
####                  "var_parms" output label to "target_var"





# Derive Adjustments to Match Conditional Mean
### cp_settings = Dataframe of settings of conditional parameters
###               *** Restricted to simulation setting under consideration
### pre_r_cond_var_parms_n_i = List of conditional pre-response variance parameter settings
###                        (e.g., those derived via MC) - MUST BE RESTRICTED TO SAMPLE SIZE
###                        [["R"]]- Conditioned on R=1
###                        [["NR"]]- Conditioned on R=0
### univ_sim_parms_n_i = Simulation parameters for a given sample size
###                      Part of output of Obtain_Univ_Sim_Params() 
### t = Time for mean of which we're taking adjustment
### path = Path which we're conditioning on for adjustment
Derive_Cond_Mean_Adj <- function(cp_settings, pre_r_cond_var_parms_n_i, 
                                 univ_sim_parms_n_i, t, path) {
  pre_r_var_str <- cp_settings[cp_settings$pathway==path, "pre_r_var_str"]
  if (cp_settings$SMART_structure[1]=="prototypical") {
    if (t %in% c(0,1)) {Output <- 0
    } else {
      dtr <- ifelse(path %in% c(1, 2), "d_1_1",
                    ifelse(path==3, "d_1_m1",
                           ifelse(path %in% c(4,5), "d_m1_1",
                                  ifelse(path==6, "d_m1_m1", NA))))
      r <- ifelse(path %in% c(1, 4), "R", "NR")
      P <- univ_sim_parms_n_i[[dtr]][["coefs"]][["P_vec"]][[r]]
      cond_errs_df <- pre_r_cond_var_parms_n_i[[r]]
      cond_errs <- 
        cond_errs_df[cond_errs_df$cond_param_setting_pre_r==pre_r_var_str, ]
      e_0 <- cond_errs$mean_e_0
      e_1 <- cond_errs$mean_e_1
      Output <- P[2]*e_0 + P[3]*e_0 + P[4]*(P[1]*e_0 + e_1) + P[5]*e_1
    }
  }
  return(Output)
}


# Calculate Mean Parameters for Prototypical SMART
### cp_settings = Dataframe of settings of conditional parameters
###               *** Restricted to simulation setting under consideration
### pre_r_marg_parms = Dataframe of marginal pre-response variance parameter settings
### pre_r_cond_var_parms_n_i = List of conditional pre-response variance parameter settings
###                        (e.g., those derived via MC) - MUST BE RESTRICTED TO SAMPLE SIZE
###                        [["R"]]- Conditioned on R=1
###                        [["NR"]]- Conditioned on R=0
### univ_sim_parms_n_i = Simulation parameters for a given sample size
###                      Part of output of Obtain_Univ_Sim_Params() 
Calculate_Mean_Parms_Prototypical <- function(cp_settings, pre_r_marg_parms,
                                              pre_r_cond_var_parms_n_i, 
                                              univ_sim_parms_n_i) {
  
  rownames(cp_settings) <- cp_settings$pathway
  var_set_1 <- cp_settings[1, "pre_r_var_str"]
  var_set_m1 <- cp_settings[4, "pre_r_var_str"]
  p_r_1 <- 
    pre_r_marg_parms[pre_r_marg_parms$cond_param_setting_pre_r==var_set_1, "p_r"]
  p_r_m1 <- 
    pre_r_marg_parms[pre_r_marg_parms$cond_param_setting_pre_r==var_set_m1, "p_r"]
  coef_mat <- matrix(data=c(
    1, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 1, 1, 0, 0, 0, 0, 0, 0,
    1, 1, -1, 0, 0, 0, 0, 0, 0,
    1, 1, 1, 1, 1, 0, 0, (1-p_r_1), (1-p_r_1),
    1, 1, 1, 1, 1, 1/(1-p_r_1), 1/(1-p_r_1), -p_r_1, -p_r_1,
    1, 1, 1, 1, 1, -1/(1-p_r_1), -1/(1-p_r_1), -p_r_1, -p_r_1,
    1, 1, -1, 1, -1, 0, 0, (1-p_r_m1), -(1-p_r_m1),
    1, 1, -1, 1, -1, 1/(1-p_r_m1), -1/(1-p_r_m1), -p_r_m1, p_r_m1,
    1, 1, -1, 1, -1, -1/(1-p_r_m1), 1/(1-p_r_m1), -p_r_m1, p_r_m1),
    nrow=9, ncol=9, byrow=TRUE)
  
  mean_vec <- matrix(data=c(
    cp_settings[cp_settings$pathway==1,"mean_0"],
    cp_settings[cp_settings$pathway==1,"mean_1"],
    cp_settings[cp_settings$pathway==4,"mean_1"],
    cp_settings[cp_settings$pathway==1,"mean_2"],
    cp_settings[cp_settings$pathway==2,"mean_2"],
    cp_settings[cp_settings$pathway==3,"mean_2"],
    cp_settings[cp_settings$pathway==4,"mean_2"],
    cp_settings[cp_settings$pathway==5,"mean_2"], 
    cp_settings[cp_settings$pathway==6,"mean_2"]), ncol=1)
  
  #Need to adjust to ensure conditional mean specifications
  adj_vec <- matrix(data=0, nrow=9, ncol=1)
  for (row in (4:9)) {
    adj_vec[row, 1] <- 
      Derive_Cond_Mean_Adj(cp_settings=cp_settings, 
                           pre_r_cond_var_parms_n_i=pre_r_cond_var_parms_n_i, 
                           univ_sim_parms_n_i=univ_sim_parms_n_i, 
                           t=2, 
                           path=row-3)
  }
  gamma_vec <- solve(coef_mat)%*%(mean_vec-adj_vec)
  Output <- NULL
  for (i in (0:6)) {Output[[paste0("gamma_",i)]] <- gamma_vec[i+1,1]}
  Output[["lambda_1"]] <- gamma_vec[8,1]
  Output[["lambda_2"]] <- gamma_vec[9,1]
  return(Output)
}

# Calculate Conditional Means for Prototypical SMARTs for a given DTR/n_i
### mean_parms = Mean parameters, output of Calculate_Mean_Parms_Prototypical() 
### univ_sim_parms_dtr_n_i = Simulation parameters for a given dtr/sample size
###                      Part of output of Obtain_Univ_Sim_Params() 
Calculate_Cond_Mean_DTR_Prototypical <- function(mean_parms, univ_sim_parms_dtr_n_i) {
  cond_var_parms_pre_r <- univ_sim_parms_dtr_n_i[["var_parms"]][["pre_r"]][["cond"]]
  a_1 <- univ_sim_parms_dtr_n_i[["dtr"]][["a_1"]]
  a_2nr <- univ_sim_parms_dtr_n_i[["dtr"]][["a_2nr"]]
  p_r <- univ_sim_parms_dtr_n_i[["var_parms"]][["pre_r"]][["marg"]][["p_r"]]
  
  Output <- NULL
  for (r_lab in c("R", "NR")) {
    r <- ifelse(r_lab=="R", 1, 0)
    delta <- NULL
    
    e_0 <- cond_var_parms_pre_r[[r_lab]][["mean_e_0"]]
    e_1 <- cond_var_parms_pre_r[[r_lab]][["mean_e_1"]]
    P <- univ_sim_parms_dtr_n_i[["coefs"]][["P_vec"]][[r_lab]]
    
    # Note that we're only concerned with differences between means for this purpose
    # so we don't need to incorporate covariates here
    delta[["delta_0"]] <- mean_parms[["gamma_0"]] + e_0
    delta[["delta_1"]] <- mean_parms[["gamma_0"]] + mean_parms[["gamma_1"]] + 
      mean_parms[["gamma_2"]]*a_1 + P[1]*e_0 + e_1
    delta[["delta_2"]] <- mean_parms[["gamma_0"]] + mean_parms[["gamma_1"]] + 
      mean_parms[["gamma_2"]]*a_1 + mean_parms[["gamma_3"]] + mean_parms[["gamma_4"]]*a_1 + 
      (mean_parms[["gamma_5"]] + mean_parms[["gamma_6"]]*a_1)*a_2nr*((1-r)/(1-p_r)) + 
      (mean_parms[["lambda_1"]] + mean_parms[["lambda_2"]]*a_1)*(r-p_r) + 
      P[2]*e_0 + P[3]*e_0 + P[4]*(P[1]*e_0 + e_1) + P[5]*e_1
    Output[[r_lab]] <- delta
  }
  return(Output)
}

# Calculate Marginal Means for Prototypical SMARTs for a given DTR/n_i
### cp_settings = Dataframe of settings of conditional parameters
###               *** Restricted to simulation setting under consideration
### univ_sim_parms_dtr_n_i = Simulation parameters for a given dtr/sample size
###                      Part of output of Obtain_Univ_Sim_Params() 
Calculate_Marg_Mean_DTR_Prototypical <- function(cp_settings, univ_sim_parms_dtr_n_i) {
  a_1 <- univ_sim_parms_dtr_n_i[["dtr"]][["a_1"]]
  a_2nr <- univ_sim_parms_dtr_n_i[["dtr"]][["a_2nr"]]
  p_r <- univ_sim_parms_dtr_n_i[["var_parms"]][["pre_r"]][["marg"]][["p_r"]]
  
  pathways_df <- Define_SMART_Pathways("prototypical")
  paths_dtr <- pathways_df[(pathways_df$a_1==a_1&pathways_df$a_2nr==a_2nr), ]
  
  means_r <- cp_settings[cp_settings$pathway==paths_dtr$path_r, 
                         c("mean_0", "mean_1", "mean_2")]
  means_nr <- cp_settings[cp_settings$pathway==paths_dtr$path_nr, 
                          c("mean_0", "mean_1", "mean_2")]
  marg_means <- p_r*means_r + (1-p_r)*means_nr
  Output <- list(
    mean_0 = marg_means$mean_0,
    mean_1 = marg_means$mean_1,
    mean_2 = marg_means$mean_2
  )
  return(Output)
}


# Calculate Marginal (Post-Response) Variances for a Protypical SMART
# (For a Single DTR/Sample Size)
### mean_parms = Mean paramters, output of Calculate_Mean_Parms_Prototypical() 
### univ_sim_parms_dtr_n_i = Simulation parameters for a given dtr/sample size
###                      Part of output of Obtain_Univ_Sim_Params() 
Calculate_Marg_Var_Post_R_DTR_Prototypical <- function(mean_parms, univ_sim_parms_dtr_n_i) {
  c_means <- Calculate_Cond_Mean_DTR_Prototypical(mean_parms=mean_parms,
                                                  univ_sim_parms_dtr_n_i=univ_sim_parms_dtr_n_i)
  c_vars <- univ_sim_parms_dtr_n_i[["var_parms"]][["post_r"]][["cond"]]
  p_r <- univ_sim_parms_dtr_n_i[["var_parms"]][["pre_r"]][["marg"]][["p_r"]]
  Output <- NULL
  for (t in c(0,1,2)) {
    if (t!=2) {
      rho_lab <- paste0("rho_", t, "2")
      phi_lab <- paste0("phi_", t, "2")
    } else {
      rho_lab <- "rho_2"
      phi_lab <- "sigma2_2"
    }
    delta_lab <- paste0("delta_", t)
    mean_term <- 
      (c_means[["NR"]][[delta_lab]] - c_means[["R"]][[delta_lab]])*
      (c_means[["NR"]][["delta_2"]] - c_means[["R"]][["delta_2"]])*
      ( p_r*(1-p_r) )
    
    Output[[rho_lab]] <- p_r*c_vars[["R"]][[rho_lab]] + 
      (1-p_r)*c_vars[["NR"]][[rho_lab]] + mean_term
    Output[[phi_lab]] <- p_r*c_vars[["R"]][[phi_lab]] + 
      (1-p_r)*c_vars[["NR"]][[phi_lab]] + mean_term
  }
  #Output[["c_means"]] <- c_means
  return(Output)
}





# Generate Response Function for a Given DTR/Sample Size
### pre_r_marg_parms = Dataframe of marginal pre-response variance parameter settings
### univ_sim_parms_dtr_n_i = Simulation parameters for a given dtr/sample size
###                      Part of output of Obtain_Univ_Sim_Params() 
Generate_Response_Function_DTR_n <- function(pre_r_marg_parms, univ_sim_parms_dtr_n_i) {
  # Use code from Pre_R_Cond_Param_Sim.R file - Uses pre_r_marg_parms as driver
  # Builds response function taking mean outcome deviations as sole input
  exp_marg_parms <- Expand_Driver_File(raw_driver_file=pre_r_marg_parms)
  pre_r_var_setting <- 
    univ_sim_parms_dtr_n_i[["var_parms"]][["pre_r"]][["marg"]][["cond_param_setting_pre_r"]]
  n_i <- univ_sim_parms_dtr_n_i[["var_parms"]][["pre_r"]][["marg"]][["n_i"]]
  # Isolate relevant driver row
  rel_driver_row <- 
    exp_marg_parms[(exp_marg_parms$cond_param_setting_pre_r==pre_r_var_setting&exp_marg_parms$sample_size==n_i),]
  pre_r_var_parms <- Process_Driver_Row_Pre_R_MC(driver_row=rel_driver_row)
  Output <- Build_Single_Input_R_Function(parms=pre_r_var_parms)
  return(Output)
}

# Calculate Mean-Dependent Parameters and Adds to Universal Simulation Parameters
## Also adds response function
### cp_settings = Dataframe of settings of conditional parameters
###               *** Restricted to simulation setting under consideration
### pre_r_marg_parms = Dataframe of marginal pre-response variance parameter settings
### univ_sim_parms = Simulation parameters outputted by Obtain_Univ_Sim_Params() 
#### Returns copy of univ_sim_parms with marginal post-response variances and mean parameters
Complete_Simulation_Params <- function(cp_settings, pre_r_marg_parms, 
                                       pre_r_cond_var_parms, univ_sim_parms) {
  temp_output <- univ_sim_parms
  SMART_structure <- cp_settings$SMART_structure[1]
  for (sample_size in names(univ_sim_parms)) {
    n_i <- as.numeric(strsplit(sample_size, "_")[[1]][2])
    univ_sim_parms_n_i <- univ_sim_parms[[sample_size]]
    pre_r_cond_var_parms_n_i <- list(
      R=pre_r_cond_var_parms[["R"]][pre_r_cond_var_parms[["R"]]$n_i==n_i, ],
      NR=pre_r_cond_var_parms[["NR"]][pre_r_cond_var_parms[["NR"]]$n_i==n_i, ]
    )
    if (SMART_structure=="prototypical") {
      mean_parms <- 
        Calculate_Mean_Parms_Prototypical(cp_settings=cp_settings, 
                                          pre_r_marg_parms=pre_r_marg_parms,
                                          pre_r_cond_var_parms_n_i=pre_r_cond_var_parms_n_i, 
                                          univ_sim_parms_n_i=univ_sim_parms_n_i)
    }
    
    for (dtr in names(univ_sim_parms_n_i)) {
      univ_sim_parms_dtr_n_i <- univ_sim_parms_n_i[[dtr]]
      temp_output[[sample_size]][[dtr]][["mean_parms"]] <- mean_parms
      names(temp_output[[sample_size]][[dtr]])[names(temp_output[[sample_size]][[dtr]])=="var_parms"] <- 
        "target_var"
      if (SMART_structure=="prototypical") {
        temp_output[[sample_size]][[dtr]][["target_var"]][["post_r"]][["marg"]] <- 
          Calculate_Marg_Var_Post_R_DTR_Prototypical(mean_parms=mean_parms,
                                                     univ_sim_parms_dtr_n_i=univ_sim_parms_dtr_n_i)
        temp_output[[sample_size]][[dtr]][["target_means"]] <- list(
          cond=Calculate_Cond_Mean_DTR_Prototypical(mean_parms=mean_parms,
                                                    univ_sim_parms_dtr_n_i=univ_sim_parms_dtr_n_i),
          marg=Calculate_Marg_Mean_DTR_Prototypical(cp_settings=cp_settings,
                                                    univ_sim_parms_dtr_n_i=univ_sim_parms_dtr_n_i)
        )
        
        temp_output[[sample_size]][[dtr]][["response_function"]] <- 
          Generate_Response_Function_DTR_n(pre_r_marg_parms=pre_r_marg_parms,
                                           univ_sim_parms_dtr_n_i=univ_sim_parms_dtr_n_i)
      }
    }
  }
  Output <- temp_output
  return(Output)
}



# Obtain All Simulation Parameters for a Given Row
### driver_parms = Output of Process_Driver_Row_Main() function
### cp_settings_df = Dataframe of settings of conditional parameters
### pre_r_marg_parms = Dataframe of marginal pre-response variance parameter settings
### pre_r_cond_var_parms = List of conditional pre-response variance parameter settings
###                        (e.g., those derived via MC)
###                        [["R"]]- Conditioned on R=1
###                        [["NR"]]- Conditioned on R=0
### post_r_var_parms = Dataframe of post-response variance parameter settings
#### Outputs list for each embedded DTR
####    var_params = Variance parameters
####        pre_r: Pre-response
####          marg: Marginal variance parameters
####          cond: Conditional variance parameters
####        post_r: Post-response
####          marg: Marginal variance parameters
####          cond: Conditional variance parameters (need mean parameters to derive marg)
####                 (R and NR)
####    coefs = Simulation coefficients
####        Sigma_01: Pre-response Error variance matrix (t=0,1)
####        Sigma_2: Post-response Error variance matrix (t=2) 
####                 (R and NR)
####        P_vec: Vectors of P coefficients
####                 (R and NR)
####    mean_parms = List of parameters for mean simulation
####    response_function = Response function
Obtain_All_Sim_Params <- function(driver_parms, cp_settings_df, pre_r_marg_parms, 
                                  pre_r_cond_var_parms, post_r_var_parms){
  #Obtain relevant conditional parameter settings
  cp_settings <- 
    cp_settings_df[cp_settings_df$cond_parm_setting==driver_parms[["cp_setting"]],]
  
  #Obtain parameters which don't depend on mean specification
  univ_sim_parms <- Obtain_Univ_Sim_Params(cp_settings=cp_settings, 
                                      n_i=driver_parms[["cluster_sizes"]],
                                      pre_r_marg_parms=pre_r_marg_parms,
                                      pre_r_cond_var_parms=pre_r_cond_var_parms, 
                                      post_r_var_parms=post_r_var_parms)
  #Augment with mean-specific parameters
  Output <- Complete_Simulation_Params(cp_settings=cp_settings, 
                                       pre_r_marg_parms=pre_r_marg_parms,
                                       pre_r_cond_var_parms=pre_r_cond_var_parms,
                                       univ_sim_parms=univ_sim_parms)
  return(Output)
}