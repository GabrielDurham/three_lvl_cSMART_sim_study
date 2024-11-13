####################################
####################################
####################################
### DERIVE SIMULATION PARAMETERS ###
####################################
####################################
####################################

#### PURPOSE: This file contains code to support the Derive_Sim_Parms() function,
####          which derives the simulation parameters for the specified simulation
####          setting. It is a cleaned-up and unified version of previous files
####          "Derive_Error_Distributions.R" and "Generate_All_Simulation_Parameters.R"


#### DATE CREATED:  03 JUL 2024
#### PROGRAMMER:    GABRIEL DURHAM (GJD)
#### EDITS:         08 JUL 2024 (GJD) - Changed name of post_r_cond_parms input
####                    to Derive_Sim_Parms() to post_r_var_parms



# Calculate Marginal and Conditional Means
### driver_parms = Output of Process_Driver_Row_Main() function
### cp_settings = Dataframe of settings of conditional parameters
###               *** Restricted to simulation setting under consideration
### pre_r_marg_parms = Dataframe of marginal pre-response variance parameter settings
### pre_r_cond_parms_n_i = List of conditional pre-response variance parameter settings
###                        (e.g., those derived via MC) - MUST BE RESTRICTED TO SAMPLE SIZE
###                        [["R"]]- Conditioned on R=1
###                        [["NR"]]- Conditioned on R=0
Calculate_Means_DTR <- function(cp_settings, dtr, pre_r_marg_parms, 
                                pre_r_cond_parms_n_i) {
  var_spec <- cp_settings$var_spec[1]
  SMART_structure_df <- 
    Define_SMART_Pathways(SMART_structure=cp_settings$SMART_structure[1])
  DTR_structure <- SMART_structure_df[SMART_structure_df$dtr==dtr, ]
  mc_parms_n_i_r <- pre_r_cond_parms_n_i[["R"]]
  mc_parms_n_i_nr <- pre_r_cond_parms_n_i[["NR"]]
  
  marg_means <- NULL
  cond_means <- NULL
  if (var_spec=="cond") {
    cp_settings_r <- cp_settings[cp_settings$pathway==DTR_structure$path_r, ]
    cp_settings_nr <- cp_settings[cp_settings$pathway==DTR_structure$path_nr, ]
    #Must have same pre-r variance specification for R/NR path
    pre_r_var_str <- cp_settings_r$pre_r_var_str
    mc_parms_r <- 
      mc_parms_n_i_r[mc_parms_n_i_r$cond_param_setting_pre_r==pre_r_var_str, ]
    mc_parms_nr <- 
      mc_parms_n_i_nr[mc_parms_n_i_r$cond_param_setting_pre_r==pre_r_var_str, ]
    p_r <- 
      pre_r_marg_parms[pre_r_marg_parms$cond_param_setting_pre_r==pre_r_var_str, "p_r"]
    # Take specified means (should be equal for pre-response)
    d_means_r <- cp_settings_r[, c("mean_0", "mean_1", "mean_2")]
    d_means_nr <- cp_settings_nr[, c("mean_0", "mean_1", "mean_2")]
    d_means_marg <- (p_r*d_means_r) + ((1-p_r)*d_means_nr)
    
    # Grab Monte-Carlo estimates for conditional means
    d_means_r$mean_0 <- d_means_r$mean_0 + mc_parms_r$mean_y_0
    d_means_r$mean_1 <- d_means_r$mean_1 + mc_parms_r$mean_y_1
    d_means_nr$mean_0 <- d_means_nr$mean_0 + mc_parms_nr$mean_y_0
    d_means_nr$mean_1 <- d_means_nr$mean_1 + mc_parms_nr$mean_y_1
    
    # Store means
    marg_means <- d_means_marg
    cond_means <- list(R=d_means_r, NR=d_means_nr)
    # Grab means if specified marginally
  } else if (var_spec %in% c("marg_r", "marg_nr")) {
    cp_settings_marg <- 
      cp_settings[(!is.na(cp_settings$dtr))&(cp_settings$dtr==dtr), ]
    pre_r_var_str <- cp_settings_marg$pre_r_var_str[1]
    mc_parms_r <- 
      mc_parms_n_i_r[mc_parms_n_i_r$cond_param_setting_pre_r==pre_r_var_str, ]
    mc_parms_nr <- 
      mc_parms_n_i_nr[mc_parms_n_i_r$cond_param_setting_pre_r==pre_r_var_str, ]
    p_r <- 
      pre_r_marg_parms[pre_r_marg_parms$cond_param_setting_pre_r==pre_r_var_str, "p_r"]
    
    #Derive details of the specified conditional distribution
    if (var_spec=="marg_r") {
      path_rs <- DTR_structure$path_r
      p_rs <- p_r
      mc_parms_rs <- mc_parms_r
      mc_parms_rns <- mc_parms_nr
    } else if (var_spec=="marg_nr") {
      path_rs <- DTR_structure$path_nr
      p_rs <- 1-p_r
      mc_parms_rs <- mc_parms_nr
      mc_parms_rns <- mc_parms_r
    }
    cp_settings_rs <- 
      cp_settings[(!is.na(cp_settings$pathway))&(cp_settings$pathway==path_rs), ]
    
    # d_means_marg and d_means_rs should match at times 0 and 1
    d_means_marg <- cp_settings_marg[, c("mean_0", "mean_1", "mean_2")]
    d_means_rs <- cp_settings_rs[, c("mean_0", "mean_1", "mean_2")]
    d_means_rns <- (d_means_marg-(p_rs*d_means_rs))/(1-p_rs)
    
    # Add MC parameter corrections
    d_means_rs$mean_0 <- d_means_rs$mean_0 + mc_parms_rs$mean_y_0
    d_means_rs$mean_1 <- d_means_rs$mean_1 + mc_parms_rs$mean_y_1
    d_means_rns$mean_0 <- d_means_rns$mean_0 + mc_parms_rns$mean_y_0
    d_means_rns$mean_1 <- d_means_rns$mean_1 + mc_parms_rns$mean_y_1
    
    marg_means <- d_means_marg
    if (var_spec=="marg_r") {
      cond_means <- list(R=d_means_rs, NR=d_means_rns)
    } else if (var_spec=="marg_nr") {
      cond_means <- list(R=d_means_rns, NR=d_means_rs)
    }
  }
  
  Output <- list(marg=marg_means, cond=cond_means)
  return(Output)
}




# Obtain Pre-R Variance Parameters
### pre_r_setting = Variance setting (cond_param_setting_pre_r) 
###                 Note - this is conditional on path up to given point
###                 (So marginal given it's pre-response)
### pre_r_marg_parms = Dataframe of marginal pre-response variance parameter settings
### pre_r_cond_var_parms_n_i = List of conditional pre-response variance parameter settings
###                        (e.g., those derived via MC) - NEEDS TO BE LIMITED
###                        TO SAMPLE SIZE UNDER CONSIDERATION
###                        [["R"]]- Conditioned on R=1
###                        [["NR"]]- Conditioned on R=0
#### Returns a list:
####    Output[["marg"]] - List of marginal parameters
####    Output[["cond"]][["R"]] - List of parameters conditioned on response=1
####    Output[["cond"]][["NR"]] - List of parameters conditioned on response=0
Obtain_Pre_R_Var_Parms <- function(pre_r_setting, pre_r_marg_parms, pre_r_cond_parms_n_i){
  
  pre_r_marg_var_str <- 
    pre_r_marg_parms[pre_r_marg_parms$cond_param_setting_pre_r==pre_r_setting,]
  pre_r_cond_var_str <- NULL
  for (r in c("R","NR")) {
    pre_r_cond_var_str[[r]] <- 
      pre_r_cond_parms_n_i[[r]][pre_r_cond_parms_n_i[[r]]$cond_param_setting_pre_r==
                                  pre_r_setting,]
  }
  
  temp_output <- NULL
  temp_output[["marg"]] <- NULL
  temp_output[["cond"]] <- NULL
  for (comp in colnames(pre_r_marg_var_str)) {
    temp_output[["marg"]][[comp]] <- pre_r_marg_var_str[[comp]]
  }
  temp_output[["marg"]][["n_i"]] <- pre_r_cond_var_str[["R"]][["n_i"]]
  
  for (r in c("R","NR")) {
    temp_output[["cond"]][[r]] <- NULL
    for (comp in colnames(pre_r_cond_var_str[[r]])) {
      temp_output[["cond"]][[r]][[comp]] <- pre_r_cond_var_str[[r]][[comp]]
    }
  }
  Output <- temp_output
  return(Output)
}




# Derive Conditional Variances from Conditional and Marginal Post-R Variances
### mean_parms = Output of Calculate_Means_DTR() function
### cond_var_r = Output of Obtain_Spec_Post_R_Var_Parms() function
###     (corresponding to conditional variances for responders)
### cond_var_nr = Output of Obtain_Spec_Post_R_Var_Parms() function
###     (corresponding to conditional variances for nonresponders)
### p_r = Probability of response in the DTR
Derive_Cond_From_CM_Var_Post_R_DTR <- function(mean_parms, marg_var, cond_var_rs,
                                               p_rs) {
  c_means <- mean_parms[["cond"]]
  Output <- NULL
  for (t in c(0,1,2)) {
    if (t!=2) {
      rho_lab <- paste0("rho_", t, "2")
      phi_lab <- paste0("phi_", t, "2")
    } else {
      rho_lab <- "rho_2"
      phi_lab <- "sigma2_2"
    }
    delta_lab <- paste0("mean_", t)
    mean_term <- 
      (c_means[["NR"]][[delta_lab]] - c_means[["R"]][[delta_lab]])*
      (c_means[["NR"]][["mean_2"]] - c_means[["R"]][["mean_2"]])*
      ( p_rs*(1-p_rs) )
    
    Output[[rho_lab]] <- (marg_var[[rho_lab]] - 
                            p_rs*cond_var_rs[[rho_lab]] - mean_term)/(1-p_rs)
    Output[[phi_lab]] <- (marg_var[[phi_lab]] - 
                            p_rs*cond_var_rs[[phi_lab]] - mean_term)/(1-p_rs)
  }
  
  return(Output)
}


# Derive Marginal Variances from Conditional Post-R Variances
### mean_parms = Output of Calculate_Means_DTR() function
### cond_var_r = Output of Obtain_Spec_Post_R_Var_Parms() function
###     (corresponding to conditional variances for responders)
### cond_var_nr = Output of Obtain_Spec_Post_R_Var_Parms() function
###     (corresponding to conditional variances for nonresponders)
### p_r = Probability of response in the DTR
Derive_Marg_From_Cond_Var_Post_R_DTR <- function(mean_parms, cond_var_r, 
                                                 cond_var_nr, p_r) {
  c_means <- mean_parms[["cond"]]
  c_vars <- list(R=cond_var_r, NR=cond_var_nr)
  Output <- NULL
  for (t in c(0,1,2)) {
    if (t!=2) {
      rho_lab <- paste0("rho_", t, "2")
      phi_lab <- paste0("phi_", t, "2")
    } else {
      rho_lab <- "rho_2"
      phi_lab <- "sigma2_2"
    }
    delta_lab <- paste0("mean_", t)
    mean_term <- 
      (c_means[["NR"]][[delta_lab]] - c_means[["R"]][[delta_lab]])*
      (c_means[["NR"]][["mean_2"]] - c_means[["R"]][["mean_2"]])*
      ( p_r*(1-p_r) )
    
    Output[[rho_lab]] <- p_r*c_vars[["R"]][[rho_lab]] + 
      (1-p_r)*c_vars[["NR"]][[rho_lab]] + mean_term
    Output[[phi_lab]] <- p_r*c_vars[["R"]][[phi_lab]] + 
      (1-p_r)*c_vars[["NR"]][[phi_lab]] + mean_term
  }
  
  return(Output)
}


# Consolidate Lists (Extract Desired Elements)
### lists = List of lists to search from
### elements_to_store = List elements to store
Consolidate_Lists <- function(lists, elements_to_store) {
  temp_out <- list(stored=NULL, other=NULL)
  for (i in (1:length(lists))) {
    temp_list <- lists[[i]]
    for (name in names(temp_list)) {
      if (name %in% elements_to_store) {
        temp_out[["stored"]][[name]] <- temp_list[[name]]
      } else {temp_out[["other"]][[name]] <- temp_list[[name]]}
    }
  }
  Output <- temp_out
  return(Output)
}


# Obtain Specified Post Response Variance Parameters
### post_r_setting = Setting of post-response variance parameters to consider
### post_r_settings_parms = Dataframe of post-response variance parameter settings
Obtain_Spec_Post_R_Var_Parms <- function(post_r_setting, post_r_settings_parms){
  post_r_cond_var_str <- 
    post_r_settings_parms[post_r_settings_parms$cond_param_setting_post_r==post_r_setting,]
  temp_output <- NULL
  for (comp in colnames(post_r_cond_var_str)) {
    temp_output[[comp]] <- post_r_cond_var_str[[comp]]
  }
  Output <- temp_output
  return(Output)
}



# Derive Variance Parameters for a Single DTR
### cp_settings = Dataframe of settings of conditional parameters
###               *** Restricted to simulation setting under consideration
### dtr = DTR label in question
### pre_r_marg_parms = Dataframe of marginal pre-response variance parameter settings
### pre_r_cond_parms_n_i = List of conditional pre-response variance parameter settings
###                        (e.g., those derived via MC) - MUST BE RESTRICTED TO SAMPLE SIZE
###                        [["R"]]- Conditioned on R=1
###                        [["NR"]]- Conditioned on R=0
### post_r_cond_parms = Dataframe of conditional post-response variance parameter settings
### mean_parms = Output of Calculate_Means_DTR()
Derive_Var_Parms_DTR <- function(cp_settings, dtr, pre_r_marg_parms, pre_r_cond_parms_n_i,
                                 post_r_cond_parms, mean_parms) {
  SMART_structure_df <- 
    Define_SMART_Pathways(SMART_structure=cp_settings$SMART_structure[1])
  DTR_structure <- SMART_structure_df[SMART_structure_df$dtr==dtr, ]
  path_r <- DTR_structure$path_r
  path_nr <- DTR_structure$path_nr
  var_spec <- cp_settings$var_spec[1]
  
  # Grab variance specifications
  if (var_spec=="cond") {
    # Pre-response variance should not depend on R/NR path
    cp_settings_r <- cp_settings[cp_settings$pathway==path_r, ]
    cp_settings_nr <- cp_settings[cp_settings$pathway==path_nr, ]
    pre_r_var_str <- cp_settings_r$pre_r_var_str
    post_r_var_str_r <- cp_settings_r$post_r_var_str
    post_r_var_str_nr <- cp_settings_nr$post_r_var_str
  } else if (var_spec %in% c("marg_r", "marg_nr")) {
    # Grab the marginal and specified conditional distribution settings
    cp_settings_marg <- 
      cp_settings[(!is.na(cp_settings$dtr))&cp_settings$dtr==dtr, ]
    if (var_spec=="marg_r") {path_rs <- path_r
    } else if (var_spec=="marg_nr") {path_rs <- path_nr}
    cp_settings_rs <- 
      cp_settings[(!is.na(cp_settings$pathway))&cp_settings$pathway==path_rs, ]
    pre_r_var_str <- cp_settings_marg$pre_r_var_str
    post_r_var_str_marg <- cp_settings_marg$post_r_var_str
    post_r_var_str_rs <- cp_settings_rs$post_r_var_str
  }
  pre_r_var_parms <- Obtain_Pre_R_Var_Parms(pre_r_setting=pre_r_var_str, 
                                            pre_r_marg_parms=pre_r_marg_parms, 
                                            pre_r_cond_parms_n_i=pre_r_cond_parms_n_i)
  
  # Grab post-response distributions
  if (var_spec=="marg_r") {p_rs <- pre_r_var_parms[["marg"]][["p_r"]]
  } else if (var_spec=="marg_nr") {p_rs <- 1-pre_r_var_parms[["marg"]][["p_r"]]}
  
  if (var_spec=="cond") {
    post_r_var_parms_r <- 
      Obtain_Spec_Post_R_Var_Parms(post_r_setting=post_r_var_str_r,
                                   post_r_settings_parms=post_r_cond_parms)
    post_r_var_parms_nr <- 
      Obtain_Spec_Post_R_Var_Parms(post_r_setting=post_r_var_str_nr,
                                   post_r_settings_parms=post_r_cond_parms)
    post_r_var_parms_marg <- 
      Derive_Marg_From_Cond_Var_Post_R_DTR(mean_parms=mean_parms, 
                                           cond_var_r=post_r_var_parms_r, 
                                           cond_var_nr=post_r_var_parms_nr, 
                                           p_r=pre_r_var_parms[["marg"]][["p_r"]])
  } else if (var_spec=="marg_r") {
    post_r_var_parms_r <- 
      Obtain_Spec_Post_R_Var_Parms(post_r_setting=post_r_var_str_rs,
                                   post_r_settings_parms=post_r_cond_parms)
    post_r_var_parms_marg <- 
      Obtain_Spec_Post_R_Var_Parms(post_r_setting=post_r_var_str_marg,
                                   post_r_settings_parms=post_r_cond_parms)
    post_r_var_parms_nr <- 
      Derive_Cond_From_CM_Var_Post_R_DTR(mean_parms=mean_parms, 
                                         marg_var=post_r_var_parms_marg, 
                                         cond_var_rs=post_r_var_parms_r,
                                         p_rs=p_rs)
    
  } else if (var_spec=="marg_nr") {
    post_r_var_parms_nr <- 
      Obtain_Spec_Post_R_Var_Parms(post_r_setting=post_r_var_str_rs,
                                   post_r_settings_parms=post_r_cond_parms)
    post_r_var_parms_marg <- 
      Obtain_Spec_Post_R_Var_Parms(post_r_setting=post_r_var_str_marg,
                                   post_r_settings_parms=post_r_cond_parms)
    post_r_var_parms_r <- 
      Derive_Cond_From_CM_Var_Post_R_DTR(mean_parms=mean_parms, 
                                         marg_var=post_r_var_parms_marg, 
                                         cond_var_rs=post_r_var_parms_nr,
                                         p_rs=p_rs)
  }
  # Consolidate Variance Parameters
  pre_r_parm_names <- c("sigma2_0", "sigma2_1", "rho_0", "rho_1", "rho_01", "phi_01")
  post_r_parm_names <- c("sigma2_2", "phi_02", "phi_12", "rho_02", "rho_12", "rho_2")
  
  
  temp_out <- NULL
  all_marg_parms <- 
    Consolidate_Lists(lists=list(pre_r_var_parms[["marg"]], post_r_var_parms_marg),
                      elements_to_store=c(pre_r_parm_names, post_r_parm_names))
  all_cond_parms_r <- 
    Consolidate_Lists(lists=list(pre_r_var_parms[["cond"]][["R"]], post_r_var_parms_r),
                      elements_to_store=c(pre_r_parm_names, post_r_parm_names))
  all_cond_parms_nr <- 
    Consolidate_Lists(lists=list(pre_r_var_parms[["cond"]][["NR"]], post_r_var_parms_nr),
                      elements_to_store=c(pre_r_parm_names, post_r_parm_names))
  
  temp_out[["marg"]] <- all_marg_parms[["stored"]]
  temp_out[["cond"]] <- list(R=all_cond_parms_r[["stored"]], 
                             NR=all_cond_parms_nr[["stored"]])
  temp_out[["other_parms"]] <- list(
    marg=all_marg_parms[["other"]],
    cond=list(R=all_cond_parms_r[["other"]],
              NR=all_cond_parms_nr[["other"]])
  )
  
  Output <- temp_out
  return(Output)
  
}



# Derive Key Distribution Details for a Given DTR
### cp_settings = Dataframe of settings of conditional parameters
###               *** Restricted to simulation setting under consideration
### dtr = DTR label in question
### pre_r_marg_parms = Dataframe of marginal pre-response variance parameter settings
### pre_r_cond_parms_n_i = List of conditional pre-response variance parameter settings
###                        (e.g., those derived via MC) - MUST BE RESTRICTED TO SAMPLE SIZE
###                        [["R"]]- Conditioned on R=1
###                        [["NR"]]- Conditioned on R=0
### post_r_cond_parms = Dataframe of conditional post-response variance parameter settings
Derive_Dist_DTR <- function(cp_settings, dtr, pre_r_marg_parms, pre_r_cond_parms_n_i,
                            post_r_cond_parms) {
  dtr_means <- 
    Calculate_Means_DTR(cp_settings=cp_settings, 
                        dtr=dtr, 
                        pre_r_marg_parms=pre_r_marg_parms,
                        pre_r_cond_parms_n_i=pre_r_cond_parms_n_i)
  dtr_var <- Derive_Var_Parms_DTR(cp_settings=cp_settings, 
                                  dtr=dtr, 
                                  pre_r_marg_parms=pre_r_marg_parms, 
                                  pre_r_cond_parms_n_i=pre_r_cond_parms_n_i,
                                  post_r_cond_parms=post_r_cond_parms, 
                                  mean_parms=dtr_means)
  Output <- list(
    marg=list(mean=dtr_means[["marg"]], var=dtr_var[["marg"]]),
    cond=list(R=list(mean=dtr_means[["cond"]][["R"]], 
                     var=dtr_var[["cond"]][["R"]]),
              NR=list(mean=dtr_means[["cond"]][["NR"]], 
                      var=dtr_var[["cond"]][["NR"]])),
    other_parms=dtr_var[["other_parms"]]
  )
  return(Output)
}


# Calculate P_1
### dtr_dist = Output of Derive_Dist_DTR() function
Calculate_P_1 <- function(dtr_dist){
  rho_01 <- dtr_dist[["marg"]][["var"]][["rho_01"]]
  rho_0 <- dtr_dist[["marg"]][["var"]][["rho_0"]]
  Output <- ifelse(rho_01==0, 0, rho_01/rho_0)
  return(Output)
}


# Make Upsilon Matrix
### dtr_dist = Output of Derive_Dist_DTR() function
### R_status = Response status ("R"/"NR")
Make_Upsilon_Matrix <- function(dtr_dist, R_status){
  #Extract Parameters
  P_1 <- Calculate_P_1(dtr_dist=dtr_dist)
  cond_parms <- dtr_dist[["other_parms"]][["cond"]][[R_status]]
  
  n_i <- cond_parms[["n_i"]]
  s2_0 <- cond_parms[["s2_0"]]
  s2_1 <- cond_parms[["s2_1"]]
  s_01 <- cond_parms[["s_01"]]
  c_0 <- cond_parms[["c_0"]]
  c_1 <- cond_parms[["c_1"]]
  c_01 <- cond_parms[["c_01"]]
  
  
  #Define matrix components
  p_11 <- s2_0
  p_12 <- (s2_0/n_i) + ((n_i-1)/n_i)*c_0
  p_13 <- P_1*s2_0 + s_01
  p_14 <- (s_01/n_i) + ((n_i-1)/n_i)*c_01
  
  p_21 <- P_1*s2_0 + s_01
  p_22 <- ((P_1*s2_0 + s_01)/n_i) + ((n_i-1)/n_i)*(P_1*c_0 + c_01)
  p_23 <- ((P_1^2)*s2_0) + 2*P_1*s_01 + s2_1
  p_24 <- ((P_1*s_01 + s2_1)/n_i) + ((n_i-1)/n_i)*(P_1*c_01 + c_1)
  
  p_31 <- c_0
  p_32 <- (s2_0/n_i) + ((n_i-1)/n_i)*c_0
  p_33 <- P_1*c_0 + c_01
  p_34 <- (s_01/n_i) + ((n_i-1)/n_i)*c_01
  
  p_41 <- P_1*c_0 + c_01
  p_42 <- ((P_1*s2_0 + s_01)/n_i) + ((n_i-1)/n_i)*(P_1*c_0 + c_01)
  p_43 <- ((P_1^2)*c_0) + 2*P_1*c_01 + c_1
  p_44 <- ((P_1*s_01 + s2_1)/n_i) + ((n_i-1)/n_i)*(P_1*c_01 + c_1)
  
  Output <- matrix(c(p_11, p_12, p_13, p_14,
                     p_21, p_22, p_23, p_24,
                     p_31, p_32, p_33, p_34,
                     p_41, p_42, p_43, p_44),
                   nrow=4, ncol=4, byrow=TRUE)
  return(Output)
}


# Derive the vector of P-coefficients for responders and nonresponders
### dtr_dist = Output of Derive_Dist_DTR() function
Derive_P_Vectors <- function(dtr_dist){
  P_vecs <- NULL
  P_1 <- Calculate_P_1(dtr_dist=dtr_dist)
  for (r in c("R","NR")) {
    Upsilon <- Make_Upsilon_Matrix(dtr_dist=dtr_dist, 
                                   R_status=r)
    cond_var_parms <- dtr_dist[["cond"]][[r]][["var"]]
    phi_02 <- cond_var_parms[["phi_02"]]
    phi_12 <- cond_var_parms[["phi_12"]]
    rho_02 <- cond_var_parms[["rho_02"]]
    rho_12 <- cond_var_parms[["rho_12"]]
    cond_parm_vec <- as.matrix(c(phi_02, phi_12, rho_02, rho_12), ncol=1)
    
    new_Ps <- solve(Upsilon) %*% cond_parm_vec
    P_vecs[[r]] <- c(P_1)
    for (row in (1:nrow(new_Ps))) {
      P_vecs[[r]] <- c(P_vecs[[r]], new_Ps[row, 1])
    }
  }
  Output <- P_vecs
  return(Output)
}


# Derive the zeta coefficient for a given response status
### dtr_dist = Output of Derive_Dist_DTR() function
### P_vecs = P vectors, output of Derive_P_Vectors() function
### R_status = Response status ("R"/"NR")
Derive_Zeta <- function(dtr_dist, P_vecs, R_status){
  cond_parms <- dtr_dist[["other_parms"]][["cond"]][[R_status]]
  P <- P_vecs[[R_status]]
  
  n_i <- cond_parms[["n_i"]]
  s2_0 <- cond_parms[["s2_0"]]
  s2_1 <- cond_parms[["s2_1"]]
  s_01 <- cond_parms[["s_01"]]
  c_0 <- cond_parms[["c_0"]]
  c_1 <- cond_parms[["c_1"]]
  c_01 <- cond_parms[["c_01"]]
  
  P1 <- P[1]
  P2 <- P[2]
  P3 <- P[3]
  P4 <- P[4]
  P5 <- P[5]
  
  #Breaking up sum into a couple different lines
  line_1 <- ( (P2 + (P3/n_i) + P4*P1)^2 )*s2_0 + 
    ( (P3/n_i)^2 )*(n_i-1)*( s2_0 + ((n_i-2)*c_0) )
  line_2 <- ( (P4+(P5/n_i))^2 )*s2_1 + 
    ( (P5/n_i)^2 )*(n_i-1)*(s2_1 + (n_i-2)*c_1)
  line_3 <- 2*(P2 + (P3/n_i) + P4*P1)*( 
    (P3/n_i)*(n_i-1)*c_0 + (P4+(P5/n_i))*s_01 + (P5/n_i)*(n_i-1)*c_01
  )
  line_4 <- 2*(P3/n_i)*(n_i-1)*(
    (P4+(P5/n_i))*c_01 + (P5/n_i)*( s_01 + (n_i-2)*c_01 )
  )
  line_5 <- 2*(P4+(P5/n_i))*(P5/n_i)*(n_i-1)*c_1
  
  Output <- line_1 + line_2 + line_3 + line_4 + line_5
  return(Output)
}


# Derive the zeta prime coefficient for a given response status
### dtr_dist = Output of Derive_Dist_DTR() function
### P_vecs = P vectors, output of Derive_P_Vectors() function
### R_status = Response status ("R"/"NR")
Derive_Zeta_Prime <- function(dtr_dist, P_vecs, R_status){
  cond_parms <- dtr_dist[["other_parms"]][["cond"]][[R_status]]
  P <- P_vecs[[R_status]]
  
  n_i <- cond_parms[["n_i"]]
  s2_0 <- cond_parms[["s2_0"]]
  s2_1 <- cond_parms[["s2_1"]]
  s_01 <- cond_parms[["s_01"]]
  c_0 <- cond_parms[["c_0"]]
  c_1 <- cond_parms[["c_1"]]
  c_01 <- cond_parms[["c_01"]]
  
  P1 <- P[1]
  P2 <- P[2]
  P3 <- P[3]
  P4 <- P[4]
  P5 <- P[5]
  
  #Breaking up sum into a couple different lines
  line_1 <- (P2 + (P3/n_i) + P4*P1)*(
    (P2 + (P3/n_i)*(n_i-1) + P4*P1)*c_0 + (P3/n_i)*s2_0
  )
  line_2 <- (P2 + (P3/n_i) + P4*P1)*(
    (P4 + (P5/n_i)*(n_i-1))*c_01 + (P5/n_i)*s_01
  )
  line_3 <- (P3/n_i)*(P2 + (P3/n_i) + P4*P1)*s2_0 + ( (P3/n_i)^2 )*(n_i-1)*c_0
  line_4 <- (P3/n_i)*(P4 + (P5/n_i))*s_01 + (P3*P5/(n_i^2))*(n_i-1)*c_01
  line_5 <- (P3/n_i)*(n_i-2)*(P2 + (2*P3/n_i) + P4*P1)*c_0 + 
    ( (P3/n_i)^2 )*(n_i-2)*( s2_0 + (n_i-3)*c_0 )
  line_6 <- (P3/n_i)*(n_i-2)*(P4+(2*P5/n_i))*c_01 + 
    (P3*P5/(n_i^2))*(n_i-2)*( s_01 + (n_i-3)*c_01 )
  line_7 <- ( P4 + (P5/n_i) )*(
    ( (P2 + (P3/n_i)*(n_i-1) + P4*P1)*c_01 + (P3/n_i)*s_01)
  )
  line_8 <- ( P4 + (P5/n_i) )*(
    (P4 + (P5/n_i)*(n_i-1))*c_1 + (P5/n_i)*s2_1
  )
  line_9 <- (P5/n_i)*(P2 + (P3/n_i) + P4*P1)*s_01 + 
    (P3*P5/(n_i^2))*(n_i-1)*c_01
  line_10 <- (P5/n_i)*(P4 + (P5/n_i))*s2_1 + 
    ((P5/n_i)^2)*(n_i-1)*c_1
  line_11 <- (P5/n_i)*(n_i-2)*(P2 + (2*P3/n_i) + P4*P1)*c_01 + 
    (P3*P5/(n_i^2))*(n_i-2)*( s_01 + (n_i-3)*c_01 )
  line_12 <- (P5/n_i)*(n_i-2)*(P4+(2*P5/n_i))*c_1 + 
    ( (P5/n_i)^2 )*(n_i-2)*( s2_1 + (n_i-3)*c_1 )
  
  
  Output <- line_1 + line_2 + line_3 + line_4 + line_5 + line_6 + 
    line_7 + line_8 + line_9 + line_10 + line_11 + line_12
  return(Output)
}





# Obtain Y2 Error Variances for a Given DTR
### dtr_dist = Output of Derive_Dist_DTR() function
### P_vecs = P vectors, output of Derive_P_Vectors() function
Obtain_Sigma_2_Mats <- function(dtr_dist, P_vecs){
  Sigma2_Matrices <- NULL
  for (r in c("R","NR")) {
    zeta <- Derive_Zeta(dtr_dist=dtr_dist, 
                        P_vecs=P_vecs, 
                        R_status=r)
    zeta_prime <- Derive_Zeta_Prime(dtr_dist=dtr_dist, 
                                    P_vecs=P_vecs, 
                                    R_status=r)
    d_term <- dtr_dist[["cond"]][[r]][["var"]][["sigma2_2"]] - zeta
    off_d_term <- dtr_dist[["cond"]][[r]][["var"]][["rho_2"]] - zeta_prime
    n_i <- dtr_dist[["other_parms"]][["cond"]][[r]][["n_i"]]
    
    Sigma2_Matrices[[r]] <- diag(x=(d_term-off_d_term), nrow=n_i) + 
      matrix(data=off_d_term, nrow=n_i, ncol=n_i)
  }
  Output <- Sigma2_Matrices
  return(Output)
}


# Obtain Pre-Response Error Variance Matrix
### dtr_dist = Output of Derive_Dist_DTR() function
Obtain_Sigma_01_Mat <- function(dtr_dist){
  P_1 <- Calculate_P_1(dtr_dist=dtr_dist)
  parms <- dtr_dist[["marg"]][["var"]]
  n_i <- dtr_dist[["other_parms"]][["marg"]][["n_i"]]
  
  V1_ni <- as.matrix(integer(n_i)+1, ncol=1)
  
  Sigma_0 <- (parms[["sigma2_0"]]-parms[["rho_0"]])*diag(n_i) + 
    parms[["rho_0"]]*V1_ni%*%t(V1_ni)
  
  Sigma_1 <- (parms[["sigma2_1"]]-parms[["rho_1"]])*diag(n_i) + 
    parms[["rho_1"]]*V1_ni%*%t(V1_ni) - (P_1^2)*Sigma_0 -
    2*P_1*(parms[["phi_01"]]-(P_1*parms[["sigma2_0"]]))*diag(n_i)
  
  off_diag <- (parms[["phi_01"]] - P_1*parms[["sigma2_0"]])*diag(n_i)
  
  top_block <- cbind(Sigma_0, off_diag)
  bottom_block <- cbind(off_diag, Sigma_1)
  
  Output <- as.matrix(rbind(top_block, bottom_block))
  return(Output)
}



# Derive Adjustments to Match Conditional Mean
### cp_settings = Dataframe of settings of conditional parameters
###               *** Restricted to simulation setting under consideration
### dtr_dist = Output of Derive_Dist_DTR() function
### t = Time for mean of which we're taking adjustment
### Response = Response of the corresponding pathway
Derive_Cond_Mean_Adj <- function(cp_settings, dtr_dist, t, response) {
  if (cp_settings$SMART_structure[1]=="prototypical") {
    if (t %in% c(0,1)) {Output <- 0
    } else {
      P <- Derive_P_Vectors(dtr_dist=dtr_dist)[[response]]
      e_0 <- dtr_dist[["other_parms"]][["cond"]][[response]][["mean_e_0"]]
      e_1 <- dtr_dist[["other_parms"]][["cond"]][[response]][["mean_e_1"]]
      Output <- P[2]*e_0 + P[3]*e_0 + P[4]*(P[1]*e_0 + e_1) + P[5]*e_1
    }
  }
  return(Output)
}


# Calculate Mean Parameters for Prototypical SMART
### cp_settings = Dataframe of settings of conditional parameters
###               *** Restricted to simulation setting under consideration
### pre_r_marg_parms = Dataframe of marginal pre-response variance parameter settings
### pre_r_cond_parms_n_i = List of conditional pre-response variance parameter settings
###                        (e.g., those derived via MC) - MUST BE RESTRICTED TO SAMPLE SIZE
###                        [["R"]]- Conditioned on R=1
###                        [["NR"]]- Conditioned on R=0
Calculate_Mean_Parms_Proto <- function(cp_settings, pre_r_marg_parms,
                                       pre_r_cond_parms_n_i, post_r_cond_parms) {
  SMART_structure <- 
    Define_SMART_Pathways(SMART_structure=cp_settings$SMART_structure[1])
  d_dists <- NULL
  for (dtr in SMART_structure$dtr) {
    d_dists[[dtr]] <- Derive_Dist_DTR(cp_settings=cp_settings, 
                                      dtr=dtr, 
                                      pre_r_marg_parms=pre_r_marg_parms, 
                                      pre_r_cond_parms_n_i=pre_r_cond_parms_n_i,
                                      post_r_cond_parms=post_r_cond_parms)
  }
  # Grab DTRs with given pathways
  d_1_1 <- 
    SMART_structure[(SMART_structure$a_1==1)&(SMART_structure$a_2nr==1), "dtr"][1]
  d_1_m1 <- 
    SMART_structure[(SMART_structure$a_1==1)&(SMART_structure$a_2nr==-1), "dtr"][1]
  d_m1_1 <- 
    SMART_structure[(SMART_structure$a_1==-1)&(SMART_structure$a_2nr==1), "dtr"][1]
  d_m1_m1 <- 
    SMART_structure[(SMART_structure$a_1==-1)&(SMART_structure$a_2nr==-1), "dtr"][1]
  p_r_1 <- d_dists[[d_1_1]][["other_parms"]][["marg"]][["p_r"]]
  p_r_m1 <- d_dists[[d_m1_1]][["other_parms"]][["marg"]][["p_r"]]
  
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
    d_dists[[d_1_1]][["marg"]][["mean"]][["mean_0"]], 
    d_dists[[d_1_1]][["marg"]][["mean"]][["mean_1"]],
    d_dists[[d_m1_1]][["marg"]][["mean"]][["mean_1"]],
    d_dists[[d_1_1]][["cond"]][["R"]][["mean"]][["mean_2"]],
    d_dists[[d_1_1]][["cond"]][["NR"]][["mean"]][["mean_2"]],
    d_dists[[d_1_m1]][["cond"]][["NR"]][["mean"]][["mean_2"]],
    d_dists[[d_m1_1]][["cond"]][["R"]][["mean"]][["mean_2"]],
    d_dists[[d_m1_1]][["cond"]][["NR"]][["mean"]][["mean_2"]],
    d_dists[[d_m1_m1]][["cond"]][["NR"]][["mean"]][["mean_2"]]), ncol=1)
  
  #Need to adjust to ensure conditional mean specifications
  adj_vec <- matrix(data=0, nrow=9, ncol=1)
  for (row in (4:9)) {
    path <- row-3
    response <- ifelse(path %in% c(1,4), "R", "NR")
    dtr <- ifelse(response=="R", 
                  SMART_structure[SMART_structure$path_r==path, "dtr"][1],
                  SMART_structure[SMART_structure$path_nr==path, "dtr"][1])
    adj_vec[row, 1] <- Derive_Cond_Mean_Adj(cp_settings=cp_settings, 
                                            dtr_dist=d_dists[[dtr]], 
                                            t=2, 
                                            response=response)
  }
  gamma_vec <- solve(coef_mat)%*%(mean_vec-adj_vec)
  
  Output <- NULL
  for (i in (0:6)) {Output[[paste0("gamma_",i)]] <- gamma_vec[i+1,1]}
  Output[["lambda_1"]] <- gamma_vec[8,1]
  Output[["lambda_2"]] <- gamma_vec[9,1]
  return(Output)
}


# Calculate All Mean Parameters for a SMART
### cp_settings = Dataframe of settings of conditional parameters
###               *** Restricted to simulation setting under consideration
### pre_r_marg_parms = Dataframe of marginal pre-response variance parameter settings
### pre_r_cond_var_parms_n_i = List of conditional pre-response variance parameter settings
###                        (e.g., those derived via MC) - MUST BE RESTRICTED TO SAMPLE SIZE
###                        [["R"]]- Conditioned on R=1
###                        [["NR"]]- Conditioned on R=0
Calculate_Mean_Parms <- function(cp_settings, pre_r_marg_parms, 
                                 pre_r_cond_parms_n_i, post_r_cond_parms) {
  SMART_structure <- cp_settings$SMART_structure[1]
  if (SMART_structure=="prototypical") {
    temp_output <- 
      Calculate_Mean_Parms_Proto(cp_settings=cp_settings, 
                                 pre_r_marg_parms=pre_r_marg_parms,
                                 pre_r_cond_parms_n_i=pre_r_cond_parms_n_i,
                                 post_r_cond_parms=post_r_cond_parms)
    
  }
  Output <- temp_output
  return(Output)
}



# Calculate the Static Parameters for Outcome Generation for an Embedded DTR
#   in a Prototypical SMART
### dtr_dist = Output of Derive_Dist_DTR() function
### mean_sim_parms = Output of Calculate_Mean_Parms() function
### P_vecs = Output of Derive_P_Vectors() function
### dtr = DTR under consideration
Calculate_Mean_Static_Helper_Parms_DTR_Proto <- function(dtr_dist, mean_sim_parms, 
                                                         P_vecs, dtr) {
  DTR_structure <- 
    Define_SMART_Pathways("prototypical")[Define_SMART_Pathways("prototypical")$dtr==dtr, ]
  a_1 <- DTR_structure$a_1
  a_2nr <- DTR_structure$a_2nr
  p_r <- dtr_dist[["other_parms"]][["marg"]][["p_r"]]
  gamma <- NULL
  lambda <- NULL
  for (i in (0:6)) {gamma[[toString(i)]] <- mean_sim_parms[[paste0("gamma_", i)]]} 
  for (i in (1:2)) {lambda[[toString(i)]] <- mean_sim_parms[[paste0("lambda_", i)]]}
  
  t_0_term <- gamma[["0"]]
  t_1_term <- (1-P_vecs[["R"]][1])*gamma[["0"]] + gamma[["1"]] + gamma[["2"]]*a_1
  t_2_term <- NULL
  for (r in c("R", "NR")) {
    P <- P_vecs[[r]]
    r_i <- ifelse(r=="R", 1, 0)
    t_2_term[[r]] <- (1 - P[2] - P[4])*gamma[["0"]] + 
      (1 - P[4])*(gamma[["1"]] + gamma[["2"]]*a_1) + 
      (gamma[["3"]] + gamma[["4"]]*a_1) + 
      #P[2]*outcomes$y_0 + 
      #P[3]*mean(outcomes$e_0) + 
      #P[4]*outcomes$y_1 + 
      #P[5]*mean(outcomes$e_1) + 
      ( (gamma[["5"]] + gamma[["6"]]*a_1)*a_2nr )*( (1-r_i)/(1-p_r) ) + 
      (lambda[["1"]] + lambda[["2"]]*a_1)*(r_i-p_r) # + 
    #e_2
  }
  Output <- list(
    t_0=t_0_term,
    t_1=t_1_term,
    t_2=t_2_term
  )
  return(Output)
}

# Calculate the Static Parameters for Outcome Generation for an Embedded DTR
#   in a Prototypical SMART
### dtr_dist = Output of Derive_Dist_DTR() function
### mean_sim_parms = Output of Calculate_Mean_Parms() function
### P_vecs = Output of Derive_P_Vectors() function
### dtr = DTR under consideration
Calculate_Mean_Static_Helper_Parms_DTR <- function(dtr_dist, mean_sim_parms, 
                                                   P_vecs, dtr, SMART_structure) {
  if (SMART_structure=="prototypical") {
    temp_output <- 
      Calculate_Mean_Static_Helper_Parms_DTR_Proto(dtr_dist=dtr_dist, 
                                                   mean_sim_parms=mean_sim_parms,
                                                   P_vecs=P_vecs, 
                                                   dtr=dtr)
    
  }
  Output <- temp_output
  return(Output)
}


# Generate Response Function for a Given DTR/Sample Size
### dtr_dist = Output of Derive_Dist_DTR() function
### pre_r_marg_parms = Dataframe of marginal pre-response variance parameter settings
Generate_Response_Function_DTR_n <- function(dtr_dist, pre_r_marg_parms) {
  # Use code from Pre_R_Cond_Param_Sim.R file - Uses pre_r_marg_parms as driver
  # Builds response function taking mean outcome deviations as sole input
  pre_r_var_setting <- 
    dtr_dist[["other_parms"]][["marg"]][["cond_param_setting_pre_r"]]
  
  # Isolate relevant driver row
  rel_driver_row <- 
    pre_r_marg_parms[(pre_r_marg_parms$cond_param_setting_pre_r==pre_r_var_setting), ]
  # Overwrite sample size to that under consideration
  rel_driver_row$sample_size <- dtr_dist[["other_parms"]][["marg"]][["n_i"]]
  
  pre_r_var_parms <- Process_Driver_Row_Pre_R_MC(driver_row=rel_driver_row)
  Output <- Build_Single_Input_R_Function(parms=pre_r_var_parms)
  return(Output)
}



# Derive DTR-Specific Data for Given DTR and Sample Size
### cp_settings = Dataframe of settings of conditional parameters
###               *** Restricted to simulation setting under consideration
### dtr = DTR label in question
### pre_r_marg_parms = Dataframe of marginal pre-response variance parameter settings
### pre_r_cond_parms_n_i = List of conditional pre-response variance parameter settings
###                        (e.g., those derived via MC) - MUST BE RESTRICTED TO SAMPLE SIZE
###                        [["R"]]- Conditioned on R=1
###                        [["NR"]]- Conditioned on R=0
### post_r_cond_parms = Dataframe of conditional post-response variance parameter settings
Derive_DTR_Data_n <- function(cp_settings, dtr, pre_r_marg_parms, 
                              pre_r_cond_parms_n_i, post_r_cond_parms) {
  dist_parms <- Derive_Dist_DTR(cp_settings=cp_settings, 
                                dtr=dtr,
                                pre_r_marg_parms=pre_r_marg_parms,
                                pre_r_cond_parms_n_i=pre_r_cond_parms_n_i,
                                post_r_cond_parms=post_r_cond_parms)
  P_vecs <- Derive_P_Vectors(dtr_dist=dist_parms)
  Sigma_01 <- Obtain_Sigma_01_Mat(dtr_dist=dist_parms)
  Sigma_2 <- Obtain_Sigma_2_Mats(dtr_dist=dist_parms, P_vecs=P_vecs)
  mean_sim_parms <- 
    Calculate_Mean_Parms(cp_settings=cp_settings, 
                         pre_r_marg_parms=pre_r_marg_parms, 
                         pre_r_cond_parms_n_i=pre_r_cond_parms_n_i,
                         post_r_cond_parms=post_r_cond_parms)
  mean_helper_parms <- 
    Calculate_Mean_Static_Helper_Parms_DTR(dtr_dist=dist_parms, 
                                           mean_sim_parms=mean_sim_parms, 
                                           P_vecs=P_vecs, 
                                           dtr=dtr, 
                                           SMART_structure=cp_settings$SMART_structure[1])
  
  
  response_function <- 
    Generate_Response_Function_DTR_n(dtr_dist=dist_parms, 
                                     pre_r_marg_parms=pre_r_marg_parms)
  SMART_structure <- 
    Define_SMART_Pathways(SMART_structure=cp_settings$SMART_structure[1])
  DTR_structure <- as.list(SMART_structure[SMART_structure$dtr==dtr, ])
  Output <- list(
    dist=dist_parms,
    sim_parms=list(
      P=P_vecs,
      Sigma_01=Sigma_01,
      Sigma_2=Sigma_2,
      mean_parms=mean_sim_parms,
      mean_helper_parms=mean_helper_parms
    ),
    response_function=response_function,
    dtr=DTR_structure
  )
  return(Output)
}


# Derive Simulation Parameters for a Given Driver Row
### driver_parms = Output of Process_Driver_Row_Main() function
### cp_settings_df = Dataframe of settings of conditional parameters
### pre_r_marg_parms = Dataframe of marginal pre-response variance parameter settings
### pre_r_cond_parms = List of conditional pre-response variance parameter settings
###                        (e.g., those derived via MC)
###                        [["R"]]- Conditioned on R=1
###                        [["NR"]]- Conditioned on R=0
### post_r_var_parms = Dataframe of post-response variance parameter settings
Derive_Sim_Parms <- function(driver_parms, cp_settings_df, pre_r_marg_parms, 
                             pre_r_cond_parms, post_r_var_parms) {
  cp_settings <- 
    cp_settings_df[cp_settings_df$cond_parm_setting==driver_parms[["cond_parm_setting"]], ]
  embedded_dtrs <- 
    unique(Define_SMART_Pathways(SMART_structure=driver_parms[["SMART_structure"]])$dtr)
  temp_output <- NULL
  for (n_i in driver_parms[["cluster_sizes"]]) {
    temp_output[[paste0("n_", n_i)]] <- NULL
    pre_r_cond_parms_n_i <- list(
      R=pre_r_cond_parms[["R"]][pre_r_cond_parms[["R"]]$n_i==n_i, ],
      NR=pre_r_cond_parms[["NR"]][pre_r_cond_parms[["NR"]]$n_i==n_i, ]
    )
    for (dtr in embedded_dtrs) {
      temp_output[[paste0("n_", n_i)]][[dtr]] <- 
        Derive_DTR_Data_n(cp_settings=cp_settings, 
                          dtr=dtr,
                          pre_r_marg_parms=pre_r_marg_parms,
                          pre_r_cond_parms_n_i=pre_r_cond_parms_n_i,
                          post_r_cond_parms=post_r_var_parms)
    }
  }
  Output <- temp_output
  return(Output)
}



