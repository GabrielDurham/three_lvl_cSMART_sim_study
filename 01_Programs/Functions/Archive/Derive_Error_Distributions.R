##################################
##################################
##################################
### DERIVE ERROR DISTRIBUTIONS ###
##################################
##################################
##################################

#### PURPOSE: This file contains code to run calculate the error distributions
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


#### DATE CREATED:  28 DEC 2023
#### PROGRAMMER:    GABRIEL DURHAM (GJD)
#### EDITS:         03 JUN 2024 (GJD) - Read in error conditional variance 
####                  parameters rather than outcome conditional variance 
####                  parameters (for Upsilon, Zeta, and Zeta prime definition)


# Obtain Pre-R Parameters
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
Obtain_Pre_R_Params <- function(pre_r_setting, pre_r_marg_parms, pre_r_cond_var_parms_n_i){
  
  pre_r_marg_var_str <- 
    pre_r_marg_parms[pre_r_marg_parms$cond_param_setting_pre_r==pre_r_setting,]
  pre_r_cond_var_str <- NULL
  for (r in c("R","NR")) {
    pre_r_cond_var_str[[r]] <- 
      pre_r_cond_var_parms_n_i[[r]][pre_r_cond_var_parms_n_i[[r]]$cond_param_setting_pre_r==
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



# Obtain Post Response Conditional Parameters
### post_r_settings = List ("R"/"NR") of desired post-response variance settings
### post_r_settings_parms = Dataframe of conditional post-response variance parameter settings
#### Returns a list ("R"/"NR") of lists of variance parameters
Obtain_Post_R_Cond_Params <- function(post_r_settings, post_r_settings_parms){
  post_r_cond_params <- NULL
  for (r in c("R", "NR")) {
    post_r_setting <- post_r_settings[[r]]
    post_r_cond_var_str <- 
      post_r_settings_parms[post_r_settings_parms$cond_param_setting_post_r==post_r_setting,]
    temp_output <- NULL
    for (comp in colnames(post_r_cond_var_str)) {
      temp_output[[comp]] <- post_r_cond_var_str[[comp]]
    }
    post_r_cond_params[[r]] <- temp_output
  }
  Output <- post_r_cond_params
  return(Output)
}  


# Calculate P_1
### pre_r_var_parms = Output of Obtain_Pre_R_Params() function
Calculate_P_1 <- function(pre_r_var_parms){
  rho_01 <- pre_r_var_parms[["marg"]][["rho_01"]]
  rho_0 <- pre_r_var_parms[["marg"]][["rho_0"]]
  Output <- ifelse(rho_01==0, 0, rho_01/rho_0)
  return(Output)
}

# Make Upsilon Matrix
### pre_r_var_parms = Output of Obtain_Pre_R_Params() function
### R_status = Response status ("R"/"NR")
Make_Upsilon_Matrix <- function(pre_r_var_parms, R_status){
  #Extract Parameters
  P_1 <- Calculate_P_1(pre_r_var_parms=pre_r_var_parms)
  cond_parms <- pre_r_var_parms[["cond"]][[R_status]]
  
  n_i <- cond_parms[["n_i"]]
  # s2_0 <- cond_parms[["sigma2_0"]]
  # s2_1 <- cond_parms[["sigma2_1"]]
  # s_01 <- cond_parms[["phi_01"]]
  # c_0 <- cond_parms[["rho_0"]]
  # c_1 <- cond_parms[["rho_1"]]
  # c_01 <- cond_parms[["rho_01"]]
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
### pre_r_var_parms = Output of Obtain_Pre_R_Params() function
### post_r_cond_params = Output of Obtain_Post_R_Cond_Params() function
Derive_P_Vectors <- function(pre_r_var_parms, post_r_cond_params){
  P_vecs <- NULL
  P_1 <- Calculate_P_1(pre_r_var_parms=pre_r_var_parms)
  for (r in c("R","NR")) {
    Upsilon <- Make_Upsilon_Matrix(pre_r_var_parms=pre_r_var_parms, 
                                   R_status=r)
    phi_02 <- post_r_cond_params[[r]][["phi_02"]]
    phi_12 <- post_r_cond_params[[r]][["phi_12"]]
    rho_02 <- post_r_cond_params[[r]][["rho_02"]]
    rho_12 <- post_r_cond_params[[r]][["rho_12"]]
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
### pre_r_var_parms = Output of Obtain_Pre_R_Params() function
### P_vecs = P vectors, output of Derive_P_Vectors() function
### R_status = Response status ("R"/"NR")
Derive_Zeta <- function(pre_r_var_parms, P_vecs, R_status){
  cond_parms <- pre_r_var_parms[["cond"]][[R_status]]
  P <- P_vecs[[R_status]]
  
  n_i <- cond_parms[["n_i"]]
  # s2_0 <- cond_parms[["sigma2_0"]]
  # s2_1 <- cond_parms[["sigma2_1"]]
  # s_01 <- cond_parms[["phi_01"]]
  # c_0 <- cond_parms[["rho_0"]]
  # c_1 <- cond_parms[["rho_1"]]
  # c_01 <- cond_parms[["rho_01"]]
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
### pre_r_var_parms = Output of Obtain_Pre_R_Params() function
### P_vecs = P vectors, output of Derive_P_Vectors() function
### R_status = Response status ("R"/"NR")
Derive_Zeta_Prime <- function(pre_r_var_parms, P_vecs, R_status){
  cond_parms <- pre_r_var_parms[["cond"]][[R_status]]
  P <- P_vecs[[R_status]]
  
  n_i <- cond_parms[["n_i"]]
  # s2_0 <- cond_parms[["sigma2_0"]]
  # s2_1 <- cond_parms[["sigma2_1"]]
  # s_01 <- cond_parms[["phi_01"]]
  # c_0 <- cond_parms[["rho_0"]]
  # c_1 <- cond_parms[["rho_1"]]
  # c_01 <- cond_parms[["rho_01"]]
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
### pre_r_var_parms = Output of Obtain_Pre_R_Params() function
### post_r_cond_params = Output of Obtain_Post_R_Cond_Params() function
### P_vecs = P vectors, output of Derive_P_Vectors() function
Obtain_Sigma_2_Mats <- function(pre_r_var_parms, post_r_cond_params, P_vecs){
  Sigma2_Matrices <- NULL
  for (r in c("R","NR")) {
    zeta <- Derive_Zeta(pre_r_var_parms=pre_r_var_parms, 
                        P_vecs=P_vecs, 
                        R_status=r)
    zeta_prime <- Derive_Zeta_Prime(pre_r_var_parms=pre_r_var_parms, 
                                    P_vecs=P_vecs, 
                                    R_status=r)
    d_term <- post_r_cond_params[[r]][["sigma2_2"]] - zeta
    off_d_term <- post_r_cond_params[[r]][["rho_2"]] - zeta_prime
    n_i <- pre_r_var_parms[["cond"]][[r]][["n_i"]]
    
    Sigma2_Matrices[[r]] <- diag(x=(d_term-off_d_term), nrow=n_i) + 
      matrix(data=off_d_term, nrow=n_i, ncol=n_i)
  }
  Output <- Sigma2_Matrices
  return(Output)
}


# Obtain Pre-Response Error Variance Matrix
### pre_r_var_parms = Output of Obtain_Pre_R_Params() function
Obtain_Sigma_01_Mat <- function(pre_r_var_parms){
  P_1 <- Calculate_P_1(pre_r_var_parms=pre_r_var_parms)
  parms <- pre_r_var_parms[["marg"]]
  
  V1_ni <- as.matrix(integer(parms[["n_i"]])+1, ncol=1)
  
  Sigma_0 <- (parms[["sigma2_0"]]-parms[["rho_0"]])*diag(parms[["n_i"]]) + 
    parms[["rho_0"]]*V1_ni%*%t(V1_ni)
  
  Sigma_1 <- (parms[["sigma2_1"]]-parms[["rho_1"]])*diag(parms[["n_i"]]) + 
    parms[["rho_1"]]*V1_ni%*%t(V1_ni) - (P_1^2)*Sigma_0 -
    2*P_1*(parms[["phi_01"]]-(P_1*parms[["sigma2_0"]]))*diag(parms[["n_i"]])
  
  off_diag <- (parms[["phi_01"]] - P_1*parms[["sigma2_0"]])*diag(parms[["n_i"]])
  
  top_block <- cbind(Sigma_0, off_diag)
  bottom_block <- cbind(off_diag, Sigma_1)
  
  Output <- as.matrix(rbind(top_block, bottom_block))
  return(Output)
}



# Obtain All Parameters for a Given DTR
## Obtain all marg/cond pre-response parameters for a embedded DTR,
## all cond post-response parameters, the error matrices, and P-vectors
### dtr_path_str = Path structure for embedded DTR
###                Row of output of Define_SMART_Pathways() function
### cp_settings = Dataframe of settings of conditional parameters
###               *** Restricted to simulation setting under consideration
### pre_r_marg_parms = Dataframe of marginal pre-response variance parameter settings
### pre_r_cond_var_parms_n_i = List of conditional pre-response variance parameter settings
###                        (e.g., those derived via MC) - NEEDS TO BE LIMITED
###                        TO SAMPLE SIZE UNDER CONSIDERATION
###                        [["R"]]- Conditioned on R=1
###                        [["NR"]]- Conditioned on R=0
### post_r_var_parms = Dataframe of post-response variance parameter settings
Obtain_Params_Single_DTR <- function(dtr_path_str, cp_settings, pre_r_marg_parms, 
                                     pre_r_cond_var_parms_n_i, post_r_var_parms){
  path <- NULL
  path[["R"]] <- dtr_path_str$path_r
  path[["NR"]] <- dtr_path_str$path_nr
  
  #Grab pre-response marginal variance setting
  #Needs to be same for R and NR pathways
  pre_r_marg_var_set <- cp_settings[cp_settings$pathway==path[["R"]], "pre_r_var_str"]
  #Bring in MC-derived conditional parameters
  pre_r_params <- Obtain_Pre_R_Params(pre_r_setting=pre_r_marg_var_set, 
                                      pre_r_marg_parms=pre_r_marg_parms, 
                                      pre_r_cond_var_parms_n_i=pre_r_cond_var_parms_n_i)
  
  #Grab post-response conditional variance settings
  post_r_cond_settings <- NULL
  for (r in c("R", "NR")) {
    post_r_cond_settings[[r]] <- cp_settings[cp_settings$pathway==path[[r]], "post_r_var_str"]
  }
  post_r_cond_params <- Obtain_Post_R_Cond_Params(post_r_settings=post_r_cond_settings,
                                                  post_r_settings_parms=post_r_var_parms)
  
  P_vecs <- Derive_P_Vectors(pre_r_var_parms=pre_r_params,
                             post_r_cond_params=post_r_cond_params)
  
  Sigma_01_Mat <- Obtain_Sigma_01_Mat(pre_r_var_parms=pre_r_params)
  
  Sigma_2_Mats <- Obtain_Sigma_2_Mats(pre_r_var_parms=pre_r_params, 
                                      post_r_cond_params=post_r_cond_params, 
                                      P_vecs=P_vecs)
  
  
  var_parms <- list(
    pre_r=pre_r_params,
    post_r=list(cond=post_r_cond_params)
  )
  coefs <- list(
    Sigma_01=Sigma_01_Mat,
    Sigma_2=Sigma_2_Mats,
    P_vec=P_vecs
  )
  Output <- list(
    var_parms=var_parms,
    coefs=coefs
  )
  
  return(Output)
}


# Obtain All Parameters for a Given Sample Size
## Obtain all marg/cond pre-response parameters for a embedded DTR,
## all cond post-response parameters, the error matrices, and P-vectors
### cp_settings = Dataframe of settings of conditional parameters
###               *** Restricted to simulation setting under consideration
### n_i = Sample size under consideration
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
####          cond: Conditional variance parameters (need mean parameters to derive marg)
####                 (R and NR)
####    coefs = Simulation coefficients
####        Sigma_01: Pre-response Error variance matrix (t=0,1)
####        Sigma_2: Post-response Error variance matrix (t=2) 
####                 (R and NR)
####        P_vec: Vectors of P coefficients
####                 (R and NR)
Obtain_Univ_Sim_Params_n_i <- function(cp_settings, n_i, pre_r_marg_parms,
                                  pre_r_cond_var_parms, post_r_var_parms){
  SMART_structure <- cp_settings$SMART_structure[1]
  SMART_pathways <- Define_SMART_Pathways(SMART_structure=SMART_structure)
  #Isolate conditional variance parameters
  pre_r_cond_var_parms_n_i <- NULL
  for (r in names(pre_r_cond_var_parms)) {
    mc_var_parms <- pre_r_cond_var_parms[[r]]
    pre_r_cond_var_parms_n_i[[r]] <- mc_var_parms[mc_var_parms$n_i==n_i,]
  }
  
  Output <- NULL
  for (d in (1:nrow(SMART_pathways))) {
    dtr_label <- SMART_pathways$dtr[d]
    Output[[dtr_label]] <- Obtain_Params_Single_DTR(dtr_path_str=
                                                      SMART_pathways[d,], 
                                                    cp_settings=
                                                      cp_settings, 
                                                    pre_r_marg_parms=
                                                      pre_r_marg_parms, 
                                                    pre_r_cond_var_parms_n_i=
                                                      pre_r_cond_var_parms_n_i, 
                                                    post_r_var_parms=
                                                      post_r_var_parms)
    Output[[dtr_label]][["dtr"]] <- list(a_1=SMART_pathways$a_1[d],
                                         a_2r=SMART_pathways$a_2r[d],
                                         a_2nr=SMART_pathways$a_2nr[d])
  }
  return(Output)
}



# Obtain All Parameters for a Given Simulation Specification
## Obtain all parameters which don't depend on mean components. I.e.,
## marg/cond pre-response parameters for a embedded DTR,
## all cond post-response parameters, the error matrices, and P-vectors
### cp_settings = Dataframe of settings of conditional parameters
###               *** Restricted to simulation setting under consideration
### n_i_vec = Vector of sample sizes to consider
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
####          cond: Conditional variance parameters (need mean parameters to derive marg)
####                 (R and NR)
####    coefs = Simulation coefficients
####        Sigma_01: Pre-response Error variance matrix (t=0,1)
####        Sigma_2: Post-response Error variance matrix (t=2) 
####                 (R and NR)
####        P_vec: Vectors of P coefficients
####                 (R and NR)
Obtain_Univ_Sim_Params <- function(cp_settings, n_i_vec, pre_r_marg_parms,
                                   pre_r_cond_var_parms, post_r_var_parms){
  Output <- NULL
  for (n in n_i_vec) {
    n_label <- paste0("n_", n)
    Output[[n_label]] <- Obtain_Univ_Sim_Params_n_i(cp_settings=cp_settings, 
                                               n_i=n, 
                                               pre_r_marg_parms=pre_r_marg_parms,
                                               pre_r_cond_var_parms=pre_r_cond_var_parms, 
                                               post_r_var_parms=post_r_var_parms)
  }
  return(Output)
}

