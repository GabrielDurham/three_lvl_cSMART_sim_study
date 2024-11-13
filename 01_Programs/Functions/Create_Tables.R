#####################
#####################
#####################
### CREATE_TABLES ###
#####################
#####################
#####################


#### PURPOSE: This file contains functions for creating tables from inference output

#### DATE CREATED:  10 SEP 2024
#### PROGRAMMER:    GABRIEL DURHAM (GJD)
#### EDITS:         


# Process Output from a Single Fit
### fit_data = Single fit output from Conduct_Inference_for_Driver_Run() output
### filter_statement = Text statement filtering rows of fit_data to relevant rows
Summarize_Single_Fit_Data <- function(fit_data, filter_statement) {
  relevant_data <- 
    eval(parse(text=paste0("fit_data %>% filter(", filter_statement, ")")))
  relevant_data$reject <- 
    ifelse(relevant_data$CI.Lower>0|relevant_data$CI.Higher<0, 1, 0)
  relevant_data$reject_correct <- 
    ifelse(relevant_data$true_value>0,
           ifelse(relevant_data$CI.Lower>0, 1, 0),
           ifelse(relevant_data$CI.Higher<0, 1, 0))
  
  Output <- data.frame(bias=mean(relevant_data$bias),
                       sd=sd(relevant_data$Estimate),
                       RMSE=sqrt(mean(relevant_data$bias^2)),
                       coverage=mean(relevant_data$coverage),
                       power=mean(relevant_data$reject),
                       ci_width=mean(relevant_data$ci_width))
  return(Output)
}


# Summarize Multiple Fit Outputs
### fit_list = List of fits, indexed by description of fits
###            fit_list[["ex"]] should be a list with elements "fit_data" and
###            [["filter_statement"]] where
###             fit_list[["ex"]][["fit_data"]] = Single fit output from 
###               Conduct_Inference_for_Driver_Run() output
###             fit_list[["ex"]][["filter_statement"]] = Text statement filtering
###               rows of fit_list[["ex"]][["fit_data"]] to relevant rows
Summarize_Multiple_Fits <- function(fit_list) {
  temp_output <- data.frame()
  for (rowname in names(fit_list)) {
    row_output <- 
      Summarize_Single_Fit_Data(fit_data=fit_list[[rowname]][["fit_data"]],
                                filter_statement=fit_list[[rowname]][["filter_statement"]])
    temp_output <- rbind(temp_output,
                         data.frame(Description=rowname,
                                    Bias=row_output$bias,
                                    SD=row_output$sd,
                                    RMSE=row_output$RMSE,
                                    Coverage=row_output$coverage,
                                    Power=row_output$power,
                                    CI_Width=row_output$ci_width))
  }
  Output <- temp_output
  return(Output)
}



# Round and Format Numeric Columns of a Table
### df = Dataframe to round/format
### n_digits = Number of digits to round/format
Round_and_Format_DF <- function(df, n_digits=3) {
  temp_output <- df
  for (col in colnames(df)) {
    if (is.numeric(df[[col]])) {
      temp_output[[col]] <- format(round(temp_output[[col]], n_digits), 
                                   nsmall=n_digits)
    }
  }
  Output <- temp_output
  return(Output)
}