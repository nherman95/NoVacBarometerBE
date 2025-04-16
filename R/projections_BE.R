########################################################################### #
# This file is part of the Stochastic Compartmental Model for SARS-COV-2 
# transmission in Belgium, conceived by members of SIMID group during the 
# COVID19 pandemic.
#
# This file can be used to run multiple simulation runs with different 
# social contact and/or vaccination uptake assumptions for Belgium or 
# one of the three regions. 
#
# Current parallelism: for different vaccine uptake file
#
# This file can be executed with command line arguments. For example:
# - terminal: Rscript R/projections_BE.R belgium 0_0_0_0 
# - terminal: Rscript R/projections_BE.R belgium 0_0_0_0 &
# - within R: source('R/projections_BE.R')
# - within R: system('Rscript R/projections_BE.R belgium 0_0_0_0 &')
#
# Copyright 2023, SIMID, University of Antwerp & Hasselt University, UNamur                                       
########################################################################### #

# Note: when "executing" an R-file, this is called line by line, which is error prone 
# when one edits the file while running a simulation (in the background). To enable 
# multiple runs of this script the same time, the code is captured in a function and 
# this function is called at the end of this script. 

# parse command line arguments
cl_args = commandArgs(trailingOnly=TRUE)
# clear workspace (except the command line arguments)
rm(list=ls()[ls()!='cl_args'])
get_arg_with_default <- function(arg_index, default_value) {
  if (length(cl_args) >= arg_index) {
    arg_value <- cl_args[arg_index]
    default_type <- typeof(default_value)
    #print(default_type)
    if (default_type == "double") {
      return(as.numeric(arg_value))
    } else if (default_type == "integer") {
      return(as.integer(arg_value))
    } else if (default_type == "logical") {
      return(as.logical(arg_value))
    } else {
      return(arg_value)
    }
  } else {
    return(default_value)
  }
}
bool_vac <- get_arg_with_default(1,FALSE)
bool_bar <- get_arg_with_default(2,TRUE)
output_tag <- get_arg_with_default(3,"test_diff")
if(bool_bar == TRUE){
  # Function to get the value of an argument with a default value
  
  # Assign values with or without default values
  or_prop <- get_arg_with_default(4, 0.7)
  red_prop <- get_arg_with_default(5, 0.4)
  new_hosp_or <- get_arg_with_default(6, 65)
  new_hosp_red <- get_arg_with_default(7, 150)
  icu_or <- get_arg_with_default(8, 300)
  icu_red <- get_arg_with_default(9, 500)
#print(or_prop)
#xprint(red_prop)
# Assigner les valeurs aux variables
gran=or_prop-red_prop
ndayslong=3
step_prop=gran/ndayslong/4
}

## load required packages (quitely)
suppressPackageStartupMessages(library(simid.rtools))

if(bool_vac==FALSE){
if(bool_bar==TRUE){
  scen_tag <- paste("barometer_orange_", or_prop, "_", new_hosp_or, "_", icu_or, "_red_", red_prop, "_", new_hosp_red, "_", icu_red, sep = "")

} else{
  scen_tag <- 'nobarometer'
}
}

# load functions and data 
source('R/main_vaccination.R')

# set output tag (= directory name)
#output_tag <- 'test_bar'

# select region?
sel_region <- 'belgium'
# sel_region <- 'brussels'
# sel_region <- 'flanders'
# sel_region <- 'wallonia'

# get parameters: default
#chains_param_files = dir('data/config/',pattern='MCMCmulti_20230810_d1213_e2_i1000_n10_p10_crit7_hubR5foronly_opti_belgium_foronly',full.names = T,recursive = T)
chains_param_files = dir('data/config/',pattern='MCMCmulti_20230810_d1213_e2_i1000_n10_p10_crit7_hubR5foronly_pessi_belgium_foronly',full.names = T,recursive = T)



# select social contact scenario?
cnt_adjust_value <-  c(0,0,0,0,1) # c(0,0,1.1,1.2,1.3)  # c(0,0,0,0,0) c(1,0,0,0,1) # 1,1.3,1.5
cnt_adjust_date  <- c('2022-09-01',  # restart june 2022
                      '2022-09-01',  # hub seasonality
                      '2022-09-01',  # increase in contacts, cfr. last wave
                      '2021-09-01',  # increase in contacts, cfr. last wave
                      '2023-07-01')  # reuse all matrices last year  -> careful, file lib_social_contacts must be adapted


# number of model parameter sets (=chains), stochastic realisations and time horizon
num_chains          <- 2
num_stochastic_real <- 1
#num_days_sim        <- 519
num_days_sim <- sim_date2day('2023-07-01')

#select vaccine uptake file(s) ----
if(bool_vac==FALSE){
vacc_files <- dir('data/uptake/uptake_belgium_novac',pattern = 'uptake_extrabooster.*csv',full.names = T)
} else {
vacc_files <- dir('data/uptake/uptake_fall23campaign',pattern = 'uptake_extrabooster.*csv',full.names = T)
vacc_files <- vacc_files[c(1)]
scen_tag <- 'vac'
or_prop <-1
red_prop <- 1
new_hosp_or <- 0
new_hosp_red <- 0
icu_or <- 0
icu_red <- 0
}

# option to adjust some specific parameters
adjusted_parameters <- data.frame()                            # default
# adjusted_parameters <- data.frame(log_VOC_omicron_init = 0)  # example

#to use a parameter file
#ve_param <- read.table('./data/vaccine_parameters_20230820HubR5_pessi.csv',sep=',',header=T)
#adjusted_parameters <- ve_param
#print("Use VE parameter file",fill=T)


# call function
projections_vaccination(output_tag=output_tag,
                        sel_region=sel_region,
                        chains_param_files=chains_param_files,
                        cnt_adjust_value=cnt_adjust_value,
                        cnt_adjust_date=cnt_adjust_date,
                        scen_tag=scen_tag,
                        num_chains=num_chains,
                        num_stochastic_real=num_stochastic_real,
                        num_days_sim=num_days_sim,
                        vacc_files=vacc_files,
                        adjusted_parameters=adjusted_parameters,
                        cl_args=list())

