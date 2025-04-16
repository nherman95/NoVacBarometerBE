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
# - terminal: Rscript R/projections_vaccination.R belgium 0_0_0_0 
# - terminal: Rscript R/projections_vaccination.R belgium 0_0_0_0 &
# - within R: source('R/projections_vaccination.R')
# - within R: system('Rscript R/projections_vaccination.R belgium 0_0_0_0 &')
#
# Copyright 2022, SIMID, University of Antwerp & Hasselt University                                        
########################################################################### #

# Note: when "excecuting" an R-file, this is called line by line, which is error prone 
# when one edits the file while running a simulation (in the background). To enable 
# multiple runs of this script the same time, the code is captured in a function and 
# this function is called at the end of this script. 

# parse command line arguments
cl_args = commandArgs(trailingOnly=TRUE)

# clear workspace (except the command line arguments)
rm(list=ls()[ls()!='cl_args'])

## load required packages (quitely)
suppressPackageStartupMessages(library(simid.rtools))

# load functions and data 
source('R/main_vaccination.R')

# set output tag (= directory name)
output_tag <- 'hub_R4'

# select region?
sel_region <- 'belgium'
# sel_region <- 'brussels'
# sel_region <- 'flanders'
# sel_region <- 'wallonia'

# get parameters: default
chains_param_files = dir('data/config/',pattern='MCMCmulti_20230419_correctedfrom_hubR4final6for_opti_belgium',full.names = T,recursive = T)
chains_param_files = dir('data/config/',pattern='MCMCmulti_20230419_d1108_e2_i1200_n10_p10_crit7_hubR4final6for_pessi_belgium_foronly',full.names = T,recursive = T)
scen_tag <- 'PESSI' #'OPTI' #'PESSI'  #'OPTVAR' #'PESSVAR'


# select social contact scenario?
cnt_adjust_value <-  c(0,0,0,0,0) # c(0,0,1.1,1.2,1.3)  # c(0,0,0,0,0) c(1,0,0,0,1) # 1,1.3,1.5
cnt_adjust_date  <- c('2022-09-01',  # restart june 2022
                      '2022-09-01',  # hub seasonality
                      '2022-09-01',  # increase in contacts, cfr. last wave
                      '2021-09-01',  # increase in contacts, cfr. last wave
                      '2023-03-15')  # reuse all matrices last year  -> careful, file lib_social_contacts must be adapted

# number of model parameter sets (=chains), stochastic realisations and time horizon
num_chains          <- 50
num_stochastic_real <- 2
num_days_sim        <- 915
num_days_sim <- sim_date2day('2024-01-15')

#select vaccine uptake file(s) ----
vacc_files <- dir('data/uptake/uptake_hubR4',pattern = 'uptake_extrabooster.*csv',full.names = T)
#vacc_files <- vacc_files[3]

# option to adjust some specific parameters
adjusted_parameters <- data.frame()                            # default
# adjusted_parameters <- data.frame(log_VOC_omicron_init = 0)  # example

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
                        cl_args=cl_args)

