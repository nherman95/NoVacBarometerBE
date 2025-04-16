########################################################################### #
# This file is part of the Stochastic Compartmental Model for SARS-COV-2 
# transmission in Belgium, conceived by members of SIMID group during the 
# COVID19 pandemic.
#
# This file can be used to estimate model parameters using MCMC
#
# This file can be excecuted with command line arguments. For example:
# - terminal: Rscript R/calibration_mcmc.R 1 out1 1000 10 2020wave1 belgium
# - terminal: Rscript R/calibration_mcmc.R 1 out1 1000 10 2020wave1 belgium &
# - within R: source('R/calibration_mcmc.R')
# - within R: system('Rscript R/calibration_mcmc.R 1 out1 1000 10 2020wave1 belgium &')
#
# Copyright 2022, SIMID, University of Antwerp & Hasselt University                                        
########################################################################### #

# Note: when "excecuting" an R-file, this is called line by line, which is error prone 
# when one edits the file while running a simulation (in the background). To enable 
# multiple runs of this script the same time, the code is captured in a function and 
# this function is called at the end of this script with the command line arguments.

# parse command line arguments
args = commandArgs(trailingOnly=TRUE)

# clear workspace (except the command line arguments)
rm(list=ls()[ls()!='args'])

# load functions and data ----
source('R/main_vaccination.R')

################################################################ #
## SETTINGS AND PREAMBLE ----
################################################################ #

# main function parameters:
# - f_task_id           the id of the current task/run to conduct one MCMC chain (0:n)
# - f_subdir            the name of the subdirectory to store the model output
# - f_n_iter            the number of MCMC iterations 
# - f_n_period          the MCMC period
# - f_param             the parameter selection
# - f_region            the region to calibrate the model for (wrt previous config files and reference data)
# - f_filename_pattern  to grasp the configuration file with initial values

# uncomment for debugging:
# f_task_id = 0; f_subdir=NA; f_n_iter = NA; f_n_period = NA; f_param = 'debug'; f_region = NA;f_filename_pattern = 'MCMCmulti_20220403'
# f_task_id = 1; f_subdir=NA; f_n_iter = NA; f_n_period = NA; f_param = 'hubR4start2'; f_region = NA;f_filename_pattern = NA
run <- function(f_task_id = NA, f_subdir=NA, 
                f_n_iter = NA, f_n_period = NA, 
                f_param = 'foronly', f_region = NA,
                f_filename_pattern = NA){
  
  
# set default run tag and output dir
run_tag <- ''                           
output_subdir <- './output/MCMC_SCM'    

# set start parameter filename pattern
#filename_pattern <- 'MCMCmulti_20230807_d1213_e64_i200_n10_p10_crit4_hubR5mort_opti_belgium_mort' # belgium
filename_pattern <- 'MCMCmulti_20230807_d1213_e64_i200_n10_p10_crit4_hubR5mort_pessi_belgium_mort' # belgium


# get parameter file(s) to start from (selection is done later)
chains_param_files = dir('data/config/',full.names = T,recursive = T) 

# set directory with vaccine uptake files (region-specific selection is done later)
# June 2022
# vaccine_uptake_dir     <- 'data/uptake/uptake_manuscript'
# vaccine_uptake_pattern <- 'vSCENmay06_uadult80_uchild40_ubooster60_jun01_belgium'
# July 2022
vaccine_uptake_dir     <- 'data/uptake/uptake_fall23campaign'
vaccine_uptake_pattern <- 'vaccine_uptake_extrabooster_vSCENbelgium_fall23campaign_0flu_biv_belgium'


# set default MCMC parameters
n_iter   <- 11
n_period <- 10
n_repeat <- 10
n_status <- 10
sel_crit <- 'crit2'   # the log-likelihood criteria from the stochastic model
s_deviance <- 0.1

# set boolean to plot model output (using initial and final parameters)
bool_plot_model_output <- TRUE

################################################################ #
# LOAD AND PRE-PROCESS ----
################################################################ #

# optional: if chain_id missing, use "1"
if(!exists('f_task_id') || is.na(f_task_id)){ 
  f_chain_id <- 40
} else{ # else increment by 1, to cope with the SLURM id starting from 0
  f_chain_id <- as.numeric(f_task_id) #+1  removed for TORQUE 
}

# optional: add sub directory
if(exists('f_subdir') && !is.na(f_subdir)){                          
  run_tag       <- f_subdir 
  output_subdir <- paste0('./output/MCMC_',f_subdir)
} 

# back up run tag (manual setup)
if(nchar(run_tag)==0){ 
  run_tag <- 'f1a' 
}  

# optional: change region
if(!exists('f_region') || is.na(f_region)){ 
  f_region <- "belgium"
}  


# add region to run_tag
run_tag <- paste(run_tag,f_region,sep='_')

# add f_param tag
run_tag <- paste(run_tag,f_param,sep='_')

# add chain id to run_tag
run_tag <- paste(run_tag,f_chain_id,sep='_')

# create small delay when running different chains
if(as.numeric(f_chain_id)){
  Sys.sleep(f_chain_id)
}

# optional: adjust filname pattern for the parameter file
if(exists('f_filename_pattern') && !is.na(f_filename_pattern)){                          
  filename_pattern   <- f_filename_pattern
} 

# reference data
be_ref_data <- get_latest_incidence_data(sel_region = f_region)

# set default time horizon
ndays_sim <- sim_date2day(max(be_ref_data$date[!is.na(be_ref_data$hospital_admissions_other)]))
  
# option to use command line MCMC paramters
if(exists('f_n_iter') && !is.na(f_n_iter)){ n_iter <- as.numeric(f_n_iter) }  
if(exists('f_n_period') && !is.na(f_n_period)){ n_period <- as.numeric(f_n_period) }  

# option to select region-specific file
chains_param_file <- ifelse(any(grepl(f_region,chains_param_files)),
                            chains_param_files[grepl(f_region,chains_param_files)],
                            chains_param_files[1])

# other option to select a specific parameter file
chains_param_file <- ifelse(any(grepl(filename_pattern,chains_param_files)),
                            chains_param_files[grepl(filename_pattern,chains_param_files)],
                            chains_param_files[1])

# read parameter values
parms_chains      = read.table(chains_param_file, sep = ",", header = T)
print(chains_param_file)
dim(parms_chains)

# Use vaccine parameters file
#ve_param <- read.table('./data/vaccine_parameters_20230820HubR5_opti.csv',sep=',',header=T)
#ve_param <- read.table('./data/vaccine_parameters_20230820HubR5_pessi.csv',sep=',',header=T)
#parms_chains[names(ve_param)] <- ve_param
#print("Warning: Use VE specific parameter file",fill=T)

# to include new CoMix waves (if new data is available)
while(identify_nb_waves(names(parms_chains),bool_comix = TRUE) < length(get_CoMix_change_date())){
    parms_chains <- include_new_wave(parms_chains,bool_comix = TRUE)
}

# select parameters set (option to vary the parameter set using the chain_id)
parms      <- unlist(parms_chains[min(as.numeric(f_chain_id),nrow(parms_chains)),])

# update parameters if needed: new parameter (names)
parms <- update2latest_model_parameter_config(parms)

# set parameter names
parms_names = names(parms)

# option to aggregate waves
opt_aggr_waves <- NA

# update region_id in param
parms['region_id'] <- get_region_id(f_region)

# # update hospital hazard ratio
# parms <- parms[!grepl('log_VOC_.*hosp',names(parms))]
# parms['log_VOC_hr_hosp'] <- 1e-15
# parms['log_VOC_delta_hr_hosp'] <- log(2.26)
# parms[paste0('log_VOC_omicron_hr_hosp_age',1:10)] <- log(c(1.0,0.89,0.67,0.57,0.54,0.42,0.32,0.42,0.49,0.49))
# parms_names <- names(parms)

# check names
if(any(is.na(names(parms)))){
  stop('ISSUE WITH PARAMETER NAMES... PLEASE CHECK FOR "NA" ')
}

################################################################ #
# SELECT MODEL PARAMETERS TO ESTIMATE ----
################################################################ #

# DEFAULT PARAMETER SET: all (except "ndays_calibration" and "h, "region_id" and mcmc meta info)
parms_names_estim <- parms_names[-which(parms_names %in% c("ndays_calibration","h","region_id") | grepl('mcmc_',parms_names))]
nwave_Comix <- identify_nb_waves(parms_names,bool_comix = TRUE)
nwave_other <- identify_nb_waves(parms_names,bool_comix = FALSE)

# option: 2020 wave1  ----
if(exists('f_param') && !is.na(f_param) && f_param == '2020wave1'){
  ndays_sim <- 120 
  
  
  # check and limit adjusted omega parameters
  flag_omega <- grepl('log_omega',parms_names)
  parms[flag_omega & parms < log(1/7)] <- log(1/2)
  parms[flag_omega & parms > log(2)]   <- log(1/2)
  parms_chains[,flag_omega]            <- log(1/2)
  
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    sel_id = 1:min(which(grepl('log_comix',parms_names))-1),
                                    tag_list = c(paste0('coef_w',1:3,'_'),
                                                 paste0('mat_w',1,'_')),
                                    ignore_list = c('log_phi1',
                                                    'log_beta1',
                                                    'log_delta3',
                                                    'log_delta4',
                                                    'log_mu')
  )
  sel_crit <- 'crit1'
}

# option: 2020 wave 2----
if(exists('f_param') && !is.na(f_param) && f_param == '2020wave2'){
  ndays_sim <- 305 
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    tag_list = c(paste0('coef_w',3:11,'_'),
                                                 paste0('mat_w',3:8,'_'))
  )
  sel_crit <- 'crit1'
}

# option: 2021 alpha  ----
if(exists('f_param') && !is.na(f_param) && f_param == '2021alpha'){
  ndays_sim <- 457
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    tag_list = c(paste0('coef_w',11:20,'_'),
                                                 'log_VOC_alpha'),
                                    ignore_list = c('hosp',
                                                    'phi1_add')
  )
  sel_crit <- 'crit2'
}

# option: 2021 delta  ----
if(exists('f_param') && !is.na(f_param) && f_param == '2021delta'){
  ndays_sim <- 563 
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    tag_list = c(paste0('coef_w',20:30,'_'),
                                                 'log_VOC_delta'),
                                    ignore_list = c('hosp',
                                                    'phi1_add')
  )
  sel_crit <- 'crit2'
}

# option: 2021 september-november  ----
if(exists('f_param') && !is.na(f_param) && f_param == '2021sept'){
  ndays_sim <- 650 
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    tag_list = c(paste0('coef_w',30:36,'_')),
                                    ignore_list = c('hosp',
                                                    'phi1_add')
  )
  sel_crit <- 'crit2'
}

# option: omicron (ba1ba2)----
if(exists('f_param') && !is.na(f_param) && grepl('omicron',f_param)){
  ndays_sim <- 730 
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    tag_list = c('log_VOC_omicron_init',
                                                 'log_VOC_omicron_transm',
                                                 'log_VOC_omicron_gamma_factor',
                                                 paste0('coef_w',36:42,'_')))
  sel_crit <- 'crit2'
}

# option: ba4ba5 ----
if(exists('f_param') && !is.na(f_param) && grepl('ba4ba5',f_param)){
  
  # select parameters
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    tag_list = c('log_VOC_ba4ba5_init',
                                                 'log_VOC_ba4ba5_transm',
                                                 paste0('coef_w',42:nwave_Comix,'_')))
  sel_crit <- 'crit2'
}

# option: omicron excl ----
if(exists('f_param') && !is.na(f_param) && grepl('omicron_excl',f_param)){
  
  # select parameters
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    tag_list = c('log_VOC_omicron_init',
                                                 'log_VOC_omicron_transm',
                                                 'log_VOC_omicron_gamma_factor'))
  sel_crit <- 'crit5'
}

# option: ba4ba5 excl ----
if(exists('f_param') && !is.na(f_param) && grepl('ba4ba5_excl',f_param)){
  
  sel_param_omicron <- parms[grepl('VOC_omicron',parms_names)]
  parms[gsub('omicron','ba4ba5',names(sel_param_omicron))] <- sel_param_omicron
  parms['VOC_ba4ba5_start'] <- 755
  parms['log_VOC_ba4ba5_transm'] <- parms['log_VOC_omicron_transm'] * 1.2
  parms_names <- names(parms)
  
  # select parameters
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    tag_list = c('log_VOC_ba4ba5_init',
                                                 'log_VOC_ba4ba5_transm'))
  sel_crit <- 'crit5'
}

# option: ba4ba5 + hospi ----
if(exists('f_param') && !is.na(f_param) && grepl('ba4ba5_hosp',f_param)){
  
  sel_param_omicron <- parms[grepl('VOC_omicron',parms_names)]
  parms[gsub('omicron','ba4ba5',names(sel_param_omicron))] <- sel_param_omicron
  parms['VOC_ba4ba5_start'] <- 755
  parms['log_VOC_ba4ba5_transm'] <- parms['log_VOC_omicron_transm'] * 1.2
  parms_names <- names(parms)
  
  # select parameters
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    tag_list = c('log_VOC_ba4ba5_init',
                                                 'log_VOC_ba4ba5_transm',
                                                 paste0('coef_w',(nwave_Comix-6):nwave_Comix,'_')))
  sel_crit <- 'crit2'
}

# option: bq1 ----
if(exists('f_param') && !is.na(f_param) && grepl('bq1',f_param)){
  
  # select parameters
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    tag_list = c('VOC_bq1',
                                                 paste0('coef_w',47:nwave_Comix,'_')))
  sel_crit <- 'crit2'
  print(c("data up to ",ndays_sim))
}


# option: all voc + for covid ----
if(exists('f_param') && !is.na(f_param) && grepl('vocandfor',f_param)){
  
  # select parameters
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    tag_list = c('init',
                                                 'transm',
                                                 'for_covid_coef'),
                                    ignore_list = c("ve_transmission"))
  sel_crit <- 'crit7'
}

# option: for covid ----
if(exists('f_param') && !is.na(f_param) && grepl('foronly',f_param)){
  
  # parms['for_covid_coef'] <- 10000
   #parms_names <- names(parms)
  # select parameters
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    tag_list = c('for_covid_coef',"ve_transmission"))
  #s_deviance <- 0.0001
  sel_crit <- 'crit7'
  
  #date specified
 # print(c("data up to ",ndays_sim))
#  ndays_sim <- sim_date2day('2023-03-14')
 # print(c("calibration up to ",ndays_sim))
}


# option: recap13 (last 13 waves)  ----
if(exists('f_param') && !is.na(f_param) && f_param == 'recap13'){
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    tag_list = paste0('coef_w',(nwave_Comix-12):nwave_Comix,'_'),
  )
  sel_crit <- 'crit1'
  #ndays_sim <- 760 
}

# option: recap7 (last 7 waves)  ----
if(exists('f_param') && !is.na(f_param) && f_param == 'recap7'){
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    tag_list = paste0('coef_w',(nwave_Comix-6):nwave_Comix,'_'),
  )
  sel_crit <- 'crit1'
  #date specified
  print(c("data up to ",ndays_sim))
  ndays_sim <- sim_date2day('2022-09-10')
  print(c("calibration up to ",ndays_sim))
}

# option: recap3 (last 3 waves)  ----
if(exists('f_param') && !is.na(f_param) && f_param == 'recap3'){
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    tag_list = paste0('coef_w',(nwave_Comix-2):nwave_Comix,'_'),
  )
  sel_crit <- 'crit1'
}

# option: initial lockdown  ----
if(exists('f_param') && !is.na(f_param) && f_param == 'lockdown'){
  ndays_sim <- 120 
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    tag_list = c(paste0('coef_w',1:3,'_'),
                                                 paste0('mat_w',1:2,'_'),
                                                 'n0'    # initial number of cases
                                    )
  )
  sel_crit <- 'crit1'
}

# option: hosp and ICU load ----
if(exists('f_param') && !is.na(f_param) && grepl('hload',f_param)){
  
  # parms['icu_fall2021_start'] <- sim_date2day('2021-10-01')
  # parms['log_fall2021_phi1_add'] <- parms['log_VOC_delta_phi1_add']
  # parms_names <- names(parms)
  # 
  # parms['icu_spring2022_start'] <- sim_date2day('2022-03-01')
  # parms['log_spring2022_phi1_add'] <- parms['log_VOC_omicron_phi1_add']
  # parms_names <- names(parms)
  
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    tag_list = c('log_delta3',
                                                 'phi1_add')

  )
  #date specified
  print(c("data up to ",sim_day2date(ndays_sim)))
  print(sim_day2date(ndays_sim))
 # ndays_sim <- sim_date2day('2022-09-10')
  print(c("calibration up to ",sim_day2date(ndays_sim)))
  print(sim_day2date(ndays_sim))
  sel_crit <- 'crit3'
}

# option: mortality ----
if(exists('f_param') && !is.na(f_param) && grepl('mort',f_param)){
  #ndays_sim <- 305 
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    tag_list = c('log_mu')
                                    
  )
 # print(c("data up to ",ndays_sim))
 # ndays_sim <- sim_date2day('2023-03-14')
 # print(c("calibration up to ",ndays_sim))
  sel_crit <- 'crit4'
}

# waning immunity ----
if(exists('f_param') && !is.na(f_param) && f_param == 'waning'){
  
  # add paramters for waning immunity
  parms['ve_waning_immunity_rate']     <- 1/(6*30)
  parms['ve_alpha_infection_waning']   <- parms['ve_alpha_infection_rna2']/2         # reduction with respect to infection with waning immunity
  parms['ve_alpha_incr_severe_waning'] <- parms['ve_alpha_incr_severe_rna2']/2       # incremental reduction within dynamic model
  parms['ve_delta_infection_waning']   <- parms['ve_delta_infection_rna2']/2   # reduction with respect to infection with waning immunity
  parms['ve_delta_incr_severe_waning'] <- parms['ve_delta_incr_severe_rna2']/2
  parms_names <- names(parms)
  
  # select CoMix wave from summer 2021
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    tag_list = c(paste0('coef_w',25:nwave_Comix,'_')))
  sel_crit <- 'crit1'
}


# hub ----
if(exists('f_param') && !is.na(f_param) && f_param == 'hubR1'){
  
  # add parameters for waning immunity after 2nd dose
  parms['ve_waning_immunity_rate']     <- 1/(3*30)
  parms['ve_alpha_infection_waning']   <- parms['ve_alpha_infection_rna2']*0.4         # reduction with respect to infection with waning immunity
  parms['ve_delta_infection_waning']   <- parms['ve_delta_infection_rna2']*0.4     # reduction with respect to infection with waning immunity
  parms['ve_omicron_infection_waning'] <- parms['ve_omicron_infection_rna2']*0.4     # reduction with respect to infection with waning immunity
  parms['ve_alpha_incr_severe_waning'] <- parms['ve_alpha_incr_severe_rna2']*0.8^(3/3)       # incremental reduction within dynamic model
  parms['ve_delta_incr_severe_waning'] <- parms['ve_delta_incr_severe_rna2']*0.8^(3/3)
  parms['ve_omicron_incr_severe_waning'] <- parms['ve_omicron_incr_severe_rna2']*0.8^(3/3)
  
  #add parameters for waning immunity after booster
  parms['ve_waning_booster_rate'] <- parms['ve_waning_immunity_rate']
  parms['ve_alpha_infection_booster_waning']   <- parms['ve_infection_booster']*0.4         # reduction with respect to infection with waning immunity
  parms['ve_delta_infection_booster_waning']   <- parms['ve_delta_infection_booster']*0.4   # reduction with respect to infection with waning immunity
  parms['ve_omicron_infection_booster_waning'] <- parms['ve_omicron_infection_booster']*0.4   # reduction with respect to infection with waning immunity
  parms['ve_alpha_incr_severe_booster_waning'] <- parms['ve_incr_severe_booster']*0.8^(3/3)       # incremental reduction within dynamic model
  parms['ve_delta_incr_severe_booster_waning'] <- parms['ve_delta_incr_severe_booster']*0.8^(3/3)
  parms['ve_omicron_incr_severe_booster_waning'] <- parms['ve_omicron_incr_severe_booster']*0.8^(3/3)


  #add parameters for waning immunity after infection
  parms['ve_waning_infection_rate'] <- parms['ve_waning_booster_rate']
  parms['ve_waning_infection_booster_rate'] <- parms['ve_waning_booster_rate']
  parms['ve_alpha_infection_reinf']     <- 0.4#parms['ve_infection_booster']*0.4
  parms['ve_delta_infection_reinf']     <- 0.4#parms['ve_delta_infection_booster']*0.4
  parms['ve_omicron_infection_reinf']   <- 0.4#parms['ve_omicron_infection_booster']*0.4
  parms['ve_alpha_infection_reinfvac']     <- 0.4#parms['ve_infection_booster']*0.4
  parms['ve_delta_infection_reinfvac']     <- 0.4#parms['ve_delta_infection_booster']*0.4
  parms['ve_omicron_infection_reinfvac']   <- 0.4#parms['ve_omicron_infection_booster']*0.4
  parms['ve_alpha_incr_severe_reinf']   <- parms['ve_incr_severe_booster']*0.8^(8/3)
  parms['ve_delta_incr_severe_reinf']   <- parms['ve_delta_incr_severe_booster']*0.8^(8/3)
  parms['ve_omicron_incr_severe_reinf'] <- parms['ve_omicron_incr_severe_booster']*0.8^(8/3)
  parms['ve_alpha_incr_severe_reinfvac']   <- parms['ve_incr_severe_booster']*0.8^(8/3)
  parms['ve_delta_incr_severe_reinfvac']   <- parms['ve_delta_incr_severe_booster']*0.8^(8/3)
  parms['ve_omicron_incr_severe_reinfvac'] <- parms['ve_omicron_incr_severe_booster']*0.8^(8/3)
  
   parms_names <- names(parms)
  
  
  # select CoMix wave from sept20
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    tag_list = c(paste0('coef_w',9:nwave_Comix,'_'),
                                                 paste0('mat_w',3:nwave_other,'_'))
  )
  sel_crit <- 'crit1'
  
  # Hub change
  print(c("data up to ",ndays_sim))
  ndays_sim <- sim_date2day('2022-07-20')
  print(c("calibration up to ",ndays_sim))
  
}

# hub R5  calibration (all restart from MCMCmulti_20230629_d1213_e40_i150_n10_p10_crit1_workingwave63_belgium_workingend) ----
# Hub R5 start (-> jan22) + set new parameters VE (now from file)
if(exists('f_param') && !is.na(f_param) && f_param == 'hubR5startnewVE'){
  
  
   # select CoMix wave from sept20 up to omicron - 
   # select CoMix wave up to jan22 (from start)
   parms_names_estim <- get_colnames(parms_names = parms_names,
                                     tag_list = c(paste0('coef_w',1:37,'_'),
                                                  paste0('mat_w',1:nwave_other,'_'),
                                                  'log_VOC_alpha_init',
                                                  'log_VOC_delta_init',
                                                  'log_VOC_alpha_transm',
                                                  'log_VOC_delta_transm')
   )
   sel_crit <- 'crit2'

   print(c("data up to ",sim_day2date(ndays_sim)))
   print(sim_day2date(ndays_sim))
   ndays_sim <- sim_date2day('2022-01-01')
   print(c("calibration up to ",sim_day2date(ndays_sim)))
   print(sim_day2date(ndays_sim))
}

# Hub R5 mid (nov21-> sept22) + differenciate new parameters VE (from file) opti and pessi
if(exists('f_param') && !is.na(f_param) && f_param == 'hubR5mid'){
  
  
  # select CoMix wave from nov21 (omicron) up to endsep22 (ba4ba5)
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    tag_list = c(paste0('coef_w',35:48,'_'),
                                                 'log_VOC_omicron_init',
                                                 'log_VOC_ba4ba5_init',
                                                 'log_VOC_omicron_transm',
                                                 'log_VOC_ba4ba5_transm',
                                                 'log_VOC_omicron_transm',
                                                 'log_VOC_ba4ba5_transm',
                                                 'log_VOC_omicron_gamma_factor',
                                                 'log_VOC_ba4ba5_gamma_factor',
                                                 've_VOC_omicron_vereduction',
                                                 've_VOC_ba4ba5_vereduction')
  )
  sel_crit <- 'crit2'
  
  
  print(c("data up to ",sim_day2date(ndays_sim)))
  print(sim_day2date(ndays_sim))
  ndays_sim <- sim_date2day('2022-10-01')
  print(c("calibration up to ",sim_day2date(ndays_sim)))
  print(sim_day2date(ndays_sim))
}

# Hub R5 end (sept22-> end) 
if(exists('f_param') && !is.na(f_param) && f_param == 'hubR5end'){
  
  
  # select CoMix wave from sept22 
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    tag_list = c(paste0('coef_w',48:nwave_Comix,'_'),
                                                 'log_VOC_bq1ba275xbb_init',
                                                 'log_VOC_bq1ba275xbb_transm',
                                                 'log_VOC_bq1ba275xbb_gamma_factor',
                                                 've_VOC_bq1ba275xbb_vereduction',
                                                 'log_VOC_xbb15_init',
                                                 'log_VOC_xbb15_transm',
                                                 'log_VOC_xbb15_gamma_factor',
                                                 've_VOC_xbb15_vereduction'))
  sel_crit <- 'crit2'
  
  
  print(c("data up to ",sim_day2date(ndays_sim)))
  print(sim_day2date(ndays_sim))
 # ndays_sim <- sim_date2day('2022-10-01')
  print(c("calibration up to ",sim_day2date(ndays_sim)))
  print(sim_day2date(ndays_sim))
}



# hub R4  calibration ----
if(exists('f_param') && !is.na(f_param) && f_param == 'hubR4final3'){
  
  # waning_inf <- 0.4
  # waning_sev <- 0.8
  # # add parameters for waning immunity after 2nd dose
  # parms['ve_waning_immunity_rate']     <- log(2)/(6*30) #median time = 6 months
  # parms['ve_waning_booster_rate'] <- parms['ve_waning_immunity_rate']
  # parms['ve_waning_infection_rate'] <- parms['ve_waning_booster_rate']
  # parms['ve_waning_infection_booster_rate'] <- parms['ve_waning_booster_rate']
  # 
  # parms[grepl('infection_waning',parms_names)] <- parms[grepl('infection_rna2',parms_names)]*waning_inf 
  # parms[grepl('incr_severe_waning',parms_names)] <- parms[grepl('incr_severe_rna2',parms_names)]*waning_sev 
  # parms[grepl('infection_booster_waning',parms_names)] <- parms[grepl('infection_booster',parms_names) & !grepl('infection_booster_waning',parms_names) & !grepl('infection_booster_rate',parms_names)]*waning_inf 
  # parms[grepl('incr_severe_booster_waning',parms_names)] <- parms[grepl('incr_severe_booster',parms_names) & !grepl('incr_severe_booster_waning',parms_names)]*waning_sev 
  # 
  # parms[grepl('infection_reinf',parms_names)] <- waning_inf
  # parms[grepl('incr_severe_reinf',parms_names)] <- waning_sev
  # 
  # 
  # parms_names <- names(parms)
  
  # # select CoMix wave from sept20
  # parms_names_estim <- get_colnames(parms_names = parms_names,
  #                                   tag_list = c(paste0('coef_w',9:47,'_'),
  #                                                paste0('mat_w',3:nwave_other,'_'),
  #                                                'init',
  #                                                'transm')
  # )
  # sel_crit <- 'crit2'
  # print(c("data up to ",sim_day2date(ndays_sim)))
  # ndays_sim <- sim_date2day('2022-09-01')
  # print(c("calibration up to ",sim_day2date(ndays_sim)))
  
  # parms_names_estim <- get_colnames(parms_names = parms_names,
  #                                   tag_list = c(paste0('coef_w',35:nwave_Comix,'_'),
  #                                                'init',
  #                                                'transm',
  #                                                'vereduction')
  # )
  # sel_crit <- 'crit2'
  # # Hub change
  # print(c("data up to ",ndays_sim))
  # ndays_sim <- sim_date2day('2023-03-14')
  # print(c("calibration up to ",ndays_sim))
  
  # parms_names_estim <- get_colnames(parms_names = parms_names,
  #                                   tag_list = c(paste0('coef_w',25:46,'_'),
  #                                                'init',
  #                                                'transm',
  #                                                'vereduction')
  # )
  # sel_crit <- 'crit2'
  # # Hub change
  # print(c("data up to ",ndays_sim))
  # ndays_sim <- sim_date2day('2022-07-01')
  # print(c("calibration up to ",ndays_sim))
  
  # parms_names_estim <- get_colnames(parms_names = parms_names,
  #                                   tag_list = c(paste0('coef_w',46:nwave_Comix,'_'),
  #                                                'init',
  #                                                'transm',
  #                                                'vereduction')
  # )
  # sel_crit <- 'crit2'
  
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    tag_list = c(paste0('coef_w',(nwave_Comix-7):nwave_Comix,'_'),
                                                                                    'xbb15_init',
                                                                                    'xbb15_transm',
                                                                                    'xbb15_vereduction')
  )
  
  sel_crit <- 'crit2'
  # Hub change
  print(c("data up to ",ndays_sim))
  ndays_sim <- sim_date2day('2023-03-14')
  print(c("calibration up to ",ndays_sim))
  
  # parms_names_estim <- get_colnames(parms_names = parms_names,
  #                                   tag_list = c(paste0('coef_w',1:nwave_Comix,'_'),
  #                                                paste0('mat_w',3:nwave_other,'_'),
  #                                                'init',
  #                                                'transm',
  #                                                'vereduction')
  # )
  # sel_crit <- 'crit2'
  # 
  # # Hub change
  # print(c("data up to ",ndays_sim))
  # ndays_sim <- sim_date2day('2023-03-14')
  # print(c("calibration up to ",ndays_sim))
  
  # select CoMix wave from sept20
  # parms_names_estim <- get_colnames(parms_names = parms_names,
  #                                   tag_list = c(paste0('coef_w',9:38,'_'),
  #                                                paste0('mat_w',4:nwave_other,'_'))
  # )
  # sel_crit <- 'crit1'
  # print(c("data up to ",sim_day2date(ndays_sim)))
  # ndays_sim <- sim_date2day('2022-01-01')
  # print(c("calibration up to ",sim_day2date(ndays_sim)))
  
}


# hub ----
# if(exists('f_param') && !is.na(f_param) && f_param == 'hubR2'){
#   
#   parms['ve_transmission'] <- 0
#   parms_names <- names(parms)
#   
#   
#   # select CoMix wave from Jan 22
#  # parms_names_estim <- get_colnames(parms_names = parms_names,
#   #                                  tag_list = c(paste0('coef_w',39:nwave_Comix,'_'),
#    #                                              'init',
#     #                                             'transm'),
#      #                               ignore_list = c("ve_transmission"))
#   
#   # select CoMix wave from sept20
#   parms_names_estim <- get_colnames(parms_names = parms_names,
#                                     tag_list = c(paste0('coef_w',1:nwave_Comix,'_'),
#                                                  paste0('mat_w',1:nwave_other,'_'),
#                                                  'init',
#                                                  'transm',
#                                                  'log_VOC_omicron_gamma_factor',
#                                                  'log_VOC_ba4ba5_gamma_factor'),
#                                     ignore_list = c("ve_transmission"))
# 
#   sel_crit <- 'crit2'
#   
#   # Hub change
#   print(c("data up to ",ndays_sim))
#   ndays_sim <- sim_date2day('2022-07-20')
#   print(c("calibration up to ",ndays_sim))
#   
# }

if(exists('f_param') && !is.na(f_param) && f_param == 'hubR2'){
  
  # #restart from R1 final
  # sel_param_omicron <- parms[grepl('VOC_omicron',parms_names)]
  # parms[gsub('omicron','ba4ba5',names(sel_param_omicron))] <- sel_param_omicron
  # parms['VOC_ba4ba5_start'] <- 755
  # parms['log_VOC_ba4ba5_transm'] <- parms['log_VOC_omicron_transm'] * 1.2
  # parms_names <- names(parms)
  # 
  # # select CoMix wave from sept 21
  # parms_names_estim <- get_colnames(parms_names = parms_names,
  #                                   tag_list = c(paste0('coef_w',30:nwave_Comix,'_'),
  #                                                'log_VOC_omicron_init',
  #                                                'log_VOC_ba4ba5_init',
  #                                                'log_VOC_omicron_transm',
  #                                                'log_VOC_ba4ba5_transm',
  #                                                'log_VOC_omicron_gamma_factor',
  #                                                'log_VOC_ba4ba5_gamma_factor'),
  #                                   ignore_list = c("ve_transmission"))
  # 
  #sel_crit <- 'crit2'
  
  parms[grepl('log_VOC_ba4ba5_hr_hosp',names(parms))] <- 0
  parms_names <- names(parms)
  # select CoMix wave from march21
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    tag_list = c(paste0('coef_w',34:nwave_Comix,'_'),
                                                 'log_VOC_ba4ba5_init',
                                                 'log_VOC_ba4ba5_transm',
                                                  'log_VOC_ba4ba5_gamma_factor'),
                ignore_list = c("ve_transmission"))

  sel_crit <- 'crit2'
  
  # Hub change
  print(c("data up to ",ndays_sim))
  ndays_sim <- sim_date2day('2022-07-23')
  print(c("calibration up to ",ndays_sim))
  
}


if(exists('f_param') && !is.na(f_param) && f_param == 'hubextra'){
  
  # select CoMix wave from fev21
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    tag_list = c(paste0('coef_w',32:nwave_Comix,'_'))
  )
  sel_crit <- 'crit1'
}

if(exists('f_param') && !is.na(f_param) && f_param == 'loadford'){
  
 # parms['for_covid_coef'] <- 20000
#  parms_names <- names(parms)
  
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    tag_list = c('for_covid_coef',
                                                'log_delta3',
                                                 'phi1_add',
                                                 'log_mu')
                                    
  )
  sel_crit <- 'crit6'
  # Hub change
  print(c("data up to ",ndays_sim))
  ndays_sim <- sim_date2day('2023-03-14')
  print(c("calibration up to ",ndays_sim))
}

# recalibration for working version
if(exists('f_param') && !is.na(f_param) && f_param == 'workingstart'){
  
  #beginning
   # select CoMix wave up to jan22 (from start)
   parms_names_estim <- get_colnames(parms_names = parms_names,
                                     tag_list = c(paste0('coef_w',1:37,'_'),
                                                  paste0('mat_w',1:nwave_other,'_'),
                                                  'log_VOC_alpha_init',
                                                  'log_VOC_delta_init',
                                                  'log_VOC_alpha_transm',
                                                  'log_VOC_delta_transm')
   )
   sel_crit <- 'crit2'
   print(c("data up to ",sim_day2date(ndays_sim)))
   ndays_sim <- sim_date2day('2022-01-01')
   print(c("calibration up to ",sim_day2date(ndays_sim)))
  
}

#try to adapt omicron parameters to immune evasion
if(exists('f_param') && !is.na(f_param) && f_param == 'workingvoc'){
  
  parms['ve_VOC_omicron_vereduction'] <- 0.5
  parms['ve_VOC_ba4ba5_vereduction'] <- 0.5
  parms['ve_VOC_bq1ba275xbb_vereduction'] <- 0.5
  parms['ve_VOC_xbb15_vereduction'] <- 0.5
  parms['log_VOC_omicron_transm'] <- 1
  parms['log_VOC_ba4ba5_transm'] <- 1
  parms['log_VOC_bq1ba275xbb_transm'] <- 1
  parms['log_VOC_xbb15_transm'] <- 1
  parms_names <- names(parms)
  
  # select parameters
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    tag_list = c('log_VOC_omicron_init',
                                                 'log_VOC_omicron_transm',
                                                 'log_VOC_omicron_gamma_factor',
                                                 've_VOC_omicron_vereduction',
                                                 'log_VOC_ba4ba5_init',
                                                 'log_VOC_ba4ba5_transm',
                                                 'log_VOC_ba4ba5_gamma_factor',
                                                 've_VOC_ba4ba5_vereduction',
                                                 'log_VOC_bq1ba275xbb_init',
                                                 'log_VOC_bq1ba275xbb_transm',
                                                 'log_VOC_bq1ba275xbb_gamma_factor',
                                                 've_VOC_bq1ba275xbb_vereduction',
                                                 'log_VOC_xbb15_init',
                                                 'log_VOC_xbb15_transm',
                                                 'log_VOC_xbb15_gamma_factor',
                                                 've_VOC_xbb15_vereduction'))
  sel_crit <- 'crit5'
  
}

#mid
# select CoMix between waves 32 and 54
if(exists('f_param') && !is.na(f_param) && f_param == 'workingmid'){
  

  
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    tag_list = c(paste0('coef_w',32:54,'_')))
  sel_crit <- 'crit1'
   print(c("data up to ",sim_day2date(ndays_sim)))
   print(sim_day2date(ndays_sim))
    ndays_sim <- sim_date2day('2023-01-01')
   print(c("calibration up to ",sim_day2date(ndays_sim)))
   print(sim_day2date(ndays_sim))
   
}

#end
# select CoMix wave up end (from wave ***) + omicron
if(exists('f_param') && !is.na(f_param) && f_param == 'workingend'){

#!!!!!!!!!!!  # note : vereduction changed to log
 # vereduction <- parms[grepl('vereduction',parms_names)]
#  vereduction[vereduction<=0] <- 0.1
  #vereduction[vereduction>1] <- 1
#  parms[grepl('vereduction',parms_names)] <- log(vereduction)
   
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    tag_list = c(paste0('coef_w',60:nwave_Comix,'_')))#,
                                                # 'log_VOC_omicron_init',
                                                # 'log_VOC_omicron_transm',
                                                # 'log_VOC_omicron_gamma_factor',
                                                # 've_VOC_omicron_vereduction',
                                                 #'log_VOC_ba4ba5_init',
                                                 #'log_VOC_ba4ba5_transm',
                                                 #'log_VOC_ba4ba5_gamma_factor',
                                                 #'ve_VOC_ba4ba5_vereduction',
                                               #  'log_VOC_bq1ba275xbb_init',
                                                # 'log_VOC_bq1ba275xbb_transm',
                                                # 'log_VOC_bq1ba275xbb_gamma_factor',
                                                # 've_VOC_bq1ba275xbb_vereduction',
                                                # 'log_VOC_xbb15_init',
                                                # 'log_VOC_xbb15_transm',
                                                # 'log_VOC_xbb15_gamma_factor',
                                                # 've_VOC_xbb15_vereduction'))
  sel_crit <- 'crit1'
 # print(c("data up to ",sim_day2date(ndays_sim)))
#  ndays_sim <- sim_date2day('2022-01-01')
 # print(c("calibration up to ",sim_day2date(ndays_sim)))
  print(c("data up to ",sim_day2date(ndays_sim)))
  print(sim_day2date(ndays_sim))
  
}

# debug using dummy set (default) ----
if(exists('f_param') && !is.na(f_param) && f_param == 'debug'){
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    sel_id = c(5,9,54,125,200),
                                    tag_list = c('log_comix_coef_w23_age')
  )
  sel_crit <- 'crit2'
  ndays_sim <- 100
}

# check parms_names_estim ----
# if f_param is not recognized and parms_names_estim did not change, end MCMC calibration
if(length(parms_names_estim) == length(parms_names)-3){
  warning("PARAMETER SELECTION TAG UNKNOWN")
  return(-1)
}

# set rng seed
tag_numeric  <- as.numeric(charToRaw(run_tag));
seed_numeric <- sum(tag_numeric*10*seq(length(tag_numeric)))
set.seed(sum(tag_numeric*1:length(tag_numeric)))

# update calibration model horizon (could be changed during parameter selection)
parms['ndays_calibration'] <- ndays_sim

# set output file path (using current date and time)
output_tag <- file.path(output_subdir,paste0(format(Sys.time(),'%Y%m%d_%H%M%S'),'_d',ndays_sim,'_e',length(parms_names_estim),'_i',n_iter,'_n',n_repeat,'_p',n_period,'_',sel_crit,'_',run_tag))
if(is_vsc()) output_tag <- gsub('./output',file.path(get_vsc_scratch_folder(),'stochastic_model/output'),output_tag)
# temp for Ceci use instead of VSC
#print("warning: output folder for CECI and not VSC")
#output_tag <- gsub('./output',file.path('/scratch/unamur/URMaC/nfranco/','stochastic_model/output'),output_tag)
if(!dir.exists(output_tag)){dir.create(output_tag,recursive = T)}
print(output_tag)

# prepare vaccine uptake
# this is done once, to save time each iteration
vacc_files_all <- dir(vaccine_uptake_dir,pattern = paste0(f_region,'.*.csv'),full.names = T)
i_vac_file     <- vacc_files_all[grepl(vaccine_uptake_pattern,vacc_files_all)]
vacc_schedule <- read.csv(i_vac_file, header = T)

print(i_vac_file)

# vaccine protection: 1st and 2nd dose
V_mat         <- get_uptake_matrix(vaccine_uptake = vacc_schedule, dose_tag = '_A_', parms = parms )
V_mat_dose2   <- get_uptake_matrix(vaccine_uptake = vacc_schedule, dose_tag = '_B_', parms = parms )
V_mat_booster <- get_uptake_matrix(vaccine_uptake = vacc_schedule, dose_tag = '_E_', parms = parms )
V_mat_2ndbooster <- get_uptake_matrix(vaccine_uptake = vacc_schedule, dose_tag = '_F_', parms = parms )
V_mat_extrabooster <- get_uptake_matrix(vaccine_uptake = vacc_schedule, dose_tag = '_X_', parms = parms )


# save parameters names
write.table(parms_names_estim,file= file.path(output_tag,'MCMC_parameter_estim.csv'),sep=',',row.names = F, col.names = F)
write.table(parms_names,file= file.path(output_tag,'MCMC_parameter_all.csv'),sep=',',row.names = F, col.names = F)

# set textfile to log model output and MCMC progress
log_file_name <- file.path(output_tag,'MCMC_logfile.txt')
cat('\noutput tag:', file = log_file_name,append = T,fill = T)
cat(output_tag, file = log_file_name,append = T,fill = T)
cat('\nparameter priors:', file = log_file_name,append = T,fill = T)
cat(chains_param_file, file = log_file_name,append = T,fill = T)
cat('\nvaccine uptake:', file = log_file_name,append = T,fill = T)
cat(i_vac_file, file = log_file_name,append = T,fill = T)
cat('\n', file = log_file_name,append = T)
cat(c('\nprevious paramer set:',f_chain_id), file = log_file_name,append = T,fill = T)
cat('\n', file = log_file_name,append = T)
cat(as.character(Sys.time()), file = log_file_name,append = T,fill = T)
cat('\n', file = log_file_name,append = T)

## Important functions
##-------------------- -
### Model 1: MCMC approach
###----------------------- -
log_likelihood <- function(parms, plots_code = "TRUE"){
  
    fit = log_likelihood_model(parms, 
                            CoMix_matrices = CoMix_mat, 
                            method = "mean",
                            plots = plots_code,
                            parms_names = parms_names, 
                            ndays_sim = ndays_sim,
                            V_mat = V_mat,
                            V_mat_dose2 = V_mat_dose2,
                            V_mat_booster = V_mat_booster,
                            V_mat_2ndbooster = V_mat_2ndbooster,
                            V_mat_extrabooster = V_mat_extrabooster,
                            be_ref_data = be_ref_data)

  return(list(crit                 = fit[[sel_crit]], 
              dev                  = -2*fit[[sel_crit]],
              new_hosp_since_intro = fit$new_hosp_since_intro,
              mean_new_hosp_icu    = fit$mean_new_hosp_icu,
              mean_new_deaths      = fit$mean_new_deaths,
              mean_new_discharged  = fit$mean_new_discharged))
}

log_prior <- function(parms){ 
  log_prior_model(parms,parms_names)
}


## Starting values ----
##---------------- -
startvalue        <- parms
startvalue[parms_names_estim] <- unlist(startvalue[parms_names_estim]) * runif(length(parms_names_estim),min = 1-s_deviance,max=1+s_deviance)

# check and limit adjusted q parameters
flag_M_coef <- grepl('log_comix_coef',parms_names) | grepl('log_add_coef_mat',parms_names) | grepl('log_coef_aggr',parms_names)
startvalue[flag_M_coef & startvalue < -5] <- -4
startvalue[flag_M_coef & startvalue > 3]  <- 2.5

# check and limit adjusted transmission parameters
flag_transm <- grepl('log_VOC.*transm',parms_names)
startvalue[flag_transm & startvalue < 0] <- 0.01
startvalue[flag_transm & startvalue > 3.5]  <- 3.4

# check and limit adjusted gamma parameters
flag_gamma <- grepl('log_VOC.*gamma_factor',parms_names)
startvalue[flag_gamma & startvalue > log(1e10)] <- log(1e9)

# account for aggregated q-param
startvalue <- adjust_q_param(startvalue,parms_names,parms_names_estim,opt_waves = opt_aggr_waves)

# startvalue['log_VOC_omicron_init'] <- 1
# startvalue['VOC_omicron_start'] <- 670

# check start values, and revert if new parameter value is not valid
scm_prior <- log_prior(startvalue)
if(length(scm_prior$invalid_param)>0){
  startvalue[scm_prior$invalid_param] <- parms[scm_prior$invalid_param]
  cat('\n\nreset invalid start values :',scm_prior$invalid_param, file = log_file_name,append = T)
  #cat('\n\nreset invalid start values (warning !!!!!!!! temporary removed !!!!!):',scm_prior$invalid_param, file = log_file_name,append = T)
}

# if still invalid values present from previous calibrations, use median of all prior values if the value is present
# note: median overcomes rounding issues 
 scm_prior <- log_prior(startvalue)
 scm_prior$invalid_param <- scm_prior$invalid_param[scm_prior$invalid_param %in% names(parms_chains)]
 if(length(scm_prior$invalid_param) == 1){
   startvalue[scm_prior$invalid_param] <- median(parms_chains[,scm_prior$invalid_param])
   cat('\n\nuse overall median for invalid start values:\n',scm_prior$invalid_param, file = log_file_name,append = T)
 } else if(length(scm_prior$invalid_param) > 1){
   startvalue[scm_prior$invalid_param] <- apply(parms_chains[,scm_prior$invalid_param],2,median)
   cat('\n\nuse overall median for invalid start values:\n',scm_prior$invalid_param, file = log_file_name,append = T)
 }

scm_prior2 <- log_prior(startvalue)
cat('\n\ninitial LL:\n',scm_prior2$log_prior_val, file = log_file_name,append = T)

# save start values
vec_chain_start     = file.path(output_tag,"MCMCstart_scm.csv")
write.table(t(startvalue), file = vec_chain_start, sep=",",row.names=F)

# save parameter info
write.table(data.frame(name= parms_names,
                       value = startvalue,
                       exp_value = exp(startvalue)),
            file= file.path(output_tag,'MCMC_parameter_values.csv'),sep=',',row.names = F, col.names = T)

## LaplacesDemon estimation (MCMC sampler)
##---------------------------------------- -
Model_CB = function(pars,Data){
  
  ## select model parameters to calibrate
  pars_estim                    <- startvalue
  pars_estim[parms_names_estim] <- pars
  
  pars_estim <- adjust_q_param(pars_estim,parms_names,parms_names_estim,opt_waves = opt_aggr_waves)

  ## log(prior densities)
  LPr = log_prior(pars_estim)
  ## log-likelihood
  LL = log_likelihood(pars_estim, plots_code = "FALSE")
  ## log-posterior
  LP = LL$crit + LPr$log_prior_val
  ## Additional parameters
  R0 = LPr$R0
  
  ## output
  modelout = list(LP      = LP, 
                  Dev     = -2*LL$crit,
                  Monitor = c(pars,LP,R0,LPr$log_prior_val),
                  yhat    = LL$mean_new_hosp_icu, 
                  parm    = pars)
  return(modelout)
}

# compile R function to speed up MCMC deamon
library(compiler)
Model_CB   <- cmpfun(Model_CB)
J          <- length(parms_names_estim)       # Total number of model parameters + residual variance
mon.names  <- c(as.parm.names(list(beta=rep(0,J))),"LP","R0","LPr")
parm.names <- as.parm.names(list(beta=rep(0,J)))
MyData     <- list(J          = J, 
                   mon.names  = mon.names, 
                   parm.names = parm.names, 
                   y          = obs_total_hosp_ext)

# check initial parameter set
if(bool_plot_model_output){
  pdf(file.path(output_tag,'MCMC_SCM_start.pdf'),10,10)
  fit <- log_likelihood(as.numeric(startvalue),plots_code = TRUE)
  dev.off()
}

print("start MCMC")
## Adaptive Metropolis-within-Gibbs algorithm (using the starting values)
##---------------------------------------------------------------------- -
ptm = proc.time()
res_AMWG <- LaplacesDemon(Model          = Model_CB, 
                          Data           = MyData, 
                          Initial.Values = unlist(startvalue[parms_names_estim]),
                          Iterations     = n_iter, 
                          Status         = n_status, 
                          LogFile        = log_file_name,
                          Thinning       = 1,
                          Algorithm      = "AMWG", 
                          Specs          = list(B=NULL, n=n_repeat, Periodicity=n_period))
print(proc.time()-ptm)/60

# add info to log file
cat(as.character(Sys.time()), file = log_file_name,append = T,fill = T)
cat('\n', file = log_file_name,append = T)

# Save Deamon
saveRDS(res_AMWG,file=file.path(output_tag,'MCMC_deamon.rds'))

# save image
#save.image(file=file.path(output_tag,'MCMC_image.RData'))  #temp removed for genius (problems)

## Save estimated parameters
##-------------------- -
chain_step1       = res_AMWG$Monitor[,1:length(parms_names_estim)]
add_chain_step1   = res_AMWG$Monitor[,length(parms_names_estim)+1]
dev_chain_step1   = res_AMWG$Monitor[,length(parms_names_estim)+2]
prior_chain_step1 = res_AMWG$Monitor[,length(parms_names_estim)+3]

## merge estimated and original parameters
chain_estim        <- data.frame(chain_step1)
names(chain_estim) <- parms_names_estim
chain_other        <- t(as.matrix(startvalue[!parms_names %in% parms_names_estim]))
chain_full         <- data.frame(chain_estim,chain_other)
dim(chain_full)

# option to account for aggregated q-param
chain_full <- adjust_q_param(chain_full,parms_names,parms_names,opt_waves = opt_aggr_waves)

# reorder columns (ordered by initial column names)
chain_full     <- chain_full[,parms_names]
chain_full_out <- as.matrix(chain_full)

## explore
if(bool_plot_model_output){
  pdf(file.path(output_tag,'MCMC_SCM_final_step.pdf'),10,10)
  fit <- log_likelihood(as.numeric(chain_full[nrow(chain_full),]),plots_code = TRUE)
  dev.off()
  
  plot(res_AMWG, 
       BurnIn=0, 
       MyData, 
       PDF=TRUE, 
       Parms=NULL, 
       FileName = file.path(output_tag,'MCMC_marginal_posterior_samples.pdf'))
  
  pdf(file.path(output_tag,'MCMC_caterpillar_plot.pdf'),10,10)
  caterpillar.plot(res_AMWG, Parms="beta")
  dev.off()
  
  cat('MCMC COMPLETE', file = log_file_name,append = T)
  cat(as.character(Sys.time()), file = log_file_name,append = T,fill = T)
  cat('\n', file = log_file_name,append = T)
}

# add MCMC info ----
chain_full$mcmc_chain_id       <- f_chain_id
chain_full$mcmc_iter_id        <- 1:nrow(chain_full)
chain_full$mcmc_ll_prior       <- prior_chain_step1
chain_full$mcmc_ll_posterior   <- add_chain_step1
chain_full$mcmc_core_version   <- get_scm_version()

# reorder columns to present MCMC info first
names_all    <- names(chain_full)
names_sorted <- c(names_all[grepl('mcmc_',names_all)],names_all[!grepl('mcmc_',names_all)])
chain_full   <- chain_full[,names_sorted]

# save chain as RDS
saveRDS(chain_full,file=file.path(output_tag,'MCMC_full.rds'))

# set filenames
vec_chain_step1       = file.path(output_tag,"MCMCchain_scm.csv")
vec_dev_chain_step1   = file.path(output_tag,"MCMCdevchain_scm.csv")

# store results
write.table(chain_full,file = vec_chain_step1, sep=",",row.names=F,col.names = T)
write.table(dev_chain_step1,file = vec_dev_chain_step1, sep=",",row.names=F,col.names = F)


}

################################################################ #
## RUN THE MCMC FUNCTION ----
################################################################ #

# load this script, and run it
print("script loaded... start running")
if(length(args)==2){
  run(args[[1]],args[[2]])
} else if(length(args)==4){
  run(args[[1]],args[[2]],args[[3]],args[[4]])
} else if(length(args)==5){
  run(args[[1]],args[[2]],args[[3]],args[[4]],args[[5]])
} else if(length(args)==6){
  run(args[[1]],args[[2]],args[[3]],args[[4]],args[[5]],args[[6]])
} else if(length(args)==7){
  run(args[[1]],args[[2]],args[[3]],args[[4]],args[[5]],args[[6]],args[[7]])
} else {
  run()
}



