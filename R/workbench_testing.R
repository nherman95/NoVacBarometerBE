########################################################################### #
# This file is part of the Stochastic Compartmental Model for SARS-COV-2 
# transmission in Belgium, conceived by members of SIMID group during the 
# COVID19 pandemic.
#
# This file can be used to run the core function with and without log-likelihood 
# calculation and comparison with previously stored results as unit-test.
#
# Copyright 2022, SIMID, University of Antwerp & Hasselt University                                        
########################################################################### #

rm(list=ls())

# load functions and data ----
source('R/main_vaccination.R')
source('R/lib_unit_testing.R')

################################################################ #
## SETTINGS AND PREAMBLE ----
################################################################ #

# set parameter file
chains_param_file <- 'data/config/MCMCmulti_20220630_d850_e70_i100_n10_p10_crit1_186158_belgium_recap7.csv'

# boolean to enable the comparison of the new output with the previously stored output (FALSE allows to store reference output)
bool_compare <- TRUE

# load file
chains_param      = read.table(chains_param_file, sep = ",", header = T)

# select one parameter set
parms      <- unlist(chains_param[1,]) #nrow(chains_param)-4
length(parms)

# store parameter names
parms_names = names(parms)

# set rng seed
set.seed(666)

# set plot configuration
plots_code <- TRUE

# set model output configuration
bool_full_output <- TRUE

# set vaccine uptake scheme
#vacc_schedule      <- read.table('data/uptake/uptake_booster/vaccine_uptake_booster_vSCENapr05_uadult80_uchild40_ubooster75_mar01_belgium.csv', sep=',',header = T)
vacc_schedule      <- read.table("data/uptake/uptake_manuscript//vaccine_uptake_booster_vSCENmay06_uadult80_uchild40_ubooster60_jun01_belgium.csv", sep=',',header = T)

################################################################ #
## PARAMETER UPDATES ----
################################################################ #

parms       <- update2latest_model_parameter_config(parms)
parms_names <- names(parms)

################################################################ #
## RUN ----
################################################################ #

# Deterministic (mean) model ----
print(paste(Sys.time(),"START DETERMINISTIC MODEL (v.x)")); ptm = proc.time()
fit_mean_scm = log_likelihood_model(parms, 
                                CoMix_matrices = CoMix_mat,
                                method = "mean",
                                plots = plots_code,
                                parms_names = names(parms),
                                vaccine_uptake = vacc_schedule)
print(proc.time()-ptm) ;
if(bool_compare) compare_output(fit_mean_scm,'mean',prev_file_names)

# Stochastic model ----
print(paste(Sys.time(),"START STOCHASTIC MODEL (v.x)")); ptm = proc.time()
fit_stochastic_scm = log_likelihood_model(parms, 
                                      CoMix_matrices = CoMix_mat,
                                      method = "stochastic",
                                      plots = plots_code,
                                      parms_names = names(parms),
                                      vaccine_uptake = vacc_schedule)
print(proc.time()-ptm)
if(bool_compare) compare_output(fit_stochastic_scm,'stochastic',prev_file_names)

## MODEL EXPLORATION ----
print(paste(Sys.time(),"EXPLORE STOCHASTIC MODEL (v.x)")); ptm = proc.time()
scm_out = run_model(parms,
                CoMix_matrices = CoMix_mat,
                method = "stochastic",
                parms_names = names(parms),
                ndays_sim = 850,
                vaccine_uptake = vacc_schedule)

print(proc.time()-ptm) ;
if(bool_compare) compare_output(scm_out,'scm_out',prev_file_names)

print(paste(Sys.time(),"EXPLORE DETERMINISTIC MODEL (v.x)")); ptm = proc.time()
scm_mean = run_model(parms,
                CoMix_matrices = CoMix_mat,
                method = "mean",
                parms_names = names(parms),
                ndays_sim = 900,
                vaccine_uptake = vacc_schedule)

print(proc.time()-ptm) ;
if(bool_compare) compare_output(scm_mean,'scm_mean',prev_file_names)

## PRIOR ----
print(paste(Sys.time(),"EXPLORE STOCHASTIC MODEL (v.x)")); ptm = proc.time()
scm_prior = log_prior_model(parms,
                            parms_names = names(parms))

if(any(grepl('invalid_param',names(scm_prior)))){
  print(scm_prior$invalid_param)
  scm_prior <- scm_prior[1:2]
}

print(proc.time()-ptm) ;
if(bool_compare) compare_output(scm_prior,'scm_prior',prev_file_names)

## EXTENDED TIME HORIZON
print(paste(Sys.time(),"EXPLORE EXTENDED TIME HORIZON (v.x)")); ptm = proc.time()
fit_extended_scm = log_likelihood_model(parms, 
                                  CoMix_matrices = CoMix_mat,
                                  method = "mean",
                                  plots = TRUE,
                                  parms_names = names(parms),
                                  ndays_sim = 850,
                                  vaccine_uptake = vacc_schedule)
print(proc.time()-ptm) ;
if(bool_compare) compare_output(fit_extended_scm,'extended',prev_file_names)

# # Regional sochastic models (w/o value testing)---- 
# print(paste(Sys.time(),"START REGIONAL STOCHASTIC MODELS (v.x)")); ptm = proc.time()
# parms_region <- parms
# parms_region['region_id'] = 2
# be_ref_data_region <- get_latest_incidence_data(sel_region = "brussels")
# fit_stochastic_bxl = log_likelihood_model(parms_region, 
#                                       CoMix_matrices = CoMix_mat,
#                                       method = "stochastic",
#                                       plots = plots_code,
#                                       parms_names = names(parms),
#                                       be_ref_data = be_ref_data_region,
#                                       vaccine_uptake = vacc_schedule)
# 
# parms_region['region_id'] = 3
# be_ref_data_region <- get_latest_incidence_data(sel_region = "flanders")
# fit_stochastic_fl = log_likelihood_model(parms_region, 
#                                       CoMix_matrices = CoMix_mat,
#                                       method = "stochastic",
#                                       plots = plots_code,
#                                       parms_names = names(parms),
#                                       be_ref_data = be_ref_data_region,
#                                       vaccine_uptake = vacc_schedule)
# print(proc.time()-ptm)

################################################################ #
## BENCHMARK ----
################################################################ #

# note: reformating the uptake once, instead of every iteration, saves 30% time! 
# decreased from 1.2s to 0.8s per run on Nov 30, 2021
# increased to 0.94s with waning immunity (default reference data)
# 1.3s with waning immunity and latest reference data (December 1st)
V_mat         <- get_uptake_matrix(vaccine_uptake = vacc_schedule, dose_tag = '_A_', parms = parms )
V_mat_dose2   <- get_uptake_matrix(vaccine_uptake = vacc_schedule, dose_tag = '_B_', parms = parms )
V_mat_booster <- get_uptake_matrix(vaccine_uptake = vacc_schedule, dose_tag = '_E_', parms = parms )

print(paste(Sys.time(),"BENCHMARK STOCHASTIC MODEL (v.x)")); ptm = proc.time()
nrun <- 10
for(i in 1:nrun){
  fit_mean_bench = log_likelihood_model(parms, 
                                     CoMix_matrices = CoMix_mat,
                                     method = "mean",
                                     plots = FALSE,
                                     parms_names = names(parms),
                                     #vaccine_uptake = vacc_schedule,
                                     #be_ref_data = be_ref_data,
                                     V_mat = V_mat,
                                     V_mat_dose2 = V_mat_dose2,
                                     V_mat_booster = V_mat_booster)
}
print((total_time = proc.time()-ptm));
print(total_time/nrun) ;
print(paste(Sys.time(),"WORKBENCH FINISHED"))


################################################################ #
# # reset reference values (optional)
################################################################ #
# rrv()






