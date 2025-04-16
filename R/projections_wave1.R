########################################################################### #
# This file is part of the Stochastic Compartmental Model for SARS-COV-2 
# transmission in Belgium, conceived by members of SIMID group during the 
# COVID19 pandemic.
#
# This file is created to run multiple COVID-19 simulations for Belgium. 
#
# This file can be executed from the command line
# - terminal: Rscript R/projections_vaccination.R  
# - terminal: Rscript R/projections_vaccination.R &
# - within R: source('R/projections_vaccination.R')
# - within R: system('Rscript R/projections_vaccination.R &')
#
# Copyright 2021, SIMID, University of Antwerp & Hasselt University                                        
########################################################################### #

# Note: when "excecuting" an R-file, this is called line by line, which is error prone 
# when one edits the file while running a simulation (in the background). To enable 
# multiple runs of this script the same time, the code is captured in a function and 
# this function is called at the end of this script. 

# clear workspace
rm(list=ls())

## load required packages (quitely)
suppressPackageStartupMessages(library(simid.rtools))

# load functions and data 
source('R/main_vaccination.R')

# define function
projections_vaccination <- function(){

################################################################ #
## SETTINGS AND PREAMBLE ----
################################################################ #

# set output tag (= directory name)
output_tag <- 'scm23_wave1'

# get parameters: default
chains_param_file <- 'data/config/MCMCmulti_20211123_wave1.csv'

# number or model runs and time horizon
num_runs            <- 5
num_days_sim        <- 250

# FYI
if(0==1){
  sim_date2day('2020-05-01')        # help function: calendar date to model day index
  sim_day2date(61)                  # help function: model day index to calendar date
  get_M_change_day()                # contact behaviour change points (day index)
  sim_day2date(get_M_change_day())  # contact behaviour change points (calendar date)
  get_CoMix_change_day()            # CoMix contact behaviour change points (calendar date)
  get_CoMix_change_date()           # CoMix contact behaviour change points (calendar date)
  get_wave_colnames(wave_id = 1,bool_comix = TRUE)  # column names for first CoMix wave
  get_wave_colnames(wave_id = 1,bool_comix = FALSE) # column names for first additional wave
}

################################################################ #
# LOAD AND PRE-PROCESS ----
################################################################ #

# select and load parameter file
print(chains_param_file)
parms_chains     = read.table(chains_param_file, sep = ",", header = T)

parms_chains <- parms_chains[,! names(parms_chains) %in% c(get_wave_colnames(3:8,bool_comix = FALSE),
                                                           get_wave_colnames(9:33,bool_comix = TRUE)
                                                           )]

# select a pareter set for each model run
# MCMC estimation provides a sorted file, on iteration and log-likelihood. As such,
# the bottom row provides the most optimal parameter set of the last iteration.
parms_chains     = parms_chains[nrow(parms_chains) - (0:num_runs),]

# vaccine uptake and region not applicable for now ----
vacc_files <- NA
sel_region <- 'belgium'

# set output directory
exp_tag    <- format(Sys.time(),format = '%Y%m%d_%H%M%S')
output_dir <- paste0('output/',output_tag,'/',exp_tag,'_d',num_days_sim,'_n',num_runs)
output_dir <- paste0(output_dir,'/')

# create output directory if not existing yet
if(!dir.exists(output_dir)){
  dir.create(output_dir,recursive = T)
}
print(output_dir)

# add adjusted parameter config to output folder
write.table(parms_chains,
            file=file.path(output_dir,basename(chains_param_file)),
            sep=',',row.names = F)


# load reference data
be_ref_data <- get_latest_incidence_data()

################################################################ #
# MAIN: CALL MODEL ----
################################################################ #

# explore contact and Q(age) parameters
pdf(file=paste0(output_dir,'next_gen.pdf'),14,7)
explore_next_gen(parms_chains)
dev.off()

# initialise summary variables
i_seq             <- round(seq(1,nrow(parms_chains),length.out = num_runs))
hosp_age_adm      <- array(dim = c(num_runs, num_days_sim, 11))
scen_cases        <- array(dim = c(num_runs, num_days_sim, 11))
mortality_age     <- array(dim = c(num_runs, num_days_sim, 11))
cumul_cases       <- array(dim = c(num_runs, num_days_sim, 11))

prev_I_asymp_sims <- array(dim = c(num_runs, num_days_sim, 10))
prev_I_mild_sims  <- array(dim = c(num_runs, num_days_sim, 10))
prev_I_severe_sims<- array(dim = c(num_runs, num_days_sim, 10))
prev_I_hosp_sims  <- array(dim = c(num_runs, num_days_sim, 10))
inc_I_hosp_sims   <- array(dim = c(num_runs, num_days_sim, 10))

mean_susceptible_full <- array(dim = c(num_runs, num_days_sim, 11))
mean_susceptible_vac_d1 <- array(dim = c(num_runs, num_days_sim, 11))
mean_susceptible_vac_d2 <- array(dim = c(num_runs, num_days_sim, 11))
mean_vaccinated_age  <- array(dim = c(num_runs, num_days_sim, 11))

scen_hosp_adm     <- matrix(NA, nrow= num_days_sim, ncol =num_runs) # for post processing
scen_Rt_infection <- matrix(NA, nrow= num_days_sim, ncol =num_runs)
scen_Rt_mild_sev  <- matrix(NA, nrow= num_days_sim, ncol =num_runs)
scen_Rt_hosp      <- matrix(NA, nrow= num_days_sim, ncol =num_runs)
scen_new_infect   <- matrix(NA, nrow= num_days_sim, ncol =num_runs)
scen_new_sympt    <- matrix(NA, nrow= num_days_sim, ncol =num_runs)
scen_mortality    <- matrix(NA, nrow= num_days_sim, ncol =num_runs)
scen_VOC_mild_alpha     <- matrix(NA, nrow= num_days_sim, ncol =num_runs)
scen_VOC_mild_delta     <- matrix(NA, nrow= num_days_sim, ncol =num_runs)

scen_hosp_load    <- matrix(NA, nrow= num_days_sim, ncol =num_runs) 
scen_ICU_load     <- matrix(NA, nrow= num_days_sim, ncol =num_runs) 

scen_hosp_exit    <- matrix(NA, nrow= num_days_sim, ncol =num_runs) 
  
nvac_hosp_adm      <- matrix(NA, nrow= num_days_sim, ncol =num_runs) 
nvac_age_hosp_adm  <- array(dim = c(num_runs, num_days_sim, 11))
new_hosp_icu_vac_d2 <- matrix(NA, nrow= num_days_sim, ncol =num_runs) 
new_hosp_icu_vac_age_d2 <- array(dim = c(num_runs, num_days_sim, 11))

mean_rna1_vaccinated    <- array(dim = c(num_runs, num_days_sim, 11)) 
mean_rna2_vaccinated    <- array(dim = c(num_runs, num_days_sim, 11)) 
mean_adeno1_vaccinated  <- array(dim = c(num_runs, num_days_sim, 11))
mean_adeno2_vaccinated  <- array(dim = c(num_runs, num_days_sim, 11))
mean_waning_vaccinated  <- array(dim = c(num_runs, num_days_sim, 11))

i <- 1
time_stamp_loop <- Sys.time()
for (i in 1:num_runs){
  set.seed(20210101 + i)
  
  smd_print_progress(i,num_runs,time_stamp_loop=time_stamp_loop)
  run_param <- unlist(parms_chains[nrow(parms_chains)-(i_seq[i]-1), ])
  
  pred_all <- run_model(parms = run_param, 
                     CoMix_matrices = CoMix_mat, 
                     method = "stochastic", 
                     ndays_sim = num_days_sim)
  
  # hospital admissions
  pred_hosp1          <- pred_all$total_new_hosp_icu
  hosp_age_adm[i,,]   <- as.matrix(pred_hosp1)       
  scen_hosp_adm[,i]   <- as.matrix(rowSums(pred_all$total_new_hosp_icu[,-1]))

  # infections
  new_cases           <- pred_all$total_new_infections
  scen_cases[i,,]     <- as.matrix(new_cases)      
  scen_new_infect[,i] <- as.matrix(rowSums(new_cases[,-1])) # new_cases
  
  # sympt cases
  new_mild_cases       <- pred_all$total_new_mild_infections
  scen_new_sympt[,i]   <- as.matrix(rowSums(new_mild_cases[,-1])) # new_cases
  
  # hospital and icu load 
  scen_hosp_load[,i]   <- as.matrix(rowSums(pred_all$hosp_icu_load[,-1]))
  scen_ICU_load[,i]    <- as.matrix(rowSums(pred_all$icu_load[,-1]))

  # prevalence
  prev_I_asymp_sims[i,,]   <- as.matrix(pred_all$prev_I_asym[,-1])
  prev_I_mild_sims[i,,]    <- as.matrix(pred_all$prev_I_mild[,-1])
  prev_I_severe_sims[i,,]  <- as.matrix(pred_all$prev_I_sev[,-1])
  prev_I_hosp_sims[i,,]    <- as.matrix(pred_all$hosp_icu_load[,-1])
  inc_I_hosp_sims[i,,]     <- as.matrix(pred_all$total_new_hosp_icu[,-1])
  
  # hospital exit
  scen_hosp_exit[,i]   <- as.matrix(rowSums(pred_all$total_hosp_exit[,-1]))
  
  # Rt
  scen_Rt_infection[,i] <- get_Rt(vect_cases = rowSums(new_cases[,-1]),
                                  vect_dates = new_cases$day)$Rt
  scen_Rt_mild_sev[,i] <- get_Rt(vect_cases = rowSums(new_mild_cases[,-1]),
                                vect_dates = new_mild_cases$day)$Rt
  scen_Rt_hosp[,i]     <- get_Rt(vect_cases = rowSums(pred_hosp1[,-1]),
                                 vect_dates = pred_hosp1[,1])$Rt
  
  # mortality
  scen_mortality[,i]   <- as.matrix(rowSums(pred_all$total_new_deaths[,-1]))
  mortality_age[i,,]   <- as.matrix(pred_all$total_new_deaths) 
  
  # VOC
  scen_VOC_mild_alpha[,i]   <- as.matrix(pred_all$p_B117)
  scen_VOC_mild_delta[,i]   <- as.matrix(pred_all$p_delta)

  # cumulative infections
  scen_cases_cum     <- apply(as.matrix(new_cases),2,cumsum)
  scen_cases_cum[,1] <- as.matrix(new_cases)[,1]
  cumul_cases[i,,]   <- scen_cases_cum
  
  # susceptible
  mean_susceptible_full[i,,]     <- as.matrix(pred_all$mean_susceptible_full)  
  mean_susceptible_vac_d1[i,,]   <- as.matrix(pred_all$mean_susceptible_vac_d1)  
  mean_susceptible_vac_d2[i,,]   <- as.matrix(pred_all$mean_susceptible_vac_d2)  
  mean_vaccinated_age[i,,]       <- as.matrix(pred_all$mean_vaccinated_age)

  # non-vaccinated hosp admissions
  nvac_hosp_adm[,i]   <- as.matrix(rowSums(pred_all$total_new_hosp_icu_nvac[,-1]))
  nvac_age_hosp_adm[i,,]   <- as.matrix(pred_all$total_new_hosp_icu_nvac)
  new_hosp_icu_vac_d2[,i] <- as.matrix(rowSums(pred_all$total_new_hosp_icu_vac_d2[,-1]))
  new_hosp_icu_vac_age_d2[i,,] <- as.matrix(pred_all$total_new_hosp_icu_vac_d2)
  
  # vaccinated
  # note !flag_rna and !flag_d1 are not working anymore with potential booster doses
  flag_rna     <- grepl('rna',names(pred_all$mean_dose_vaccinated))
  flag_adeno   <- grepl('adeno',names(pred_all$mean_dose_vaccinated))  
  flag_d1      <- grepl('d1',names(pred_all$mean_dose_vaccinated))
  flag_d2      <- grepl('d2',names(pred_all$mean_dose_vaccinated))
  flag_booster <- grepl('booster',names(pred_all$mean_dose_vaccinated))
  names(pred_all$mean_dose_vaccinated)[flag_rna & flag_d1]
  mean_rna1_vaccinated[i,,]    <- as.matrix(pred_all$mean_dose_vaccinated[,c(1,which(flag_rna & flag_d1))])
  mean_rna2_vaccinated[i,,]    <- as.matrix(pred_all$mean_dose_vaccinated[,c(1,which(flag_rna & flag_d2))])
  mean_adeno1_vaccinated[i,,]  <- as.matrix(pred_all$mean_dose_vaccinated[,c(1,which(flag_adeno & flag_d1))])
  mean_adeno2_vaccinated[i,,]  <- as.matrix(pred_all$mean_dose_vaccinated[,c(1,which(flag_adeno & flag_d2))])
  mean_waning_vaccinated[i,,]  <- as.matrix(pred_all$mean_dose_vaccinated[,c(1,which(flag_booster))])
}
 
saveRDS(hosp_age_adm,file=paste0(output_dir,'scenario1.rds'))
saveRDS(hosp_age_adm,file=paste0(output_dir,'hosp_age_adm.rds'))
saveRDS(scen_cases,file=paste0(output_dir,'scen_cases.rds'))
saveRDS(cumul_cases,file=paste0(output_dir,'cumul_cases.rds'))
saveRDS(mean_susceptible_full,file=paste0(output_dir,'mean_susceptible_full.rds'))
saveRDS(mean_susceptible_vac_d1,file=paste0(output_dir,'mean_susceptible_vac_d1.rds'))
saveRDS(mean_susceptible_vac_d2,file=paste0(output_dir,'mean_susceptible_vac_d2.rds'))
saveRDS(mean_vaccinated_age,file=paste0(output_dir,'mean_vaccinated_age.rds'))

saveRDS(mean_rna1_vaccinated,file=paste0(output_dir,'mean_rna1_vaccinated.rds'))
saveRDS(mean_rna2_vaccinated,file=paste0(output_dir,'mean_rna2_vaccinated.rds'))
saveRDS(mean_adeno1_vaccinated,file=paste0(output_dir,'mean_adeno1_vaccinated.rds'))
saveRDS(mean_adeno2_vaccinated,file=paste0(output_dir,'mean_adeno2_vaccinated.rds'))
saveRDS(mean_waning_vaccinated,file=paste0(output_dir,'mean_waning_vaccinated.rds'))

saveRDS(scen_new_infect,file=paste0(output_dir,'new_infect.rds'))
saveRDS(scen_new_sympt,file=paste0(output_dir,'new_sympt.rds'))
saveRDS(scen_hosp_load,file=paste0(output_dir,'hosp_load.rds'))
saveRDS(scen_Rt_infection,file=paste0(output_dir,'scen_Rt_infection.rds'))
saveRDS(scen_Rt_mild_sev,file=paste0(output_dir,'scen_Rt_mild_sev.rds'))
saveRDS(scen_Rt_hosp,file=paste0(output_dir,'scen_Rt_hosp.rds'))
saveRDS(scen_hosp_adm,file=paste0(output_dir,'hosp_adm.rds'))
saveRDS(scen_ICU_load,file=paste0(output_dir,'icu_load.rds'))
saveRDS(scen_hosp_exit,file=paste0(output_dir,'hosp_exit.rds'))
saveRDS(scen_mortality,file=paste0(output_dir,'scen_mortality.rds'))
saveRDS(mortality_age,file=paste0(output_dir,'mortality_age.rds'))
saveRDS(scen_VOC_mild_alpha,file=paste0(output_dir,'scen_VOC_mild_alpha.rds'))
saveRDS(scen_VOC_mild_delta,file=paste0(output_dir,'scen_VOC_mild_delta.rds'))

# burden of disease
saveRDS(prev_I_asymp_sims,file=paste0(output_dir,'prev_I_asymp_sims.rds'))
saveRDS(prev_I_mild_sims,file=paste0(output_dir,'prev_I_mild_sims.rds'))
saveRDS(prev_I_severe_sims,file=paste0(output_dir,'prev_I_severe_sims.rds'))
saveRDS(prev_I_hosp_sims,file=paste0(output_dir,'prev_I_hosp_sims.rds'))
saveRDS(inc_I_hosp_sims,file=paste0(output_dir,'inc_I_hosp_sims.rds'))

vac_hosp_adm    <- scen_hosp_adm - nvac_hosp_adm
saveRDS(nvac_hosp_adm,file=paste0(output_dir,'nvac_hosp_adm.rds'))
saveRDS(vac_hosp_adm,file=paste0(output_dir,'vac_hosp_adm.rds'))

vac_age_hosp_adm    <- hosp_age_adm - nvac_age_hosp_adm
vac_age_hosp_adm[,,1] <- hosp_age_adm[,,1]
saveRDS(vac_age_hosp_adm,file=paste0(output_dir,'vac_age_hosp_adm.rds'))
saveRDS(nvac_age_hosp_adm,file=paste0(output_dir,'nvac_age_hosp_adm.rds'))

saveRDS(new_hosp_icu_vac_d2,file=paste0(output_dir,'new_hosp_adm_vac_d2.rds'))
saveRDS(new_hosp_icu_vac_age_d2,file=paste0(output_dir,'new_hosp_adm_vac_age_d2.rds'))


## plot results ----
pdf(file=paste0(output_dir,'projections.pdf'))

x_axis_scm_day = c(0,num_days_sim)
x_axis_vaccine = c(350,num_days_sim)

# set default plot color and pch
col_lib <- data.frame(tags = c('xxxx',
                               output_tag),
                      col  = c('black','darkred'),
                      label = c(paste0('Reported (',max(be_ref_data$date),')'),
                                paste('Scenario')),
                      lwd = c(NA,2),
                      pch = c(1,NA))

output_files <- dir(output_dir,full.names = T)
multi_plot_incidence_time(output_files        = output_files,
                          col_lib             = col_lib,
                          be_ref_data         = be_ref_data,
                          x_axis_scm_day      = x_axis_scm_day)
dev.off()

# # optional: plot age-specific results
# pdf(file=paste0(output_dir,'scenario_cases_mean_age.pdf'))
# multi_plot_incidence_age_time(output_files   = output_files,
#                               col_lib        = col_lib,
#                               x_axis_scm_day = x_axis_scm_day,
#                               x_axis_vaccine = x_axis_vaccine)
#
#dev.off() # close PDF stream
  
}

# call function
projections_vaccination()
