########################################################################### #
# This file is part of the Stochastic Compartmental Model for SARS-COV-2 
# transmission in Belgium, conceived by members of SIMID group during the 
# COVID19 pandemic.
#
# This file contains the main function to run multiple simulation runs with 
# different social contact and/or vaccination uptake assumptions for Belgium 
# or one of the three regions. 
#
# Current parallelism: for different vaccine uptake file
#
# Copyright 2022, SIMID, University of Antwerp & Hasselt University                                        
########################################################################### #

# Note: when "excecuting" an R-file, this is called line by line, which is error prone 
# when one edits the file while running a simulation (in the background). To enable 
# multiple runs of this script the same time, the code is captured in a function and 
# this function is called at the end of this script. 

# define function
projections_vaccination <- function(output_tag,
                                    sel_region,
                                    chains_param_files,
                                    cnt_adjust_value,
                                    cnt_adjust_date,
                                    scen_tag,
                                    num_chains,
                                    num_stochastic_real,
                                    num_days_sim,
                                    vacc_files,
                                    adjusted_parameters = NULL,
                                    cl_args){

################################################################ #
## SETTINGS AND PREAMBLE ----
################################################################ #

# use function arguments
  
# use default y axis?
default_y_axis <- TRUE

################################################################ #
# LOAD AND PRE-PROCESS ----
################################################################ #

# adjust parameters based on optional command line arguments
if(length(cl_args)>=1){ # region
  sel_region <- cl_args[1]
  }
 
 
if(length(cl_args)>=2){ # social/risk behaviour
  cnt_adjust_str <- cl_args[2]
  cnt_adjust_value <- as.numeric(unlist(strsplit(cnt_adjust_str,'_')))
  
  if(length(cnt_adjust_value) != length(cnt_adjust_date)){
    warning('GIVEN CONTACT ADJUST VALES ARE NOT IN LINE WITH CONTACT ADJUST DATES, ADJUSTED DATE LIST',immediate. = T)
    cnt_adjust_date <- cnt_adjust_date[1:length(cnt_adjust_value)]
  }
}

# if(length(cnt_adjust_value)>=2 & cnt_adjust_value[2]>0){
#   scen_tag      <- 'sep20'
#   num_days_sim  <- 320
#   vacc_files <- vacc_files[2]
# }
# 
# if(length(cnt_adjust_date)>=4 & cnt_adjust_value[4]>0){
#   scen_tag      <- 'sep21'
#   vacc_files <- vacc_files[2]
# }
# 
# if(length(cnt_adjust_date)>=5 & cnt_adjust_value[5]>0){
#   scen_tag      <- 'sep21'
#   vacc_files <- vacc_files[2]
# }

# include contact scenario info into output tag
cnt_adjust_tag <- paste(round(cnt_adjust_value*100),collapse='_')

# convert cnt_adjust dates into simulation day
cnt_adjust_day   <- sim_date2day(cnt_adjust_date)

# select and load parameter file
chains_param_file <- chains_param_files[grepl(sel_region,chains_param_files) ]

if(length(chains_param_file) == 0){
  chains_param_file <- chains_param_files[1] 
}

print(chains_param_file)
parms_chains     = read.table(chains_param_file, sep = ",", header = T)

# select the parameter sets and make a copy for each stochastic realisations
parms_chains  <- parms_chains[rep(1:num_chains,num_stochastic_real),]
num_runs      <- nrow(parms_chains)

# option to adjust some specific parameters
if(!is.null(adjusted_parameters) && length(adjusted_parameters)>0){
  for(i_param in 1:length(adjusted_parameters)){
    if(names(adjusted_parameters)[i_param] %in% names (parms_chains)){
      parms_chains[names(adjusted_parameters)[i_param]] <- adjusted_parameters[i_param]
    } else{
      warning('Issue in function projections_vaccination: "adjusted_parameters" contains a parameter that is not present in "parms_chains"')
    }
  }
}

# set parameter names
parms_names <- names(parms_chains)
if(!exists('sel_region')){
  sel_region <- get_region(unique(parms_chains[,'region_id']))
} else{
  parms_chains[,'region_id'] <- get_region_id(sel_region)
}

# set output directory
exp_dir <- gsub('data/config','',chains_param_file)
#exp_dir <- basename(dirname(dirname(exp_dir)))
exp_dir <- basename((dirname(exp_dir)))
output_dir <- paste0('output/',output_tag,'/',exp_dir,'/',scen_tag,'_',num_days_sim,'_cnt',cnt_adjust_tag,'_n',num_runs,'_s',num_stochastic_real,'_',sel_region)
output_dir <- paste0(output_dir,'/')
output_vac_dir <- paste0('output/',output_tag,'/',exp_dir,'/','vac','_',num_days_sim,'_cnt',cnt_adjust_tag,'_n',num_runs,'_s',num_stochastic_real,'_',sel_region)
output_vac_dir <- paste0(output_vac_dir,'/')
# create output directory if not existing yet
if(!dir.exists(output_dir)){
  dir.create(output_dir,recursive = T)
}
print(output_dir)

# aadd adjusted contact behaviour to the model parameter matrix
for(i_cp in 1:length(cnt_adjust_value)){
  if(cnt_adjust_value[i_cp]!=0){
    parms_chains[,paste0('cnt_adjust_value',i_cp)] <- cnt_adjust_value[i_cp]
    parms_chains[,paste0('cnt_adjust_day',i_cp)]   <- cnt_adjust_day[i_cp]
  }
}

# add adjusted parameter config to output folder
write.table(parms_chains,
            file=file.path(output_dir,basename(chains_param_file)),
            sep=',',row.names = F)


# load reference data
if(sel_region != 'belgium'){
  be_ref_data <- get_latest_incidence_data(sel_region = sel_region)
} else{
  be_ref_data <- get_latest_incidence_data(enable_data_reuse=T)
}

## select and load vaccine uptake
vacc_files <- vacc_files[(grepl(sel_region,vacc_files) )]
print(vacc_files)
vacc_refs <- dir('data/uptake/uptake_fall23campaign',pattern = 'uptake_extrabooster.*csv',full.names = T)
vacc_refs <- vacc_refs[c(1)]
vacc_file_tag <- tolower(gsub('\\.csv','',gsub('.*Output_v','',vacc_files)))
vacc_file_tag <- basename(vacc_file_tag) #tmp fix when using different uptake files
vacc_refs_tag <- tolower(gsub('\\.csv','',gsub('.*Output_v','',vacc_refs)))
vacc_refs_tag <- basename(vacc_refs_tag) #tmp fix when using different uptake files

# set default plot color and pch
col_lib <- data.frame(tags = c('xxxx',
                               cnt_adjust_tag,
                               'mild_alpha',
                               'mild_delta'),
                      col  = c('black','darkred','darkblue','darkgreen'),
                      label = c(paste0('Reported (',max(be_ref_data$date),')'),
                                paste('Scenario:',cnt_adjust_tag),
                                'Alpha+Beta+Gamma',
                                'Delta'),
                      lwd = c(NA,2,2,2),
                      pch = c(1,NA,NA,NA))

# set colors for post-processing figures
vacc_colors    <- vacc_file_tag
is_invalid_col <- !vacc_colors %in% colors()
vacc_colors[is_invalid_col] <- (1:sum(is_invalid_col))+1

################################################################ #
# MAIN: CALL MODEL ----
################################################################ #

# explore contact and Q(age) parameters
pdf(file=paste0(output_dir,'next_gen_',gsub('\\.','p',cnt_adjust_tag),'.pdf'),14,7)
explore_next_gen(parms_chains)
dev.off()

# # start parallel nodes ----
#smd_start_cluster(num_proc = min(length(vacc_files),10),
#                  timeout  = 600) # 10min

# run ----
i_vac <- 1
# foreach(i_vac = 1:length(vacc_files),
#         .export = c('par_nodes_info'),
#         .packages = c('EpiEstim',
#                       'simid.rtools',
#                       'zoo',
#                       'openxlsx',
#                       'scales',
#                       'RColorBrewer')
#         )   #%dopar%
{
  # temp: to cope with the parallel environment
  source('R/main_vaccination.R')

  i_vac_file        <- vacc_files[i_vac]
  vacc_schedule     <- read.csv(i_vac_file, header = T)
  vacc_schedule_tag <- vacc_file_tag[i_vac]
  vacc_ref_tag.     <- vacc_refs_tag[1]
  col_scen          <- vacc_colors[i_vac]
  
  # fix for color_"tag"
  col_scen          <- unlist(strsplit(col_scen,'_'))[1]
  
  # initialise summary variables
  i_seq             <- round(seq(1,nrow(parms_chains),length.out = num_runs))
  hosp_age_adm      <- array(dim = c(num_runs, num_days_sim, 11))
  hosp_age_adforcovid      <- array(dim = c(num_runs, num_days_sim, 11))
  scen_cases        <- array(dim = c(num_runs, num_days_sim, 11))
  scen_mild_cases   <- array(dim = c(num_runs, num_days_sim, 11))
  mortality_age     <- array(dim = c(num_runs, num_days_sim, 11))
  cumul_cases       <- array(dim = c(num_runs, num_days_sim, 11))
  cumul_deaths       <- array(dim = c(num_runs, num_days_sim, 11))
  scen_numvac        <- array(dim = c(num_runs, num_days_sim, 11))
  
  prev_I_asymp_sims <- array(dim = c(num_runs, num_days_sim, 10))
  prev_I_mild_sims  <- array(dim = c(num_runs, num_days_sim, 10))
  prev_I_severe_sims<- array(dim = c(num_runs, num_days_sim, 10))
  prev_I_hosp_sims  <- array(dim = c(num_runs, num_days_sim, 10))
  prev_I_icu_sims   <- array(dim = c(num_runs, num_days_sim, 10)) #new
  inc_I_hosp_sims   <- array(dim = c(num_runs, num_days_sim, 10))
  inc_I_icu_sims    <- array(dim = c(num_runs, num_days_sim, 10)) #new
  inc_I_hosp_for_covid_sims   <- array(dim = c(num_runs, num_days_sim, 10))
  
  mean_susceptible_full <- array(dim = c(num_runs, num_days_sim, 11))
  mean_susceptible_initpluswaning <- array(dim = c(num_runs, num_days_sim, 11))
  mean_susceptible_vac_d1 <- array(dim = c(num_runs, num_days_sim, 11))
  mean_susceptible_vac_d2 <- array(dim = c(num_runs, num_days_sim, 11))
  mean_vaccinated_age  <- array(dim = c(num_runs, num_days_sim, 11))

  scen_hosp_adm     <- matrix(NA, nrow= num_days_sim, ncol =num_runs) # for post processing
  scen_hosp_adforcovid     <- matrix(NA, nrow= num_days_sim, ncol =num_runs) # for post processing
  scen_Rt_infection <- matrix(NA, nrow= num_days_sim, ncol =num_runs)
  scen_Rt_mild_sev  <- matrix(NA, nrow= num_days_sim, ncol =num_runs)
  scen_Rt_hosp      <- matrix(NA, nrow= num_days_sim, ncol =num_runs)
  scen_new_infect   <- matrix(NA, nrow= num_days_sim, ncol =num_runs)
  scen_new_reinfect <- matrix(NA, nrow= num_days_sim, ncol =num_runs)
  scen_new_sympt    <- matrix(NA, nrow= num_days_sim, ncol =num_runs)
  scen_mortality    <- matrix(NA, nrow= num_days_sim, ncol =num_runs)
  scen_VOC_mild_alpha     <- matrix(NA, nrow= num_days_sim, ncol =num_runs)
  scen_VOC_mild_delta     <- matrix(NA, nrow= num_days_sim, ncol =num_runs)
  scen_VOC_mild_omicron   <- matrix(NA, nrow= num_days_sim, ncol =num_runs) #BA1BA2
  scen_VOC_mild_ba4ba5    <- matrix(NA, nrow= num_days_sim, ncol =num_runs)
  scen_VOC_mild_bq1ba275xbb    <- matrix(NA, nrow= num_days_sim, ncol =num_runs)
  scen_VOC_mild_xbb15    <- matrix(NA, nrow= num_days_sim, ncol =num_runs) #XBB15 include XBB1.9
  
  scen_hosp_load    <- matrix(NA, nrow= num_days_sim, ncol =num_runs) 
  scen_ICU_load     <- matrix(NA, nrow= num_days_sim, ncol =num_runs) 
  scen_ICU_adm      <- matrix(NA, nrow= num_days_sim, ncol =num_runs) 

  scen_hosp_exit    <- matrix(NA, nrow= num_days_sim, ncol =num_runs) 
  social_coef    <- matrix(NA, nrow= num_days_sim, ncol =num_runs) 
    
  nvac_hosp_adm      <- matrix(NA, nrow= num_days_sim, ncol =num_runs) 
  nvac_age_hosp_adm  <- array(dim = c(num_runs, num_days_sim, 11))
  new_hosp_icu_vac_d2 <- matrix(NA, nrow= num_days_sim, ncol =num_runs) 
  new_hosp_icu_vac_age_d2 <- array(dim = c(num_runs, num_days_sim, 11))
  new_hosp_icu_vac_booster<- matrix(NA, nrow= num_days_sim, ncol =num_runs) 
  
  mean_rna1_vaccinated    <- array(dim = c(num_runs, num_days_sim, 11)) 
  mean_rna2_vaccinated    <- array(dim = c(num_runs, num_days_sim, 11)) 
  mean_adeno1_vaccinated  <- array(dim = c(num_runs, num_days_sim, 11))
  mean_adeno2_vaccinated  <- array(dim = c(num_runs, num_days_sim, 11))
  mean_waning_vaccinated  <- array(dim = c(num_runs, num_days_sim, 11))
    
  i <- 1
  time_stamp_loop <- Sys.time()
  for (i in 1:num_runs){
    set.seed(20210101 + i)
    print("RUN")
    print(i)
    #smd_print_progress(i,num_runs,time_stamp_loop=time_stamp_loop,par_nodes_info = par_nodes_info)
    run_param <- unlist(parms_chains[nrow(parms_chains)-(i_seq[i]-1), ])
    
    pred_all <- run_model(parms = run_param, 
                       CoMix_matrices = CoMix_mat, 
                       method = "stochastic", 
                       ndays_sim = num_days_sim,
                       vaccine_uptake = vacc_schedule,
                       cnt_change_pt = cnt_adjust_day, cnt_change = cnt_adjust_value)
    
    # hospital admissions
    pred_hosp1          <- pred_all$total_new_hosp_icu
    hosp_age_adm[i,,]   <- as.matrix(pred_hosp1)       
    scen_hosp_adm[,i]   <- as.matrix(rowSums(pred_all$total_new_hosp_icu[,-1]))
    
    
    # hospital admissions for covid (with estimation/correction)
    pred_hosp1_for_covid       <- pred_all$total_new_hosp_icu_for_covid
    hosp_age_adforcovid[i,,]   <- as.matrix(round(pred_hosp1_for_covid))       
    scen_hosp_adforcovid[,i]   <- as.matrix(round(rowSums(pred_all$total_new_hosp_icu_for_covid[,-1])))

    # infections
    new_cases           <- pred_all$total_new_infections
    scen_cases[i,,]     <- as.matrix(new_cases)      
    scen_new_infect[,i] <- as.matrix(rowSums(new_cases[,-1])) # new_cases
    
    # reinfections
    new_reinfect          <- pred_all$total_new_reinfections
    scen_new_reinfect[,i] <- as.matrix(rowSums(new_reinfect[,-1])) # new reinfections
    
    # numver of vaccines
    new_numvac           <- pred_all$total_new_numvac
    scen_numvac[i,,]     <- as.matrix(new_numvac)      
    #scen_new_infect[,i] <- as.matrix(rowSums(new_cases[,-1])) # new_cases
    
    # mild cases
    new_mild_cases        <- pred_all$total_new_mild_infections
    scen_mild_cases[i,,]  <- as.matrix(new_mild_cases)  
    scen_new_sympt[,i]    <- as.matrix(rowSums(new_mild_cases[,-1])) # new_cases
    
    # hospital and icu load 
    scen_hosp_load[,i]   <- as.matrix(rowSums(pred_all$hosp_icu_load[,-1]))
    scen_ICU_load[,i]    <- as.matrix(rowSums(pred_all$icu_load[,-1]))
    scen_ICU_adm[,i]     <- as.matrix(rowSums(pred_all$total_new_icu[,-1]))

    # prevalence and incidence (DALY)
    prev_I_asymp_sims[i,,]   <- as.matrix(pred_all$prev_I_asym[,-1])
    prev_I_mild_sims[i,,]    <- as.matrix(pred_all$prev_I_mild[,-1])
    prev_I_severe_sims[i,,]  <- as.matrix(pred_all$prev_I_sev[,-1])
    prev_I_hosp_sims[i,,]    <- as.matrix(pred_all$hosp_icu_load[,-1])
    prev_I_icu_sims[i,,]     <- as.matrix(pred_all$icu_load[,-1])
    inc_I_hosp_sims[i,,]     <- as.matrix(pred_all$total_new_hosp_icu[,-1])
    inc_I_icu_sims[i,,]      <- as.matrix(pred_all$total_new_icu[,-1])
    
    inc_I_hosp_for_covid_sims[i,,]     <- as.matrix(pred_all$total_new_hosp_icu_for_covid[,-1])
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
    
    # make sure that "VOC_ba4ba5_start" is part of run_param
    if(!"VOC_ba4ba5_start" %in% names(run_param)){
      run_param['VOC_ba4ba5_start'] <- num_days_sim
    }
  
    # make sure that "VOC_bq1ba275xb_start" is part of run_param
    if(!"VOC_bq1ba275xbb_start" %in% names(run_param)){
        run_param['VOC_bq1ba275xbb_start'] <- sim_date2day('2022-08-15')
      }

    # make sure that "VOC_xbb15_start" is part of run_param
    if(!"VOC_xbb15_start" %in% names(run_param)){
      run_param['VOC_xbb15_start'] <- sim_date2day('2023-01-01')
    }

    # VOC
    d_alpha   <- min(run_param['VOC_alpha_start'],num_days_sim):min(run_param['VOC_omicron_start'],num_days_sim)
    d_delta   <- min(run_param['VOC_delta_start'],num_days_sim):min(run_param['VOC_ba4ba5_start'],num_days_sim)
    d_omicron <- min(run_param['VOC_omicron_start'],num_days_sim):min(run_param['VOC_bq1ba275xbb_start'],num_days_sim)
    d_ba4ba5  <- min(run_param['VOC_ba4ba5_start'],num_days_sim):min(run_param['VOC_xbb15_start'],num_days_sim)
    d_bq1ba275xbb  <- min(run_param['VOC_bq1ba275xbb_start'],num_days_sim):num_days_sim
    d_xbb15  <- min(run_param['VOC_xbb15_start']+1,num_days_sim):num_days_sim
    scen_VOC_mild_alpha[d_alpha,i]     <- as.matrix(pred_all$p_VOC)[d_alpha]
    scen_VOC_mild_delta[d_delta,i]     <- as.matrix(1-pred_all$p_VOC)[d_delta]
    scen_VOC_mild_omicron[d_omicron,i] <- as.matrix(pred_all$p_VOC)[d_omicron]
    scen_VOC_mild_ba4ba5[d_ba4ba5,i]   <- as.matrix(1-pred_all$p_VOC)[d_ba4ba5]
    scen_VOC_mild_bq1ba275xbb[d_bq1ba275xbb,i]   <-  as.matrix(pred_all$p_VOC)[d_bq1ba275xbb]
    scen_VOC_mild_xbb15[d_xbb15,i]   <- as.matrix(1-pred_all$p_VOC)[d_xbb15]
    
    # cumulative infections
    scen_cases_cum     <- apply(as.matrix(new_cases),2,cumsum)
    scen_cases_cum[,1] <- as.matrix(new_cases)[,1]
    cumul_cases[i,,]   <- scen_cases_cum
    
    # cumulative deaths
    scen_deaths_cum     <- apply(as.matrix(pred_all$total_new_deaths),2,cumsum)
    scen_deaths_cum[,1] <- as.matrix(pred_all$total_new_deaths)[,1]
    cumul_deaths[i,,]   <- scen_deaths_cum
  
    
    # susceptible
    mean_susceptible_full[i,,]     <- as.matrix(pred_all$mean_susceptible_full)  
    mean_susceptible_initpluswaning[i,,]     <- as.matrix(pred_all$mean_susceptible_initpluswaning)  
    mean_susceptible_vac_d1[i,,]   <- as.matrix(pred_all$mean_susceptible_vac_d1)  
    mean_susceptible_vac_d2[i,,]   <- as.matrix(pred_all$mean_susceptible_vac_d2)  
    mean_vaccinated_age[i,,]       <- as.matrix(pred_all$mean_vaccinated_age)
  
    # non-vaccinated hosp admissions
    nvac_hosp_adm[,i]   <- as.matrix(rowSums(pred_all$total_new_hosp_icu_nvac[,-1]))
    nvac_age_hosp_adm[i,,]   <- as.matrix(pred_all$total_new_hosp_icu_nvac)
    new_hosp_icu_vac_d2[,i] <- as.matrix(rowSums(pred_all$total_new_hosp_icu_vac_d2[,-1]))
    new_hosp_icu_vac_age_d2[i,,] <- as.matrix(pred_all$total_new_hosp_icu_vac_d2)
    new_hosp_icu_vac_booster[,i] <- as.matrix(rowSums(pred_all$total_new_hosp_icu_vac_booster[,-1]))
    
    # vaccinated
    # note !flag_rna and !flag_d1 are not working anymore with potential booster doses
    flag_rna     <- grepl('rna',names(pred_all$mean_dose_vaccinated))
    flag_adeno   <- grepl('adeno',names(pred_all$mean_dose_vaccinated))  
    flag_d1      <- grepl('d1',names(pred_all$mean_dose_vaccinated))
    flag_d2      <- grepl('d2',names(pred_all$mean_dose_vaccinated))
    flag_booster <- grepl('booster',names(pred_all$mean_dose_vaccinated)) & ! grepl('waning',names(pred_all$mean_dose_vaccinated))
    #names(pred_all$mean_dose_vaccinated)[flag_rna & flag_d1]
    mean_rna1_vaccinated[i,,]    <- as.matrix(pred_all$mean_dose_vaccinated[,c(1,which(flag_rna & flag_d1))])
    mean_rna2_vaccinated[i,,]    <- as.matrix(pred_all$mean_dose_vaccinated[,c(1,which(flag_rna & flag_d2))])
    mean_adeno1_vaccinated[i,,]  <- as.matrix(pred_all$mean_dose_vaccinated[,c(1,which(flag_adeno & flag_d1))])
    mean_adeno2_vaccinated[i,,]  <- as.matrix(pred_all$mean_dose_vaccinated[,c(1,which(flag_adeno & flag_d2))])
    mean_waning_vaccinated[i,,]  <- as.matrix(pred_all$mean_dose_vaccinated[,c(1,which(flag_booster))])
    social_coef[,i] <- as.matrix(pred_all$mean_bar_coef[,-1])
  }
  # include contact scenario info into output tag - readded for hub
  cnt_adjust_tagvac <- paste0(cnt_adjust_tag,'_','vac')
  cnt_adjust_tag <- paste0(cnt_adjust_tag,'_',scen_tag)
  
  if(bool_vac == FALSE){
  scen_hosp_adm_vac <- readRDS(paste0(output_vac_dir,'hosp_adm_incr',gsub('\\.','p',cnt_adjust_tagvac),'_',vacc_ref_tag.,'.rds'))
  delta_hosp_new <- scen_hosp_adm-scen_hosp_adm_vac
  saveRDS(delta_hosp_new,file=paste0(output_dir,'delta_hosp_new_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  saveRDS(apply(delta_hosp_new,2,cumsum),file=paste0(output_dir,'delta_hosp_new_cumul',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
 
  scen_hosp_load_vac <- readRDS(paste0(output_vac_dir,'hosp_load_incr',gsub('\\.','p',cnt_adjust_tagvac),'_',vacc_ref_tag.,'.rds'))
  delta_hosp_load <- scen_hosp_load-scen_hosp_load_vac
  saveRDS(delta_hosp_load,file=paste0(output_dir,'delta_hosp_load',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  
  scen_new_infect_vac <- readRDS(paste0(output_vac_dir,'scen_new_infect',gsub('\\.','p',cnt_adjust_tagvac),'_',vacc_ref_tag.,'.rds'))
  delta_new_infect <- scen_new_infect - scen_new_infect_vac
  saveRDS(delta_new_infect,file=paste0(output_dir,'delta_new_infect',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  
  cumul_cases_vac <- readRDS(paste0(output_vac_dir,'cumul_cases_incr',gsub('\\.','p',cnt_adjust_tagvac),'_',vacc_ref_tag.,'.rds'))
  cum_aggr <- t(apply(cumul_cases[,,2:11], c(1,2), sum))
  delta_cumul_cases<- cum_aggr - cumul_cases_vac
  saveRDS(delta_cumul_cases,file=paste0(output_dir,'delta_cumul_cases',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  
  cumul_deaths_vac <- readRDS(paste0(output_vac_dir,'cumul_deaths_incr',gsub('\\.','p',cnt_adjust_tagvac),'_',vacc_ref_tag.,'.rds'))
  cum_aggr <- t(apply(cumul_deaths[,,2:11], c(1,2), sum))
  delta_cumul_deaths<- cum_aggr - cumul_deaths_vac
  saveRDS(delta_cumul_deaths,file=paste0(output_dir,'delta_cumul_deaths',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  
  scen_Rt_infection_vac <- readRDS(paste0(output_vac_dir,'scen_rt_infection_incr',gsub('\\.','p',cnt_adjust_tagvac),'_',vacc_ref_tag.,'.rds'))
  delta_Rt_infection <- scen_Rt_infection - scen_Rt_infection_vac
  saveRDS(delta_Rt_infection,file=paste0(output_dir,'delta_Rt_infection',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  
  scen_Rt_hosp_vac <- readRDS(paste0(output_vac_dir,'scen_Rt_hosp_incr',gsub('\\.','p',cnt_adjust_tagvac),'_',vacc_ref_tag.,'.rds'))
  delta_Rt_hosp <- scen_Rt_hosp - scen_Rt_hosp_vac
  saveRDS(delta_Rt_hosp,file=paste0(output_dir,'delta_Rt_hosp',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  
  scen_icu_load_vac <- readRDS(paste0(output_vac_dir,'icu_load_incr',gsub('\\.','p',cnt_adjust_tagvac),'_',vacc_ref_tag.,'.rds'))
  delta_icu_load <- scen_ICU_load - scen_icu_load_vac
  saveRDS(delta_icu_load,file=paste0(output_dir,'delta_icu_load',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  
  scen_mortality_vac <- readRDS(paste0(output_vac_dir,'scen_mortality',gsub('\\.','p',cnt_adjust_tagvac),'_',vacc_ref_tag.,'.rds'))
  delta_mortality <- scen_mortality - scen_mortality_vac
  saveRDS(delta_mortality,file=paste0(output_dir,'delta_mortality',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  

  }
 # saveRDS(pred_all$norm_asy,file=paste0(output_dir,'norm_asy',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
#  saveRDS(pred_all$norm_sy,file=paste0(output_dir,'norm_sy',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  saveRDS(hosp_age_adm,file=paste0(output_dir,'scenario1_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  saveRDS(hosp_age_adm,file=paste0(output_dir,'hosp_age_adm_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  saveRDS(hosp_age_adforcovid,file=paste0(output_dir,'hosp_age_adforcovid_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  saveRDS(scen_cases,file=paste0(output_dir,'scen_cases_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  saveRDS(scen_mild_cases,file=paste0(output_dir,'scen_mild_cases_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  saveRDS(t(apply(cumul_cases[,,2:11], c(1,2), sum)),file=paste0(output_dir,'cumul_cases_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  saveRDS(t(apply(cumul_deaths[,,2:11], c(1,2), sum)),file=paste0(output_dir,'cumul_deaths_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  saveRDS(mean_susceptible_full,file=paste0(output_dir,'mean_susceptible_full_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  saveRDS(mean_susceptible_initpluswaning,file=paste0(output_dir,'mean_susceptible_initpluswaning_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  #saveRDS(mean_susceptible_vac_d1,file=paste0(output_dir,'mean_susceptible_vac_d1_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  #saveRDS(mean_susceptible_vac_d2,file=paste0(output_dir,'mean_susceptible_vac_d2_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  #saveRDS(mean_vaccinated_age,file=paste0(output_dir,'mean_vaccinated_age_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  #saveRDS(scen_numvac,file=paste0(output_dir,'scen_numvac',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  
  #saveRDS(mean_rna1_vaccinated,file=paste0(output_dir,'mean_rna1_vaccinated_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  #saveRDS(mean_rna2_vaccinated,file=paste0(output_dir,'mean_rna2_vaccinated_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  #saveRDS(mean_adeno1_vaccinated,file=paste0(output_dir,'mean_adeno1_vaccinated_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  #saveRDS(mean_adeno2_vaccinated,file=paste0(output_dir,'mean_adeno2_vaccinated_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  #saveRDS(mean_waning_vaccinated,file=paste0(output_dir,'mean_waning_vaccinated_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  
  saveRDS(scen_new_infect,file=paste0(output_dir,'scen_new_infect',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  #saveRDS(scen_new_reinfect,file=paste0(output_dir,'new_reinfect',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  #saveRDS(scen_new_sympt,file=paste0(output_dir,'new_sympt',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  saveRDS(scen_hosp_load,file=paste0(output_dir,'hosp_load_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  saveRDS(scen_Rt_infection,file=paste0(output_dir,'scen_Rt_infection_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  saveRDS(scen_Rt_mild_sev,file=paste0(output_dir,'scen_Rt_mild_sev_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  saveRDS(scen_Rt_hosp,file=paste0(output_dir,'scen_Rt_hosp_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  saveRDS(scen_hosp_adm,file=paste0(output_dir,'hosp_adm_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  saveRDS(apply(scen_hosp_adm,2,cumsum),file=paste0(output_dir,'hosp_adm_cum',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  #saveRDS(scen_hosp_adforcovid,file=paste0(output_dir,'hosp_adforcovid_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  saveRDS(scen_ICU_load,file=paste0(output_dir,'icu_load_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  #saveRDS(scen_ICU_adm,file=paste0(output_dir,'icu_adm_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  #saveRDS(scen_hosp_exit,file=paste0(output_dir,'hosp_exit_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  saveRDS(scen_mortality,file=paste0(output_dir,'scen_mortality',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  saveRDS(mortality_age,file=paste0(output_dir,'mortality_age_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  #saveRDS(scen_VOC_mild_alpha,file=paste0(output_dir,'scen_VOC_mild_alpha_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  #saveRDS(scen_VOC_mild_delta,file=paste0(output_dir,'scen_VOC_mild_delta_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  #saveRDS(scen_VOC_mild_omicron,file=paste0(output_dir,'scen_VOC_mild_omicron_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  #saveRDS(scen_VOC_mild_ba4ba5,file=paste0(output_dir,'scen_VOC_mild_ba4ba5_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  #saveRDS(scen_VOC_mild_bq1ba275xbb,file=paste0(output_dir,'scen_VOC_mild_bq1ba275xbb_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  #saveRDS(scen_VOC_mild_xbb15,file=paste0(output_dir,'scen_VOC_mild_xbb15_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  saveRDS(inc_I_icu_sims,file=paste0(output_dir,'inc_I_icu_sims_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  
  # burden of disease figures
  saveRDS(prev_I_asymp_sims,file=paste0(output_dir,'prev_I_asymp_sims_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  saveRDS(prev_I_mild_sims,file=paste0(output_dir,'prev_I_mild_sims_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  saveRDS(prev_I_severe_sims,file=paste0(output_dir,'prev_I_severe_sims_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  saveRDS(prev_I_hosp_sims,file=paste0(output_dir,'prev_I_hosp_sims_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  saveRDS(prev_I_icu_sims,file=paste0(output_dir,'prev_I_icu_sims_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  saveRDS(inc_I_hosp_sims,file=paste0(output_dir,'inc_I_hosp_sims_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  saveRDS(inc_I_icu_sims,file=paste0(output_dir,'inc_I_icu_sims_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  saveRDS(inc_I_hosp_for_covid_sims,file=paste0(output_dir,'inc_I_hosp_for_covid_sims_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  
  #vac_hosp_adm    <- scen_hosp_adm - nvac_hosp_adm
  #saveRDS(nvac_hosp_adm,file=paste0(output_dir,'nvac_hosp_adm_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  #saveRDS(vac_hosp_adm,file=paste0(output_dir,'vac_hosp_adm_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  
  #vac_age_hosp_adm    <- hosp_age_adm - nvac_age_hosp_adm
  #vac_age_hosp_adm[,,1] <- hosp_age_adm[,,1]
  #saveRDS(vac_age_hosp_adm,file=paste0(output_dir,'vac_age_hosp_adm_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  #saveRDS(nvac_age_hosp_adm,file=paste0(output_dir,'nvac_age_hosp_adm_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  
  #saveRDS(new_hosp_icu_vac_d2,file=paste0(output_dir,'new_hosp_adm_vac_d2_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  #saveRDS(new_hosp_icu_vac_age_d2,file=paste0(output_dir,'new_hosp_adm_vac_age_d2_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  #saveRDS(new_hosp_icu_vac_booster,file=paste0(output_dir,'new_hosp_adm_vac_booster_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  saveRDS(social_coef,file=paste0(output_dir,'social_coef',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  
  cum_death_bar <- rowSums(cumul_deaths[, dim(cumul_deaths)[2], 2:11 ])
  if(bool_vac == FALSE){
  delta_death <- delta_cumul_deaths[dim(delta_cumul_deaths)[1],]
  delta_cases <- colSums(delta_new_infect)
  } else {
    delta_death <- matrix(0,1,5)
    delta_cases <- matrix(0,1,5)  
  }
  
  num_days_critic_hosp <- apply(scen_hosp_load, 2, function(col) sum(col > 8000))
  num_days_critic_ICU <- apply(scen_ICU_load, 2, function(col) sum(col > 1500))
  num_days_overload_hosp <- apply(scen_hosp_load, 2, function(col) sum(col > 10000))
  num_days_overload_ICU <- apply(scen_ICU_load, 2, function(col) sum(col > 2000))
  if (bool_bar == TRUE)
  {
  ndays_or <- apply(social_coef, 2, function(col) sum(col < 1 & col >= or_prop))
  ndays_red <- apply(social_coef, 2, function(col) sum(col < or_prop & col >= red_prop))
  restrict_score <- num_days_sim-colSums(social_coef)
  }
  else
  {
  ndays_or <- matrix(0,1,5)
  ndays_red <- matrix(0,1,5)
  restrict_score <- matrix(0,1,5)
  restrict_score <- matrix(0,1,5)
  or_prop <- 1.0
  red_prop <- 1.0
  new_hosp_or <- 0.0
  new_hosp_red <- 0.0
  icu_or <- 0.0
  icu_red <- 0.0
  }
  bar_summary <- list(
    cum_death_bar = cum_death_bar,
    delta_death = delta_death,
    delta_cases = delta_cases,
    num_days_overload_hosp = num_days_overload_hosp,
    num_days_overload_ICU = num_days_overload_ICU,
    num_days_critic_hosp = num_days_critic_hosp,
    num_days_critic_ICU = num_days_critic_ICU,
    ndays_or = ndays_or,
    ndays_red = ndays_red,
    restrict_score = restrict_score
  )
  
  
  
  # Write the median values to a text file
  write(paste("Parameters,",paste(c(or_prop,red_prop,new_hosp_or, new_hosp_red, icu_or, icu_red), collapse = ",")), file = paste0(output_dir,'bar_summary',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.csv'))
  write(paste("Indicator,",paste(names(bar_summary), collapse = ",")), file = paste0(output_dir,'bar_summary',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.csv'), append = TRUE)
  write(paste("Median,",paste(sapply(bar_summary, median), collapse = ",")), file = paste0(output_dir,'bar_summary',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.csv'), append = TRUE)
  
  # Write the confidence interval bounds to a text file
  write(paste("Lower CI,", paste(sapply(bar_summary, function(x) quantile(x, 0.025)), collapse = ",")), file =paste0(output_dir,'bar_summary',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.csv'), append = TRUE)
  write(paste("Upper CI,", paste(sapply(bar_summary, function(x) quantile(x, 0.975)), collapse = ",")), file =paste0(output_dir,'bar_summary',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.csv'), append = TRUE)
  
  # Write the summary statistics to a text file
  #write.table(bar_summary_df, file = paste0(output_dir,'bar_summary',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.txt'), sep = "\t", row.names = FALSE)
  
  #write.table(paste(c(or_prop, red_prop, new_hosp_or, new_hosp_red, icu_or, icu_red, bar_summary_med), collapse = " "), file = paste0(output_dir,'bar_summary',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.txt'), row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  saveRDS(bar_summary,file=paste0(output_dir,'bar_summary',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds'))
  
  
  
    ## plot results ----
  pdf(file=paste0(output_dir,'scenario_cases_mean_incr',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.pdf'),9,6)
  par(mar=c(4.5,5,1,1))
  scm_scenario_cp <- sim_day2date(cnt_adjust_day[cnt_adjust_value>0]) + rep(0:6,each=sum(cnt_adjust_value>0))
  #temp 1 month - to change
  scm_scenario_cp <- sim_day2date(cnt_adjust_day[cnt_adjust_value>0]) + rep(0:30,each=sum(cnt_adjust_value>0))
  x_axis_scm_day = c(0,num_days_sim) #c(0,num_days_sim)
  x_axis_vaccine = c(0,num_days_sim)
  
  output_files <- dir(output_dir,full.names = T,pattern = vacc_schedule_tag)
  multi_plot_incidence_time(output_files       = output_files,
                            bool_polygon_multi = T,
                            col_lib            = col_lib,
                            be_ref_data        = be_ref_data,
                            default_y_axis     = default_y_axis,
                            x_axis_scm_day      = x_axis_scm_day,
                            scm_scenario_cp     = scm_scenario_cp)
  dev.off()
  
  # # optional: plot age-specific results
  # pdf(file=paste0(output_dir,'scenario_cases_mean_age',gsub('\\.','p',cnt_adjust_tag),'_',vacc_schedule_tag,'.pdf'))
  # multi_plot_incidence_age_time(output_files   = output_files,
  #                               col_lib        = col_lib,
  #                               x_axis_scm_day  = x_axis_scm_day,
  #                               x_axis_vaccine = x_axis_vaccine)
  # 
  # dev.off()
  
} # end pallel foreach
#smd_stop_cluster()
}
