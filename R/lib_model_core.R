########################################################################### #
# This file is part of the Stochastic Compartmental Model for SARS-COV-2 
# transmission in Belgium, conceived by members of SIMID group during the 
# COVID19 pandemic.
#
# This file contains the core function of the stochastic model, in combination
# with the function to calculate the loglikelihood and log prior distributions. 
# This file also contains a function to plot rudimental model output.
#
# Copyright 2022, SIMID, University of Antwerp & Hasselt University                                        
########################################################################### #

if(0==1){ # debug parameters
  parms <- run_param 
  parms_names <- names(parms)
  ndays_sim <- 740;
  vaccine_uptake <- vacc_schedule # vacc_schedule_uptake
  CoMix_matrices = CoMix_mat; CoMix = "full"; plots = "FALSE"
  method = "mean"; # method = "stochastic"; 
  cnt_change_pt = c(0,0,0); cnt_change = c(0,0,0)
}
run_model <- function(parms, CoMix_matrices, method = "mean",
                      parms_names = NA, ndays_sim, 
                      cnt_change_pt = 0, cnt_change = 0,
                      bool_full_output=TRUE, 
                      vaccine_uptake = NA,
                      V_mat, V_mat_dose2,
                      V_mat_booster,V_mat_2ndbooster,
                      V_mat_extrabooster){  
  
  # (re-)store parameter names if they are not given, or the parameter matrix don't contain them in e.g. the MCMC deamon
  if(any(is.na(parms_names))){
    parms_names <- names(parms)
  } else {
    names(parms) <- parms_names
  }
  
  ## make sure that the model parameters are in line with latest requirements
  # as such, update colnames and add default values if needed
  ##-------------------------------------------------------- -
  parms <- update2latest_model_parameter_config(parms)
  parms_names <- names(parms)
  
  ## General model parameters -----
  ##-------------------------------------------------------- -
  gamma  = exp(parms['log_gamma']);   # average length of the latency period - 2 days
  theta  = exp(parms['log_theta']);   # infectious period at pre-symptomatic stage - 3.5 days
  delta1 = exp(parms['log_delta1']);  # infectious period at asymptomatic stage - 3.5 days
  delta2 = exp(parms['log_delta2']);  # infectious period at (mild) symptomatic stage - 3.5 days
  
  phi0 = expit(parms[grepl('log_phi0',parms_names)]);  # proportion of symptomatically infected with mild symptoms
  p = p_vec;                                           # proportion of asymptomatic cases
  
  omega = exp(parms[grepl('log_omega',parms_names)]);          # waiting time between symptom onset and hospitalization
  phi1  = expit(parms[grepl('log_phi1_age',parms_names)] * parms['log_phi1_add']);  # proportion of severly infected with regular hospitalization
  q = exp(parms['log_q']);                # proportionality constant
  
  mu_param = parms[grepl('mu_sev',parms_names)]
  mu_sev = c(0,expit(mu_param[c(1,1:8)]));   # age-specific hospital fatality ratio - probability of dying
  delta3 = exp(parms['log_delta3']);         # recovery period in severely infected/length of stay - 7 days
  
  # update mortality from September 2020
  mortality_changepoint = parms[grepl('mortality_changepoint',parms_names)]
  mu2_param = parms[grepl('mu2_sev',parms_names)]
  mu2_sev   = c(0,expit(mu2_param[c(1,1:8)])); # age-specific hospital fatality ratio - probability of dying
  delta3_wave2 = exp(parms['log_delta3_wave2']);    # recovery rate from hospital (from sept 2021)
  
  ## Model parameters related to the social contact behaviour (and CoMix matrices) -----
  ##---------------------------------------------------------------------------------- -
  nCoMix <- sum(grepl('log_comix_coef',parms_names)) / 10  # number of comix related column names, divided by 10 age groups
  nOther <- sum(grepl('log_add_coef',parms_names)) / 10 # number of additional matrix related column names, divided by 10 age groups
  
  CoMix_coef_mat = matrix(exp(parms[grepl('log_comix_coef',parms_names)]), nrow = nCoMix, ncol = 10, byrow = T);
  Add_coef_mat   = matrix(exp(parms[grepl('log_add_coef',parms_names)]), nrow = nOther, ncol = 10, byrow = T);
  
  ## Model initialization: population & region  ---- 
  ##------------------------- -
  n0 = exp(parms['log_n0'])
  imported_cases = round(rel_freq_cc*n0*(1/(1-p_vec)),0);     
  
  sel_region  <- get_region(parms['region_id'])
  cohort.size <- get_regional_pop(region=sel_region)
  
  ## Transmissibility parameters  ---- 
  ##-------------------------------- -
  f = 0.51;   # relative infectiousness of asymptomatic vs. symptomatic cases
  
  # combine parameters
  q_asy = f * q
  q_sy  = q 
  
  ## Stochastic model parameters  ---- 
  ##-------------------------------- -
  h = parms['h'];             # resolution of the binomial chains
  
  ## Data augmentation  ---- 
  ##---------------------- -
  times = seq(0, ndays_sim-h, h)
  times_day <- floor(times)
  
  ## Initialization VOC (generic)  ---- 
  ##----------------------- -
  #warning, from 2022-08-01 new Voc can have immune evasion
  start_date_voc_immune_evasion <- '2021-11-01' #2022-08-01'
  
  VOC_start  = parms[grepl('VOC_.*_start',parms_names)]
  
  #fix extrabooster ve as delta-like
 # parms['ve_VOC_omicron_infection_extrabooster']         <- parms['ve_VOC_delta_infection_booster']
  #parms['ve_VOC_ba4ba5_infection_extrabooster']         <- parms['ve_VOC_delta_infection_booster']
  
  parms_names <- names(parms)
  
  # add one new voc---
  start_date_new_hypothetic_voc <- '2123-05-01'
  newdate <- start_date_new_hypothetic_voc
  numnewvoc <- 0
     VOC_start <- c(VOC_start, sim_date2day(newdate) )
     names(VOC_start) <- c(names(head(VOC_start,-1)),paste0("VOC_newvoc",(numnewvoc +1),"_start"))
     newdata <- as.data.frame(t(c(paste0("VOC_newvoc",(numnewvoc +1)),newdate,newdate,(2-numnewvoc%%2),sim_date2day(newdate),sim_date2day(newdate))))
      names(newdata) <- names(global_lib_voc)
      global_lib_voc <- rbind(global_lib_voc, newdata)
      
      sel_param <- parms[grepl('VOC_xbb15',parms_names)]
      parms[gsub('xbb15',paste0("newvoc",(numnewvoc +1)),names(sel_param))] <- sel_param
      parms_names <- names(parms)
      
      # reduction in immunity against infection (except vs itself)
        # newvereduction <- 0.70
         #parms["ve_VOC_newvoc1_vereduction"] <- parms["ve_VOC_xbb15_vereduction"]*newvereduction  # to change with exp/log
         #parms_names <- names(parms)

        parms[paste0("log_VOC_newvoc1_init")]  <- 0
        parms[paste0("VOC_newvoc1_init")]  <- 0
        parms[paste0("VOC_newvoc1_start")] <- 1e10 # default: no introduction
     
     # adding repetitive new voc---
  # numnewvoc <- 0
  # 
  # while(sim_date2day(as.Date(start_date_voc_immune_evasion)%m+% months(90*numnewvoc)) %in% times_day){
  #   newdate <- as.Date(start_date_voc_immune_evasion)%m+% months(90*numnewvoc)
  #   VOC_start <- c(VOC_start, sim_date2day(newdate) )
  #   names(VOC_start) <- c(names(head(VOC_start,-1)),paste0("VOC_newvoc",(numnewvoc +1),"_start"))
  #   newdata <- as.data.frame(t(c(paste0("VOC_newvoc",(numnewvoc +1)),newdate,newdate,(2-numnewvoc%%2),sim_date2day(newdate),sim_date2day(newdate))))
  #   names(newdata) <- names(global_lib_voc)
  #   global_lib_voc <- rbind(global_lib_voc, newdata)
  #   
  #   sel_param <- parms[grepl('VOC_ba4ba5',parms_names)]
  #   parms[gsub('ba4ba5',paste0("newvoc",(numnewvoc +1)),names(sel_param))] <- sel_param
  #   parms_names <- names(parms)
  #   
  #   # reduction in immunity against infection (except vs itself)
  #   vereduction <- 0.70
  #   if(numnewvoc==0){ # bq1
  #     vereduction <- exp(parms["log_VOC_bq1_vereduction"])
  #   }
  #   sel_param <- parms[grepl(paste0("ve_VOC_newvoc",(numnewvoc +1),"_infection") ,parms_names) & !grepl('_reinf',parms_names)]
  #   parms[grepl(paste0("ve_VOC_newvoc",(numnewvoc +1),"_infection"),parms_names) & !grepl('_reinf',parms_names)] <- sel_param*(vereduction)
  #   parms_names <- names(parms)
  #   sel_param <- parms[grepl(paste0("ve_VOC_newvoc",(numnewvoc +1),"_infection") ,parms_names) & grepl('_reinf',parms_names)]
  #   parms[paste0(names(sel_param),"_oldvoc")]<- sel_param*(vereduction)
  #   parms_names <- names(parms)
  #   
  #   # reduction in immunity against severe (delta as base)
  #   vesevreduction <- 1
  #   if(numnewvoc==0){ # bq1
  #     vesevreduction <- exp(parms["log_VOC_bq1_vesevreduction"])
  #   }
  #   sel_param <- parms[grepl(paste0("ve_VOC_","delta","_incr_severe"),parms_names)]
  #   parms[grepl(paste0("ve_VOC_newvoc",(numnewvoc +1),"_incr_severe"),parms_names)] <- sel_param*(vesevreduction) #^(numnewvoc +1)
  #   parms_names <- names(parms)
  #   #print( parms[grepl(paste0("ve_VOC_newvoc",(numnewvoc +1),"_incr_severe"),parms_names)])
  #   #print( parms[grepl(paste0("ve_VOC_delta_incr_severe_extrabooster"),parms_names)])
  #   #print( parms[grepl(paste0("ve_VOC_ba4ba5_incr_severe_extrabooster"),parms_names)])
  #   
  #   parms[paste0("VOC_newvoc",(numnewvoc +1),'_init')]  <- 0
  #   parms[paste0("VOC_newvoc",(numnewvoc +1),'_start')] <- 1e10 # default: no introduction
  #   
  #   numnewvoc <- numnewvoc +1
  # }  
  
  #newvoc 1 = bq1
 # parms["log_VOC_newvoc1_init"]  <- parms["log_VOC_bq1_init"]
  #parms["VOC_newvoc1_start"] <-  parms["VOC_bq1_start"]
  
  #---
  #print(parms[grepl(paste0("ba4ba5_infection"),parms_names)]) # & !grepl('_reinf',parms_names)
  #print(parms[grepl(paste0("newvoc1_infection"),parms_names)]) # & !grepl('_reinf',parms_names)
  #print(parms[grepl(paste0("newvoc2_infection"),parms_names)]) # & !grepl('_reinf',parms_names)
  

  
  VOC_name   = gsub('_start','',names(VOC_start))
  
  #browser()
  
  if(!all(VOC_name %in% global_lib_voc$name)){
    stop("ERROR: UNKNOWN VOC NAME IN PARAMETER FILE")
  }
  
  if(any(duplicated(VOC_start))){
    VOC_start[duplicated(VOC_start)] <- VOC_start[duplicated(VOC_start)]+1
    warning("WARNING: duplicated VOC_start value, increased '",names(VOC_start)[duplicated(VOC_start)],"' by one day")
  }
  
  # keep track of current VOC-related adjustment of the hospital admission probability
  log_VOC_XXX_hosp <- NA
  
  # VOC-specific transmission factor ----
  transm_fctr  = 1                    # start with factor 1 (= 2020 wild type)
  transm_fctr_VOC    = 1                    # start with factor 1 (= 2020 wild type)
  
  # updated ICU/hosp proportion from fall 2021 and spring 2022
  phi1_fall2021      = expit(parms[grepl('log_phi1_age',parms_names)] * parms['log_fall2021_phi1_add']); # extra: proportion of severly infected with regular hospitalization
  icu_fall2021_start = parms['icu_fall2021_start']
  phi1_spring2022      = expit(parms[grepl('log_phi1_age',parms_names)] * parms['log_spring2022_phi1_add']); # extra: proportion of severly infected with regular hospitalization
  icu_spring2022_start = parms['icu_spring2022_start']
  
  ## Regional simulation?  ---- 
  ##----------------------- -
  pop_be     <- get_regional_pop(region = 'belgium')
  pop_model  <- get_regional_pop(region = sel_region)
  pop_factor   <- sum(pop_model)/sum(pop_be)
  if(pop_factor <1){
    imported_cases  <- round(imported_cases * pop_factor)
  }  
  
  ## Vaccination parameters  ---- 
  ##--------------------------- -
  
  # start with equal VE against infection for wild type and alpha variant
  ve_infection <- get_ve_infection(parms,ve_tag='ve_VOC_alpha')
  ve_voc_infection <- ve_infection
  
  # protection against infectiousness/transmission
  ve_transmission <- 0*parms['ve_transmission']  #temp put to zero
  
  # vaccine related waning immunity
  ve_waning_immunity_rate <- parms['ve_waning_immunity_rate'] * h
  ve_waning_booster_rate  <- parms['ve_waning_booster_rate'] * h
  ve_waning_extrabooster_rate  <- parms['ve_waning_booster_rate'] * h
  
  # infection related waning immunity
  ve_waning_infection_rate <- parms['ve_waning_infection_rate'] * h
  ve_waning_infection_booster_rate <- parms['ve_waning_infection_booster_rate'] * h
  
  # check vaccine uptake: dose 1 and 2
  if(any(!is.na(vaccine_uptake))){
    V_mat         <- get_uptake_matrix(vaccine_uptake = vaccine_uptake, dose_tag = '_A_', parms = parms )
    V_mat_dose2   <- get_uptake_matrix(vaccine_uptake = vaccine_uptake, dose_tag = '_B_', parms = parms )
    V_mat_booster <- get_uptake_matrix(vaccine_uptake = vaccine_uptake, dose_tag = '_E_', parms = parms )
    V_mat_2ndbooster <- get_uptake_matrix(vaccine_uptake = vaccine_uptake, dose_tag = '_F_', parms = parms )
    V_mat_extrabooster <- get_uptake_matrix(vaccine_uptake = vaccine_uptake, dose_tag = '_X_', parms = parms )
  } else if(missing('V_mat') || missing('V_mat_dose2')){
    V_mat <- matrix(0,nrow=1,ncol=20)
    V_mat_dose2 <- matrix(0,nrow=1,ncol=20)
  }
  # check booster vaccine uptake matrix
  if(missing('V_mat_booster')){
    V_mat_booster <- V_mat*0
  }
  # check 2ndbooster vaccine uptake matrix
  if(missing('V_mat_2ndbooster')){
    V_mat_2ndbooster <- V_mat*0
  }
  # check extra booster vaccine uptake matrix
  if(missing('V_mat_extrabooster')){
    V_mat_extrabooster <- V_mat*0
  }
  
  # specify time steps with vaccine uptake
  bool_n_uptake                             <- rowSums(V_mat) > 0
  bool_n_uptake[rowSums(V_mat_dose2) > 0]   <- TRUE
  bool_n_uptake[rowSums(V_mat_booster) > 0] <- TRUE
  bool_n_uptake[rowSums(V_mat_2ndbooster) > 0] <- TRUE
  bool_n_uptake[rowSums(V_mat_extrabooster) > 0] <- TRUE
  #temp, add vacc dates for buffer vaccination -> to check if this is still needed (only if more vaccinated people than available in waining)
  bool_n_uptake[(4*sim_date2day('2022-09-12')):length(bool_n_uptake)] <- TRUE
  if(length(bool_n_uptake)<length(times)){
    bool_n_uptake <- c(bool_n_uptake,rep(FALSE,length(times) - length(bool_n_uptake)))
  }
  

  ## Hospital probability ----
  # phi0 # age-specific probability NOT to be hospitalized
  get_phi0_vac <- function(phi0,ve_hosp){
    return(phi0 + (1-phi0)*ve_hosp)
  }
  
  ## Sensitivity function (serial serological data)  ---- 
  ##--------------------------------------------------- -
  sens_fun <- function(t){
    return(expit(-8.397 + 0.928*t))
  }    
  
  ## Transmission function (mean or stochastic)  ---- 
  ##--------------------------------------------------- -
  if (method == "mean"){   #floor added to to with stochastic one - temp removed
    get_new_transitions <- function(size_vec, new_prob_vec){
      return(size_vec*new_prob_vec)
    }
    bool_stochastic <- FALSE
  } else if (method == "stochastic"){  #floor added to avoid some NA  - temp removed
    get_new_transitions <- function(size_vec, new_prob_vec){
      if(any(size_vec%%1!=0)){browser()}
      output <- rbinom(length(size_vec), size = size_vec, prob = new_prob_vec)
      if(any(is.na(output))){browser()} 
      return(output)
    }
    bool_stochastic <- TRUE
  } else{
    warning(paste("UNKNOWN method:",method))
  }
  
  ## Initialization of the compartments  ---- 
  ##--------------------------------------- -
  
  ######################################################## #
  #  Description of the compartments (tmp_vec):            #
  #                                                        #
  #  susceptible - S: columns 1-10                         #
  #  exposed - E: columns 11-20                            #
  #  pre-symptomatic - I_presym: columns 21-30             #
  #  asymptomatic - I_asym: columns 31-40                  #
  #  mildly infected - I_mild: columns 41-50               #
  #  severely infected - I_sev: columns 51-60              #
  #  hospitalized - I_hosp: columns 61-70                  #
  #  icu admitted - I_icu: columns 71-80                   #
  #  deaths - D: columns 81-90                             #
  #  recovered - R: columns 91-100                         #
  #                                                        #
  #  Additional computations:                              #
  #  new hospitalized total: columns 101-110               #
  #  asymptotic and mildly infected: columns 111-120       #
  #  recovered from hospital only: columns 121-130 
  #
  #  exposed (VOC) - Evoc: columns 131-140                 #
  #  pre-symptomatic (VOC) - Ivoc_presym: columns 141-150  #
  #  asymptomatic (VOC) - Ivoc_asym: columns 151-160       #
  #  mildly infected (VOC) - Ivoc_mild: columns 161-170    #
  #  severely infected (VOC) - Ivoc_sev: columns 171-180   #
  #  hospitalized (VOC) - Ivoc_hosp: columns 181-190       #
  #  icu admitted (VOC) - Ivoc_icu: columns 191-200        #
  #  deaths (VOC) - Dvoc: columns: 201-210                 #
  #  recovered (VOC) - Rvoc: columns: 211-220              #
  #                                                        #
  #  Additional computations:                              #
  #  new hospitalized total (VOC): columns 221-230         #
  #  asymptotic and mildly infected (VOC): columns 231-240 #
  #  recovered from hospital only (VOC): columns 241-250   #
  #                                                        #
  #  Additional computations:                              #
  #  new mild infections                                   #
  #  new mild infections (VOC)                             #
  #  new recoveries after mild infection                   #
  #  new recoveries after mild infection (VOC)             #
  #  ICU admissions                                        #
  #  ICU admissions (VOC)                                  #
  ######################################################## #      
  
  icol_tmp_vec <- data.frame(S        = 1:10,
                             E        = 11:20,
                             I_presym = 21:30,
                             I_asym   = 31:40,
                             I_mild   = 41:50,               
                             I_sev    = 51:60,              
                             I_hosp   = 61:70,                  
                             I_icu    = 71:80,                   
                             D        = 81:90,  
                             R        = 91:100, 
                             
                             #  Additional computations          
                             new_hosp_total     = 101:110,
                             asym_mild_infected = 111:120,
                             hosp_recovered     = 121:130,
                             
                             # VOC
                             Evoc        = 131:140,  
                             Ivoc_presym = 141:150, 
                             Ivoc_asym   = 151:160,  
                             Ivoc_mild   = 161:170,  
                             Ivoc_sev    = 171:180,  
                             Ivoc_hosp   = 181:190, 
                             Ivoc_icu    = 191:200,  
                             Dvoc        = 201:210,    
                             Rvoc        = 211:220,   
                             
                             #  Additional computations (VOC)  
                             new_hosp_total_voc     = 221:230,
                             asym_mild_infected_voc = 231:240,
                             hosp_recovered_voc     = 241:250,
                             
                             # to register the mild infections
                             mild_infected     = 251:260,
                             mild_infected_voc = 261:270,
                             recov_mild_infection = 271:280,
                             recov_mild_infection_voc = 281:290,
                             
                             # to register the ICU admissions
                             new_icu = 291:300,
                             new_icu_voc = 301:310
  )
  
  # set flags to specify the ODE and log parameters  
  icol_flag_ode <- c(rep(1,100),rep(0,30),rep(1,90),rep(0,90))
  
  ## ATTACH ----
  names(icol_tmp_vec) <- paste0('c_',names(icol_tmp_vec))
  attach(icol_tmp_vec)
  
  ##################################################################################### #
  #  nfull: keep track of wild type and VOC compartments over time                      
  #  nfull_vac_rna: vaccinated with mRNA                                               
  #  nfull_vac_adeno: vaccinated with adeno-based vaccine   
  #  nfull_vac_waning: waning immunity of vaccine-induced immunity
  #  nfull_vac_booster: vaccinated with booster dose
  #  nfull_vac_booster_waning: waning immunity of booster dose
  #  nfull_vac_reinf: reinfection path (without vaccinations)
  #  nfull_vac_reinf_oldvoc: reinfection path (without vaccinations) from previous voc
  #  nfull_vac_reinfvac: reinfection path (for vaccinated)
  #  nfull_vac_reinfvac_oldvoc: reinfection path (for vaccinated) from previous voc
  #  nfull_vac_extrabooster1: vaccinated with extrabooster dose (high ve vs voc 1)
  #  nfull_vac_extrabooster2: vaccinated with extrabooster dose (high ve vs voc 2)
  ##################################################################################### #
  nfull <- as.list(data.frame(matrix(0, nrow = max(icol_tmp_vec), ncol = length(times))));
  nfull_vac_rna1         <- nfull;
  nfull_vac_rna2         <- nfull;
  nfull_vac_adeno1       <- nfull;
  nfull_vac_adeno2       <- nfull;
  nfull_vac_waning       <- nfull;
  nfull_vac_booster      <- nfull;
  nfull_vac_booster_waning <- nfull;
  nfull_vac_reinf        <- nfull;
  nfull_vac_reinf_oldvoc        <- nfull;
  nfull_vac_reinfvac     <- nfull;
  nfull_vac_reinfvac_oldvoc      <- nfull;
  nfull_vac_extrabooster1      <- nfull;
  nfull_vac_extrabooster2      <- nfull;
  
  # log re/infections
  nsubset_infection   <- as.list(data.frame(matrix(0, nrow = 1*length(c_R), ncol = length(times))))
  nsubset_reinfection <- as.list(data.frame(matrix(0, nrow = 1*length(c_R), ncol = length(times))))
  # count number of vaccine
  nsubset_numvac <- as.list(data.frame(matrix(0, nrow = 1*length(c_R), ncol = length(times))))
  
  ##################################################################################### #
  #  initialise: nfull                                                                  #
  ##################################################################################### #
  nfull[[1]][c(c_S,c_E)] <- c(cohort.size - imported_cases, imported_cases);
  
  # vaccine uptake states
  # c_vac <- c(c_S,c_E,c_I_presym,c_I_asym,c_R,c_Evoc,c_Ivoc_presym,c_Ivoc_asym,c_Rvoc)
  c_vac <- c(c_S,c_R,c_Rvoc)
  # reinfection states
  c_reinfection <- c(c_R,c_Rvoc)
  #voc/nonvoc
  c_nonvoc <- c(c_E,c_I_presym,c_I_asym,c_I_mild,c_I_sev,c_I_hosp,c_I_icu,c_D,c_R)
  c_voc <- c(c_Evoc,c_Ivoc_presym,c_Ivoc_asym,c_Ivoc_mild,c_Ivoc_sev,c_Ivoc_hosp,c_Ivoc_icu,c_Dvoc,c_Rvoc)
  
  
  ## Transition probabilities  ---- 
  ##----------------------------- -
  col_vec = c(c_S,              #1:10,
              c_E,  #11:20
              
              # Ipresym to I_asym en I_mild
              c_I_presym,       #21:30,
              c_I_presym,       #21:30,
              
              # Severe infections 
              c_I_mild,       #41:50
              
              # Hospital and ICU admission
              c_I_sev,          #51:60,
              c_I_sev,         #51:60,
              
              # Mortality
              c_I_hosp,
              c_I_icu,          #61:80,
              
              # Recovery
              c_I_asym, 
              c_I_mild,         #31:50,
              c_I_hosp,  
              c_I_icu,          #61:80, # to "regular" hospital compartment
              
              # idem for VoC
              c_S,                    #1:10,
              c_Evoc,c_Ivoc_presym,   #131:150,
              c_Ivoc_presym,          #141:150,
              c_Ivoc_mild,c_Ivoc_sev, #161:180,
              c_Ivoc_sev,c_Ivoc_hosp,c_Ivoc_icu,#171:200,
              c_Ivoc_asym, c_Ivoc_mild, #151:170
              c_Ivoc_hosp, c_Ivoc_icu    #181:200
  );
  
  ### Transitions for nfull
  ###-------------------- -
  prob_vec <- c(
    #Enew = (1 - exp(-h*lambda)),                  #
    I_presym_new =  rep((1 - exp(-h*gamma)),10),
    
    I_asym_new   =  (1 - exp(-h*p*theta)),
    I_mild_new   =  (1 - exp(-h*(1-p)*theta)), 
    
    I_sev_new    =  (1 - exp(-h*(1-phi0)*delta2)), 
    
    I_hosp_new   =  (1 - exp(-h*phi1*omega)),
    I_icu_new    =  (1 - exp(-h*(1 - phi1)*omega)),
    
    D_hosp_new   =  (1 - exp(-h*delta3*mu_sev)),
    D_icu_new    =  (1 - exp(-h*delta3*mu_sev)),
    
    R_asym_new   =  rep(1 - exp(-h*delta1), each = 10),
    R_mild_new   =  (1 - exp(-h*phi0*delta2)), 
    R_hosp_new   =  (1 - exp(-h*delta3*(1-mu_sev))),
    R_icu_new    =  (1 - exp(-h*delta3*(1-mu_sev)))
  )
  
  ### Transitions with vaccination (nfull_vac_rna, nfull_vac_adeno)
  get_prob_vec_vacc <- function(prob_vec,sel_parms,h,sel_phi0,delta2,ve_tag){
    
    # changed hospital admission probability
    prob_vec_vacc <- matrix(rep(prob_vec,14),ncol=length(prob_vec),byrow = T)
    prob_vec_vacc[1,grepl('I_sev_new',names(prob_vec))]   <- 1 - exp(-h*(1-sel_phi0)*delta2)  # no vaccine
    prob_vec_vacc[1,grepl('R_mild_new',names(prob_vec))]  <- 1 - exp(-h* sel_phi0   *delta2)  # no vaccine
    prob_vec_vacc[2,grepl('I_sev_new',names(prob_vec))]   <- 1 - exp(-h*(1-get_phi0_vac(sel_phi0,sel_parms[paste0(ve_tag,'_incr_severe_rna1')]))*delta2)
    prob_vec_vacc[2,grepl('R_mild_new',names(prob_vec))]  <- 1 - exp(-h*   get_phi0_vac(sel_phi0,sel_parms[paste0(ve_tag,'_incr_severe_rna1')]) *delta2)
    prob_vec_vacc[3,grepl('I_sev_new',names(prob_vec))]   <- 1 - exp(-h*(1-get_phi0_vac(sel_phi0,sel_parms[paste0(ve_tag,'_incr_severe_rna2')]))*delta2)
    prob_vec_vacc[3,grepl('R_mild_new',names(prob_vec))]  <- 1 - exp(-h*   get_phi0_vac(sel_phi0,sel_parms[paste0(ve_tag,'_incr_severe_rna2')]) *delta2)
    prob_vec_vacc[4,grepl('I_sev_new',names(prob_vec))]   <- 1 - exp(-h*(1-get_phi0_vac(sel_phi0,sel_parms[paste0(ve_tag,'_incr_severe_adeno1')]))*delta2)
    prob_vec_vacc[4,grepl('R_mild_new',names(prob_vec))]  <- 1 - exp(-h*   get_phi0_vac(sel_phi0,sel_parms[paste0(ve_tag,'_incr_severe_adeno1')]) *delta2)   
    prob_vec_vacc[5,grepl('I_sev_new',names(prob_vec))]   <- 1 - exp(-h*(1-get_phi0_vac(sel_phi0,sel_parms[paste0(ve_tag,'_incr_severe_adeno2')]))*delta2)
    prob_vec_vacc[5,grepl('R_mild_new',names(prob_vec))]  <- 1 - exp(-h*   get_phi0_vac(sel_phi0,sel_parms[paste0(ve_tag,'_incr_severe_adeno2')]) *delta2)
    prob_vec_vacc[6,grepl('I_sev_new',names(prob_vec))]   <- 1 - exp(-h*(1-get_phi0_vac(sel_phi0,sel_parms[paste0(ve_tag,'_incr_severe_waning')]))*delta2)
    prob_vec_vacc[6,grepl('R_mild_new',names(prob_vec))]  <- 1 - exp(-h*   get_phi0_vac(sel_phi0,sel_parms[paste0(ve_tag,'_incr_severe_waning')]) *delta2)
    prob_vec_vacc[7,grepl('I_sev_new',names(prob_vec))]   <- 1 - exp(-h*(1-get_phi0_vac(sel_phi0,sel_parms[paste0(ve_tag,'_incr_severe_booster')]))*delta2)
    prob_vec_vacc[7,grepl('R_mild_new',names(prob_vec))]  <- 1 - exp(-h*   get_phi0_vac(sel_phi0,sel_parms[paste0(ve_tag,'_incr_severe_booster')]) *delta2)
    prob_vec_vacc[8,grepl('I_sev_new',names(prob_vec))]   <- 1 - exp(-h*(1-get_phi0_vac(sel_phi0,sel_parms[paste0(ve_tag,'_incr_severe_booster_waning')]))*delta2)
    prob_vec_vacc[8,grepl('R_mild_new',names(prob_vec))]  <- 1 - exp(-h*   get_phi0_vac(sel_phi0,sel_parms[paste0(ve_tag,'_incr_severe_booster_waning')]) *delta2)
    prob_vec_vacc[9,grepl('I_sev_new',names(prob_vec))]   <- 1 - exp(-h*(1-get_phi0_vac(sel_phi0,sel_parms[paste0(ve_tag,'_incr_severe_reinf')]))*delta2)
    prob_vec_vacc[9,grepl('R_mild_new',names(prob_vec))]  <- 1 - exp(-h*   get_phi0_vac(sel_phi0,sel_parms[paste0(ve_tag,'_incr_severe_reinf')]) *delta2)
    prob_vec_vacc[10,grepl('I_sev_new',names(prob_vec))]   <- 1 - exp(-h*(1-get_phi0_vac(sel_phi0,sel_parms[paste0(ve_tag,'_incr_severe_reinf')]))*delta2)  #oldvoc
    prob_vec_vacc[10,grepl('R_mild_new',names(prob_vec))]  <- 1 - exp(-h*   get_phi0_vac(sel_phi0,sel_parms[paste0(ve_tag,'_incr_severe_reinf')]) *delta2)
    prob_vec_vacc[11,grepl('I_sev_new',names(prob_vec))]   <- 1 - exp(-h*(1-get_phi0_vac(sel_phi0,sel_parms[paste0(ve_tag,'_incr_severe_reinfvac')]))*delta2)
    prob_vec_vacc[11,grepl('R_mild_new',names(prob_vec))]  <- 1 - exp(-h*   get_phi0_vac(sel_phi0,sel_parms[paste0(ve_tag,'_incr_severe_reinfvac')]) *delta2)
    prob_vec_vacc[12,grepl('I_sev_new',names(prob_vec))]   <- 1 - exp(-h*(1-get_phi0_vac(sel_phi0,sel_parms[paste0(ve_tag,'_incr_severe_reinfvac')]))*delta2)  #oldvoc
    prob_vec_vacc[12,grepl('R_mild_new',names(prob_vec))]  <- 1 - exp(-h*   get_phi0_vac(sel_phi0,sel_parms[paste0(ve_tag,'_incr_severe_reinfvac')]) *delta2)
    prob_vec_vacc[13,grepl('I_sev_new',names(prob_vec))]   <- 1 - exp(-h*(1-get_phi0_vac(sel_phi0,sel_parms[paste0(ve_tag,'_incr_severe_extrabooster')]))*delta2)
    prob_vec_vacc[13,grepl('R_mild_new',names(prob_vec))]  <- 1 - exp(-h*   get_phi0_vac(sel_phi0,sel_parms[paste0(ve_tag,'_incr_severe_extrabooster')]) *delta2)
    prob_vec_vacc[14,grepl('I_sev_new',names(prob_vec))]   <- 1 - exp(-h*(1-get_phi0_vac(sel_phi0,sel_parms[paste0(ve_tag,'_incr_severe_extrabooster')]))*delta2)  #to adapt
    prob_vec_vacc[14,grepl('R_mild_new',names(prob_vec))]  <- 1 - exp(-h*   get_phi0_vac(sel_phi0,sel_parms[paste0(ve_tag,'_incr_severe_extrabooster')]) *delta2)
    #add correction if incr_severe is greater than 1, should be corrected before - temp
    prob_vec_vacc         <- pmax(prob_vec_vacc,0)
    # return
    return(prob_vec_vacc)
  }
  ###------------------------- -
  # adjust hospital admission probability to vaccine-related protection 
  prob_vec_vacc     <- get_prob_vec_vacc(prob_vec,parms,h,phi0,delta2,ve_tag='ve_VOC_alpha')
  prob_vec_vacc_voc <- prob_vec_vacc
  
  
  # set number of vaccine-related types
  nrow_vec_vacc <- nrow(prob_vec_vacc)
  
  ############################################################### #
  #  Description of transitions:                                  #
  #                                                               #
  #  S_t+h = S_t - E_t+h                                          # 
  #        = tmp[1:10] - new_trans[1:10]                          #
  #                                                               #
  #  with new_trans including                                     #
  #     Enew: 1-10          (1 - exp(-h*lambda))                  #
  #     I_presym_new: 11-20 (1 - exp(-h*gamma))                  #
  #     I_asym_new: 21-30   (1 - exp(-h*p*theta))                 #                    
  #     I_mild_new: 31-40   (1 - exp(-h*(1-p)*theta))             #                     
  #     I_sev_new: 41-50    (1 - exp(-h*(1-phi0)*delta2)          #  
  #     I_hosp_new: 51-60   (1 - exp(-h*phi1*omega))              #
  #     I_icu_new: 61-70    (1 - exp(-h*(1-phi1)*omega))              #
  #     D_hosp_new: 71-80   (1 - exp(-h*mu_sev*delta3))             #
  #     D_icu_new: 81-90    (1 - exp(-h*mu_sev*delta3))             #
  #     R_asym_new: 91-100  (1 - exp(-h*delta1))                  #
  #     R_mild_new: 101-110 (1 - exp(-h*phi0*delta2))             #
  #     R_hosp_new: 111-120 (1 - exp(-h*(1-mu_sev)*delta3))       #    
  #     R_icu_new: 121-130  (1 - exp(-h*(1-mu_sev)*delta3))       #
  #                                                               #
  #     Evoc_new: 131-140         (1 - exp(-h*lambda*transm_fctr_VOC)) #
  #     Ivoc_presym_new: 141-150  (1 - exp(-h*gamma))            #
  #     Ivoc_asym_new: 151-160    (1 - exp(-h*p*theta))           #
  #     Ivoc_mild_new: 161-170    (1 - exp(-h*(1-p)*theta))       #                     
  #     Ivoc_sev_new: 171-180     (1 - exp(-h*(1-phi0)*delta2)    #  
  #     Ivoc_hosp_new: 181-190    (1 - exp(-h*phi1*omega))        #
  #     Ivoc_icu_new: 191-200     (1 - exp(-h*(1-phi1)*omega))        #
  #     Dvoc_hosp_new: 201-210    (1 - exp(-h*mu_sev*delta3))       #
  #     Dvoc_icu_new: 211-220     (1 - exp(-h*mu_sev*delta3))       #
  #     Rvoc_asym_new: 221-230    (1 - exp(-h*delta1))            #
  #     Rvoc_mild_new: 231-240    (1 - exp(-h*phi0*delta2))       #
  #     Rvoc_hosp_new: 241-250    (1 - exp(-h*(1-mu_sev)*delta3)) #    
  #     Rvoc_icu_new: 251-260     (1 - exp(-h*(1-mu_sev)*delta3)) #
  ############################################################### #
  
  t_names <- c(names(icol_tmp_vec)[2:8],
               paste(names(icol_tmp_vec)[9],c('hosp','icu'),sep='_'),
               paste(names(icol_tmp_vec)[10],c('asym','mild','hosp','icu'),sep='_'))
  t_names <- gsub('c_','t_',t_names)
  t_names <- c(t_names,paste0(substr(t_names,0,3),'voc',substr(t_names,4,20)))
  icol_transitions <- data.frame(matrix(1:(length(t_names)*10),nrow=10,byrow = F))
  names(icol_transitions) <- t_names
  attach(icol_transitions)
  
  min_transition_vec1 = 1 + c(t_E, #1:10,
                              t_I_presym, #11:20,
                              t_I_asym, #21:30,
                              t_R_asym, #91:100,
                              t_I_sev, #41:50,
                              t_I_hosp, #51:60,
                              t_D_hosp, #71:80,
                              t_D_icu, #81:90,
                              rep(0,50), # R, D and additional computations (no exit) 
                              
                              t_Ivoc_presym, #141:150,
                              t_Ivoc_asym, #151:160,
                              t_Rvoc_asym, #221:230,
                              t_Ivoc_sev, #171:180,
                              t_Ivoc_hosp, #181:190,
                              t_Dvoc_hosp, #201:210,
                              t_Dvoc_icu, #211:220,
                              rep(0,50), # R, D and additional computations (no exit) 
                              rep(0,60)) # to register the mild infections and icu admissions
  
  min_transition_vec2 = 1 + c(t_Evoc, #131:140,
                              rep(0,10),
                              t_I_mild, #31:40,
                              rep(0,10),
                              t_R_mild, #101:110,
                              t_I_icu,  #61:70,
                              t_R_hosp, #111:120,
                              t_R_icu,  #121:130,
                              rep(0,50),
                              rep(0,10),
                              t_Ivoc_mild, #161:170,
                              rep(0,10),
                              t_Rvoc_mild, #231:240,
                              t_Ivoc_icu,  #191:200,
                              t_Rvoc_hosp, #241:250,
                              t_Rvoc_icu,  #251:260,
                              rep(0,50),
                              rep(0,60)) # to register the mild infections and icu admissions
  
  plus_transition_vec1 = 1 + c(rep(0,10),
                               1:70,
                               t_D_hosp,# rep(0,10),   # no transition from I_hosp to Death
                               t_R_asym, #91:100,
                               t_I_hosp, #51:60,    # additional computations
                               t_I_asym, #21:30,    # additional computations
                               t_R_hosp, #111:120,  # additional computations
                               131:200,
                               t_Dvoc_hosp, # rep(0,10), # no transition from I_hosp to Death
                               t_Rvoc_asym, #221:230,
                               t_Ivoc_hosp, #181:190, # additional computations
                               t_Ivoc_asym, #151:160, # additional computations
                               t_Rvoc_hosp, #241:250) # additional computations
                               t_I_mild,    # to register the mild infections
                               t_Ivoc_mild, # to register the mild infections
                               t_R_mild,    # to register the recovery from mild infections
                               t_Rvoc_mild, # to register the recovery from mild infections
                               t_I_icu,     # to register new ICU admissions
                               t_Ivoc_icu) # to register new ICU_voc admissions
  
  plus_transition_vec2 = 1 + c(rep(0,10),
                               rep(0,70),
                               t_D_icu,  #81:90, 
                               t_R_mild, #101:110,
                               t_I_icu,  #61:70,   # additional computations
                               t_I_mild, #31:40,   # additional computations
                               t_R_icu, #rep(0,10), #t_R_icu,  #121:130, # ICU is transered to I_hosp now, not accounted anymore as "hospital recovery"
                               rep(0,70),
                               t_Dvoc_icu, #211:220,
                               t_Rvoc_mild, #231:240,
                               t_Ivoc_icu, #191:200,  # additional computations
                               t_Ivoc_mild, #161:170, # additional computations
                               t_Rvoc_icu, #rep(0,10), #  t_Rvoc_icu, #251:260) # this should go
                               rep(0,60))   # to register the mild infections and ICU admissions
  
  plus_transition_vec3 = 1 + c(rep(0,90),
                               t_R_hosp, #111:120,
                               rep(0,30),
                               rep(0,80),
                               t_Rvoc_hosp, #241:250,
                               rep(0,30),
                               rep(0,60))   # to register the mild infections and icu admissions
  
  detach(icol_transitions)
  
  # Initialise population details for time step "t"
  tmp             = nfull[[1]];
  tmp_vac_rna1    = nfull_vac_rna1[[1]];
  tmp_vac_rna2    = nfull_vac_rna2[[1]];
  tmp_vac_adeno1  = nfull_vac_adeno1[[1]];
  tmp_vac_adeno2  = nfull_vac_adeno2[[1]];
  tmp_vac_waning  = nfull_vac_waning[[1]];
  tmp_vac_booster = nfull_vac_booster[[1]];
  tmp_vac_booster_waning = nfull_vac_booster_waning[[1]];
  tmp_vac_reinf   = nfull_vac_reinf[[1]];
  tmp_vac_reinf_oldvoc   = nfull_vac_reinf_oldvoc[[1]];
  tmp_vac_reinfvac = nfull_vac_reinfvac[[1]];
  tmp_vac_reinfvac_oldvoc = nfull_vac_reinfvac_oldvoc[[1]];
  tmp_vac_extrabooster1 = nfull_vac_extrabooster1[[1]];
  tmp_vac_extrabooster2 = nfull_vac_extrabooster2[[1]];
  
  # CoMix matrices and social contact behavior  ----
  #----------------------------------------------- -
  M_list_all <- get_susceptiblity_matrices(CoMix_matrices, 
                                           CoMix_coef_mat, 
                                           Add_coef_mat,
                                           cnt_change_pt,cnt_change)
  
  # adapt to regional simulation?
  if(pop_factor <1){
    sel_items <- which(grepl('sy',names(M_list_all)))
    for(i in sel_items){
      M_list_all[[i]] <- M_list_all[[i]] * 1/pop_factor
    }
  }
  
  # start matrices
  M_asy <- M_list_all$M2010_asy
  M_sy  <- M_list_all$M2010_sy
  M_asy_bar <- M_list_all$M2010_asy
  M_sy_bar  <- M_list_all$M2010_sy
  bar_flag <-rep("yellow",length(times))
  bar_coef <-rep(1,length(times))
  # run parameters
  bool_vaccine    <- FALSE
  bool_voc        <- FALSE
  
  new_transitions            <- rep(0,max(icol_transitions)+1) #rep(0,261)
  new_transitions_vac_rna1   <- new_transitions
  new_transitions_vac_rna2   <- new_transitions
  new_transitions_vac_adeno1 <- new_transitions
  new_transitions_vac_adeno2 <- new_transitions
  new_transitions_vac_waning <- new_transitions
  new_transitions_vac_booster<- new_transitions
  new_transitions_vac_booster_waning <- new_transitions
  new_transitions_vac_reinf  <- new_transitions
  new_transitions_vac_reinf_oldvoc  <- new_transitions
  new_transitions_vac_reinfvac <- new_transitions
  new_transitions_vac_reinfvac_oldvoc <- new_transitions
  new_transitions_vac_extrabooster1<- new_transitions
  new_transitions_vac_extrabooster2<- new_transitions
  
  testrealvacin <- rep(0,10)  
  testvacin     <- rep(0,10)
  testvacout    <- rep(0,10)
  
  total_vac_rna1_remain    <- rep(0,length(c_S))
  total_vac_rna2_remain    <- total_vac_rna1_remain
  total_vac_adeno1_remain  <- total_vac_rna1_remain
  total_vac_adeno2_remain  <- total_vac_rna1_remain
  total_vac_booster_remain <- total_vac_rna1_remain
  total_vac_2ndbooster_remain  <- rep(0,10)
  total_vac_extrabooster_remain  <- rep(0,10)
  norm_sy <- rep(0,length(times))
  norm_asy <- rep(0,length(times))
  step_index <- 1
  # for compatibility
  sel_VOC <- "VOC_delta"
  # note: 
  # - each time step starts from the situation at nfull[[step_index]]
  # - the results of a time step are stored at nfull[[step_index+1]]
  # - so last "step_index" is restricted to "length(times)-1"
  ## LOOP OVER TIME STEPS ----
  for (step_index in 1:(length(times)-1)){
    #if(step_index == 1250) {cat('BREAK!'); break}
    
    calendar_time = times[step_index]
   # print(calendar_time)
   # print(sim_day2date(calendar_time))
   # if(calendar_time == 974) {cat('BREAK!'); browser() }
    #if(calendar_time == VOC_start[1]) {cat('BREAK!'); break }
    
    
    #update contact matrix?
    cnt_out <- get_contact_details(calendar_time,step_index,
                                   M_list_all,
                                   # beta0, beta1,
                                   M_asy,M_sy,
                                   cnt_change_pt,cnt_change,nfull,nfull_vac_reinf,nfull_vac_reinf_oldvoc)
    M_asy <- cnt_out$M_asy
    M_sy  <- cnt_out$M_sy
    M_sy_bar <- M_sy
    M_asy_bar <- M_asy
    #if (calendar_time > sim_date2day("2024-01-02")) browser()
    if(bool_bar==TRUE){
    if(calendar_time >sim_date2day("2021-07-01")){
    #print(step_index)
         icu_beds=sum(nfull[[step_index]][c_I_icu]+nfull[[step_index]][c_Ivoc_icu]+nfull_vac_reinf[[step_index]][c_I_icu]+nfull_vac_reinf[[step_index]][c_Ivoc_icu]+nfull_vac_reinf_oldvoc[[step_index]][c_I_icu]+nfull_vac_reinf_oldvoc[[step_index]][c_Ivoc_icu])
         new_hospi=0
         begin=step_index-3
         for (i in begin:step_index){
         new_hospi=new_hospi+sum(nfull[[i]][c_new_hosp_total]+nfull[[i]][c_new_hosp_total_voc]+nfull_vac_reinf[[i]][c_new_hosp_total]+nfull_vac_reinf[[i]][c_new_hosp_total_voc]+nfull_vac_reinf_oldvoc[[i]][c_new_hosp_total]+nfull_vac_reinf_oldvoc[[i]][c_new_hosp_total_voc])
         }  
         if((new_hospi >= new_hosp_or && new_hospi < new_hosp_red) || (icu_beds>= icu_or &icu_beds < icu_red)) {
              bar_flag[step_index] <- "orange"
              if(bar_coef[step_index-1]>or_prop) bar_coef[step_index] <- pmax(bar_coef[step_index-1]-step_prop,or_prop)
              else bar_coef[step_index] <- pmin(bar_coef[step_index-1]+step_prop,or_prop)
              }
         else if(new_hospi>= new_hosp_red || icu_beds >= icu_red ){
              bar_flag[step_index] <-"red"
              bar_coef[step_index] <- pmax(bar_coef[step_index-1]-step_prop,red_prop)
         }
         else bar_coef[step_index] <- pmin(bar_coef[step_index-1]+step_prop,1)
    }
    M_sy_bar <- bar_coef[step_index]*M_sy
    M_asy_bar <-bar_coef[step_index]*M_asy
    }
    # if(calendar_time >sim_date2day("2022-01-01") && calendar_time < sim_date2day("2022-04-01")){
    #   #print(calendar_time)
    #   icu_beds=sum(nfull[[step_index]][c_I_icu]+nfull[[step_index]][c_Ivoc_icu]+nfull_vac_reinf[[step_index]][c_I_icu]+nfull_vac_reinf[[step_index]][c_Ivoc_icu]+nfull_vac_reinf_oldvoc[[step_index]][c_I_icu]+nfull_vac_reinf_oldvoc[[step_index]][c_Ivoc_icu])
    #   new_hospi=sum(nfull[[step_index]][c_new_hosp_total]+nfull[[step_index]][c_new_hosp_total_voc]+nfull_vac_reinf[[step_index]][c_new_hosp_total]+nfull_vac_reinf[[step_index]][c_new_hosp_total_voc]+nfull_vac_reinf_oldvoc[[step_index]][c_new_hosp_total]+nfull_vac_reinf_oldvoc[[step_index]][c_new_hosp_total_voc])
    #   if((new_hospi > 64 && new_hospi <150) || (icu_beds>300 &icu_beds<=500)) {
    #     print("orange")
    #     M_sy=0.9*M_sy
    #     M_asy=0.9*M_asy
    #   }
    #   else if(new_hospi>149 || icu_beds > 500 ){
    #     print("red")
    #     M_sy=0.7*M_sy
    #     M_asy=0.7*M_asy
    #   }
    # }
    # if(calendar_time >sim_date2day("2022-01-01") && calendar_time < sim_date2day("2022-04-01")){
    #   #print(calendar_time)
    #   icu_beds=sum(nfull[[step_index]][c_I_icu]+nfull[[step_index]][c_Ivoc_icu]+nfull_vac_reinf[[step_index]][c_I_icu]+nfull_vac_reinf[[step_index]][c_Ivoc_icu]+nfull_vac_reinf_oldvoc[[step_index]][c_I_icu]+nfull_vac_reinf_oldvoc[[step_index]][c_Ivoc_icu])
    #   new_hospi=sum(nfull[[step_index]][c_new_hosp_total]+nfull[[step_index]][c_new_hosp_total_voc]+nfull_vac_reinf[[step_index]][c_new_hosp_total]+nfull_vac_reinf[[step_index]][c_new_hosp_total_voc]+nfull_vac_reinf_oldvoc[[step_index]][c_new_hosp_total]+nfull_vac_reinf_oldvoc[[step_index]][c_new_hosp_total_voc])
    #   if((new_hospi > 64 && new_hospi <150) || (icu_beds>300 &icu_beds<=500)) {
    #     print("orange")
    #     M_sy=0.9*M_sy
    #     M_asy=0.9*M_asy
    #   }
    #   else if(new_hospi>149 || icu_beds > 500 ){
    #     print("red")
    #     M_sy=0.7*M_sy
    #     M_asy=0.7*M_asy
    #   }
    # }
    # if(calendar_time >sim_date2day("2022-04-01") && calendar_time < sim_date2day("2022-09-01")){
    #   #print(calendar_time)
    #   icu_beds=sum(nfull[[step_index]][c_I_icu]+nfull[[step_index]][c_Ivoc_icu]+nfull_vac_reinf[[step_index]][c_I_icu]+nfull_vac_reinf[[step_index]][c_Ivoc_icu]+nfull_vac_reinf_oldvoc[[step_index]][c_I_icu]+nfull_vac_reinf_oldvoc[[step_index]][c_Ivoc_icu])
    #   new_hospi=sum(nfull[[step_index]][c_new_hosp_total]+nfull[[step_index]][c_new_hosp_total_voc]+nfull_vac_reinf[[step_index]][c_new_hosp_total]+nfull_vac_reinf[[step_index]][c_new_hosp_total_voc]+nfull_vac_reinf_oldvoc[[step_index]][c_new_hosp_total]+nfull_vac_reinf_oldvoc[[step_index]][c_new_hosp_total_voc])
    #   if((new_hospi > 64 && new_hospi <150) || (icu_beds>300 &icu_beds<=500)) {
    #     print("orange")
    #     M_sy=0.99*M_sy
    #     M_asy=0.99*M_asy
    #   }
    #   else if(new_hospi>149 || icu_beds > 500 ){
    #     print("red")
    #     M_sy=0.85*M_sy
    #     M_asy=0.85*M_asy
    #   }
    # }
    # if(calendar_time >=sim_date2day("2022-09-01") && calendar_time < sim_date2day("2023-05-01")){
    #   #print(calendar_time)
    #   icu_beds=sum(nfull[[step_index]][c_I_icu]+nfull[[step_index]][c_Ivoc_icu]+nfull_vac_reinf[[step_index]][c_I_icu]+nfull_vac_reinf[[step_index]][c_Ivoc_icu]+nfull_vac_reinf_oldvoc[[step_index]][c_I_icu]+nfull_vac_reinf_oldvoc[[step_index]][c_Ivoc_icu])
    #   new_hospi=sum(nfull[[step_index]][c_new_hosp_total]+nfull[[step_index]][c_new_hosp_total_voc]+nfull_vac_reinf[[step_index]][c_new_hosp_total]+nfull_vac_reinf[[step_index]][c_new_hosp_total_voc]+nfull_vac_reinf_oldvoc[[step_index]][c_new_hosp_total]+nfull_vac_reinf_oldvoc[[step_index]][c_new_hosp_total_voc])
    #   if((new_hospi > 64 && new_hospi <150) || (icu_beds>300 &icu_beds<=500)) {
    #     print("orange")
    #     M_sy=0.7*M_sy
    #     M_asy=0.7*M_asy
    #   }
    #   else if(new_hospi>149 || icu_beds > 500 ){
    #     print("red")
    #     M_sy=0.5*M_sy
    #     M_asy=0.5*M_asy
    #   }
    # }
    # if(step_index==300){
    #   spot
    # }
    #print(calendar_time)
    # conf_date1 <- "2021-07-01"
    # if(calendar_time>=sim_date2day(conf_date1))
    #  {
    #    M_asy <- M_asy*(0.8+0.2*min(c((calendar_time-sim_date2day(conf_date1))/75,1)))
    #    M_sy <- M_sy*(0.5+0.2*min(c((calendar_time-sim_date2day(conf_date1))/75,1)))
    # }
    # conf_date2 <- "2021-10-01"
    # if(calendar_time>=sim_date2day(conf_date2))
    # {
    #   M_asy <- M_asy*(0.5+0.5*min(c((calendar_time-sim_date2day(conf_date2))/100,1)))
    #   M_sy <- M_sy*0.7
    # }
    # 
    # conf_date3 <- "2022-01-15"
    # if(calendar_time>=sim_date2day(conf_date3))
    # {
    #   M_asy <- M_asy*(0.7+0.3*min(c((calendar_time-sim_date2day(conf_date3))/75,1)))
    #   M_sy <- M_sy*0.7
    # }
    # 
    # conf_date4 <- "2022-05-15"
    # if(calendar_time>=sim_date2day(conf_date4))
    # {
    #   M_asy <- M_asy*(0.5+0.5*min(c((calendar_time-sim_date2day(conf_date4))/400,1)))
    #   M_sy <- M_sy*0.7
    # }
    # 
    # conf_date5 <- "2023-09-15"
    # if(calendar_time>=sim_date2day(conf_date5))
    # {
    #   M_asy <- M_asy*(0.3+0.7*min(c((calendar_time-sim_date2day(conf_date5))/400,1)))
    #   M_sy <- M_sy*0.7
    # }


    
    
    
    betas_asy = q_asy * M_asy_bar
    betas_sy  = q_sy * M_sy_bar
    norm_asy[step_index]=norm(M_asy_bar,type="F")
    norm_sy[step_index]=norm(M_sy_bar,type="F")
    #print(betas_asy)
    #print(betas_sy)

    ## Update Mortality and Recovery      ----
    ##----------------------------- -
    if(calendar_time == mortality_changepoint){
      
      prob_vec_vacc[,grepl('D_icu_new',names(prob_vec))]               =
        prob_vec_vacc[,grepl('D_hosp_new',names(prob_vec))]              = matrix(rep(1 - exp(-h*delta3_wave2*mu2_sev),nrow_vec_vacc),ncol=10,byrow = T)
      
      prob_vec_vacc[,grepl('R_icu_new',names(prob_vec))]                = 
        prob_vec_vacc[,grepl('R_hosp_new',names(prob_vec))]               = matrix(rep(1 - exp(-h*delta3_wave2*(1-mu2_sev)),nrow_vec_vacc),ncol=10,byrow = T)
    }
    
    #  if(calendar_time==1300){browser()}
    
    ## Update VOC_delta & VOC_omicron          ----
    ##----------------------------- -
    if(any(calendar_time == VOC_start)){
      
      #print(calendar_time)
      # update optimisation boolean
      bool_voc <- TRUE;
      
      ## 1. DEFINE VOC-SPECIFIC VALUES
      
      # set VOC name and parameters
      sel_VOC              <- VOC_name[VOC_start == calendar_time]
      sel_VOC_param        <- parms[grepl(sel_VOC,parms_names)]
      names(sel_VOC_param) <- gsub(paste0(sel_VOC,'_'),'',names(sel_VOC_param))
      names(sel_VOC_param) <- gsub(paste0('_',sel_VOC),'',names(sel_VOC_param))
      sel_VOC_param
      
      # VOC parameters  
      VOC_init        <- exp(sel_VOC_param['log_init']) * pop_factor # population factor for regional simulation
      VOC_hr_hosp     <- exp(sel_VOC_param[grepl('log_hr_hosp',names(sel_VOC_param))]);   # hospital hazard ratio vs previous VOC
      log_VOC_XXX_hosp<- get_log_VOC_hosp(phi0,delta2,h, odds_ratio = VOC_hr_hosp, log_VOC_previous = log_VOC_XXX_hosp)
      VOC_hosp        <- expit(log_VOC_XXX_hosp)
      VOC_phi0        <- phi0 - (phi0 * VOC_hosp)
      VOC_phi1        <- expit(parms[grepl('log_phi1_age',parms_names)] * sel_VOC_param['log_phi1_add']);  # proportion of severly infected with regular hospitalization
      VOC_delta3      <- exp(sel_VOC_param['log_delta3']);    # recovery rate from hospital (from delta VOC)
      
      log_VOC_mu_param <- sel_VOC_param[grepl('log_mu_sev',names(sel_VOC_param))]
      VOC_mu_sev       <- c(0,expit(log_VOC_mu_param[c(1,1:8)])); # age-specific hospital fatality ratio - probability of dying
      
      VOC_gamma        <- exp(parms['log_gamma']) * ifelse('log_gamma_factor' %in% names(sel_VOC_param), exp(sel_VOC_param['log_gamma_factor']), 1)
      
      # update severity, with respect to hospital admission
      VOC_prob_vec_vacc <- get_prob_vec_vacc(prob_vec,sel_VOC_param,h,VOC_phi0,delta2,ve_tag='ve')
      
      # update severity, with respect to probability to be admitted to ICU
      VOC_prob_vec_vacc[,grepl('I_hosp_new',names(prob_vec))]   <- matrix(rep(1 - exp(-h*  VOC_phi1  *omega),nrow_vec_vacc),ncol=10,byrow = T)
      VOC_prob_vec_vacc[,grepl('I_icu_new',names(prob_vec))]    <- matrix(rep(1 - exp(-h*(1-VOC_phi1)*omega),nrow_vec_vacc),ncol=10,byrow = T)
      
      # update hospital length of stay
      VOC_prob_vec_vacc[,grepl('D_icu_new',names(prob_vec))]               =
        VOC_prob_vec_vacc[,grepl('D_hosp_new',names(prob_vec))]              = matrix(rep(1 - exp(-h*VOC_delta3*VOC_mu_sev),nrow_vec_vacc),ncol=10,byrow = T)
      
      VOC_prob_vec_vacc[,grepl('R_icu_new',names(prob_vec))]               = 
        VOC_prob_vec_vacc[,grepl('R_hosp_new',names(prob_vec))]              = matrix(rep(1 - exp(-h*VOC_delta3*(1-VOC_mu_sev)),nrow_vec_vacc),ncol=10,byrow = T)
      
      # update latency period (optional)
      VOC_prob_vec_vacc[,grepl('I_presym_new',names(prob_vec))]              = (1 - exp(-h*VOC_gamma))
      
      ## 2. REGISTER VOC-SPECIFIC VALUES in STRAIN 1 OR 2
      if(global_lib_voc$model_strain[global_lib_voc$name == sel_VOC] == 1){
        # switch transmission potential
        transm_fctr <- exp(sel_VOC_param['log_transm']);   # additional VoC transmissibility
        
        # switch VE against infection for wild type strain
        ve_infection <- get_ve_infection(sel_VOC_param,ve_tag='ve')
        
        #change immunity vs old_current voc
        if(calendar_time >= sim_date2day(as.Date(start_date_voc_immune_evasion))){
          ve_infection$reinf     <- ve_infection$reinf
          ve_voc_infection$reinf     <- ve_infection$reinf
          ve_infection$reinf_oldvoc     <- ve_infection$reinf_oldvoc
          ve_voc_infection$reinf_oldvoc     <- ve_infection$reinf
         # print(sel_VOC)
         # print(c( ve_infection$reinf,ve_voc_infection$reinf,ve_infection$reinf_oldvoc,ve_voc_infection$reinf_oldvoc))
         # print("transfert of all old voc in non current voc (1)")
          tmp[c_voc] <- tmp[c_voc] + tmp[c_nonvoc]
          tmp[c_nonvoc] <- 0
          tmp_vac_rna1[c_voc] <- tmp_vac_rna1[c_voc] + tmp_vac_rna1[c_nonvoc]
          tmp_vac_rna1[c_nonvoc] <- 0
          tmp_vac_rna2[c_voc] <- tmp_vac_rna2[c_voc] + tmp_vac_rna2[c_nonvoc]
          tmp_vac_rna2[c_nonvoc] <- 0
          tmp_vac_adeno1[c_voc] <- tmp_vac_adeno1[c_voc] + tmp_vac_adeno1[c_nonvoc]
          tmp_vac_adeno1[c_nonvoc] <- 0
          tmp_vac_adeno2[c_voc] <- tmp_vac_adeno2[c_voc] + tmp_vac_adeno2[c_nonvoc]
          tmp_vac_adeno2[c_nonvoc] <- 0
          tmp_vac_waning[c_voc] <- tmp_vac_waning[c_voc] + tmp_vac_waning[c_nonvoc]
          tmp_vac_waning[c_nonvoc] <- 0
          tmp_vac_booster[c_voc] <- tmp_vac_booster[c_voc] + tmp_vac_booster[c_nonvoc]
          tmp_vac_booster[c_nonvoc] <- 0
          tmp_vac_booster_waning[c_voc] <- tmp_vac_booster_waning[c_voc] + tmp_vac_booster_waning[c_nonvoc]
          tmp_vac_booster_waning[c_nonvoc] <- 0
          tmp_vac_reinf[c_voc] <- tmp_vac_reinf[c_voc] + tmp_vac_reinf[c_nonvoc]
          tmp_vac_reinf[c_nonvoc] <- 0
          tmp_vac_reinf_oldvoc[c_voc] <- tmp_vac_reinf_oldvoc[c_voc] + tmp_vac_reinf_oldvoc[c_nonvoc]
          tmp_vac_reinf_oldvoc[c_nonvoc] <- 0
          tmp_vac_reinfvac[c_voc] <- tmp_vac_reinfvac[c_voc] + tmp_vac_reinfvac[c_nonvoc]
          tmp_vac_reinfvac[c_nonvoc] <- 0
          tmp_vac_reinfvac_oldvoc[c_voc] <- tmp_vac_reinfvac_oldvoc[c_voc] + tmp_vac_reinfvac_oldvoc[c_nonvoc]
          tmp_vac_reinfvac_oldvoc[c_nonvoc] <- 0
          tmp_vac_extrabooster1[c_voc] <- tmp_vac_extrabooster1[c_voc] + tmp_vac_extrabooster1[c_nonvoc]
          tmp_vac_extrabooster1[c_nonvoc] <- 0
          tmp_vac_extrabooster2[c_voc] <- tmp_vac_extrabooster2[c_voc] + tmp_vac_extrabooster2[c_nonvoc]
          tmp_vac_extrabooster2[c_nonvoc] <- 0
          
        }

        # introduce cases (done later for newvoc -> changed)
        tmp[c_E] <- round(VOC_init/10, 0)
        
        # update transmission rates/probabilies
        prob_vec_vacc <- VOC_prob_vec_vacc
        
      } else {
        transm_fctr_VOC <- exp(sel_VOC_param['log_transm']);   # additional VoC transmissibility
        
        # switch VE against infection for voc strain
        ve_voc_infection <- get_ve_infection(sel_VOC_param,ve_tag='ve')
        
        # #change immunity vs old_current voc
        if(calendar_time >= sim_date2day(as.Date(start_date_voc_immune_evasion))){
          ve_infection$reinf     <- ve_voc_infection$reinf 
          ve_voc_infection$reinf     <- ve_voc_infection$reinf 
          ve_infection$reinf_oldvoc     <- ve_voc_infection$reinf 
          ve_voc_infection$reinf_oldvoc     <- ve_voc_infection$reinf_oldvoc 
         # print(sel_VOC)
         # print(c( ve_infection$reinf,ve_voc_infection$reinf,ve_infection$reinf_oldvoc,ve_voc_infection$reinf_oldvoc))
         # print("transfert of all old voc in non current voc (2)")
          tmp[c_nonvoc] <- tmp[c_nonvoc] + tmp[c_voc]
          tmp[c_voc] <- 0
          tmp_vac_rna1[c_nonvoc] <- tmp_vac_rna1[c_nonvoc] + tmp_vac_rna1[c_voc]
          tmp_vac_rna1[c_voc] <- 0
          tmp_vac_rna2[c_nonvoc] <- tmp_vac_rna2[c_nonvoc] + tmp_vac_rna2[c_voc]
          tmp_vac_rna2[c_voc] <- 0
          tmp_vac_adeno1[c_nonvoc] <- tmp_vac_adeno1[c_nonvoc] + tmp_vac_adeno1[c_voc]
          tmp_vac_adeno1[c_voc] <- 0
          tmp_vac_adeno2[c_nonvoc] <- tmp_vac_adeno2[c_nonvoc] + tmp_vac_adeno2[c_voc]
          tmp_vac_adeno2[c_voc] <- 0
          tmp_vac_waning[c_nonvoc] <- tmp_vac_waning[c_nonvoc] + tmp_vac_waning[c_voc]
          tmp_vac_waning[c_voc] <- 0
          tmp_vac_booster[c_nonvoc] <- tmp_vac_booster[c_nonvoc] + tmp_vac_booster[c_voc]
          tmp_vac_booster[c_voc] <- 0
          tmp_vac_booster_waning[c_nonvoc] <- tmp_vac_booster_waning[c_nonvoc] + tmp_vac_booster_waning[c_voc]
          tmp_vac_booster_waning[c_voc] <- 0
          tmp_vac_reinf[c_nonvoc] <- tmp_vac_reinf[c_nonvoc] + tmp_vac_reinf[c_voc]
          tmp_vac_reinf[c_voc] <- 0
          tmp_vac_reinf_oldvoc[c_nonvoc] <- tmp_vac_reinf_oldvoc[c_nonvoc] + tmp_vac_reinf_oldvoc[c_voc]
          tmp_vac_reinf_oldvoc[c_voc] <- 0
          tmp_vac_reinfvac[c_nonvoc] <- tmp_vac_reinfvac[c_nonvoc] + tmp_vac_reinfvac[c_voc]
          tmp_vac_reinfvac[c_voc] <- 0
          tmp_vac_reinfvac_oldvoc[c_nonvoc] <- tmp_vac_reinfvac_oldvoc[c_nonvoc] + tmp_vac_reinfvac_oldvoc[c_voc]
          tmp_vac_reinfvac_oldvoc[c_voc] <- 0
          tmp_vac_extrabooster1[c_nonvoc] <- tmp_vac_extrabooster1[c_nonvoc] + tmp_vac_extrabooster1[c_voc]
          tmp_vac_extrabooster1[c_voc] <- 0
          tmp_vac_extrabooster2[c_nonvoc] <- tmp_vac_extrabooster2[c_nonvoc] + tmp_vac_extrabooster2[c_voc]
          tmp_vac_extrabooster2[c_voc] <- 0
        }
        
        # introduce cases (done later for newvoc -> changed)
        tmp[c_Evoc] = round(VOC_init/10, 0)
        
        # update transmission rates/probabilies
        prob_vec_vacc_voc <- VOC_prob_vec_vacc
      }
      
      
      #general update severe ve for dedicated vaccine
      if(calendar_time >= sim_date2day(as.Date(start_date_voc_immune_evasion))){
        prob_vec_vacc[13,grepl('I_sev_new',names(prob_vec))]   <- 1 - exp(-h*(1-get_phi0_vac(VOC_phi0,parms[grepl(paste0("ve_VOC_delta_incr_severe_extrabooster"),parms_names)]))*delta2)
        prob_vec_vacc[13,grepl('R_mild_new',names(prob_vec))]  <- 1 - exp(-h*   get_phi0_vac(VOC_phi0,parms[grepl(paste0("ve_VOC_delta_incr_severe_extrabooster"),parms_names)]) *delta2)
        prob_vec_vacc_voc[14,grepl('I_sev_new',names(prob_vec))]   <- 1 - exp(-h*(1-get_phi0_vac(VOC_phi0,parms[grepl(paste0("ve_VOC_delta_incr_severe_extrabooster"),parms_names)]))*delta2)
        prob_vec_vacc_voc[14,grepl('R_mild_new',names(prob_vec))]  <- 1 - exp(-h*   get_phi0_vac(VOC_phi0,parms[grepl(paste0("ve_VOC_delta_incr_severe_extrabooster"),parms_names)]) *delta2)
      }
      
    }
    
    ## Updated ICU ratio for fall 2021       ----
    ##----------------------------- -
    if(calendar_time == icu_fall2021_start){
      # update severity, with respect to probability to be admitted to ICU
      prob_vec_vacc[,grepl('I_hosp_new',names(prob_vec))]   <- matrix(rep(1 - exp(-h*  phi1_fall2021  *omega),nrow_vec_vacc),ncol=10,byrow = T)
      prob_vec_vacc[,grepl('I_icu_new',names(prob_vec))]    <- matrix(rep(1 - exp(-h*(1-phi1_fall2021)*omega),nrow_vec_vacc),ncol=10,byrow = T)
    }
    
    ## Updated ICU ratio for spring 2022       ----
    ##----------------------------- -
    if(calendar_time == icu_spring2022_start){
      # update severity, with respect to probability to be admitted to ICU
      prob_vec_vacc_voc[,grepl('I_hosp_new',names(prob_vec))]   <- matrix(rep(1 - exp(-h*  phi1_spring2022  *omega),nrow_vec_vacc),ncol=10,byrow = T)
      prob_vec_vacc_voc[,grepl('I_icu_new',names(prob_vec))]    <- matrix(rep(1 - exp(-h*(1-phi1_spring2022)*omega),nrow_vec_vacc),ncol=10,byrow = T)
    }
    
    ## Force of infection (see definition p*(k) in manuscript); reduction in transmissibility  ---- 
    ##------------------------------------------------------------------------------------------------------- -
    if(bool_voc | bool_vaccine){  # temp : previously infect with reduced infectiousness or not ???
      tmp_aggr <- tmp + tmp_vac_reinf + tmp_vac_reinf_oldvoc + (tmp_vac_rna1 + tmp_vac_rna2 + tmp_vac_adeno1 + tmp_vac_adeno2 + tmp_vac_waning + tmp_vac_booster + tmp_vac_booster_waning + tmp_vac_reinfvac + tmp_vac_reinfvac_oldvoc + tmp_vac_extrabooster1+ tmp_vac_extrabooster2)* (1-ve_transmission)
      
      
      lambda_t = transm_fctr*(betas_asy %*% (tmp_aggr[c_I_presym] + tmp_aggr[c_I_asym]) +
                                betas_sy  %*% (tmp_aggr[c_I_mild] + tmp_aggr[c_I_sev]));
      
      lambda_t_voc = transm_fctr_VOC*(betas_asy %*% (tmp_aggr[c_Ivoc_presym] + tmp_aggr[c_Ivoc_asym]) + 
                                        betas_sy  %*% (tmp_aggr[c_Ivoc_mild] + tmp_aggr[c_Ivoc_sev]));      
    } else{
      tmp_aggr <- tmp + tmp_vac_reinf + tmp_vac_reinf_oldvoc
      
      lambda_t = betas_asy %*% (tmp_aggr[c_I_presym] + tmp_aggr[c_I_asym]) +
        betas_sy  %*% (tmp_aggr[c_I_mild] + tmp_aggr[c_I_sev]);
      
      lambda_t_voc <- lambda_t*0
    }
    #print(lambda_t)
    #print(lambda_t_voc)
    # perform a recurent operation once
    neg_h_lambda_t     <- -h*lambda_t
    neg_h_lambda_t_voc <- -h*lambda_t_voc
    
    ## Chain binomial sequences  ---- 
    ##----------------------------- -
    # infections
    new_prob_vec    = c(1 - exp(neg_h_lambda_t), prob_vec_vacc[1,],
                        1 - exp(neg_h_lambda_t_voc), prob_vec_vacc_voc[1,])
    new_transitions[-1] = get_new_transitions(tmp[col_vec],new_prob_vec)
    
    # reinfections
    new_prob_reinf = c(1 - exp((1-ve_infection$reinf)*neg_h_lambda_t),     prob_vec_vacc[9,],      
                       1 - exp((1-ve_voc_infection$reinf)*neg_h_lambda_t_voc), prob_vec_vacc_voc[9,])
    new_transitions_vac_reinf[-1] = get_new_transitions(tmp_vac_reinf[col_vec],new_prob_reinf)
    new_prob_reinf_oldvoc = c(1 - exp((1-ve_infection$reinf_oldvoc)*neg_h_lambda_t),     prob_vec_vacc[10,],      
                              1 - exp((1-ve_voc_infection$reinf_oldvoc)*neg_h_lambda_t_voc), prob_vec_vacc_voc[10,])
    new_transitions_vac_reinf_oldvoc[-1] = get_new_transitions(tmp_vac_reinf_oldvoc[col_vec],new_prob_reinf_oldvoc)
    
    if(bool_vaccine){
      new_prob_rna1 = c(1 - exp((1-ve_infection$rna1)*neg_h_lambda_t),      prob_vec_vacc[2,],      # mRNA dose 1
                        1 - exp((1-ve_voc_infection$rna1)*neg_h_lambda_t_voc),  prob_vec_vacc_voc[2,])
      new_prob_rna2 = c(1 - exp((1-ve_infection$rna2)*neg_h_lambda_t),      prob_vec_vacc[3,],      # mRNA dose 2
                        1 - exp((1-ve_voc_infection$rna2)*neg_h_lambda_t_voc),  prob_vec_vacc_voc[3,])
      new_prob_adeno1 = c(1 - exp((1-ve_infection$adeno1)*neg_h_lambda_t),     prob_vec_vacc[4,],      # adeno dose 1
                          1 - exp((1-ve_voc_infection$adeno1)*neg_h_lambda_t_voc), prob_vec_vacc_voc[4,])
      new_prob_adeno2 = c(1 - exp((1-ve_infection$adeno2)*neg_h_lambda_t),     prob_vec_vacc[5,],      # adeno dose 2
                          1 - exp((1-ve_voc_infection$adeno2)*neg_h_lambda_t_voc), prob_vec_vacc_voc[5,])
      new_prob_waning = c(1 - exp((1-ve_infection$waning)*neg_h_lambda_t),     prob_vec_vacc[6,],      # waning immunity
                          1 - exp((1-ve_voc_infection$waning)*neg_h_lambda_t_voc), prob_vec_vacc_voc[6,])
      new_prob_booster= c(1 - exp((1-ve_infection$booster)*neg_h_lambda_t),     prob_vec_vacc[7,],      # booster dose
                          1 - exp((1-ve_voc_infection$booster)*neg_h_lambda_t_voc), prob_vec_vacc_voc[7,])
      new_prob_booster_waning = c(1 - exp((1-ve_infection$booster_waning)*neg_h_lambda_t),     prob_vec_vacc[8,],      # booster dose with waning immunity
                                  1 - exp((1-ve_voc_infection$booster_waning)*neg_h_lambda_t_voc), prob_vec_vacc_voc[8,])
      new_prob_reinfvac = c(1 - exp((1-ve_infection$reinf)*neg_h_lambda_t),     prob_vec_vacc[11,],      # reinfection paths after vaccination
                            1 - exp((1-ve_voc_infection$reinf)*neg_h_lambda_t_voc), prob_vec_vacc_voc[11,]) #note: using reinf instead of reinfvac since identical for the moment
      new_prob_reinfvac_oldvoc = c(1 - exp((1-ve_infection$reinf_oldvoc)*neg_h_lambda_t),     prob_vec_vacc[12,],      # reinfection paths after vaccination
                                   1 - exp((1-ve_voc_infection$reinf_oldvoc)*neg_h_lambda_t_voc), prob_vec_vacc_voc[12,])
      new_prob_extrabooster1= c(1 - exp((1-ve_infection$extrabooster)*neg_h_lambda_t),     prob_vec_vacc[13,],      # extrabooster dose
                                1 - exp((1-ve_voc_infection$extrabooster)*neg_h_lambda_t_voc), prob_vec_vacc_voc[13,])
      new_prob_extrabooster2= c(1 - exp((1-ve_infection$extrabooster)*neg_h_lambda_t),     prob_vec_vacc[14,],      # extrabooster dose to adapt
                                1 - exp((1-ve_voc_infection$extrabooster)*neg_h_lambda_t_voc), prob_vec_vacc_voc[14,])
      if(calendar_time > sim_date2day(as.Date(start_date_voc_immune_evasion))){
        new_prob_extrabooster1= c(1 - exp((1-parms['ve_VOC_ba4ba5_infection_extrabooster'])*neg_h_lambda_t),     prob_vec_vacc[13,],      # extrabooster dose
                                  1 - exp((1-ve_voc_infection$extrabooster)*neg_h_lambda_t_voc), prob_vec_vacc_voc[13,])
        new_prob_extrabooster2= c(1 - exp((1-ve_infection$extrabooster)*neg_h_lambda_t),     prob_vec_vacc[14,],      # extrabooster dose to adapt
                                  1 - exp((1-parms['ve_VOC_ba4ba5_infection_extrabooster'])*neg_h_lambda_t_voc), prob_vec_vacc_voc[14,])
      }
      
      
      new_transitions_vac_rna1[-1]   = get_new_transitions(tmp_vac_rna1[col_vec],new_prob_rna1) 
      new_transitions_vac_rna2[-1]   = get_new_transitions(tmp_vac_rna2[col_vec],new_prob_rna2) 
      new_transitions_vac_adeno1[-1] = get_new_transitions(tmp_vac_adeno1[col_vec],new_prob_adeno1)
      new_transitions_vac_adeno2[-1] = get_new_transitions(tmp_vac_adeno2[col_vec],new_prob_adeno2)
      new_transitions_vac_waning[-1] = get_new_transitions(tmp_vac_waning[col_vec],new_prob_waning)
      new_transitions_vac_booster[-1] = get_new_transitions(tmp_vac_booster[col_vec],new_prob_booster)
      new_transitions_vac_booster_waning[-1] = get_new_transitions(tmp_vac_booster_waning[col_vec],new_prob_booster_waning)
      new_transitions_vac_reinfvac[-1] = get_new_transitions(tmp_vac_reinfvac[col_vec],new_prob_reinfvac)
      new_transitions_vac_reinfvac_oldvoc[-1] = get_new_transitions(tmp_vac_reinfvac_oldvoc[col_vec],new_prob_reinfvac_oldvoc)
      new_transitions_vac_extrabooster1[-1] = get_new_transitions(tmp_vac_extrabooster1[col_vec],new_prob_extrabooster1)
      new_transitions_vac_extrabooster2[-1] = get_new_transitions(tmp_vac_extrabooster2[col_vec],new_prob_extrabooster2)
    }  
    
    
    ### 1. Flows between disease states  ---- 
    ###------------------------------------ -
    tmp = (tmp * icol_flag_ode - 
             new_transitions[min_transition_vec1] - 
             new_transitions[min_transition_vec2] + 
             new_transitions[plus_transition_vec1] + 
             new_transitions[plus_transition_vec2] + 
             new_transitions[plus_transition_vec3]);
    
    tmp_vac_reinf = (tmp_vac_reinf  * icol_flag_ode -
                       new_transitions_vac_reinf[min_transition_vec1] -
                       new_transitions_vac_reinf[min_transition_vec2] +
                       new_transitions_vac_reinf[plus_transition_vec1] +
                       new_transitions_vac_reinf[plus_transition_vec2] +
                       new_transitions_vac_reinf[plus_transition_vec3]);
    
    tmp_vac_reinf_oldvoc = (tmp_vac_reinf_oldvoc  * icol_flag_ode -
                              new_transitions_vac_reinf_oldvoc[min_transition_vec1] -
                              new_transitions_vac_reinf_oldvoc[min_transition_vec2] +
                              new_transitions_vac_reinf_oldvoc[plus_transition_vec1] +
                              new_transitions_vac_reinf_oldvoc[plus_transition_vec2] +
                              new_transitions_vac_reinf_oldvoc[plus_transition_vec3]);
    
    if(bool_vaccine){
      tmp_vac_rna1 = (tmp_vac_rna1 * icol_flag_ode - 
                        new_transitions_vac_rna1[min_transition_vec1] - 
                        new_transitions_vac_rna1[min_transition_vec2] + 
                        new_transitions_vac_rna1[plus_transition_vec1] + 
                        new_transitions_vac_rna1[plus_transition_vec2] + 
                        new_transitions_vac_rna1[plus_transition_vec3]);
      
      tmp_vac_rna2 = (tmp_vac_rna2 * icol_flag_ode - 
                        new_transitions_vac_rna2[min_transition_vec1] - 
                        new_transitions_vac_rna2[min_transition_vec2] + 
                        new_transitions_vac_rna2[plus_transition_vec1] + 
                        new_transitions_vac_rna2[plus_transition_vec2] + 
                        new_transitions_vac_rna2[plus_transition_vec3]);
      
      tmp_vac_adeno1 = (tmp_vac_adeno1  * icol_flag_ode - 
                          new_transitions_vac_adeno1[min_transition_vec1] - 
                          new_transitions_vac_adeno1[min_transition_vec2] + 
                          new_transitions_vac_adeno1[plus_transition_vec1] + 
                          new_transitions_vac_adeno1[plus_transition_vec2] + 
                          new_transitions_vac_adeno1[plus_transition_vec3]);
      
      tmp_vac_adeno2 = (tmp_vac_adeno2  * icol_flag_ode - 
                          new_transitions_vac_adeno2[min_transition_vec1] - 
                          new_transitions_vac_adeno2[min_transition_vec2] + 
                          new_transitions_vac_adeno2[plus_transition_vec1] + 
                          new_transitions_vac_adeno2[plus_transition_vec2] + 
                          new_transitions_vac_adeno2[plus_transition_vec3]);
      
      tmp_vac_waning = (tmp_vac_waning  * icol_flag_ode - 
                          new_transitions_vac_waning[min_transition_vec1] - 
                          new_transitions_vac_waning[min_transition_vec2] + 
                          new_transitions_vac_waning[plus_transition_vec1] + 
                          new_transitions_vac_waning[plus_transition_vec2] + 
                          new_transitions_vac_waning[plus_transition_vec3]);
      
      tmp_vac_booster = (tmp_vac_booster  * icol_flag_ode - 
                           new_transitions_vac_booster[min_transition_vec1] - 
                           new_transitions_vac_booster[min_transition_vec2] + 
                           new_transitions_vac_booster[plus_transition_vec1] + 
                           new_transitions_vac_booster[plus_transition_vec2] + 
                           new_transitions_vac_booster[plus_transition_vec3]);
      
      tmp_vac_booster_waning = (tmp_vac_booster_waning  * icol_flag_ode -
                                  new_transitions_vac_booster_waning[min_transition_vec1] -
                                  new_transitions_vac_booster_waning[min_transition_vec2] +
                                  new_transitions_vac_booster_waning[plus_transition_vec1] +
                                  new_transitions_vac_booster_waning[plus_transition_vec2] +
                                  new_transitions_vac_booster_waning[plus_transition_vec3]);
      
      
      tmp_vac_reinfvac = (tmp_vac_reinfvac  * icol_flag_ode -
                            new_transitions_vac_reinfvac[min_transition_vec1] -
                            new_transitions_vac_reinfvac[min_transition_vec2] +
                            new_transitions_vac_reinfvac[plus_transition_vec1] +
                            new_transitions_vac_reinfvac[plus_transition_vec2] +
                            new_transitions_vac_reinfvac[plus_transition_vec3]);
      
      tmp_vac_reinfvac_oldvoc = (tmp_vac_reinfvac_oldvoc  * icol_flag_ode -
                                   new_transitions_vac_reinfvac_oldvoc[min_transition_vec1] -
                                   new_transitions_vac_reinfvac_oldvoc[min_transition_vec2] +
                                   new_transitions_vac_reinfvac_oldvoc[plus_transition_vec1] +
                                   new_transitions_vac_reinfvac_oldvoc[plus_transition_vec2] +
                                   new_transitions_vac_reinfvac_oldvoc[plus_transition_vec3]);
      
      tmp_vac_extrabooster1 = (tmp_vac_extrabooster1  * icol_flag_ode - 
                                 new_transitions_vac_extrabooster1[min_transition_vec1] - 
                                 new_transitions_vac_extrabooster1[min_transition_vec2] + 
                                 new_transitions_vac_extrabooster1[plus_transition_vec1] + 
                                 new_transitions_vac_extrabooster1[plus_transition_vec2] + 
                                 new_transitions_vac_extrabooster1[plus_transition_vec3]);
      
      tmp_vac_extrabooster2 = (tmp_vac_extrabooster2  * icol_flag_ode - 
                                 new_transitions_vac_extrabooster2[min_transition_vec1] - 
                                 new_transitions_vac_extrabooster2[min_transition_vec2] + 
                                 new_transitions_vac_extrabooster2[plus_transition_vec1] + 
                                 new_transitions_vac_extrabooster2[plus_transition_vec2] + 
                                 new_transitions_vac_extrabooster2[plus_transition_vec3]);
      
      tmp_vac_rna1         <- pmax(tmp_vac_rna1,0)
      tmp_vac_rna2         <- pmax(tmp_vac_rna2,0)
      tmp_vac_adeno1       <- pmax(tmp_vac_adeno1,0)
      tmp_vac_adeno2       <- pmax(tmp_vac_adeno2,0)
      tmp_vac_waning       <- pmax(tmp_vac_waning,0)
      tmp_vac_booster      <- pmax(tmp_vac_booster,0)
      tmp_vac_booster_waning <- pmax(tmp_vac_booster_waning,0)
      tmp_vac_reinfvac     <- pmax(tmp_vac_reinfvac,0)
      tmp_vac_reinfvac_oldvoc     <- pmax(tmp_vac_reinfvac_oldvoc,0)
      tmp_vac_extrabooster1      <- pmax(tmp_vac_extrabooster1,0)
      tmp_vac_extrabooster2      <- pmax(tmp_vac_extrabooster2,0)
      
    }
    
    tmp           <- pmax(tmp,0)
    tmp_vac_reinf <- pmax(tmp_vac_reinf,0)
    tmp_vac_reinf_oldvoc <- pmax(tmp_vac_reinf_oldvoc,0)
    
    # count re/infections (before updating S compartments)    
    nsubset_infection[[step_index + 1]] <- nfull[[step_index]][c_S] - tmp[c_S]
    nsubset_reinfection[[step_index + 1]] <- nfull_vac_reinf[[step_index]][c_S] - tmp_vac_reinf[c_S] + nfull_vac_reinf_oldvoc[[step_index]][c_S] - tmp_vac_reinf_oldvoc[c_S]
    
    ### newvoc as part of current infections (introduction) ---- check
     if(calendar_time >= sim_date2day(as.Date(start_date_new_hypothetic_voc))){
       if(global_lib_voc$model_strain[global_lib_voc$name == sel_VOC] == 1){
         currentinf <- tmp[c_E]+tmp[c_Evoc]
         if(sum(tmp[c_E]) < sum(currentinf)*0.001){
           tmp[c_E] <- pmax(currentinf*0.001, 1)  #new
          # tmp[c_Evoc] <- pmax(currentinf - tmp[c_E],0)
         }
       } else {
         currentinf <- tmp[c_E]+tmp[c_Evoc]
         if(sum(tmp[c_Evoc]) < sum(currentinf)*0.001){
           tmp[c_Evoc] = pmax(currentinf*0.001, 1)  #new
           #tmp[c_E] <- pmax(currentinf - tmp[c_Evoc],0)
         }
       }
     }
    
    
    ### 2. Vaccination events: uptake and waning immunity  ---- 
    ###-------------------------- -
    # always true once there is any vaccine uptake to account for waning immunity
    if(bool_vaccine | bool_n_uptake[step_index]){
      
      # set vaccine boolean to TRUE
      bool_vaccine <- TRUE
      
      # count of re/infection (before updating S compartments)
      nsubset_infection[[step_index + 1]] <- nsubset_infection[[step_index + 1]] +
        nfull_vac_rna1[[step_index]][c_S] - tmp_vac_rna1[c_S] +
        nfull_vac_rna2[[step_index]][c_S] - tmp_vac_rna2[c_S] +
        nfull_vac_adeno1[[step_index]][c_S] - tmp_vac_adeno1[c_S] +
        nfull_vac_adeno2[[step_index]][c_S] - tmp_vac_adeno2[c_S] +
        nfull_vac_waning[[step_index]][c_S] - tmp_vac_waning[c_S] +
        nfull_vac_booster[[step_index]][c_S] - tmp_vac_booster[c_S] +
        nfull_vac_booster_waning [[step_index]][c_S] - tmp_vac_booster_waning [c_S];
      nsubset_reinfection[[step_index + 1]] <- nsubset_reinfection[[step_index + 1]] +
        nfull_vac_reinfvac[[step_index]][c_S] - tmp_vac_reinfvac[c_S] +
        nfull_vac_reinfvac_oldvoc[[step_index]][c_S] - tmp_vac_reinfvac_oldvoc[c_S] +
        nfull_vac_extrabooster1[[step_index]][c_S] - tmp_vac_extrabooster1[c_S] +
        nfull_vac_extrabooster2[[step_index]][c_S] - tmp_vac_extrabooster2[c_S];
      
      coef_vac <- 1
      # if(step_index%%4==0){coef_vac <- 4} else {coef_vac <- 0}  #if vaccine per day
      total_vac_adeno1 = V_mat[step_index, 1:10]*coef_vac;
      total_vac_rna1   = V_mat[step_index, 11:20]*coef_vac;
      total_vac_adeno2 = V_mat_dose2[step_index, 1:10]*coef_vac;
      total_vac_rna2   = V_mat_dose2[step_index, 11:20]*coef_vac;
      total_vac_booster = V_mat_booster[step_index, 1:10]*coef_vac;
      total_vac_2ndbooster = V_mat_2ndbooster[step_index, 1:10]*coef_vac;
      total_vac_extrabooster = V_mat_extrabooster[step_index, 1:10]*coef_vac;
      
      # create temporary reduced vectors
      tmp_c_vac            <- tmp[c_vac]
      tmp_vac_adeno1_c_vac <- tmp_vac_adeno1[c_vac]
      tmp_vac_rna1_c_vac   <- tmp_vac_rna1[c_vac]
      tmp_vac_adeno2_c_vac <- tmp_vac_adeno2[c_vac]
      tmp_vac_rna2_c_vac   <- tmp_vac_rna2[c_vac]
      tmp_vac_waning_c_vac <- tmp_vac_waning[c_vac]
      tmp_vac_booster_c_vac<- tmp_vac_booster[c_vac]
      tmp_vac_booster_waning_c_vac<- tmp_vac_booster_waning[c_vac]
      tmp_vac_reinf_c_vac<- tmp_vac_reinf[c_vac]
      tmp_vac_reinf_oldvoc_c_vac<- tmp_vac_reinf_oldvoc[c_vac]
      tmp_vac_reinfvac_c_vac<- tmp_vac_reinfvac[c_vac]
      tmp_vac_reinfvac_oldvoc_c_vac<- tmp_vac_reinfvac_oldvoc[c_vac]
      tmp_vac_extrabooster1_c_vac<- tmp_vac_extrabooster1[c_vac]
      tmp_vac_extrabooster2_c_vac<- tmp_vac_extrabooster2[c_vac]
      
      # total_alive_nvac       = rowSums(matrix(tmp_c_vac,nrow=10))
      # total_alive_vac_adeno1 = rowSums(matrix(tmp_vac_adeno1_c_vac,nrow=10))
      # total_alive_vac_rna1   = rowSums(matrix(tmp_vac_rna1_c_vac,nrow=10))
      # 
      # # calculate vaccination factor   #to check !!!!
      # nvac must account for individuals both in tmp_c_vac and tmp_vac_reinf_c_vac
      factor_alive_nvac        = 1/(rowSums(matrix(tmp_c_vac,nrow=10))+rowSums(matrix(tmp_vac_reinf_c_vac,nrow=10))+rowSums(matrix(tmp_vac_reinf_oldvoc_c_vac,nrow=10)))
      factor_alive_vac_adeno1  = 1/rowSums(matrix(tmp_vac_adeno1_c_vac,nrow=10))
      factor_alive_vac_rna1    = 1/rowSums(matrix(tmp_vac_rna1_c_vac,nrow=10))
      factor_alive_vac_waning  = 1/rowSums(matrix(tmp_vac_waning_c_vac,nrow=10)+matrix(tmp_vac_adeno2_c_vac,nrow=10)+matrix(tmp_vac_rna2_c_vac,nrow=10))
      factor_alive_vac_for2ndbooster  = 1/rowSums(matrix(tmp_vac_booster_waning_c_vac,nrow=10)+matrix(tmp_vac_reinfvac_c_vac,nrow=10)+matrix(tmp_vac_reinfvac_oldvoc_c_vac,nrow=10))
      factor_alive_vac_forextrabooster  = 1/rowSums(matrix(tmp_vac_booster_waning_c_vac,nrow=10)+matrix(tmp_vac_reinfvac_c_vac,nrow=10)+matrix(tmp_vac_reinfvac_oldvoc_c_vac,nrow=10))
     #temporary adaptation for hubR5 to allow extrabooster to nonvaccinated people
      if(calendar_time >= sim_date2day(as.Date('2023-10-01'))){
        #extraboostertoall <- 1
      } 
      if(exists("extraboostertoall")){
        print("Warning : extraboostertoall",fill=T)
        factor_alive_vac_forextrabooster  = 1/rowSums(matrix(tmp_c_vac,nrow=10)+matrix(tmp_vac_reinf_c_vac,nrow=10)+matrix(tmp_vac_reinf_oldvoc_c_vac,nrow=10)+matrix(tmp_vac_booster_waning_c_vac,nrow=10)+matrix(tmp_vac_reinfvac_c_vac,nrow=10)+matrix(tmp_vac_reinfvac_oldvoc_c_vac,nrow=10))
      }
      
      # calculate vaccination factor
      factor_alive_nvac[factor_alive_nvac==Inf] <- 0    
      factor_alive_vac_adeno1[factor_alive_vac_adeno1==Inf] <- 0
      factor_alive_vac_rna1[factor_alive_vac_rna1==Inf] <- 0
      factor_alive_vac_waning[factor_alive_vac_waning==Inf] <- 0
      factor_alive_vac_for2ndbooster[factor_alive_vac_for2ndbooster==Inf] <- 0
      factor_alive_vac_forextrabooster[factor_alive_vac_forextrabooster==Inf] <- 0
      
      #### mRNA-based uptake  ----  
      ####---------------------------------- -
      
      new_vac_rna1 <- floor(tmp_c_vac * ((total_vac_rna1+total_vac_rna1_remain) * factor_alive_nvac))
      new_vac_rna1_reinf  <- floor(tmp_vac_reinf_c_vac * ((total_vac_rna1+total_vac_rna1_remain) * factor_alive_nvac))
      new_vac_rna1_reinf_oldvoc  <- floor(tmp_vac_reinf_oldvoc_c_vac * ((total_vac_rna1+total_vac_rna1_remain) * factor_alive_nvac))
      total_vac_rna1_remain =  total_vac_rna1_remain + total_vac_rna1 - (rowSums(matrix(new_vac_rna1,nrow=10))+rowSums(matrix(new_vac_rna1_reinf,nrow=10))+rowSums(matrix(new_vac_rna1_reinf_oldvoc,nrow=10)))
      new_vac_rna2  <- floor(tmp_vac_rna1_c_vac * ((total_vac_rna2+total_vac_rna2_remain) * factor_alive_vac_rna1))
      total_vac_rna2_remain =  total_vac_rna2_remain + total_vac_rna2 - rowSums(matrix(new_vac_rna2,nrow=10))
      
      #### Adeno-based uptake   ----  
      ####----------------------------------- -
      new_vac_adeno1 <- floor(tmp_c_vac * ((total_vac_adeno1+total_vac_adeno1_remain) * factor_alive_nvac))
      new_vac_adeno1_reinf  <- floor(tmp_vac_reinf_c_vac * ((total_vac_adeno1+total_vac_adeno1_remain) * factor_alive_nvac))
      new_vac_adeno1_reinf_oldvoc  <- floor(tmp_vac_reinf_oldvoc_c_vac * ((total_vac_adeno1+total_vac_adeno1_remain) * factor_alive_nvac))
      total_vac_adeno1_remain =  total_vac_adeno1_remain + total_vac_adeno1 - (rowSums(matrix(new_vac_adeno1,nrow=10))+rowSums(matrix(new_vac_adeno1_reinf,nrow=10))+rowSums(matrix(new_vac_adeno1_reinf_oldvoc,nrow=10)))
      new_vac_adeno2  <- floor(tmp_vac_adeno1_c_vac * ((total_vac_adeno2+total_vac_adeno2_remain) * factor_alive_vac_adeno1))
      total_vac_adeno2_remain =  total_vac_adeno2_remain + total_vac_adeno2 - rowSums(matrix(new_vac_adeno2,nrow=10))
      
      #### booster uptake   ----  
      ####----------------------------------- -
      new_vac_booster_w    <- floor(tmp_vac_waning_c_vac * ((total_vac_booster+total_vac_booster_remain) * factor_alive_vac_waning))
      new_vac_booster_ad    <- floor(tmp_vac_adeno2_c_vac * ((total_vac_booster+total_vac_booster_remain) * factor_alive_vac_waning))
      new_vac_booster_rna    <- floor(tmp_vac_rna2_c_vac * ((total_vac_booster+total_vac_booster_remain) * factor_alive_vac_waning))
      total_vac_booster_remain =  total_vac_booster_remain + total_vac_booster - rowSums(matrix(new_vac_booster_w+new_vac_booster_ad+new_vac_booster_rna,nrow=10))
      
      
      #### Waning immunity and booster dose uptake  ----  using get_new_transitions function.
      ####----------------------------------- -
      new_waning_rna2    <- get_new_transitions(tmp_vac_rna2_c_vac,ve_waning_immunity_rate)
      new_waning_adeno2  <- get_new_transitions(tmp_vac_adeno2_c_vac,ve_waning_immunity_rate)
      new_waning_booster <- get_new_transitions(tmp_vac_booster_c_vac,ve_waning_booster_rate)  
      new_waning_extrabooster1 <- get_new_transitions(tmp_vac_extrabooster1_c_vac,ve_waning_booster_rate)
      new_waning_extrabooster2 <- get_new_transitions(tmp_vac_extrabooster2_c_vac,ve_waning_booster_rate)
      
      #### 2ndbooster uptake   ----  
      ####----------------------------------- -
      new_vac_2ndbooster_w    <- floor(tmp_vac_booster_waning_c_vac * ((total_vac_2ndbooster+total_vac_2ndbooster_remain) * factor_alive_vac_for2ndbooster))
      new_vac_2ndbooster_reinf    <- floor(tmp_vac_reinfvac_c_vac * ((total_vac_2ndbooster+total_vac_2ndbooster_remain) * factor_alive_vac_for2ndbooster))
      new_vac_2ndbooster_reinf_oldvoc    <- floor(tmp_vac_reinfvac_oldvoc_c_vac * ((total_vac_2ndbooster+total_vac_2ndbooster_remain) * factor_alive_vac_for2ndbooster))
      
      # buffer
      tmp_vac_booster_waning[c_vac] <- tmp_vac_booster_waning_c_vac + new_waning_booster - new_vac_2ndbooster_w  #warning, this is only to check, not updated (done later with extrabooster)
      tmp_vac_reinfvac[c_vac] <-  tmp_vac_reinfvac_c_vac - new_vac_2ndbooster_reinf
      tmp_vac_reinfvac_oldvoc[c_vac] <-  tmp_vac_reinfvac_oldvoc_c_vac - new_vac_2ndbooster_reinf_oldvoc
      overvacc <- which(tmp_vac_booster_waning[c_vac]<0)
      new_vac_2ndbooster_w[overvacc] <- new_vac_2ndbooster_w[overvacc] + tmp_vac_booster_waning[c_vac][overvacc]
      tmp_vac_booster_waning[c_vac][overvacc]  <- 0
      overvacc2 <- which(tmp_vac_reinfvac[c_vac]<0)
      new_vac_2ndbooster_reinf[overvacc2] <- new_vac_2ndbooster_reinf[overvacc2] + tmp_vac_reinfvac[c_vac][overvacc2]
      tmp_vac_reinfvac[c_vac][overvacc2]  <- 0
      overvacc3 <- which(tmp_vac_reinfvac_oldvoc[c_vac]<0)
      new_vac_2ndbooster_reinf_oldvoc[overvacc3] <- new_vac_2ndbooster_reinf_oldvoc[overvacc3] + tmp_vac_reinfvac_oldvoc[c_vac][overvacc3]
      tmp_vac_reinfvac_oldvoc[c_vac][overvacc3]  <- 0
      
      total_vac_2ndbooster_remain =  total_vac_2ndbooster_remain + total_vac_2ndbooster - rowSums(matrix(new_vac_2ndbooster_w+new_vac_2ndbooster_reinf+new_vac_2ndbooster_reinf_oldvoc,nrow=10))
      # print("total2ndremain:")
      # print(total_vac_2ndbooster_remain)
      # 
    
      
      
      #### Update vaccine compartments   #to check
      tmp[c_vac]            <- tmp_c_vac - new_vac_rna1 - new_vac_adeno1
      tmp_vac_reinf[c_vac]  <- tmp_vac_reinf_c_vac - new_vac_rna1_reinf - new_vac_adeno1_reinf
      tmp_vac_reinf_oldvoc[c_vac]  <- tmp_vac_reinf_oldvoc_c_vac - new_vac_rna1_reinf_oldvoc - new_vac_adeno1_reinf_oldvoc
      tmp_vac_rna1[c_vac]   <- tmp_vac_rna1_c_vac + new_vac_rna1 + new_vac_rna1_reinf + new_vac_rna1_reinf_oldvoc      - new_vac_rna2
      tmp_vac_adeno1[c_vac] <- tmp_vac_adeno1_c_vac + new_vac_adeno1 + new_vac_adeno1_reinf + new_vac_adeno1_reinf_oldvoc - new_vac_adeno2
      
      tmp_vac_rna2[c_vac]    <- tmp_vac_rna2_c_vac + new_vac_rna2 - new_vac_booster_rna - new_waning_rna2 
      tmp_vac_adeno2[c_vac]  <- tmp_vac_adeno2_c_vac + new_vac_adeno2 - new_vac_booster_ad - new_waning_adeno2
      
      tmp_vac_waning[c_vac]  <- tmp_vac_waning_c_vac + new_waning_rna2 + new_waning_adeno2 - new_vac_booster_w
      tmp_vac_booster[c_vac] <- tmp_vac_booster_c_vac + new_vac_booster_w + new_vac_booster_ad + new_vac_booster_rna - new_waning_booster + new_vac_2ndbooster_w + new_vac_2ndbooster_reinf + new_vac_2ndbooster_reinf_oldvoc
      
      #### extrabooster uptake   ----    to check  - separated to correct overvaccination (with buffer)   to check !!!!
      ####----------------------------------- -
      new_vac_extrabooster_w    <- floor(tmp_vac_booster_waning_c_vac * ((total_vac_extrabooster+total_vac_extrabooster_remain) * factor_alive_vac_forextrabooster))
      new_vac_extrabooster_reinf    <- floor(tmp_vac_reinfvac_c_vac * ((total_vac_extrabooster+total_vac_extrabooster_remain) * factor_alive_vac_forextrabooster))
      new_vac_extrabooster_reinf_oldvoc    <- floor(tmp_vac_reinfvac_oldvoc_c_vac * ((total_vac_extrabooster+total_vac_extrabooster_remain) * factor_alive_vac_forextrabooster))
      
      
      tmp_vac_booster_waning[c_vac] <- tmp_vac_booster_waning_c_vac + new_waning_booster - new_vac_2ndbooster_w - new_vac_extrabooster_w
      tmp_vac_reinfvac[c_vac] <-  tmp_vac_reinfvac_c_vac - new_vac_2ndbooster_reinf - new_vac_extrabooster_reinf
      tmp_vac_reinfvac_oldvoc[c_vac] <-  tmp_vac_reinfvac_oldvoc_c_vac - new_vac_2ndbooster_reinf_oldvoc - new_vac_extrabooster_reinf_oldvoc
      overvacc <- which(tmp_vac_booster_waning[c_vac]<0)
      new_vac_extrabooster_w[overvacc] <- new_vac_extrabooster_w[overvacc] + tmp_vac_booster_waning[c_vac][overvacc]
      tmp_vac_booster_waning[c_vac][overvacc]  <- 0
      overvacc2 <- which(tmp_vac_reinfvac[c_vac]<0)
      new_vac_extrabooster_reinf[overvacc2] <- new_vac_extrabooster_reinf[overvacc2] + tmp_vac_reinfvac[c_vac][overvacc2]
      tmp_vac_reinfvac[c_vac][overvacc2]  <- 0
      overvacc3 <- which(tmp_vac_reinfvac_oldvoc[c_vac]<0)
      new_vac_extrabooster_reinf_oldvoc[overvacc3] <- new_vac_extrabooster_reinf_oldvoc[overvacc3] + tmp_vac_reinfvac_oldvoc[c_vac][overvacc3]
      tmp_vac_reinfvac_oldvoc[c_vac][overvacc3]  <- 0
      
      if(exists("extraboostertoall")){
        new_vac_extrabooster_1nvac    <- floor(tmp_c_vac * ((total_vac_extrabooster+total_vac_extrabooster_remain) * factor_alive_vac_forextrabooster))
        new_vac_extrabooster_2nvacreinf    <- floor(tmp_vac_reinf_c_vac * ((total_vac_extrabooster+total_vac_extrabooster_remain) * factor_alive_vac_forextrabooster))
        new_vac_extrabooster_3nvacreinfoldvoc    <- floor(tmp_vac_reinf_oldvoc_c_vac * ((total_vac_extrabooster+total_vac_extrabooster_remain) * factor_alive_vac_forextrabooster))
        tmp[c_vac]            <- tmp[c_vac]  - new_vac_extrabooster_1nvac
        overvacc <- which(tmp[c_vac]<0)
        new_vac_extrabooster_1nvac[overvacc] <- new_vac_extrabooster_1nvac[overvacc] + tmp[c_vac][overvacc]
        tmp[c_vac][overvacc]  <- 0
        tmp_vac_reinf[c_vac]  <- tmp_vac_reinf[c_vac] - new_vac_extrabooster_2nvacreinf
        overvacc <- which(tmp_vac_reinf[c_vac]<0)
        new_vac_extrabooster_2nvacreinf[overvacc] <- new_vac_extrabooster_2nvacreinf[overvacc] + tmp_vac_reinf[c_vac][overvacc]
        tmp_vac_reinf[c_vac][overvacc]  <- 0
        tmp_vac_reinf_oldvoc[c_vac]  <- tmp_vac_reinf_oldvoc[c_vac] - new_vac_extrabooster_3nvacreinfoldvoc
        overvacc <- which(tmp_vac_reinf_oldvoc[c_vac]<0)
        new_vac_extrabooster_3nvacreinfoldvoc[overvacc] <- new_vac_extrabooster_3nvacreinfoldvoc[overvacc] + tmp_vac_reinf_oldvoc[c_vac][overvacc]
        tmp_vac_reinf_oldvoc[c_vac][overvacc]  <- 0
        #put all in new_vac_extrabooster_w
        new_vac_extrabooster_w <- new_vac_extrabooster_w + new_vac_extrabooster_1nvac + new_vac_extrabooster_2nvacreinf + new_vac_extrabooster_3nvacreinfoldvoc
      }
      

      # to check if adapted vaccine : warning, adapted dedicated bivalent vaccine is after 6 months  -> temp, to adapt later
     # if(calendar_time >= sim_date2day(as.Date(start_date_voc_immune_evasion)) & min((calendar_time - VOC_start)[calendar_time - VOC_start > 0]) > 6*30){
      #always dedicated vaccine for extrabooster according to uptake
      #(note: test on start_date_voc_immune_evasion are rests from previous scenarios from hub, they should be no extrabooster before start_date_voc_immune_evasion in the reality)
      if(calendar_time >= sim_date2day(as.Date(start_date_voc_immune_evasion)) ){
          if(calendar_time >= sim_date2day(as.Date(start_date_voc_immune_evasion)) & global_lib_voc$model_strain[global_lib_voc$name == sel_VOC] == 1){
          tmp_vac_extrabooster1[c_vac] <- tmp_vac_extrabooster1_c_vac + new_vac_extrabooster_w + new_vac_extrabooster_reinf  + new_vac_extrabooster_reinf_oldvoc
        } else {
          tmp_vac_extrabooster2[c_vac] <- tmp_vac_extrabooster2_c_vac + new_vac_extrabooster_w + new_vac_extrabooster_reinf  + new_vac_extrabooster_reinf_oldvoc
        }
      } else {
        if(calendar_time < sim_date2day(as.Date(start_date_voc_immune_evasion)) | global_lib_voc$model_strain[global_lib_voc$name == sel_VOC] == 2){
          tmp_vac_extrabooster1[c_vac] <- tmp_vac_extrabooster1_c_vac + new_vac_extrabooster_w + new_vac_extrabooster_reinf  + new_vac_extrabooster_reinf_oldvoc
        } else {
          tmp_vac_extrabooster2[c_vac] <- tmp_vac_extrabooster2_c_vac + new_vac_extrabooster_w + new_vac_extrabooster_reinf  + new_vac_extrabooster_reinf_oldvoc
        }
      }
      
      total_vac_extrabooster_remain =  (total_vac_extrabooster_remain + total_vac_extrabooster - rowSums(matrix(new_vac_extrabooster_w+new_vac_extrabooster_reinf+new_vac_extrabooster_reinf_oldvoc,nrow=10)))
      
         # count only from extraboost
      nsubset_numvac[[step_index + 1]] = rowSums(matrix(new_vac_extrabooster_w+new_vac_extrabooster_reinf+new_vac_extrabooster_reinf_oldvoc,nrow=10))
                   
      #waning for extrabooster -> currently waning into tmp_vac_reinfvac_c_vac or tmp_vac_reinfvac_oldvoc_c_vac to keep track of the protection vs immune evasion for each variant
       tmp_vac_extrabooster1[c_vac] = tmp_vac_extrabooster1[c_vac] - new_waning_extrabooster1
       tmp_vac_extrabooster2[c_vac] = tmp_vac_extrabooster2[c_vac] - new_waning_extrabooster2
       if(global_lib_voc$model_strain[global_lib_voc$name == sel_VOC] == 1){
         tmp_vac_reinfvac[c_vac] = tmp_vac_reinfvac[c_vac] + new_waning_extrabooster1
         tmp_vac_reinfvac_oldvoc[c_vac] = tmp_vac_reinfvac_oldvoc[c_vac] + new_waning_extrabooster2
         } else {
           tmp_vac_reinfvac[c_vac] = tmp_vac_reinfvac[c_vac] + new_waning_extrabooster2
           tmp_vac_reinfvac_oldvoc[c_vac] = tmp_vac_reinfvac_oldvoc[c_vac] + new_waning_extrabooster1
           }
                                   
      if(any(tmp_vac_booster_waning[c_vac]<0)){
        browser()
      }
      
      #check negative doses due to vaccination, and correct for overflow if needed
      if(any(tmp[c_vac]<0)){
        flag_overfow <- tmp[c_vac] < 0
        tmp_vac_rna2[c_vac[flag_overfow]] <- tmp_vac_rna2[c_vac[flag_overfow]] + tmp[c_vac[flag_overfow]]
        tmp[c_vac[flag_overfow]] <- 0
      }
      
      if(any(tmp_vac_waning[c_vac]<0)) {
        flag_overfow <- tmp_vac_waning[c_vac] < 0
        tmp_vac_adeno2[c_vac[flag_overfow]] <- tmp_vac_adeno2[c_vac[flag_overfow]] +  tmp_vac_waning[c_vac[flag_overfow]] 
        tmp_vac_waning[c_vac[flag_overfow]] <- 0
      }
      
      if(any(tmp_vac_adeno2[c_vac]<0)){
        flag_overfow <- tmp_vac_adeno2[c_vac] < 0
        tmp_vac_rna2[c_vac[flag_overfow]] <- tmp_vac_rna2[c_vac[flag_overfow]] + tmp_vac_adeno2[c_vac[flag_overfow]]
        tmp_vac_adeno2[c_vac[flag_overfow]] <- 0
      }
      
      if(any(tmp_vac_rna2[c_vac]<0)){
        flag_overfow <- tmp_vac_rna2[c_vac] < 0
        tmp_vac_booster[c_vac[flag_overfow]] <- tmp_vac_booster[c_vac[flag_overfow]] + tmp_vac_rna2[c_vac[flag_overfow]]
        tmp_vac_rna2[c_vac[flag_overfow]] <- 0
      }
      
      ### Re-infections ----
      # waning into reinfection path - computed after vac waning to avoid negativity  
      new_waning_reinf_adeno1    <- (tmp_vac_adeno1[c(c_R,c_Rvoc)] * 0)#ve_waning_infection_booster_rate)  
      new_waning_reinf_rna1    <- (tmp_vac_rna1[c(c_R,c_Rvoc)] * 0)#ve_waning_infection_booster_rate)  
      new_waning_reinf_adeno2    <- get_new_transitions(tmp_vac_adeno2[c(c_R,c_Rvoc)],ve_waning_infection_booster_rate)  
      new_waning_reinf_rna2    <- get_new_transitions(tmp_vac_rna2[c(c_R,c_Rvoc)],ve_waning_infection_booster_rate)  
      new_waning_reinf_waning    <- get_new_transitions(tmp_vac_waning[c(c_R,c_Rvoc)],ve_waning_infection_booster_rate)  
      new_waning_reinf_booster    <- get_new_transitions(tmp_vac_booster[c(c_R,c_Rvoc)],ve_waning_infection_booster_rate)  
      new_waning_reinf_booster_waning    <-get_new_transitions(tmp_vac_booster_waning[c(c_R,c_Rvoc)],ve_waning_infection_booster_rate)  
      new_waning_reinf_extrabooster1    <- get_new_transitions(tmp_vac_extrabooster1[c(c_R,c_Rvoc)],ve_waning_infection_booster_rate)   #to check
      new_waning_reinf_extrabooster2    <- get_new_transitions(tmp_vac_extrabooster2[c(c_R,c_Rvoc)],ve_waning_infection_booster_rate)   #to check
      
      # loop inside reinfection path   
      new_waning_reinf_suppreinfvac    <- get_new_transitions(tmp_vac_reinfvac[c(c_R,c_Rvoc)],ve_waning_infection_rate)  
      new_waning_reinf_suppreinfvac_oldvoc    <- get_new_transitions(tmp_vac_reinfvac_oldvoc[c(c_R,c_Rvoc)],ve_waning_infection_rate)  
      
      
      ##update reinfection   
      tmp_vac_adeno1[c(c_R,c_Rvoc)] <- tmp_vac_adeno1[c(c_R,c_Rvoc)] - new_waning_reinf_adeno1 
      tmp_vac_rna1[c(c_R,c_Rvoc)] <- tmp_vac_rna1[c(c_R,c_Rvoc)] - new_waning_reinf_rna1
      tmp_vac_adeno2[c(c_R,c_Rvoc)] <- tmp_vac_adeno2[c(c_R,c_Rvoc)] - new_waning_reinf_adeno2
      tmp_vac_rna2[c(c_R,c_Rvoc)] <- tmp_vac_rna2[c(c_R,c_Rvoc)] - new_waning_reinf_rna2
      tmp_vac_waning[c(c_R,c_Rvoc)] <- tmp_vac_waning[c(c_R,c_Rvoc)] - new_waning_reinf_waning
      tmp_vac_booster[c(c_R,c_Rvoc)] <- tmp_vac_booster[c(c_R,c_Rvoc)] - new_waning_reinf_booster
      tmp_vac_booster_waning[c(c_R,c_Rvoc)] <- tmp_vac_booster_waning[c(c_R,c_Rvoc)] - new_waning_reinf_booster_waning
      tmp_vac_extrabooster1[c(c_R,c_Rvoc)] <- tmp_vac_extrabooster1[c(c_R,c_Rvoc)] - new_waning_reinf_extrabooster1
      tmp_vac_extrabooster2[c(c_R,c_Rvoc)] <- tmp_vac_extrabooster2[c(c_R,c_Rvoc)] - new_waning_reinf_extrabooster2
      tmp_vac_reinfvac[c(c_R,c_Rvoc)] <- tmp_vac_reinfvac[c(c_R,c_Rvoc)] - new_waning_reinf_suppreinfvac
      tmp_vac_reinfvac_oldvoc[c(c_R,c_Rvoc)] <- tmp_vac_reinfvac_oldvoc[c(c_R,c_Rvoc)] - new_waning_reinf_suppreinfvac_oldvoc
      
      if(global_lib_voc$model_strain[global_lib_voc$name == sel_VOC] == 1){
        tmp_vac_reinfvac[c(c_S)] <- tmp_vac_reinfvac[c(c_S)] + new_waning_reinf_adeno1[1:10] +
          new_waning_reinf_rna1[1:10] +
          new_waning_reinf_adeno2[1:10] + 
          new_waning_reinf_rna2[1:10] + 
          new_waning_reinf_waning[1:10] +
          new_waning_reinf_booster[1:10]  +
          new_waning_reinf_booster_waning[1:10] +
          new_waning_reinf_extrabooster1[1:10]  + 
          new_waning_reinf_extrabooster2[1:10]  + 
          new_waning_reinf_suppreinfvac[1:10] +
          new_waning_reinf_suppreinfvac_oldvoc[1:10] 
        tmp_vac_reinfvac_oldvoc[c(c_S)] <- tmp_vac_reinfvac_oldvoc[c(c_S)] + new_waning_reinf_adeno1[11:20] +  #_oldvoc
          new_waning_reinf_rna1[11:20] +
          new_waning_reinf_adeno2[11:20] + 
          new_waning_reinf_rna2[11:20] + 
          new_waning_reinf_waning[11:20] +
          new_waning_reinf_booster[11:20]  +
          new_waning_reinf_booster_waning[11:20] +
          new_waning_reinf_extrabooster1[11:20]  + 
          new_waning_reinf_extrabooster2[11:20]  + 
          new_waning_reinf_suppreinfvac[11:20] +
          new_waning_reinf_suppreinfvac_oldvoc[11:20] 
      } else {
        tmp_vac_reinfvac_oldvoc[c(c_S)] <- tmp_vac_reinfvac_oldvoc[c(c_S)] + new_waning_reinf_adeno1[1:10] + #_oldvoc
          new_waning_reinf_rna1[1:10] +
          new_waning_reinf_adeno2[1:10] + 
          new_waning_reinf_rna2[1:10] + 
          new_waning_reinf_waning[1:10] +
          new_waning_reinf_booster[1:10]  +
          new_waning_reinf_booster_waning[1:10] +
          new_waning_reinf_extrabooster1[1:10]  + 
          new_waning_reinf_extrabooster2[1:10]  + 
          new_waning_reinf_suppreinfvac[1:10] +
          new_waning_reinf_suppreinfvac_oldvoc[1:10] 
        tmp_vac_reinfvac[c(c_S)] <- tmp_vac_reinfvac[c(c_S)] + new_waning_reinf_adeno1[11:20] +
          new_waning_reinf_rna1[11:20] +
          new_waning_reinf_adeno2[11:20] + 
          new_waning_reinf_rna2[11:20] + 
          new_waning_reinf_waning[11:20] +
          new_waning_reinf_booster[11:20]  +
          new_waning_reinf_booster_waning[11:20] +
          new_waning_reinf_extrabooster1[11:20]  + 
          new_waning_reinf_extrabooster2[11:20]  + 
          new_waning_reinf_suppreinfvac[11:20] +
          new_waning_reinf_suppreinfvac_oldvoc[11:20] 
      }
      
      
    }
    
    # Waning immunity after infection ----  
    new_waning_reinf    <- get_new_transitions(tmp[c(c_R,c_Rvoc)],ve_waning_infection_rate)  
    tmp[c(c_R,c_Rvoc)] <- tmp[c(c_R,c_Rvoc)] - new_waning_reinf 
    new_waning_reinf_suppreinf    <- get_new_transitions(tmp_vac_reinf[c(c_R,c_Rvoc)],ve_waning_infection_booster_rate)
    new_waning_reinf_suppreinf_oldvoc    <- get_new_transitions(tmp_vac_reinf_oldvoc[c(c_R,c_Rvoc)],ve_waning_infection_booster_rate)
    tmp_vac_reinf[c(c_R,c_Rvoc)] <- tmp_vac_reinf[c(c_R,c_Rvoc)] - new_waning_reinf_suppreinf
    tmp_vac_reinf_oldvoc[c(c_R,c_Rvoc)] <- tmp_vac_reinf_oldvoc[c(c_R,c_Rvoc)] - new_waning_reinf_suppreinf_oldvoc
    if(global_lib_voc$model_strain[global_lib_voc$name == sel_VOC] == 1){
      tmp_vac_reinf[c(c_S)] <- tmp_vac_reinf[c(c_S)] + new_waning_reinf[1:10]
      tmp_vac_reinf_oldvoc[c(c_S)] <- tmp_vac_reinf_oldvoc[c(c_S)] + new_waning_reinf[11:20] # #_oldvoc
      tmp_vac_reinf[c(c_S)] <- tmp_vac_reinf[c(c_S)] + new_waning_reinf_suppreinf[1:10]  + new_waning_reinf_suppreinf_oldvoc[1:10] 
      tmp_vac_reinf_oldvoc[c(c_S)] <- tmp_vac_reinf_oldvoc[c(c_S)] + new_waning_reinf_suppreinf[11:20] + new_waning_reinf_suppreinf_oldvoc[11:20]
    } else {
      tmp_vac_reinf[c(c_S)] <- tmp_vac_reinf[c(c_S)] + new_waning_reinf[11:20]
      tmp_vac_reinf_oldvoc[c(c_S)] <- tmp_vac_reinf_oldvoc[c(c_S)] + new_waning_reinf[1:10] # #_oldvoc
      tmp_vac_reinf[c(c_S)] <- tmp_vac_reinf[c(c_S)] + new_waning_reinf_suppreinf[11:20]  + new_waning_reinf_suppreinf_oldvoc[11:20] 
      tmp_vac_reinf_oldvoc[c(c_S)] <- tmp_vac_reinf_oldvoc[c(c_S)] + new_waning_reinf_suppreinf[1:10] + new_waning_reinf_suppreinf_oldvoc[1:10]
    }
    
    
    
    # account for stochastic issues to end up with negative results
    # and if "uptake(a) > alive unvaccinated individuals(a)"
    if(bool_stochastic | bool_vaccine){
      tmp             <- pmax(tmp,0)
      tmp_vac_reinf   <- pmax(tmp_vac_reinf,0)
      tmp_vac_reinf_oldvoc   <- pmax(tmp_vac_reinf_oldvoc,0)
    }
    
    
    ## Update lists  ---- 
    ##----------------- -
    nfull[[step_index + 1]] = tmp;
    nfull_vac_reinf[[step_index + 1]]= tmp_vac_reinf;
    nfull_vac_reinf_oldvoc[[step_index + 1]]= tmp_vac_reinf_oldvoc;
    
    # Same for vaccine compartments
    if(bool_vaccine){
      tmp_vac_rna1         <- pmax(tmp_vac_rna1,0)
      tmp_vac_rna2         <- pmax(tmp_vac_rna2,0)
      tmp_vac_adeno1       <- pmax(tmp_vac_adeno1,0)
      tmp_vac_adeno2       <- pmax(tmp_vac_adeno2,0)
      tmp_vac_waning       <- pmax(tmp_vac_waning,0)
      tmp_vac_booster      <- pmax(tmp_vac_booster,0)
      tmp_vac_booster_waning <- pmax(tmp_vac_booster_waning,0)
      tmp_vac_reinfvac <- pmax(tmp_vac_reinfvac,0)
      tmp_vac_reinfvac_oldvoc <- pmax(tmp_vac_reinfvac_oldvoc,0)
      tmp_vac_extrabooster1      <- pmax(tmp_vac_extrabooster1,0)
      tmp_vac_extrabooster2      <- pmax(tmp_vac_extrabooster2,0)
      
      nfull_vac_rna1[[step_index + 1]] = tmp_vac_rna1;
      nfull_vac_rna2[[step_index + 1]] = tmp_vac_rna2;
      nfull_vac_adeno1[[step_index + 1]] = tmp_vac_adeno1;
      nfull_vac_adeno2[[step_index + 1]] = tmp_vac_adeno2;
      nfull_vac_waning[[step_index + 1]] = tmp_vac_waning;
      nfull_vac_booster[[step_index + 1]]= tmp_vac_booster;
      nfull_vac_booster_waning[[step_index + 1]]= tmp_vac_booster_waning;
      nfull_vac_reinfvac[[step_index + 1]]= tmp_vac_reinfvac;
      nfull_vac_reinfvac_oldvoc[[step_index + 1]]= tmp_vac_reinfvac_oldvoc;
      nfull_vac_extrabooster1[[step_index + 1]]= tmp_vac_extrabooster1;
      nfull_vac_extrabooster2[[step_index + 1]]= tmp_vac_extrabooster2;
    }
  }
  

  ## Model output  ---- 
  ##----------------- - 
  sim_results             = do.call(rbind, nfull)
  sim_results_vac_rna1    = do.call(rbind, nfull_vac_rna1)
  sim_results_vac_rna2    = do.call(rbind, nfull_vac_rna2)
  sim_results_vac_adeno1  = do.call(rbind, nfull_vac_adeno1)
  sim_results_vac_adeno2  = do.call(rbind, nfull_vac_adeno2)
  sim_results_vac_waning  = do.call(rbind, nfull_vac_waning)
  sim_results_vac_booster = do.call(rbind, nfull_vac_booster)
  sim_results_vac_booster_waning = do.call(rbind, nfull_vac_booster_waning)
  sim_results_vac_reinf = do.call(rbind, nfull_vac_reinf)
  sim_results_vac_reinf_oldvoc = do.call(rbind, nfull_vac_reinf_oldvoc)
  sim_results_vac_reinfvac = do.call(rbind, nfull_vac_reinfvac)
  sim_results_vac_reinfvac_oldvoc = do.call(rbind, nfull_vac_reinfvac_oldvoc)
  sim_results_vac_extrabooster1 = do.call(rbind, nfull_vac_extrabooster1)
  sim_results_vac_extrabooster2 = do.call(rbind, nfull_vac_extrabooster2)
  sim_results_vac_extrabooster = sim_results_vac_extrabooster1 + sim_results_vac_extrabooster2
  
  sim_results_infection = do.call(rbind,nsubset_infection)
  sim_results_reinfection = do.call(rbind,nsubset_reinfection)
  sim_results_numvac = do.call(rbind,nsubset_numvac)
  
  ## Data (log)likelihood  ---- 
  ##------------------------- -
  ## aggregate by compartiment type
  sim_aggr <- sim_results + sim_results_vac_rna1 + sim_results_vac_rna2 + sim_results_vac_adeno1 + sim_results_vac_adeno2 + sim_results_vac_waning + sim_results_vac_booster + sim_results_vac_booster_waning + sim_results_vac_reinf + sim_results_vac_reinf_oldvoc + sim_results_vac_reinfvac + sim_results_vac_reinfvac_oldvoc + sim_results_vac_extrabooster
  
  # aggregate by burden of disease type
  new_hosp_icu <- data.frame(cases = sim_aggr[,c_new_hosp_total] + sim_aggr[,c_new_hosp_total_voc],
                             day = times_day)
  total_new_hosp_icu <- aggregate(new_hosp_icu[,1:10], by=list(new_hosp_icu$day), sum)
  
  new_deaths <- data.frame(cases = rbind(0,diff(sim_aggr[,c_D] + sim_aggr[,c_Dvoc])),
                           day = times_day)
  total_new_deaths = aggregate(. ~ day, data = new_deaths, sum)
  
  nI_sev_sims <- data.frame(cases = sim_aggr[,c_I_sev] + sim_aggr[,c_Ivoc_sev],
                            day = times_day)
  total_nI_sev_sims = aggregate(nI_sev_sims[,1:10], by=list(nI_sev_sims$day), sum)
  
  # calculate the sum of the hospital and ICU prevalence (h repetitions!)
  nHosp_ICU_sims <- data.frame(cases = sim_aggr[,c_I_hosp] + sim_aggr[,c_I_icu] +
                                 sim_aggr[,c_Ivoc_hosp] + sim_aggr[,c_Ivoc_icu],
                               day = times_day)
  total_nHosp_ICU_sims = aggregate(nHosp_ICU_sims[,1:10], by=list(nHosp_ICU_sims$day), sum)
  
  # calculate the sum of the VOC-related hospital and ICU prevalence (h repetitions!)
  nHosp_ICU_voc_sims <- data.frame(cases = sim_aggr[,c_Ivoc_hosp] + sim_aggr[,c_Ivoc_icu],
                                   day = times_day)
  total_nHosp_ICU_voc_sims = aggregate(nHosp_ICU_voc_sims[,1:10], by=list(nHosp_ICU_voc_sims$day), sum)
  
  # calculate the sum of the ICU prevalence (h repetitions!) WITHOUTH ICU
  nICU_sims <- data.frame(cases = sim_aggr[,c_I_icu] + sim_aggr[,c_Ivoc_icu], 
                          day = times_day)
  total_nICU_sims = aggregate(. ~ day, data = nICU_sims, sum)
  
  # calculate the sum of the hospital prevalence (h repetitions!) WITHOUTH ICU
  nHosp_sims <- data.frame(cases = sim_aggr[,c_I_hosp] + sim_aggr[,c_Ivoc_hosp], 
                           day = times_day)
  total_nHosp_sims = aggregate(. ~ day, data = nHosp_sims, sum)
  
  prob_sero1 = (sens_fun(seq(30,0,-h))%*%sim_results[which(times <= 30),c_asym_mild_infected])/cohort.size;
  prob_sero2 = (sens_fun(seq(51,0,-h))%*%sim_results[which(times <= 51),c_asym_mild_infected])/cohort.size;
  
  n1_mild_sims = data.frame(cases = sim_aggr[,c_Ivoc_mild],
                            day   = times_day)
  total_n1_mild_sims = apply(aggregate(n1_mild_sims, by=list(n1_mild_sims$day), sum)[2:11],1,sum)
  
  n2_mild_sims = data.frame(cases = sim_aggr[,c_I_mild] + sim_aggr[,c_Ivoc_mild],
                            day   = times_day)
  total_n2_mild_sims = apply(aggregate(n2_mild_sims, by=list(n2_mild_sims$day), sum)[2:11],1,sum)
  
  p_VOC = total_n1_mild_sims/total_n2_mild_sims;
  
  # hospital recoveries
  hosp_recoveries <- data.frame(cases = sim_aggr[,c_hosp_recovered] + sim_aggr[,c_hosp_recovered_voc],
                                day = times_day)
  total_hosp_recoveries <- aggregate(hosp_recoveries[,1:10], by=list(hosp_recoveries$day), sum)
  
  # hospital exit
  total_hosp_exit <- total_hosp_recoveries + total_new_deaths
  
  # fraction in ICU
  nI_ICU_sims <- data.frame(cases = sim_aggr[,c_I_icu] + sim_aggr[,c_Ivoc_icu],
                            day = times_day)
  total_nI_ICU_sims = aggregate(nI_ICU_sims[,1:10], by=list(nI_ICU_sims$day), sum)
  p_ICU <- rowSums(total_nI_ICU_sims[-1]) / rowSums(total_nHosp_ICU_sims[,-1])
  
  # prevalence of (a)symptomatic cases
  sim_I_asym <- data.frame(cases = sim_aggr[,c_I_asym] + sim_aggr[,c_Ivoc_asym] +
                             sim_aggr[,c_I_presym] + sim_aggr[,c_Ivoc_presym],
                           day = times_day)
  prev_I_asym = aggregate(sim_I_asym[,1:10], by=list(sim_I_asym$day), mean)
  
  prev_I_sev = aggregate(nI_sev_sims[,1:10], by=list(nI_sev_sims$day), sum)
  
  sim_I_mild <- data.frame(cases = sim_aggr[,c_I_mild] + sim_aggr[,c_Ivoc_mild],
                           day = times_day)
  prev_I_mild = aggregate(sim_I_mild[,1:10], by=list(sim_I_mild$day), mean)
  
  # hospital admissions for covid (with estimation/correction)
  pred_hosp1_for_covid <- as.matrix(total_new_hosp_icu)
  pop_data_be  <- get_regional_pop("belgium")
  pop60 <- pop_data_be[6] + pop_data_be[7] + pop_data_be[8] + pop_data_be[9] + pop_data_be[10]
  esti_prev60 <-  as.matrix(rowSums(prev_I_asym[,7:11])) 
  esti_prev60 <- matrix(esti_prev60, ncol=10, nrow=length(esti_prev60))
  pred_hosp1_for_covid[,-1] <- pred_hosp1_for_covid[,-1] - (pred_hosp1_for_covid[,-1]/rowSums(pred_hosp1_for_covid[,-1]))*(parms['for_covid_coef'])*esti_prev60/pop60
  pred_hosp1_for_covid[pred_hosp1_for_covid < 0] <- 0
  pred_hosp1_for_covid[is.na(pred_hosp1_for_covid)] <- 0
  total_new_hosp_icu_for_covid <- pred_hosp1_for_covid
  
  
  # additional statistics
  if(!bool_full_output){
    
    # DETACH ----
    detach(icol_tmp_vec)
    
    ## Return
    return(list(total_nI_sev_sims = total_nI_sev_sims,
                total_new_hosp_icu = total_new_hosp_icu,
                total_new_deaths   = total_new_deaths,
                total_nHosp_ICU_sims = total_nHosp_ICU_sims,
                total_nHosp_ICU_voc_sims = total_nHosp_ICU_voc_sims,
                total_nHosp_sims = total_nHosp_sims,
                total_nICU_sims = total_nICU_sims,
                prob_sero1 = prob_sero1,
                prob_sero2 = prob_sero2,
                p_VOC     = p_VOC,
                total_new_hosp_icu_for_covid = total_new_hosp_icu_for_covid,
                total_hosp_recoveries = total_hosp_recoveries,
                total_hosp_exit = total_hosp_exit,
                p_ICU=p_ICU))
    
  } else {
    
    # compute new infection accounting for reinfection changes in S compartments - round to avoid negative numbers due to numerical errors
    new_infections <- data.frame(cases = rbind(0,round(sim_results_infection[-1,]+sim_results_reinfection[-1,])),
                                 day = times_day)
    total_new_infections = aggregate(. ~ day, data = new_infections, sum)
    
    new_reinfections <- data.frame(cases = rbind(0,round(sim_results_reinfection[-1,])),
                                   day = times_day)
    total_new_reinfections <- aggregate(. ~ day, data = new_reinfections, sum)
    
    new_numvac <- data.frame(cases = rbind(0,round(sim_results_numvac[-1,])),
                                   day = times_day)
    total_new_numvac <- aggregate(. ~ day, data = new_numvac, sum)
    
    hosp_icu_load = aggregate(nHosp_ICU_sims[,1:10], by=list(nHosp_ICU_sims$day), mean) # mean prevelance per hour
    icu_load      = aggregate(nICU_sims[,1:10], by=list(nICU_sims$day), mean) # mean prevelance per hour
    
    new_mild_infections <- data.frame(cases = sim_aggr[,c_mild_infected] + sim_aggr[,c_mild_infected_voc],
                                      day = times_day)
    total_new_mild_infections = aggregate(. ~ day, data = new_mild_infections, sum)
    
    recov_mild_infection <- data.frame(cases = sim_aggr[,c_recov_mild_infection] + sim_aggr[,c_recov_mild_infection_voc],
                                       day = times_day)
    total_recov_mild_infection = aggregate(. ~ day, data = recov_mild_infection, sum)
    
    new_icu <- data.frame(cases = sim_aggr[,c_new_icu] + sim_aggr[,c_new_icu_voc],
                          day = times_day)
    total_new_icu <- aggregate(new_icu[,1:10], by=list(new_icu$day), sum)
    
    data_bar_coef <- data.frame(cases = bar_coef,
                                 day = times_day)
    mean_bar_coef <- aggregate(. ~ day, data = data_bar_coef, mean)
    # non-vaccinated hospital admissions
    new_hosp_icu_nvac   <- data.frame(cases = sim_results[,c_new_hosp_total] + sim_results[,c_new_hosp_total_voc]+
                                        sim_results_vac_reinf[,c_new_hosp_total] + sim_results_vac_reinf_oldvoc[,c_new_hosp_total]+
                                        sim_results_vac_reinf[,c_new_hosp_total_voc] + sim_results_vac_reinf_oldvoc[,c_new_hosp_total_voc],
                                      day = times_day)
    total_new_hosp_icu_nvac <- aggregate(new_hosp_icu_nvac[,1:10], by=list(new_hosp_icu_nvac$day), sum)
    
    # vaccinated 2doses and booster hospital admissions
    new_hosp_icu_vac_d2   <- data.frame(cases = sim_results_vac_adeno2[,c_new_hosp_total] + sim_results_vac_adeno2[,c_new_hosp_total_voc] +
                                          sim_results_vac_rna2[,c_new_hosp_total] + sim_results_vac_rna2[,c_new_hosp_total_voc] +
                                          sim_results_vac_waning[,c_new_hosp_total] + sim_results_vac_waning[,c_new_hosp_total_voc]+ 
                                          sim_results_vac_booster[,c_new_hosp_total] + sim_results_vac_booster[,c_new_hosp_total_voc] +
                                          sim_results_vac_booster_waning[,c_new_hosp_total] + sim_results_vac_booster_waning[,c_new_hosp_total_voc] +
                                          sim_results_vac_reinfvac[,c_new_hosp_total] + sim_results_vac_reinfvac[,c_new_hosp_total_voc] +
                                          sim_results_vac_reinfvac_oldvoc[,c_new_hosp_total] + sim_results_vac_reinfvac_oldvoc[,c_new_hosp_total_voc] +
                                          sim_results_vac_extrabooster[,c_new_hosp_total] + sim_results_vac_extrabooster[,c_new_hosp_total_voc] ,
                                        day = times_day)
    total_new_hosp_icu_vac_d2 <- aggregate(new_hosp_icu_vac_d2[,1:10], by=list(new_hosp_icu_vac_d2$day), sum)
    
    # vaccinated with booster hospital admissions
    #  new_hosp_icu_vac_booster   <- data.frame(cases = sim_results_vac_booster[,c_new_hosp_total] + sim_results_vac_booster[,c_new_hosp_total_voc] +
    #                                                   sim_results_vac_booster_waning[,c_new_hosp_total] + sim_results_vac_booster_waning[,c_new_hosp_total_voc] +
    #                                                sim_results_vac_reinfvac[,c_new_hosp_total] + sim_results_vac_reinfvac[,c_new_hosp_total_voc],
    #                                    day = times_day)
    #TODO temporary changed to reinfection admissions for checking !!!!!!!
   # new_hosp_icu_vac_booster   <- data.frame(cases = sim_results_vac_reinf[,c_new_hosp_total] + sim_results_vac_reinf[,c_new_hosp_total_voc] +
    #                                           sim_results_vac_reinf_oldvoc[,c_new_hosp_total] + sim_results_vac_reinf_oldvoc[,c_new_hosp_total_voc] +
     #                                          sim_results_vac_reinfvac[,c_new_hosp_total] + sim_results_vac_reinfvac[,c_new_hosp_total_voc]+
      #                                         sim_results_vac_reinfvac_oldvoc[,c_new_hosp_total] + sim_results_vac_reinfvac_oldvoc[,c_new_hosp_total_voc],
       #                                      day = times_day)
    #temp test specific vaccination compartment
    new_hosp_icu_vac_booster   <- data.frame(cases =sim_results_vac_extrabooster[,c_new_hosp_total] + sim_results_vac_extrabooster[,c_new_hosp_total_voc] ,
                                                                                     day = times_day)
    total_new_hosp_icu_vac_booster <- aggregate(new_hosp_icu_vac_booster[,1:10], by=list(new_hosp_icu_vac_booster$day), sum)
    
    # no vaccine
    susceptible_full <- data.frame(cases = sim_aggr[,c_S],#sim_results[,c_S],
                                   day = times_day)
    mean_susceptible_full = aggregate(. ~ day, data = susceptible_full, mean)
    
    susceptible_vac_d1 <- data.frame(cases = sim_results_vac_rna1[,c_S] + sim_results_vac_adeno1[,c_S],
                                     day = times_day)
    mean_susceptible_vac_d1 = aggregate(. ~ day, data = susceptible_vac_d1, mean)
    
    susceptible_vac_d2 <- data.frame(cases = sim_results_vac_rna2[,c_S] + sim_results_vac_adeno2[,c_S] +sim_results_vac_waning[,c_S] +   sim_results_vac_booster[,c_S] + sim_results_vac_booster_waning[,c_S] +  sim_results_vac_reinfvac[,c_S] +  sim_results_vac_reinfvac_oldvoc[,c_S]+   sim_results_vac_extrabooster[,c_S],
                                     day = times_day)
    mean_susceptible_vac_d2 = aggregate(. ~ day, data = susceptible_vac_d2, mean)
    
    susceptible_initpluswaning <- data.frame(cases = sim_results[,c_S]+  sim_results_vac_reinf[,c_S] +  sim_results_vac_reinf_oldvoc[,c_S] +sim_results_vac_waning[,c_S] + sim_results_vac_booster_waning[,c_S] +  sim_results_vac_reinfvac[,c_S] +  sim_results_vac_reinfvac_oldvoc[,c_S],
                                     day = times_day)
    mean_susceptible_initpluswaning = aggregate(. ~ day, data = susceptible_initpluswaning, mean)
    
    total_vaccinated <- data.frame(cases = sim_results_vac_rna1 + sim_results_vac_rna2 + 
                                     sim_results_vac_adeno1 + sim_results_vac_adeno2 +
                                     sim_results_vac_waning + sim_results_vac_booster +
                                     sim_results_vac_booster_waning + sim_results_vac_reinfvac + sim_results_vac_reinfvac_oldvoc + sim_results_vac_extrabooster,
                                   day = times_day)
    
    # select ODE compartments (by excluding the columns to keep track of incidences)
    total_vaccinated <- total_vaccinated[,c(icol_flag_ode,1)==1]
    
    # aggregate
    mean_vaccinated = aggregate(. ~ day, data = total_vaccinated, mean)
    mean_vaccinated_age <- matrix(NA,nrow=nrow(mean_vaccinated),ncol=10)
    i_age <- 1
    for(i_age in 1:10){
      mean_vaccinated_age[,i_age] <-  rowSums(mean_vaccinated[,seq(i_age,ncol(mean_vaccinated)-1,10)])
    }
    tail(mean_vaccinated_age)
    mean_vaccinated_age <- data.frame(day=mean_vaccinated$day,
                                      cases=mean_vaccinated_age)
    
    # vaccinated: age, dose type
    # total_vaccinated <- data.frame(cases = sim_results_vac_rna1 + sim_results_vac_rna2 + 
    #                                  sim_results_vac_adeno1 + sim_results_vac_adeno2,
    #                                day = times_day)
    # mean_vaccinated = aggregate(. ~ day, data = total_vaccinated, mean)
    num_age_rna1           <- matrix(NA,nrow=nrow(sim_results_vac_rna1),ncol=10)
    num_age_rna2           <- num_age_rna1
    num_age_adeno1         <- num_age_rna1
    num_age_adeno2         <- num_age_rna1
    num_age_waning         <- num_age_rna1
    num_age_booster        <- num_age_rna1
    num_age_booster_waning <- num_age_rna1 
    num_age_reinf          <- num_age_rna1
    num_age_reinf_oldvoc          <- num_age_rna1
    num_age_reinfvac       <- num_age_rna1
    num_age_reinfvac_oldvoc       <- num_age_rna1
    num_age_extrabooster        <- num_age_rna1
    i_age <- 1
    
    age_ode_columns <- which(icol_flag_ode==1)
    
    for(i_age in 1:10){
      
      # select age-specific ODE compartments
      age_columns <- seq(i_age,ncol(sim_results_vac_rna1)-1,10)
      age_columns <- age_columns[age_columns %in% age_ode_columns]
      
      num_age_rna1[,i_age]           <- rowSums(sim_results_vac_rna1[,age_columns])
      num_age_rna2[,i_age]           <- rowSums(sim_results_vac_rna2[,age_columns])
      num_age_adeno1[,i_age]         <- rowSums(sim_results_vac_adeno1[,age_columns])
      num_age_adeno2[,i_age]         <- rowSums(sim_results_vac_adeno2[,age_columns])
      num_age_waning[,i_age]         <- rowSums(sim_results_vac_waning[,age_columns])
      num_age_booster[,i_age]        <- rowSums(sim_results_vac_booster[,age_columns])
      num_age_booster_waning[,i_age] <- rowSums(sim_results_vac_booster_waning[,age_columns])
      num_age_reinf[,i_age]          <- rowSums(sim_results_vac_reinf[,age_columns])
      num_age_reinf_oldvoc[,i_age]          <- rowSums(sim_results_vac_reinf_oldvoc[,age_columns])
      num_age_reinfvac[,i_age]       <- rowSums(sim_results_vac_reinfvac[,age_columns])
      num_age_reinfvac_oldvoc[,i_age]       <- rowSums(sim_results_vac_reinfvac_oldvoc[,age_columns])
      num_age_extrabooster[,i_age]        <- rowSums(sim_results_vac_extrabooster[,age_columns])
    }
    total_dose_vaccinated <- data.frame(day=times_day,
                                        rna_d1_age = num_age_rna1,
                                        rna_d2_age = num_age_rna2,
                                        adeno_d1_age = num_age_adeno1,
                                        adeno_d2_age = num_age_adeno2,
                                        waning_age   = num_age_waning,
                                        booster_age  = num_age_booster,
                                        booster_waning_age  = num_age_booster_waning,
                                        #extrabooster_age  = num_age_extrabooster,  #temp, to check
                                        reinfvac_age  = num_age_reinfvac + num_age_reinfvac_oldvoc)
    mean_dose_vaccinated = aggregate(. ~ day, data = total_dose_vaccinated, mean)
    names(mean_dose_vaccinated) <- gsub('\\.','',names(mean_dose_vaccinated))
    
    # DETACH ----
    detach(icol_tmp_vec)
    
    ## Return
    return(list(sim_results = sim_results,
                sim_results_vac_rna1 = sim_results_vac_rna1,
                sim_results_vac_rna2 = sim_results_vac_rna2,
                sim_results_vac_adeno1 = sim_results_vac_adeno1,
                sim_results_vac_adeno2 = sim_results_vac_adeno2,
                sim_results_vac_waning = sim_results_vac_waning,
                sim_results_vac_booster= sim_results_vac_booster,
                sim_results_vac_booster_waning= sim_results_vac_booster_waning,
                sim_results_vac_reinf= sim_results_vac_reinf + sim_results_vac_reinf_oldvoc,
                sim_results_vac_reinfvac= sim_results_vac_reinfvac + sim_results_vac_reinfvac_oldvoc,
                sim_results_vac_extrabooster= sim_results_vac_extrabooster,
                
                total_nI_sev_sims = total_nI_sev_sims,
                total_new_hosp_icu = total_new_hosp_icu,
                total_new_icu      = total_new_icu,
                total_new_deaths   = total_new_deaths,
                total_nHosp_ICU_sims = total_nHosp_ICU_sims,
                total_nHosp_ICU_voc_sims = total_nHosp_ICU_voc_sims,
                total_nHosp_sims   = total_nHosp_sims,
                total_nICU_sims = total_nICU_sims,
                prob_sero1 = prob_sero1,
                prob_sero2 = prob_sero2,
                p_VOC      = p_VOC,
                total_hosp_recoveries = total_hosp_recoveries,
                total_hosp_exit = total_hosp_exit,
                p_ICU = p_ICU,
                
                total_new_infections = total_new_infections,
                total_new_mild_infections = total_new_mild_infections,
                total_recov_mild_infection = total_recov_mild_infection,
                hosp_icu_load = hosp_icu_load,
                icu_load = icu_load,
                
                total_new_reinfections = total_new_reinfections,
                total_new_numvac = total_new_numvac,
                
                prev_I_asym=prev_I_asym,
                prev_I_mild=prev_I_mild,
                prev_I_sev=prev_I_sev,
                
                total_new_hosp_icu_nvac = total_new_hosp_icu_nvac,
                total_new_hosp_icu_vac_d2 = total_new_hosp_icu_vac_d2,
                total_new_hosp_icu_vac_booster = total_new_hosp_icu_vac_booster,
                total_new_hosp_icu_for_covid = total_new_hosp_icu_for_covid,
                
                mean_susceptible_full = mean_susceptible_full,
                mean_susceptible_initpluswaning = mean_susceptible_initpluswaning,
                mean_susceptible_vac_d1 = mean_susceptible_vac_d1,
                mean_susceptible_vac_d2 = mean_susceptible_vac_d2,
                mean_vaccinated_age = mean_vaccinated_age,
                mean_dose_vaccinated = mean_dose_vaccinated,norm_asy=norm_asy,norm_sy=norm_sy,mean_bar_coef=mean_bar_coef)
    )
  } 
}  

log_likelihood_model <- function(parms, CoMix_matrices, method = "mean", plots = "FALSE",
                                 parms_names, ndays_sim = NA, #vaccine_uptake = NA,
                                 vaccine_uptake = NA,
                                 V_mat, V_mat_dose2, V_mat_booster, V_mat_2ndbooster,
                                 V_mat_extrabooster,
                                 be_ref_data = NULL){
  
  # run the dynamic model ----
  if(is.na(ndays_sim)) { ndays_sim = length(obs_total_hosp_ext)};
  model_out <- run_model(parms=parms, 
                         CoMix_matrices=CoMix_matrices, 
                         method=method, 
                         parms_names=parms_names,
                         ndays_sim=ndays_sim,
                         bool_full_output=FALSE,
                         vaccine_uptake = vaccine_uptake,
                         V_mat = V_mat,
                         V_mat_dose2 = V_mat_dose2,
                         V_mat_booster = V_mat_booster,
                         V_mat_2ndbooster = V_mat_2ndbooster,
                         V_mat_extrabooster = V_mat_extrabooster)
  
  if(!is.null(be_ref_data)){
    obs_total_hosp_ext   <- be_ref_data$hospital_admissions + be_ref_data$hospital_admissions_other
    obs_forcovid_hosp    <- be_ref_data$hospital_admissions 
    
    # age_dist_hosp_mat_ext <- get_regional_hospital_age_distr(be_ref_data$region[1])
    #age_dist_hosp_mat_ext <- get_hospital_age_distr(bool_use_default = FALSE)
    age_dist_hosp_mat_ext <- t(be_ref_data[,grepl('prop_hospital_admission_age',names(be_ref_data))])
    missing_hosp_age_data <- length(obs_total_hosp_ext) - ncol(age_dist_hosp_mat_ext)
    if(missing_hosp_age_data>0){
      age_dist_hosp_mat_ext <- cbind(age_dist_hosp_mat_ext,
                                     age_dist_hosp_mat_ext[,rep(ncol(age_dist_hosp_mat_ext),missing_hosp_age_data)])   
    }
    
    obs_age_dist_mort <- be_ref_data[,grepl('covid19_deaths_age',names(be_ref_data))]
    obs_hosp_load_ext <- be_ref_data$hospital_load
    obs_icu_load_ext  <- be_ref_data$icu_load
    obs_hosp_exit_ext <- be_ref_data$hospital_exit
  } 
  
  # re-capture some parameters 
  h        = parms[parms_names == 'h'];  
  sel_region = get_region(parms[parms_names == 'region_id'])
  
  # parse model output
  total_nI_sev_sims      = model_out[['total_nI_sev_sims']]
  total_new_hosp_icu     = model_out[['total_new_hosp_icu']]
  total_new_hosp_icu_for_covid     = model_out[['total_new_hosp_icu_for_covid']]
  total_new_deaths       = model_out[['total_new_deaths']]
  total_nHosp_ICU_sims   = model_out[['total_nHosp_ICU_sims']]       # required for Binomial distribution
  total_nHosp_ICU_voc_sims = model_out[['total_nHosp_ICU_voc_sims']] # required for Binomial distribution
  total_nHosp_sims       = model_out[['total_nHosp_sims']]           # required for Binomial distribution
  total_nICU_sims        = model_out[['total_nICU_sims']]            # required for Binomial distribution
  
  total_nHosp_ICU_load = rowSums(total_nHosp_ICU_sims[,-1])*h  # = average hospital load / day
  total_nICU_load    = rowSums(total_nICU_sims[,-1])*h  # = average ICU load / day
  prob_sero1         = model_out[['prob_sero1']]
  prob_sero2         = model_out[['prob_sero2']]
  p_VOC             = model_out[['p_VOC']]
  
  total_hosp_recoveries = rowSums(model_out[['total_hosp_recoveries']][,-1])
  total_hosp_exit       = rowSums(model_out[['total_hosp_exit']][,-1])
  p_ICU                 = model_out[['p_ICU']]
  
  omega    = exp(parms[grepl('log_omega',parms_names)]); 
  delta3   = exp(parms[grepl('log_delta3',parms_names)]); 
  delta3_wave2 = exp(parms['log_delta3_wave2']);    # recovery rate from hospital (from sept 2021)
  ### 1. Data sources (hospitalizations, mortality, serology)  ---- 
  ###------------------------------------------------------------ -
  
  cohort.size <- get_regional_pop(sel_region)
  
  # ## REFERENCE DATA
  # be_ref_data <- get_latest_incidence_data()
  
  #### 1.1. Daily number of new hospitalizations by age group (from 1/3 onwards)
  ####------------------------------------------------------------------------ -
  y1 = obs_total_hosp_ext;
  y1_mat = round_tot(y1*t(age_dist_hosp_mat_ext),0)
  if(nrow(y1_mat)>ndays_sim){
    y1_mat <- y1_mat[1:ndays_sim,]
  }
  
  #### 1.2. Daily number of new deaths by age group (from 1/3 onwards) 
  ####-------------------------------------------------------------- -
  y2_mat = obs_age_dist_mort;
  y2 = rowSums(y2_mat)
  if(nrow(y2_mat)>ndays_sim){
    y2_mat <- y2_mat[1:ndays_sim,]
  }
  
  #### 1.3. Serial serological survey results based on residual samples 
  ####--------------------------------------------------------------- -
  y3 = sero_data_all$igg_cat_pos[sero_data_all$cround == 1];
  n3 = sero_data_all$igg_cat_total[sero_data_all$cround == 1];
  y4 = sero_data_all$igg_cat_pos[sero_data_all$cround == 2];
  n4 = sero_data_all$igg_cat_total[sero_data_all$cround == 2];
  
  #### 1.4. Prevalence of Variants of Concern: VOC strain in the model (alpha, omicron,bq1)
  ####----------------------------------------- -
  y5 = db_voc_raw$voc_aggr + db_voc_raw$voc_omicron + db_voc_raw$voc_bq1ba275xbb
  n5 = db_voc_raw$n_sequenced
  y5_date <- db_voc_raw$date
  y5_days <- sim_date2day(y5_date)

  
  # make daily
  y5_days <- sim_date2day(approx_sequenced(y5_date,y5)$x)
  y5      <- approx_sequenced(y5_date,y5)$y
  n5      <- approx_sequenced(y5_date,n5)$y
  
  # y5 = db_voc_raw$sftg_omicron
  # n5 = db_voc_raw$sftg_total
  # y5_days   = sim_date2day(db_voc_raw$date)
  # adjust nday_sim
  if(max(y5_days)>ndays_sim){
    y5     <- y5[y5_days<=ndays_sim]
    n5     <- n5[y5_days<=ndays_sim]
    y5_days <- y5_days[y5_days<=ndays_sim] 
  }
  
  ####----------------------------------------- -
  y6_exit <- obs_hosp_exit_ext
  if(length(y6_exit)>ndays_sim){
    y6_exit <- y6_exit[1:ndays_sim]
  }
  
  y6_load <- obs_hosp_load_ext
  if(length(y6_load)>ndays_sim){
    y6_load <- y6_load[1:ndays_sim]
  }
  
  #### 1.6. Hospital and ICU load
  ####----------------------------------------- -
  y7_icu  <- obs_icu_load_ext  
  n7_hosp <- obs_hosp_load_ext 
  y7_days <- which(!is.na(y7_icu))
  if(length(y7_icu)>ndays_sim){
    y7_icu  <- y7_icu[1:ndays_sim]
    n7_hosp <- n7_hosp[1:ndays_sim]
    y7_days <- y7_days[y7_days<=ndays_sim] 
  }
  
  #### 1.1. Daily number of new hospitalizations for covid
  ####------------------------------------------------------------------------ -
  obs_forcovid_hosp <- obs_forcovid_hosp[1:ndays_sim]
  
  ### 2. Loglikelihood contributions  ---- 
  ###----------------------------------- -
  
  #### 2.1. Binomial distribution for the daily number of new hospitalizations 
  ####---------------------------------------------------------------------- -
  ll_hosp_icu_mat <- matrix(nrow = nrow(y1_mat), ncol = 10)
  for (j in 1:10){
    ll_hosp_icu_mat[,j] = dbinom(y1_mat[,j], size = round(unlist(total_nI_sev_sims[1:nrow(y1_mat),j+1]),0), 
                                 prob = 1 - exp(-h*omega[j]), log = T)  
  }
  
  ll_hosp_icu_mat = ifelse(ll_hosp_icu_mat == -Inf, -100000, ll_hosp_icu_mat)
  
  #### 2.2. Weighed Least Squares for the daily number of new deaths (assuming equal death rates for hosp/ICU)
  ####------------------------------------------------------------------------------------------------------ -
  #note: transition rates change over time with adjusted LoS, mortality, severity...
  
  ls_hosp_mort  = sum((y2_mat - total_new_deaths[1:nrow(y2_mat),-1])^2,na.rm=T)
  wls_hosp_mort = sum((y2_mat[,-1] - total_new_deaths[1:nrow(y2_mat),-(1:2)])^2 / total_new_deaths[1:nrow(y2_mat),-(1:2)],na.rm=T)
  
  ll_mort_mat    <- matrix(nrow = nrow(y2_mat), ncol = 10)
  j <- 2
  for (j in 2:10){
    if(sum(y2_mat[,j],na.rm=T)/(ndays_sim/7)<1){
      rollmean_k <- 30
    } else{
      rollmean_k <- 7
    }
    n_mort    <- round(rollsum(unlist(total_new_deaths[1:nrow(y2_mat),j+1]),k=rollmean_k,na.rm=T,fill = NA),0)
    hosp_load <- round(rollsum(unlist(total_nHosp_ICU_sims[1:nrow(y2_mat),j+1]),k=rollmean_k,na.rm=T,fill = NA),0)
    p_mort    <- n_mort / hosp_load
    
    ll_mort_mat[,j] <- dbinom(x = rollsum(y2_mat[,j],k=rollmean_k,na.rm=T,fill = NA),
                              size = hosp_load, 
                              prob = p_mort, 
                              log = T)
  }
  
  ll_mort_mat = ifelse(ll_mort_mat == -Inf, -100, ll_mort_mat)
  # colSums(ll_mort_mat,na.rm=T)
  # sum(ll_mort_mat,na.rm=T)
  
  #### 2.3. Binomial distribution for the serial serological data 
  ####--------------------------------------------------------- -
  
  ll_sero1 = dbinom(y3, size = n3, prob = prob_sero1, log = T)
  ll_sero1 = ifelse(ll_sero1 == -Inf, -100000, ll_sero1)
  ll_sero2 = dbinom(y4, size = n4, prob = prob_sero2, log = T)
  ll_sero2 = ifelse(ll_sero2 == -Inf, -100000, ll_sero2)
  
  
  #### 2.4. Binomial distribution for the VOC-strain prevalence (start 1/1/2021)
  ####------------------------------------------------------------------ -
  ll_VOC = dbinom(y5, size = n5, prob = p_VOC[y5_days], log = T);
  ll_VOC = ll_VOC[ll_VOC != -Inf]
  
  #### 2.4.a Binomial distribution for the hospital exit
  ####------------------------------------------------------------------ -
  ll_hosp_exit = dbinom(y6_exit, size = round(rowSums(total_nHosp_ICU_sims[1:length(y6_exit),-1]),0), 
                        prob = 1 - exp(-h*delta3), log = T)
  ll_hosp_exit = ifelse(ll_hosp_exit == -Inf, -100000, ll_hosp_exit)
  
  #### 2.4.b Weighed Least Squares for the hospital load
  ####------------------------------------------------------------------ -
  ls_hosp_load  = sum((y6_load - total_nHosp_ICU_load[1:length(y6_load)])^2,na.rm=T)
  wls_hosp_load = sum((y6_load - total_nHosp_ICU_load[1:length(y6_load)])^2 / total_nHosp_ICU_load[1:length(y6_load)],na.rm=T)
  
  #### 2.5. Binomial distribution for the fraction ICU load / hospital load
  ####------------------------------------------------------------------ -
  ll_icu_fraction = dbinom(y7_icu[y7_days], size = n7_hosp[y7_days], prob = p_ICU[y7_days], log = T)
  ll_icu_fraction = ifelse(ll_icu_fraction == -Inf, -100000, ll_icu_fraction)
  
  #### 2.4.b Weighed Least Squares for the ICU load
  ####------------------------------------------------------------------ -
  wls_icu_load = sum((y7_icu[y7_days] - total_nICU_load[y7_days])^2 / total_nICU_load[y7_days],na.rm=T)
  
  
  #### Weighed Least Squares for "for covid" coef
  ####------------------------------------------------------------------ -
  startlenght = length(obs_forcovid_hosp)-90
  total_new_hosp_icu_for_covid = rowSums(total_new_hosp_icu_for_covid[,-1])
  ls_forcovid  = sum((obs_forcovid_hosp - total_new_hosp_icu_for_covid[1:length(obs_forcovid_hosp)])^2,na.rm=T)
  wls_forcovid = sum((obs_forcovid_hosp[startlenght:length(obs_forcovid_hosp)] - total_new_hosp_icu_for_covid[startlenght:length(obs_forcovid_hosp)])^2 / total_new_hosp_icu_for_covid[startlenght:length(obs_forcovid_hosp)],na.rm=T)
  
  ## Graphical exploration  ---- 
  ##-------------------------- - 
  if(plots == "TRUE"){
    plot_log_likelihood_model(total_nI_sev_sims,total_new_hosp_icu,
                              total_new_deaths,total_nHosp_ICU_load,
                              total_nICU_load,
                              total_hosp_exit,
                              y1_mat,y2_mat,
                              prob_sero1,prob_sero2,y3,n3,y4,n4,
                              y5,n5,y5_days,p_VOC,n7_hosp,y7_icu,
                              y6_exit,y6_days)
  }
  
  ## Loglikelihood function value  ---- 
  ##--------------------------------- -
  w1 = (cohort.size/n3);         # poststratification weights (inverse probability weighing)
  w2 = (cohort.size/n4);         # poststratification weights (inverse probability weighing)
  
  denom = length(ll_hosp_icu_mat)
  length(ll_sero1) + length(ll_sero2) 
  
  w1_imp = 1/(length(ll_hosp_icu_mat)/denom)                                # importance weights
  w3_imp = 1/(length(ll_sero1)/denom)                                       # importance weights
  w4_imp = 1/(length(ll_sero2)/denom)                                       # importance weights
  
  ### Weighted likelihood function: inverse probability and importance weighting
  ###----------------------------------------------------------- -
  sumll1 = w1_imp*sum(ll_hosp_icu_mat,na.rm=T) + 
    w3_imp*sum((w1/sum(w1))*ll_sero1,na.rm=T) + 
    w4_imp*sum((w2/sum(w2))*ll_sero2,na.rm=T);
  
  ### Weighted likelihood function: hosp+sero+VOC
  ###-------------------------------------------------------------------------- -
  sumll2 = sumll1 + w1_imp*sum(ll_VOC,na.rm=T); 
  
  ### Weighted likelihood function: hospital and ICU load
  ###-------------------------------------------------------------------------- -
  sumll3 = -wls_hosp_load - wls_icu_load
  
  ### WLS function: mortality
  ###-------------------------------------------------------------------------- -
  #sumll4 = -wls_hosp_mort 
  sumll4 = sum(ll_mort_mat,na.rm=T)
  
  ### Weighted likelihood function: VOC
  ###-------------------------------------------------------------------------- -
  sumll5 = sum(ll_VOC,na.rm=T)
  
  ### Weighted likelihood function: for covid + load + mortality
  ###-------------------------------------------------------------------------- -
  sumll6 = -wls_forcovid + sumll3 + sumll4
  
  ### Weighted likelihood function: for covid
  ###-------------------------------------------------------------------------- -
  sumll7 = -wls_forcovid 
  
  return(list(crit1 = sumll1, 
              crit2 = sumll2, 
              crit3 = sumll3, 
              crit4 = sumll4, 
              crit5 = sumll5,
              crit6 = sumll6,
              crit7 = sumll7,
              dev1 = -2*sumll1,
              dev2 = -2*sumll2,
              dev3 = -2*sumll3,
              dev4 = -2*sumll4,
              dev5 = -2*sumll5,
              dev6 = -2*sumll6,
              dev7 = -2*sumll7
  ))
}

## Original graphical exploration  ---- 
##-------------------------- - 
plot_log_likelihood_model <- function(total_nI_sev_sims,total_new_hosp_icu,
                                      total_new_deaths,total_nHosp_ICU_load,
                                      total_nICU_load,
                                      total_hosp_exit,
                                      y1_mat,y2_mat,
                                      prob_sero1,prob_sero2,y3,n3,y4,n4,
                                      y5,n5,y5_days,p_VOC,n7_hosp,y7_icu,
                                      y6_exit,y6_days,
                                      bool_seroprev=FALSE){
  
  # check figure margins, and abort if too large (or figure is to narrow)
  if(any(dev.size()<4)){
    warning('Issue with dev.size, figure are too narrow,, skip plotting')
    return(NULL)
  }
  
  if(!is_vsc()) par(mfrow = c(3,3))
  plot(1:nrow(y1_mat), apply(y1_mat, 1, sum),cex=0.5,
       ylab = "Number of new hospitalizations", xlab = "Days since introduction",
       #ylim = c(0,1000),
       ylim = c(0,max(rowSums(y1_mat)*1.1,na.rm = T)),
       xlim =range(0,nrow(y1_mat),nrow(total_new_hosp_icu)))
  lines(apply(total_new_hosp_icu[-1,-1],1,sum), col = 2, lwd = 2)
  
  plot(total_new_hosp_icu[-1,3], col = "orange", lwd = 2, lty = 1,
       ylab = "Number of new hospitalizations", xlab = "Days since introduction",
       ylim = c(0,max(y1_mat*1.1,na.rm=T)), type = "l")
  points(1:nrow(y1_mat), y1_mat[,2], col = "orange",cex=0.5)
  lines(1:nrow(y1_mat), total_new_hosp_icu[2:(nrow(y1_mat)+1),5], col = "purple", lwd = 2, lty = 2)
  points(1:nrow(y1_mat), y1_mat[,4], col = "purple",cex=0.5)
  lines(1:nrow(y1_mat), total_new_hosp_icu[2:(nrow(y1_mat)+1),7], col = "red", lwd = 2, lty = 3)
  points(1:nrow(y1_mat), y1_mat[,6], col = "red",cex=0.5)
  lines(1:nrow(y1_mat), total_new_hosp_icu[2:(nrow(y1_mat)+1),9], col = "brown", lwd = 2, lty = 4)
  points(1:nrow(y1_mat), y1_mat[,8], col = "brown",cex=0.5)
  lines(1:nrow(y1_mat), total_new_hosp_icu[2:(nrow(y1_mat)+1),11], col = "pink", lwd = 2, lty = 5)
  points(1:nrow(y1_mat), y1_mat[,10], col = "pink",cex=0.5)
  
  legend('topleft',
         paste(seq(2,10,2)),
         pch=16,
         ncol=5,
         cex=0.7,
         col = c('orange','purple','red','brown','pink'))
  
  plot(1:nrow(y2_mat), apply(y2_mat,1,sum),
       ylab = "Number of new deaths", xlab = "Days since introduction",
       ylim = c(0,150),
       xlim = range(0,nrow(total_new_deaths)),
       cex=0.5)
  lines(apply(total_new_deaths[-1,-1],1,sum),
        col = 2, lwd = 2)
  
  plot(total_new_deaths[-1,3], col = "orange", lwd = 2, lty = 1,
       ylab = "Number of new deaths", xlab = "Days since introduction",
       ylim = c(0,60), type = "l")
  points(1:nrow(y2_mat), y2_mat[,2], col = "orange",cex=0.5)
  lines(1:nrow(y2_mat), total_new_deaths[1:nrow(y2_mat),5], col = "purple", lwd = 2, lty = 2)
  points(1:nrow(y2_mat), y2_mat[,4], col = "purple",cex=0.5)
  lines(1:nrow(y2_mat), total_new_deaths[1:nrow(y2_mat),7], col = "red", lwd = 2, lty = 3)
  points(1:nrow(y2_mat), y2_mat[,6], col = "red",cex=0.5)
  lines(1:nrow(y2_mat), total_new_deaths[1:nrow(y2_mat),9], col = "brown", lwd = 2, lty = 4)
  points(1:nrow(y2_mat), y2_mat[,8], col = "brown",cex=0.5)
  lines(1:nrow(y2_mat), total_new_deaths[1:nrow(y2_mat),11], col = "pink", lwd = 2, lty = 5)
  points(1:nrow(y2_mat), y2_mat[,10], col = "pink",cex=0.5)
  
  legend('topleft',
         paste(seq(2,10,2)),
         pch=16,
         ncol=5,
         cex=0.7,
         col = c('orange','purple','red','brown','pink'))
  
  # VOC
  plot(y5_days,
       y5/n5, 
       xlim = range(0,length(p_VOC)),
       ylim = 0:1,
       las=1,
       ylab =  "Prevalence of VOC-strain", xlab = "Days since introduction",cex=0.5)
  grid()
  lines(p_VOC, col = "red", lwd = 2)
  
  ## HOSPITAL LOAD
  plot(n7_hosp,
       ylab = "Hospital load", xlab = "Days since introduction",
       ylim = c(0,max(n7_hosp*1.1,na.rm=T)),
       xlim = range(0,length(total_nHosp_ICU_load)),
       cex=0.5)
  lines(total_nHosp_ICU_load, col = "red", lwd = 2)
  
  ## ICU LOAD
  plot(y7_icu,
       ylab = "ICU load", xlab = "Days since introduction",
       ylim = c(0,max(y7_icu*1.1,na.rm=T)),
       xlim = range(0,length(total_nICU_load)),
       cex=0.5)
  lines(total_nICU_load, col = "red", lwd = 2)
  
  ## EXIT
  plot(y6_exit,
       ylab = "Hospital exit", xlab = "Days since introduction")
  lines(rollmean(y6_exit,k=7,align='center',fill=0),col=5,lwd=3)
  lines(total_hosp_exit,col='red',lwd=2)
  
  ## SEROPREVALENCE (optional)
  #if(bool_seroprev){
  plot(1:length(prob_sero1),y3/n3, ylim = c(0,0.2),
       xlab = "Age category", ylab = "Seroprevalence",pch=16)
  lines(1:length(prob_sero1), prob_sero1, col = 1, lwd = 2)
  
  points(1:length(prob_sero2),y4/n4,col=4,pch=16)
  lines(1:length(prob_sero2), prob_sero2, col = 4, lwd = 2)
  legend('topleft',
         c('March 2020',
           'April 2020'),
         pch=1,
         lwd=2,
         cex=0.7,
         col=c(1,4),
         ncol=2)
  #}
  
}



## Log prior distributions
##------------------------ -
log_prior_model <- function(parms,parms_names){
  # make sure that the parameter names are set
  names(parms[1:length(parms_names)]) <- parms_names
  
  gamma  = exp(parms['log_gamma']);     # average length of the latency period - 2 days
  theta  = exp(parms['log_theta']);     # infectious period at pre-symptomatic stage - 3.5 days
  delta1 = exp(parms['log_delta1']);    # infectious period at asymptomatic stage - 3.5 days
  delta2 = exp(parms['log_delta2']);    # infectious period at (mild) symptomatic stage - 3.5 days
  
  phi0 = expit(parms[grepl('log_phi0',parms_names)]); # proportion of symptomatically infected with mild symptoms
  p    = p_vec;                                       # proportion of asymptomatic cases
  
  omega = exp(parms[grepl('log_omega',parms_names)]);              # waiting time between symptom onset and hospitalization
  phi1  = expit(parms[grepl('log_phi1_age',parms_names)]);         # proportion of severly infected with regular hospitalization
  
  q = exp(parms['log_q']);                # proportionality constant
  
  # all mu parameters
  mu_param = parms[grepl('mu.*_sev',parms_names)]
  mu_sev = expit(mu_param);          # age-specific hospital fatality ratio - probability of dying
  delta3 = exp(parms[grepl('log_delta3',parms_names)]);  # recovery period in severely infected/length of stay - 7 days
  
  mortality_changepoint = parms[grepl('mortality_changepoint',parms_names)] #FYI: not part of estimation process
  
  alpha       = exp(parms[grepl('log_comix_coef',parms_names) | grepl('log_add_coef',parms_names)]);  
  
  n0 = exp(parms["log_n0"]) # imported cases in 2020
  
  voc_start = parms[grepl('VOC.*_start',parms_names)]
  
  voc_init    = exp(parms[grepl('log.*_init',parms_names)]) #exp(parms['log_VOC_alpha_init']);
  voc_transm  = exp(parms[grepl('log.*_transm',parms_names)]) #exp(parms['log_VOC_alpha_transm']);    # additional transmissibility
  voc_hr_hosp = exp(parms[grepl('log.*hr_hosp',parms_names)]) # hospitalization admission hazard ratio
  
  phi1_add  = parms[grepl('phi1_add',parms_names)]  # additional factor for regular hospitalization (instead of ICU)
  
  voc_gamma_factor = exp(parms[grepl('log.*_VOC_.*_gamma_factor',parms_names)]); # adjustment factor for the latency period with Omicron
  
  log_gamma_prior    <- dnorm(gamma, mean = 1/2, sd = 0.05, log = T);
  log_theta_prior    <- dnorm(theta, mean = 1/3.5, sd = 0.05, log = T);
  log_delta1_prior   <- dnorm(delta1, mean = 1/3.5, sd = 0.05, log = T);
  log_delta2_prior   <- dnorm(delta2, mean = 1/7, sd = 0.05, log = T);
  
  log_phi0_prior     <- dunif(phi0, 0, 1, log = T);
  log_omega_prior    <- dunif(omega, 1/7 , 1/0.5, log= T)
  log_phi1_prior     <- dunif(phi1, 0, 1, log = T);
  
  log_mu_prior       <- dunif(mu_sev, 0, 1, log = T);
  log_delta3_prior   <- dunif(delta3, 0, 0.2, log = T);
  
  log_phi1_add_prior   <- dunif(phi1_add, 0, 3, log = T);
  
  log_alpha_prior    <- dunif(alpha, exp(-5), exp(3), log = T);
  
  log_n0_prior       <- dunif(n0, 0, 5000, log = T);
  
  log_voc_start        <- dunif(voc_start, global_lib_voc$start_ll[1], global_lib_voc$start_ul[1], log = T);
  for(i_voc in 1:length(voc_start)){
    log_voc_start[i_voc]       <- dunif(voc_start[i_voc], global_lib_voc$start_ll[i_voc], global_lib_voc$start_ul[i_voc], log = T);
  }
  
  log_voc_init_prior   <- dunif(voc_init, 20, 10000, log = T);
  log_voc_transm_prior <- dunif(voc_transm, 0.001, 11, log = T);
  log_voc_hr_hosp_prior<- dunif(voc_hr_hosp, 0, 5, log = T);
  log_voc_gamma      <- dunif(voc_gamma_factor, 1,1e10, log = T)
  
  fparam    <- 0.51;
  betas_asy <- fparam * q * CoMix_mat$C_CoMix2010_asy;
  betas_sy  <- q * CoMix_mat$C_CoMix2010_sy;
  
  p_asy <- p_vec;        
  phi1  <- phi1_vec; 
  D1    <- 1/theta; 
  D2    <- 1/delta1;
  D3    <- 1/delta2;
  D4    <- 1/omega;
  
  # R0 estimation is based on national statistics
  cohort.size <- get_regional_pop(region = "belgium")
  
  R1 <- D1*(betas_asy * cohort.size) + 
    p_asy*D2*(betas_asy * cohort.size) + 
    (1-p_asy)*D3*(betas_sy * cohort.size) +
    (1-p_asy)*(1-phi0)*D4*(betas_sy * cohort.size);
  
  R0            <- Re((eigen(R1)$values[1]));
  
  log_q_prior   <- dnorm(R0, mean = 3, sd = 0.01, log = T)
  
  log_prior_val <-  c(log_gamma_prior, 
                      log_theta_prior, 
                      log_delta1_prior, 
                      log_delta2_prior,
                      log_q_prior, 
                      log_phi0_prior, 
                      log_omega_prior,
                      log_mu_prior,
                      log_delta3_prior, 
                      log_alpha_prior,
                      
                      log_n0_prior,
                      
                      log_voc_start,
                      
                      log_voc_init_prior,
                      log_voc_transm_prior, 
                      log_voc_hr_hosp_prior,
                      log_phi1_add_prior,
                      
                      log_voc_gamma)
  
  ll_penalty <- -500000
  log_prior_val <- ifelse(log_prior_val == -Inf, ll_penalty, log_prior_val)
  
  invalid_param <- NULL
  # include check  
  if(any(log_prior_val == ll_penalty)){
    #print(which(log_prior_val == -1e5))    
    invalid_param <- names(log_prior_val)[log_prior_val == ll_penalty]
  }
  
  log_prior_val <- sum(log_prior_val)
  
  return(list(log_prior_val = log_prior_val, R0val = R0, invalid_param = invalid_param))
}


