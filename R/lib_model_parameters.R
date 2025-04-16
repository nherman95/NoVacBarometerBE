########################################################################### #
# This file is part of the Stochastic Compartmental Model for SARS-COV-2 
# transmission in Belgium, conceived by members of SIMID group during the 
# COVID19 pandemic.
#
# This file is used to obtain, modify and explore model parameters.
#
# Copyright 2021, SIMID, University of Antwerp & Hasselt University                                        
########################################################################### #

# load library
library(RColorBrewer)

update2latest_model_parameter_config <- function(parms){
  

  parms_names <- names(parms)
  
  #temporary change to ba4ba5 hr parameters
 # parms[grepl('log_VOC_ba4ba5_hr_hosp',names(parms))] <- 0
  #parms_names <- names(parms)
  #print("adjust log_VOC_ba4ba5_hr_hosp to 0",fill=T)

    if(!any(grepl('VOC_alpha',names(parms)))){
    sel_voc_alpha <- grepl('VOC_',parms_names) & !grepl('delta_',parms_names) & !grepl('omicron_',parms_names)
    parms_names[sel_voc_alpha] <- gsub('VOC_','VOC_alpha_',parms_names[sel_voc_alpha])
    names(parms) <- parms_names
    print("adjust parameter names to VOC_alpha",fill=T)
    
    sel_voc_delta3 <- grepl('delta3_',parms_names) & !grepl('wave2',parms_names)
    parms_names[sel_voc_delta3] <- gsub('delta3_','delta3_VOC_',parms_names[sel_voc_delta3])
    names(parms) <- parms_names
    print("adjust parameter names to delta3_VOC_",fill=T)
    
    sel_voc_mu <- grepl('mu_',parms_names) & !grepl('mu_age',parms_names) & !grepl('mu_sev',parms_names)
    parms_names[sel_voc_mu] <- gsub('voc','alpha',parms_names[sel_voc_mu])
    parms_names[sel_voc_mu] <- gsub('mu_','mu_VOC_',parms_names[sel_voc_mu])
    names(parms) <- parms_names
    print("adjust parameter names to mu_VOC_",fill=T)
    
    sel_voc_ve <- grepl('ve_',parms_names) & !grepl('delta',parms_names) & !grepl('omicron',parms_names) & !grepl('rate',parms_names) & !grepl('transmission',parms_names)
    parms_names[sel_voc_ve] <- gsub('ve_','ve_VOC_alpha_',parms_names[sel_voc_ve])
    parms_names <- gsub('ve_delta','ve_VOC_delta',parms_names)
    parms_names <- gsub('ve_omicron','ve_VOC_omicron',parms_names)
    names(parms) <- parms_names
    print("adjust parameter names to ve_VOC_",fill=T)
    
    }
  
  
  if(!any(grepl('VOC_ba4ba5',names(parms)))){
    sel_param_omicron <- parms[grepl('VOC_omicron',parms_names)]
    parms[gsub('omicron','ba4ba5',names(sel_param_omicron))] <- sel_param_omicron
    parms_names <- names(parms)
    
    parms['VOC_ba4ba5_init']  <- 0
    parms['VOC_ba4ba5_start'] <- 1e10 #760  # default: no introduction
    #parms['log_VOC_ba4ba5_transm'] <- parms['log_VOC_omicron_transm'] * 1.3
    print("include dummy parameter names for VOC_ba4ba5",fill=T)
  }
  if(!any(grepl('VOC_bq1ba275xbb_incr',names(parms)))){
    sel_param_omicron <- parms[grepl('VOC_ba4ba5',parms_names)]
    parms[gsub('ba4ba5','bq1ba275xbb',names(sel_param_omicron))] <- sel_param_omicron
    parms_names <- names(parms)
    parms['ve_VOC_bq1ba275xbb_vereduction'] <- exp(parms['log_VOC_bq1_vereduction'])
    parms['VOC_bq1ba275xbb_init'] <- parms['log_VOC_bq1_init'] 
    parms['VOC_bq1ba275xbb_start'] <- sim_date2day('2022-08-15')
    vesevreduction <- exp(parms["log_VOC_bq1_vesevreduction"])
       sel_param <- parms[grepl(paste0("ve_VOC_","delta","_incr_severe"),parms_names)]
    parms[grepl(paste0("ve_VOC_","bq1ba275xbb","_incr_severe"),parms_names)] <- sel_param*(vesevreduction) 
    parms_names <- names(parms)
    print("include dummy parameter names for VOC_bq1ba275xbb",fill=T)
  }
  
  if(!any(grepl('VOC_xbb15',names(parms)))){
    sel_param_omicron <- parms[grepl('log_mu_VOC_bq1ba275xbb',parms_names)]
    parms[gsub('bq1ba275xbb','xbb15',names(sel_param_omicron))] <- sel_param_omicron
    print("!!!Reuse mortality parameter from bq1ba275xbb for VOC_xbb15",fill=T)
  }
  
  if(!any(grepl('ve_VOC_ba4ba5_vereduction',names(parms)))){

    parms['ve_VOC_omicron_vereduction'] <- 0.5
    parms['ve_VOC_ba4ba5_vereduction'] <- 0.5
    parms['ve_VOC_bq1ba275xbb_vereduction'] <- 0.5
    parms['ve_VOC_xbb15_vereduction'] <- 0.5
    parms['log_VOC_omicron_transm'] <- 1
    parms['log_VOC_ba4ba5_transm'] <- 1
    parms['log_VOC_bq1ba275xbb_transm'] <- 1
    parms['log_VOC_xbb15_transm'] <- 1
    parms_names <- names(parms)
    print("introduce vereduction for omicron and ba4ba5 instead of transm, reset bq1ba275xbb/xbb15",fill=T)
  }
    

    if(any(grepl('for_covid_coef',names(parms))) && parms['for_covid_coef']<0){
      parms['for_covid_coef'] <- (10000)
      parms_names <- names(parms)
      print("default 'for_covid_coef parameter",fill=T)
    }
    
  if(!any(grepl('for_covid_coef',names(parms)))){
    parms['for_covid_coef'] <- (20000)
    parms_names <- names(parms)
    print("default 'for_covid_coef parameter",fill=T)
  }
  
  if(!any(grepl('ve_.*_booster_waning',names(parms)))){
    #with waning booster similar to 2nd dose - infection only !!!  #temp
    parms['ve_waning_booster_rate'] <- 1/70 #parms['ve_waning_immunity_rate']
    parms['ve_VOC_alpha_infection_booster_waning']         <- parms['ve_VOC_alpha_infection_booster']        # reduction with respect to infection with waning immunity
    parms['ve_VOC_alpha_incr_severe_booster_waning'] <- parms['ve_VOC_alpha_incr_severe_booster']       # incremental reduction within dynamic model
    parms['ve_VOC_delta_infection_booster_waning']   <- parms['ve_VOC_delta_infection_booster']   # reduction with respect to infection with waning immunity
    parms['ve_VOC_delta_incr_severe_booster_waning'] <- parms['ve_VOC_delta_incr_severe_booster']
    parms['ve_VOC_omicron_infection_booster_waning']   <- parms['ve_VOC_omicron_infection_booster']   # reduction with respect to infection with waning immunity
    parms['ve_VOC_omicron_incr_severe_booster_waning'] <- parms['ve_VOC_omicron_incr_severe_booster']
    parms['ve_VOC_alpha_incr_severe_waning']       <- parms['ve_VOC_alpha_incr_severe_rna2']       # incremental reduction within dynamic model
    parms['ve_VOC_delta_incr_severe_waning'] <- parms['ve_VOC_delta_incr_severe_rna2']
    parms['ve_VOC_omicron_incr_severe_waning'] <- parms['ve_VOC_omicron_incr_severe_rna2']
    parms_names <- names(parms)
    print("custom ve_booster_waning parameters",fill=T)
  }
  
  if(!any(grepl('ve_.*_reinfvac',names(parms)))){
    parms['ve_waning_infection_rate'] <- 0 #1/730 #parms['ve_waning_booster_rate']
    parms['ve_VOC_alpha_infection_reinf']     <- parms['ve_VOC_alpha_infection_waning']
    parms['ve_VOC_alpha_incr_severe_reinf']   <- parms['ve_VOC_alpha_incr_severe_waning']
    parms['ve_VOC_delta_infection_reinf']     <- parms['ve_VOC_delta_infection_waning']
    parms['ve_VOC_delta_incr_severe_reinf']   <- parms['ve_VOC_delta_incr_severe_waning']
    parms['ve_VOC_omicron_infection_reinf']   <- parms['ve_VOC_omicron_infection_waning']
    parms['ve_VOC_omicron_incr_severe_reinf'] <- parms['ve_VOC_omicron_incr_severe_booster']
    
    parms['ve_waning_infection_booster_rate'] <- 0 # parms['ve_waning_booster_rate']
    parms['ve_VOC_alpha_infection_reinfvac']     <- parms['ve_VOC_alpha_infection_booster']
    parms['ve_VOC_alpha_incr_severe_reinfvac']   <- parms['ve_VOC_alpha_incr_severe_booster']   
    parms['ve_VOC_delta_infection_reinfvac']     <- parms['ve_VOC_delta_infection_booster']
    parms['ve_VOC_delta_incr_severe_reinfvac']   <- parms['ve_VOC_delta_incr_severe_booster']
    parms['ve_VOC_omicron_infection_reinfvac']   <- parms['ve_VOC_omicron_infection_booster']
    parms['ve_VOC_omicron_incr_severe_reinfvac'] <- parms['ve_VOC_omicron_incr_severe_booster']
    parms_names <- names(parms)
    print("default reinf parameters",fill=T)
  }
  
  if(!any(grepl('log_.*_hr_hosp',names(parms)))){
    parms['log_VOC_alpha_hr_hosp'] <- 1e-15
    parms['log_VOC_delta_hr_hosp'] <- log(2.26)
    parms[paste0('log_VOC_omicron_hr_hosp_age',1:10)] <- log(c(1.0,0.89,0.67,0.57,0.54,0.42,0.32,0.42,0.49,0.49))
    parms_names <- names(parms)
    print("default hr parameters",fill=T)
  }
  
  if(!any(grepl('log_spring2022_phi1_add',names(parms)))){
    parms['icu_spring2022_start']    <- sim_date2day('2022-03-20')
    parms['log_spring2022_phi1_add'] <- parms['log_VOC_omicron_phi1_add']
    parms_names <- names(parms)
    print("default spring22 parameters",fill=T)
  }
  
  if(!any(grepl('mu_VOC_alpha_sev',names(parms)))){
    colnames_mu2 <- names(parms)[grepl('mu2_sev',names(parms))]
    parms[gsub('mu2_sev','mu_VOC_alpha_sev',colnames_mu2)] <- parms[colnames_mu2]
    parms[gsub('mu2_sev','mu_VOC_delta_sev',colnames_mu2)] <- parms[colnames_mu2]
    parms_names <- names(parms)
    print("default VOC_alpha parameters",fill=T)
  }
  
  if(!any(grepl('extrabooster',names(parms)))){
     parms['ve_VOC_alpha_infection_extrabooster']         <- parms['ve_VOC_alpha_infection_booster']  
     parms['ve_VOC_delta_infection_extrabooster']         <- parms['ve_VOC_delta_infection_booster']  
     parms['ve_VOC_omicron_infection_extrabooster']         <- parms['ve_VOC_omicron_infection_booster']  
     parms['ve_VOC_ba4ba5_infection_extrabooster']         <- parms['ve_VOC_ba4ba5_infection_booster']  
     parms['ve_VOC_alpha_incr_severe_extrabooster']         <- parms['ve_VOC_alpha_incr_severe_booster']  
     parms['ve_VOC_delta_incr_severe_extrabooster']         <- parms['ve_VOC_delta_incr_severe_booster']  
     parms['ve_VOC_omicron_incr_severe_extrabooster']         <- parms['ve_VOC_omicron_incr_severe_booster']  
     parms['ve_VOC_ba4ba5_incr_severe_extrabooster']         <- parms['ve_VOC_ba4ba5_incr_severe_booster'] 
     parms_names <- names(parms)
    # print("default 'extrabooster' parameter",fill=T)
  }
  
  #hub scenarios concerning extrabooster protection (present here because it's easier for the moment...)
 # parms['ve_VOC_omicron_incr_severe_extrabooster']         <- parms['ve_VOC_omicron_incr_severe_booster']  
  #parms['ve_VOC_ba4ba5_incr_severe_extrabooster']         <- parms['ve_VOC_ba4ba5_incr_severe_booster'] 
  
  #-20pc
#  parms['ve_VOC_omicron_infection_extrabooster']         <- parms['ve_VOC_omicron_infection_booster']*0.8
 # parms['ve_VOC_ba4ba5_infection_extrabooster']         <- parms['ve_VOC_ba4ba5_infection_booster']*0.8

  #delta 
   #parms['ve_VOC_omicron_infection_extrabooster']         <- parms['ve_VOC_delta_infection_booster']
  #parms['ve_VOC_ba4ba5_infection_extrabooster']         <- parms['ve_VOC_delta_infection_booster']

   #parms_names <- names(parms)
   #print("hub specific 'extrabooster' parameter",fill=T)
  
  return(parms)
}


explore_next_gen <- function(parms_chains,x_lim=NULL){

  # set number of comix and additional waves
  num_age_groups <- 10
  nwave_Comix  <- identify_nb_waves(names(parms_chains),bool_comix = TRUE)
  nwave_other  <- identify_nb_waves(names(parms_chains),bool_comix = FALSE)

  # get q parameters
  q <- unique(exp(parms_chains['log_q'])); # proportionality constant
  
  M_day      <- get_M_change_day()
  M_date     <- sim_day2date(M_day)  
  CoMix_day  <- get_M_change_day()[is_M_CoMix()]
  CoMix_date <- sim_day2date(CoMix_day)
  add_date   <- M_date[!M_date %in% CoMix_date]
  q_date     <- c(CoMix_date[1:nwave_Comix],add_date[-1])
  
  i_param <- 1
  parm_names <- names(parms_chains)
  age_meta <- data.frame(tag = get_age_groups(),
                         # col = brewer.pal(10,'Set3'),
                         col = rainbow(num_age_groups),
                         #col = 1:10,
                         pch = 1:num_age_groups)
  
  ## ACCROSS CHAINS (polygon)
  i_age <- 1
  par(mfrow=c(3,4))
  for(i_age in 1:10){
    parms_age <- exp(parms_chains[,parm_names[grepl(paste0('log_.*_coef_.*_age',i_age),parm_names)]])
    parms_age <- parms_age[,gsub('log_.*_coef_.*_age','',names(parms_age)) == i_age]
    names(parms_age)
    dim(parms_age)
    
    q_summary <- data.frame(q_mean = apply(parms_age,2,mean),
                         # q_975 = apply(parms_age,2,quantile,0.975),
                         # q_025 = apply(parms_age,2,quantile,0.025),
                         q_LL = apply(parms_age,2,min),
                         q_UL = apply(parms_age,2,max))
    y_max <- 20
    plot(q_date,
         q_summary$q_mean,
         col=0,
         ylim = c(0,y_max),
         xaxt='n',
         main = age_meta$tag[i_age],
         ylab='coeff (age)',
         xlab='date')
    grid(ny=NULL,nx=NA)
    add_polygon(scenario_data = q_summary,
                scenario_dates = q_date,
                col=2)
    add_date_axis(q_date)
    
    #add_change_points(M_date)
    wave_id_sel <- seq(0,length(CoMix_date),5) 
    wave_id_sel[1] <- 1
    abline(v=CoMix_date[wave_id_sel],
           lty=2,
           col=8)
    text(wave_id_sel,
         x=CoMix_date[wave_id_sel],
         y=y_max,
         col=1,
         cex = 0.75)
    
    # add horizontal line with q = 1
    abline(h=1,lty=3)
  }
  
  ## Q (polygon)
  i_age <- 1
  par(mfrow=c(3,4))
  for(i_age in 1:10){
    parms_age <- exp(parms_chains[,parm_names[grepl(paste0('log_.*_coef_.*_age',i_age),parm_names)]])
    parms_age <- parms_age[,gsub('log_.*_coef_.*_age','',names(parms_age)) == i_age]
    
    CoMix_q <- unlist(q) * apply(parms_age,2,rep,1,each=nrow(q))
    
    q_summary <- data.frame(q_mean = apply(CoMix_q,2,mean),
                         # q_975 = apply(CoMix_q,2,quantile,0.975),
                         # q_025 = apply(CoMix_q,2,quantile,0.025),
                         q_LL = apply(CoMix_q,2,min),
                         q_UL = apply(CoMix_q,2,max))
    
    y_max <- 2
    plot(q_date,
         q_summary$q_mean,
         col=0,
         ylim = c(0,y_max),
         xaxt='n',
         main = age_meta$tag[i_age],
         ylab='q (age)',
         xlab='date')
    grid(ny=NULL,nx=NA)
    add_polygon(scenario_data = q_summary,
                scenario_dates = q_date,
                col=2)
    add_date_axis(M_date[-1])
    
    #add_change_points(M_date)
    wave_id_sel <- seq(0,length(CoMix_date),5) 
    wave_id_sel[1] <- 1
    abline(v=CoMix_date[wave_id_sel],
           lty=2,
           col=8)
    text(wave_id_sel,
         x=CoMix_date[wave_id_sel],
         y=y_max,
         col=1,
         cex = 0.75)
    
    # add horizontal line with q = 1
    #abline(h=1,lty=3)
  }
}

include_new_wave <- function(f_parms,bool_comix=TRUE){
  
  col_tag <- ifelse(bool_comix,
                    'log_comix_coef_w',
                    'log_add_coef_mat_w')
  
  nwave_Comix      <- identify_nb_waves(names(f_parms),bool_comix)
  i_w_age10        <- which(grepl(paste0(col_tag,nwave_Comix,'_age10'),names(f_parms)))
  i_new            <- c(1:i_w_age10,         # param 1 till wave X-1
                        i_w_age10- (9:0),    # new wave
                        (i_w_age10+1):length(f_parms)) # final param
  
  if(is.null(ncol(f_parms))){
    f_parms <- f_parms[i_new]
  } else{
    f_parms <- f_parms[,i_new]
  }
  
  i_col_prev <- i_w_age10- (9:0)
  i_col_new  <- i_w_age10+(1:10)
  names(f_parms)[i_col_new] <- gsub(nwave_Comix,nwave_Comix+1,names(f_parms[i_col_prev]))
  
  print(paste("WARNING: INCLUDED PARAMETER FOR WAVE",nwave_Comix+1))
  
  return(f_parms)
}

remove_last_wave <- function(f_parms,bool_comix=TRUE){
  
  col_tag <- ifelse(bool_comix,
                    'log_comix_coef_w',
                    'log_add_coef_mat_w')
  
  nwave_Comix      <- identify_nb_waves(names(f_parms),bool_comix)
  i_new        <- which(!grepl(paste0(col_tag,nwave_Comix,'_age'),names(f_parms)))

  if(is.null(ncol(f_parms))){
    f_parms <- f_parms[i_new]
  } else{
    f_parms <- f_parms[,i_new]
  }
  
  print(paste("WARNING: REMOVED PARAMETERS FOR",ifelse(bool_comix,"COMIX","ADDITIONAL"),"WAVE",nwave_Comix))
  
  return(f_parms)
}

get_wave_colnames <- function(wave_id,bool_comix=TRUE){
  
  if(bool_comix){
    return(paste0(rep(paste0('log_comix_coef_w',wave_id),each=10),'_age',1:10))
  }
  else {
    return(paste0(rep(paste0('log_add_coef_mat_',wave_id),each=10),'_age',1:10))
  }
}

identify_nb_waves <- function(parms_names, bool_comix = TRUE){
  
  col_tag <- ifelse(bool_comix,
                    'log_comix_coef_w',
                    'log_add_coef_mat_w')
  
  nb_age_groups    <- sum(grepl(paste0(col_tag,1,'_'),parms_names))
  nb_columns <- sum(grepl(col_tag,parms_names))
  nb_waves   <- nb_columns / nb_age_groups
  
  return(nb_waves)
}

identify_total_nb_waves <- function(parms_names){
  
  return(identify_nb_waves(parms_names,bool_comix = TRUE) + 
           identify_nb_waves(parms_names,bool_comix = FALSE))
}

get_colnames <- function(parms_names,
                         tag_list = NA,
                         ignore_list=NULL,
                         sel_id = NA){
  
  bool_selection <- matrix(FALSE,length(parms_names))
  if(!any(is.na(sel_id))){
    bool_selection[sel_id] <- TRUE
  }
  if(!any(is.na(tag_list))){
    for(i_tag in tag_list){
      bool_selection[grepl(i_tag,parms_names)] <- TRUE
    }
  }
  if(!any(is.na(ignore_list))){
    for(i_tag in ignore_list){
      bool_selection[grepl(i_tag,parms_names)] <- FALSE
    }  
  }

  return(parms_names[bool_selection])
}

get_ve_severe<- function(ve_infection,ve_incr_severe){
  return(1-(1-ve_infection)*1-ve_incr_severe)  #warning : false... but never used
}
get_ve_incr_severe <- function(ve_infection,ve_severe){
  return(1-(1-ve_severe)/(1-ve_infection))
}

# ve_tag <- 've'
get_ve_infection <- function(parms,ve_tag){
  
  vereduction <- 1
  if(any(grepl(paste0(ve_tag,'_vereduction'),names(parms)))){
    vereduction <- exp(parms[paste0(ve_tag,'_vereduction')])  # now in log (name to be changed in the future)
    vereduction <- min(1,vereduction)
    vereduction <- max(0,vereduction)
  }

  ve_infection <- data.frame(adeno1    = parms[paste0(ve_tag,'_infection_adeno1')]) *vereduction  # reduction with respect to infection after 1st adeno-based dose
  ve_infection$adeno2    <- parms[paste0(ve_tag,'_infection_adeno2')]*vereduction   # reduction with respect to infection after 2nd adeno-based dose
  ve_infection$rna1      <- parms[paste0(ve_tag,'_infection_rna1')]*vereduction   # reduction with respect to infection after 1st mRNA dose 
  ve_infection$rna2      <- parms[paste0(ve_tag,'_infection_rna2')]*vereduction   # reduction with respect to infection after 2nd mRNA dose 
  ve_infection$waning    <- parms[paste0(ve_tag,'_infection_waning')]*vereduction   # reduction with respect to infection with waning immunity
  ve_infection$booster   <- parms[paste0(ve_tag,'_infection_booster')]*vereduction   # reduction with respect to infection with booster dose
  ve_infection$booster_waning    <- parms[paste0(ve_tag,'_infection_booster_waning')]*vereduction   # reduction with respect to infection with booster dose with waning immunity
  ve_infection$reinf     <- parms[paste0(ve_tag,'_infection_reinf')]   #reinfection path
  ve_infection$reinf_oldvoc     <- parms[paste0(ve_tag,'_infection_reinf')]*vereduction   #reinfection path (same initialy for compatibility)
  ve_infection$reinfvac  <- parms[paste0(ve_tag,'_infection_reinfvac')]   #reinfection path after vaccination
  ve_infection$reinfvac_oldvoc  <- parms[paste0(ve_tag,'_infection_reinfvac')]*vereduction   #reinfection path after vaccination (same initialy for compatibility)
  ve_infection$extrabooster    <- parms[paste0(ve_tag,'_infection_extrabooster')]*vereduction   # reduction with respect to infection with extra booster dose
  
  return(ve_infection)
}

save_vaccine_parameters <- function(){
  
  # VE references
  # Bernal 2021, NEJM : VE infection
  # Stowe2021 (pre-print): VE hospital admission
  # https://www.nature.com/articles/d41586-021-02054-z: VE transmission
  # Bernard 2021: pre-print on VOC omicron in England.

  #opti ECDC Hub R5
  waning <- 0.7
  waning_severe <- 1
  #pessi ECDC Hub R5 
  #waning <- 0.4
  #waning_severe <- 0.8
  
  vaccine_param <- data.frame(
    ve_VOC_alpha_infection_adeno1    = 0.49,  # reduction with respect to infection after 1st adeno-based dose
    ve_VOC_alpha_infection_adeno2    = 0.74,  # reduction with respect to infection after 2nd adeno-based dose
    ve_VOC_alpha_severe_adeno1       = 0.76,     # reduction with respect to severe infection with hospital admission after 1st adeno-based dose
    ve_VOC_alpha_severe_adeno2       = 0.86,     # reduction with respect to severe infection with hospital admission after 2nd adeno-based dose
    
    ve_VOC_alpha_infection_rna1      = 0.48,  # reduction with respect to infection after 1st mRNA dose #Bernal2021
    ve_VOC_alpha_infection_rna2      = 0.94,  # reduction with respect to infection after 2nd mRNA dose #Bernal2021
    ve_VOC_alpha_severe_rna1         = 0.83,     # reduction with respect to severe infection with hospital admission after 1st mRNA dose
    ve_VOC_alpha_severe_rna2         = 0.95,     # reduction with respect to severe infection with hospital admission after 2nd mRNA dose
    
    ve_VOC_alpha_infection_waning   = 0.627 ,  # [Andrews 2022] reduction with respect to infection after 2nd (mRNA) dose with waning immunity
    ve_VOC_alpha_severe_waning      = 0.917 ,  # [Andrews 2022] reduction with respect to severe infection with hospital admission 2nd (mRNA) dose with waning immunity

  
    ve_VOC_alpha_infection_booster   = 0.94 ,  # [copy of ve_infection_rna2] reduction with respect to infection after mRNA booster dose
    ve_VOC_alpha_severe_booster      = 0.95 ,  # [copy of ve_severe_rna2] reduction with respect to severe infection with hospital admission after mRNA booster dose

    ve_VOC_alpha_infection_booster_waning   = 0.889 , # [Andrews 2022] reduction with respect to infection after mRNA booster dose
     ve_VOC_alpha_severe_booster_waning      = 0.917 ,  # [copy of ve_severe_rna2] reduction with respect to severe infection with hospital admission after mRNA booster dose

      
    # protection delta variant
    ve_VOC_delta_infection_adeno1 = 0.429,       # 0.30 , # reduction with respect to infection after 1st adeno-based dose
    ve_VOC_delta_infection_adeno2 = 0.828,       # 0.67 ,  # reduction with respect to infection after 2nd adeno-based dose
    ve_VOC_delta_severe_adeno1    = 0.952 * 4/5, # 0.71 ,  # reduction with respect to severe infection with hospital admission after 1st adeno-based dose
    ve_VOC_delta_severe_adeno2    = 0.952,       # 0.92 ,  # reduction with respect to severe infection with hospital admission after 2nd adeno-based dose
    
    ve_VOC_delta_infection_rna1 = 0.723,         # 0.36 ,  # reduction with respect to infection after 1st mRNA dose
    ve_VOC_delta_infection_rna2 = 0.909,         # 0.88 ,  # reduction with respect to infection after 2nd mRNA dose
    ve_VOC_delta_severe_rna1    = 0.987 * 4/5,   # 0.94 ,  # reduction with respect to severe infection with hospital admission after 1st mRNA dose
    ve_VOC_delta_severe_rna2    = 0.987,          # 0.96 ,  # reduction with respect to severe infection with hospital admission after 2nd mRNA dose
    
    ve_VOC_delta_infection_waning   = 0.627 ,  # [Andrews 2022] reduction with respect to infection after 2nd (mRNA) dose with waning immunity
    ve_VOC_delta_severe_waning      = 0.917 ,  # [Andrews 2022] reduction with respect to severe infection with hospital admission 2nd (mRNA) dose with waning immunity
    
    ve_VOC_delta_infection_booster = 0.951, # 0.88,# 0.959 ,  # reduction with respect to infection after mRNA booster dose
    ve_VOC_delta_severe_booster    = 0.987, #0.996 ,  # reduction with respect to severe infection with hospital admission after mRNA booster dose
    
    ve_VOC_delta_infection_booster_waning   = 0.889 , # [Andrews 2022] reduction with respect to infection after mRNA booster dose
    ve_VOC_delta_severe_booster_waning      = 0.917 ,  # [copy of ve_severe_rna2] reduction with respect to severe infection with hospital admission after mRNA booster dose
    
    # protection against omicron variant (all omicron-like vocs, others are copy with different immune evasion)
    ve_VOC_omicron_infection_adeno1 = 0.177, # 0.129 , # reduction with respect to infection after 1st adeno-based dose
    ve_VOC_omicron_infection_adeno2 = 0.489, #0.190 ,  # reduction with respect to infection after 2nd adeno-based dose
    ve_VOC_omicron_severe_adeno1    = 0.81 * 4/5, # 0.497 ,  # reduction with respect to severe infection with hospital admission after 1st adeno-based dose
    ve_VOC_omicron_severe_adeno2    = 0.81,       # 0.600 ,  # reduction with respect to severe infection with hospital admission after 2nd adeno-based dose
    
    ve_VOC_omicron_infection_rna1 = 0.315,      # 0.187 ,  # reduction with respect to infection after 1st mRNA dose
    ve_VOC_omicron_infection_rna2 = 0.655,      # 0.241 ,  # reduction with respect to infection after 2nd mRNA dose
    ve_VOC_omicron_severe_rna1    = 0.81 * 4/5, # 0.596 ,  # reduction with respect to severe infection with hospital admission after 1st mRNA dose
    ve_VOC_omicron_severe_rna2    = 0.81,       # 0.668 ,  # reduction with respect to severe infection with hospital admission after 2nd mRNA dose
    
   # ve_VOC_omicron_infection_waning   = 0.088 ,  # [Andrews 2022] reduction with respect to infection after 2nd (mRNA) dose with waning immunity
    #ve_VOC_omicron_severe_waning      = 0.57 ,  # [CDC] reduction with respect to severe infection with hospital admission 2nd (mRNA) dose with waning immunity
   ve_VOC_omicron_infection_waning = 0.655 * waning,
    ve_VOC_omicron_severe_waning  = 0.81 * waning_severe,
     
    ve_VOC_omicron_infection_booster = 0.672 ,  # reduction with respect to infection after mRNA booster dose
    ve_VOC_omicron_severe_booster    = 0.90 ,  # reduction with respect to severe infection with hospital admission after mRNA booster dose
    
   # ve_VOC_omicron_infection_booster_waning   = 0.457 , # [Andrews 2022] reduction with respect to infection after mRNA booster dose
   # ve_VOC_omicron_severe_booster_waning      = 0.90 * 0.90 ,  # [copy of ve_severe_rna2] reduction with respect to severe infection with hospital admission after mRNA booster dose
    ve_VOC_omicron_infection_booster_waning   = 0.672 * waning,
    ve_VOC_omicron_severe_booster_waning      = 0.90 * waning_severe,
      
    # protection against infectiousness/transmission
    ve_transmission = 0,# 0.45,
    
    # start protection
    delay_protection_rna1   = 21,
    delay_protection_rna2   = 7,
    delay_protection_adeno1 = 21,
    delay_protection_adeno2 = 7
    
   # ve_waning_immunity_rate      = 1/180,  # 6 months
    #ve_waning_booster_rate       = 1/70    # 10 weeks
  )
  
  vaccine_param$ve_VOC_alpha_incr_severe_adeno1  = get_ve_incr_severe(vaccine_param$ve_VOC_alpha_infection_adeno1,vaccine_param$ve_VOC_alpha_severe_adeno1) # incremental reduction within dynamic model
  vaccine_param$ve_VOC_alpha_incr_severe_adeno2  = get_ve_incr_severe(vaccine_param$ve_VOC_alpha_infection_adeno2,vaccine_param$ve_VOC_alpha_severe_adeno2) # incremental reduction within dynamic model
  vaccine_param$ve_VOC_alpha_incr_severe_rna1    = get_ve_incr_severe(vaccine_param$ve_VOC_alpha_infection_rna1,vaccine_param$ve_VOC_alpha_severe_rna1) # incremental reduction within dynamic model
  vaccine_param$ve_VOC_alpha_incr_severe_rna2    = get_ve_incr_severe(vaccine_param$ve_VOC_alpha_infection_rna2,vaccine_param$ve_VOC_alpha_severe_rna2) # incremental reduction within dynamic model
  vaccine_param$ve_VOC_alpha_incr_severe_waning = get_ve_incr_severe(vaccine_param$ve_VOC_alpha_infection_waning,vaccine_param$ve_VOC_alpha_severe_waning) # incremental reduction within dynamic model
  vaccine_param$ve_VOC_alpha_incr_severe_booster = get_ve_incr_severe(vaccine_param$ve_VOC_alpha_infection_booster,vaccine_param$ve_VOC_alpha_severe_booster) # incremental reduction within dynamic model
  vaccine_param$ve_VOC_alpha_incr_severe_booster_waning = get_ve_incr_severe(vaccine_param$ve_VOC_alpha_infection_booster_waning,vaccine_param$ve_VOC_alpha_severe_booster_waning) # incremental reduction within dynamic model
  
  vaccine_param$ve_VOC_delta_incr_severe_adeno1  = get_ve_incr_severe(vaccine_param$ve_VOC_delta_infection_adeno1,vaccine_param$ve_VOC_delta_severe_adeno1) # incremental reduction within dynamic model
  vaccine_param$ve_VOC_delta_incr_severe_adeno2  = get_ve_incr_severe(vaccine_param$ve_VOC_delta_infection_adeno2,vaccine_param$ve_VOC_delta_severe_adeno2) # incremental reduction within dynamic model
  vaccine_param$ve_VOC_delta_incr_severe_rna1    = get_ve_incr_severe(vaccine_param$ve_VOC_delta_infection_rna1,vaccine_param$ve_VOC_delta_severe_rna1) # incremental reduction within dynamic model
  vaccine_param$ve_VOC_delta_incr_severe_rna2    = get_ve_incr_severe(vaccine_param$ve_VOC_delta_infection_rna2,vaccine_param$ve_VOC_delta_severe_rna2) # incremental reduction within dynamic model
  vaccine_param$ve_VOC_delta_incr_severe_waning  = get_ve_incr_severe(vaccine_param$ve_VOC_delta_infection_waning,vaccine_param$ve_VOC_delta_severe_waning) # incremental reduction within dynamic model
  vaccine_param$ve_VOC_delta_incr_severe_booster = get_ve_incr_severe(vaccine_param$ve_VOC_delta_infection_booster,vaccine_param$ve_VOC_delta_severe_booster) # incremental reduction within dynamic model
  vaccine_param$ve_VOC_delta_incr_severe_booster_waning = get_ve_incr_severe(vaccine_param$ve_VOC_delta_infection_booster_waning,vaccine_param$ve_VOC_delta_severe_booster_waning) # incremental reduction within dynamic model
  
  vaccine_param$ve_VOC_omicron_incr_severe_adeno1  = get_ve_incr_severe(vaccine_param$ve_VOC_omicron_infection_adeno1,vaccine_param$ve_VOC_omicron_severe_adeno1) # incremental reduction within dynamic model
  vaccine_param$ve_VOC_omicron_incr_severe_adeno2  = get_ve_incr_severe(vaccine_param$ve_VOC_omicron_infection_adeno2,vaccine_param$ve_VOC_omicron_severe_adeno2) # incremental reduction within dynamic model
  vaccine_param$ve_VOC_omicron_incr_severe_rna1    = get_ve_incr_severe(vaccine_param$ve_VOC_omicron_infection_rna1,vaccine_param$ve_VOC_omicron_severe_rna1) # incremental reduction within dynamic model
  vaccine_param$ve_VOC_omicron_incr_severe_rna2    = get_ve_incr_severe(vaccine_param$ve_VOC_omicron_infection_rna2,vaccine_param$ve_VOC_omicron_severe_rna2) # incremental reduction within dynamic model
  vaccine_param$ve_VOC_omicron_incr_severe_waning = get_ve_incr_severe(vaccine_param$ve_VOC_omicron_infection_waning,vaccine_param$ve_VOC_omicron_severe_waning) # incremental reduction within dynamic model
  vaccine_param$ve_VOC_omicron_incr_severe_booster = get_ve_incr_severe(vaccine_param$ve_VOC_omicron_infection_booster,vaccine_param$ve_VOC_omicron_severe_booster) # incremental reduction within dynamic model
  vaccine_param$ve_VOC_omicron_incr_severe_booster_waning = get_ve_incr_severe(vaccine_param$ve_VOC_omicron_infection_booster_waning,vaccine_param$ve_VOC_omicron_severe_booster_waning) # incremental reduction within dynamic model
  
  #preview
  flag_infection <- grepl('infection',names(vaccine_param))
  flag_severe    <- grepl('severe',names(vaccine_param)) & !grepl('incr',names(vaccine_param))
  flag_incr      <- grepl('incr',names(vaccine_param))
  
  flag_delta   <- grepl('delta',names(vaccine_param))
  flag_omicron <- grepl('omicron',names(vaccine_param))
  flag_wt      <- !(flag_delta |  flag_omicron)
  
  vaccine_param_unlist <- unlist(vaccine_param)  
  
  print(cbind(alpha_infection = vaccine_param_unlist[flag_infection & flag_wt],
              alpha_severe = vaccine_param_unlist[flag_severe & flag_wt],
              
              delta_infection = vaccine_param_unlist[flag_infection & flag_delta],
              delta_severe  = vaccine_param_unlist[flag_severe & flag_delta],
              
              omicron_infection = vaccine_param_unlist[flag_infection & flag_omicron],
              omicron_severe = vaccine_param_unlist[flag_severe & flag_omicron])
  )
  
  
  #extrabooster (all similar parameters for the moment)
  parms_names <- names(vaccine_param)
  sel_params <- vaccine_param[grepl('_infection_booster',parms_names)]
  vaccine_param[gsub('_infection_booster','_infection_extrabooster',names(sel_params))] <- sel_params
  
  #additional parameters to save
  vaccine_param$ve_waning_immunity_rate     <- log(2)/(6*30) #median time = 6 months
  vaccine_param$ve_waning_booster_rate <- log(2)/(6*30) #median time = 6 months
  vaccine_param$ve_waning_infection_rate <- log(2)/(6*30) #median time = 6 months
  vaccine_param$ve_waning_infection_booster_rate <- log(2)/(6*30) #median time = 6 months
  
  #other omicron variant
  parms_names <- names(vaccine_param)
  sel_params <- vaccine_param[grepl('ve_VOC_omicron',parms_names)]
  vaccine_param[gsub('ve_VOC_omicron','ve_VOC_ba4ba5',names(sel_params))] <- sel_params
  vaccine_param[gsub('ve_VOC_omicron','ve_VOC_bq1ba275xbb',names(sel_params))] <- sel_params
  vaccine_param[gsub('ve_VOC_omicron','ve_VOC_xbb15',names(sel_params))] <- sel_params

  #waning natural and hybrid immunity (similar as waning 2 doses/booster)
  parms_names <- names(vaccine_param)
  sel_params <- vaccine_param[grepl('infection_waning',parms_names)]
  vaccine_param[gsub('infection_waning','infection_reinf',names(sel_params))] <- sel_params
  parms_names <- names(vaccine_param)
  sel_params <- vaccine_param[grepl('incr_severe_waning',parms_names)]
  vaccine_param[gsub('incr_severe_waning','incr_severe_reinf',names(sel_params))] <- sel_params
  parms_names <- names(vaccine_param)
  sel_params <- vaccine_param[grepl('infection_booster_waning',parms_names)]
  vaccine_param[gsub('infection_booster_waning','infection_reinfvac',names(sel_params))] <- sel_params
  parms_names <- names(vaccine_param)
  sel_params <- vaccine_param[grepl('incr_severe_booster_waning',parms_names)]
  vaccine_param[gsub('incr_severe_booster_waning','incr_severe_reinfvac',names(sel_params))] <- sel_params
  

  write.table(vaccine_param,
             file='output/vaccine_parameters_20230820HubR5_opti.csv',
             sep=',',
             col.names=T,
             row.names=F)
  

}
#save_vaccine_parameters()

# get_colnames(parms_names = parms_names, 
#              tag_list = c(paste0('coef_w',1:6,'_'),
#                paste0('mat_w',2:5,'_'),
#                'nu2',
#                'VOC'),
#              ignore_list ='VOC_delta'
#              )


explore_param <- function(chains_param_file){
  
  # read file  
  parms_chains     = read.table(chains_param_file, sep = ",", header = T)
  
  # set parameter names
  parms_names <- names(parms_chains)
  
  # plot all parameters
  pdf('output/param_all.pdf',9,9)
  par(mfrow=c(5,5))
  for(i_param in 1:ncol(parms_chains)){
    
    boxplot(parms_chains[,i_param],
            main = parms_names[i_param])
  }
  dev.off()
  
  # identify all parameters that are fixed
  parms_min <- apply(parms_chains,2,min)
  parms_max <- apply(parms_chains,2,max)
  print(which(parms_min == parms_max))
}


reuse_parameter_sets <- function(filename_param_chain_base,
                                 filename_param_chain_source = NA,
                                 parms_to_adopt_list,
                                 rows_to_adopt,
                                 merge_tag = NA){
  
  # load parameters
  param_chain_base <- read.table(filename_param_chain_base,sep = ',',header = T)
  dim(param_chain_base)
  
  # load parameters to adopt
  if(!is.na(filename_param_chain_source)){
    param_chain_source <- read.table(filename_param_chain_source,sep = ',',header = T)
  } else {
    param_chain_source <- param_chain_base
  }
  
  # get parameter names
  parms_names <- names(param_chain_base)
  
  parms_to_adopt <- get_colnames(parms_names = parms_names,
                                tag_list = parms_to_adopt_list)
  
  # select all columns of the given parameters
  col_parms_to_adopt <- vector(length=length(parms_names))
  for(i_parm in parms_to_adopt){
    col_parms_to_adopt <- col_parms_to_adopt | grepl(i_parm,parms_names)
  }
  table(col_parms_to_adopt)
  
  row_params_to_adopt <- nrow(param_chain_source) - rows_to_adopt
  param_chain_source[row_params_to_adopt,col_parms_to_adopt]
  
  table(is.na(param_chain_source[row_params_to_adopt,col_parms_to_adopt]))
  
  # adopt
  i_param <- parms_names[col_parms_to_adopt][1]
  for(i_param in parms_names[col_parms_to_adopt]){
    param_chain_base[,i_param] <- sample(param_chain_source[row_params_to_adopt,i_param],nrow(param_chain_source),replace = T)
  }
  
  table(is.na(param_chain_base[,col_parms_to_adopt]))
  
  param_chain_base[row_params_to_adopt,'VOC_omicron_start']
  table(param_chain_base$VOC_omicron_start)
  #param_chain_base$VOC_omicron_start <- sample(624:628)
  
  # get new filename and save parameter values
  if(is.na(merge_tag)){
    merge_tag <- paste(parms_to_adopt,collapse='_')
  }
  
  
  filename_merge <- gsub('\\.',paste0('_',merge_tag,'.'),filename_param_chain_base)
  filename_merge <- gsub(dirname(filename_merge),paste0(dirname(filename_merge),'_',merge_tag),filename_merge)
  dir.create(dirname(filename_merge),showWarnings = F)
  write.table(param_chain_base,filename_merge,sep = ',',row.names = F,col.names = T)
  print(filename_merge)
  
  # return filename
  return(filename_merge)
}

explore_disease_history_parameters <- function(chains_param_file,n_chains=NA,sel_chains = NA){
  
  # load parameter file
  cat(fill = T)
  cat(chains_param_file,fill=T)
  parms_chains     = read.table(chains_param_file, sep = ",", header = T)
  
  if(all(!is.na(c(n_chains,sel_chains)))){
    warning('SET n_chains OR sel_chains')
    return(NULL)
  }
  
  if(!is.na(n_chains)){
    parms_chains <- parms_chains[1:n_chains,]
  }
  
  if(any(!is.na(sel_chains))){
    if(('mcmc_iter_id' %in% names(parms_chains))){
      parms_chains <- parms_chains[(parms_chains$mcmc_iter_id == max(parms_chains$mcmc_iter_id) &
                                      parms_chains$mcmc_chain_id %in% sel_chains),]
    } else{
      warning('MISSING: mcmc_iter_id')
      return(NULL)
    }

  }
  
  # get transition parameters
  gamma  = unlist(exp(parms_chains['log_gamma']));   # average length of the latency period - 2 days
  theta  = unlist(exp(parms_chains['log_theta']));   # infectious period at pre-symptomatic stage - 3.5 days
  delta1 = unlist(exp(parms_chains['log_delta1']));  # infectious period at asymptomatic stage - 3.5 days
  delta2 = unlist(exp(parms_chains['log_delta2']));  # infectious period at (mild) symptomatic stage - 3.5 days
  omega = colMeans(exp(parms_chains[grepl('log_omega',names(parms_chains))]));  # waiting time between symptom onset and hospitalization
  
  cat('exposed period:',mean(1/gamma),fill = T) # average length of the latency period
  cat('pre-sympt infectious period:',mean(1/theta),fill = T) # infectious period at pre-symptomatic stage
  cat('asympt infectious  period:',mean(1/delta1),fill = T) # infectious period at asymptomatic stage
  cat('mild sympt infectious period:',mean(1/delta2),fill = T) # infectious period at (mild) symptomatic stage
  cat('severe sympt infectious period before hosp admission:',(1/omega)) # waiting time between severe symptoms and hospitalization
  cat(fill = T)
  
  omicron_gamma_factor <- colMeans(exp(parms_chains['log_VOC_omicron_gamma_factor'])); 
  ba4ba5_gamma_factor <- colMeans(exp(parms_chains['log_VOC_ba4ba5_gamma_factor'])); 
  bq1ba275xbb_gamma_factor <- colMeans(exp(parms_chains['log_VOC_bq1ba275xbb_gamma_factor'])); 
  xbb15_gamma_factor <- colMeans(exp(parms_chains['log_VOC_xbb15_gamma_factor'])); 
  # omicron_transm       <- colMeans(exp(parms_chains['log_VOC_omicron_transm'])); 
  # delta_transm         <- colMeans(exp(parms_chains['log_VOC_delta_transm'])); 
  # alpha_transm         <- colMeans(exp(parms_chains['log_VOC_transm'])); 
  #print(omicron_gamma_factor)
  
  # # check VOC increased transmission
  alpha_transm     <- 100*mean(exp(parms_chains[,'log_VOC_alpha_transm'])-1);       # additional transmissibility
  alpha_transm_CrI <-100*quantile(exp(parms_chains[,'log_VOC_alpha_transm'])-1,c(0.025,0.975));       # additional transmissibility
  delta_transm <- 100*mean(exp(parms_chains[,'log_VOC_delta_transm'])/exp(parms_chains[,'log_VOC_alpha_transm'])-1);       # additional transmissibility
  delta_transm_CrI <- 100*quantile(exp(parms_chains[,'log_VOC_delta_transm'])/exp(parms_chains[,'log_VOC_alpha_transm'])-1,c(0.025,0.975));       # additional transmissibility
  omicron_transm <- 100*mean(exp(parms_chains[,'log_VOC_omicron_transm'])/exp(parms_chains[,'log_VOC_delta_transm'])-1);       # additional transmissibility
  omicron_transm_CrI <- 100*quantile(exp(parms_chains[,'log_VOC_omicron_transm'])/exp(parms_chains[,'log_VOC_delta_transm'])-1,c(0.025,0.975));       # additional transmissibility
  ba4ba5_transm <- 100*mean(exp(parms_chains[,'log_VOC_ba4ba5_transm'])/exp(parms_chains[,'log_VOC_omicron_transm'])-1);       # additional transmissibility
  ba4ba5_transm_CrI <- 100*quantile(exp(parms_chains[,'log_VOC_ba4ba5_transm'])/exp(parms_chains[,'log_VOC_omicron_transm'])-1,c(0.025,0.975),na.rm=T);       # additional transmissibility
  bq1ba275xbb_transm <- 100*mean(exp(parms_chains[,'log_VOC_bq1ba275xbb_transm'])/exp(parms_chains[,'log_VOC_ba4ba5_transm'])-1);       # additional transmissibility
  bq1ba275xbb_transm_CrI <- 100*quantile(exp(parms_chains[,'log_VOC_bq1ba275xbb_transm'])/exp(parms_chains[,'log_VOC_ba4ba5_transm'])-1,c(0.025,0.975),na.rm=T);       # additional transmissibility
  xbb15_transm <- 100*mean(exp(parms_chains[,'log_VOC_xbb15_transm'])/exp(parms_chains[,'log_VOC_bq1ba275xbb_transm'])-1);       # additional transmissibility
  xbb15_transm_CrI <- 100*quantile(exp(parms_chains[,'log_VOC_xbb15_transm'])/exp(parms_chains[,'log_VOC_bq1ba275xbb_transm'])-1,c(0.025,0.975),na.rm=T);       # additional transmissibility
  
  omicron_vereduction <- 100*mean(exp(parms_chains[,'ve_VOC_omicron_vereduction']));       
  omicron_vereduction_CrI <- 100*quantile(exp(parms_chains[,'ve_VOC_omicron_vereduction']),c(0.025,0.975),na.rm=T);      
  ba4ba5_vereduction <- 100*mean(exp(parms_chains[,'ve_VOC_ba4ba5_vereduction']));       
  ba4ba5_vereduction_CrI <- 100*quantile(exp(parms_chains[,'ve_VOC_ba4ba5_vereduction']),c(0.025,0.975),na.rm=T);       
  bq1ba275xbb_vereduction <- 100*mean(exp(parms_chains[,'ve_VOC_bq1ba275xbb_vereduction']));      
  bq1ba275xbb_vereduction_CrI <- 100*quantile(exp(parms_chains[,'ve_VOC_bq1ba275xbb_vereduction']),c(0.025,0.975),na.rm=T);       
  xbb15_vereduction <- 100*mean(exp(parms_chains[,'ve_VOC_xbb15_vereduction']));      
  xbb15_vereduction_CrI <- 100*quantile(exp(parms_chains[,'ve_VOC_xbb15_vereduction']),c(0.025,0.975),na.rm=T);       
  
  omicron_start        <- colMeans(parms_chains['VOC_omicron_start'])
  omicron_start_range  <- range(parms_chains['VOC_omicron_start'])
  omicron_init         <- colMeans(exp(parms_chains['log_VOC_omicron_init']))
  omicron_init_range   <- range(exp(parms_chains['log_VOC_omicron_init']))
  
  ba4ba5_start        <- colMeans(parms_chains['VOC_ba4ba5_start'])
  ba4ba5_start_range  <- range(parms_chains['VOC_ba4ba5_start'])
  ba4ba5_init        <- colMeans(exp(parms_chains['log_VOC_ba4ba5_init']))
  ba4ba5_init_range  <- range(exp(parms_chains['log_VOC_ba4ba5_init']))
  
  cat(fill = T)
  cat('Alpha VOC transmission advantage (wrt Wuhan):\t',alpha_transm,'[',alpha_transm_CrI,']',fill = T)
  cat('Delta VOC transmission advantage (wrt Alpha):\t',delta_transm,'[',delta_transm_CrI,']',fill = T)
  cat('Omicron VOC transmission advantage (wrt Delta):\t',omicron_transm,'[',omicron_transm_CrI,']',fill = T)
  cat('BA4BA5 VOC transmission advantage (wrt Omicron):',ba4ba5_transm,'[',ba4ba5_transm_CrI,']',fill = T)
  cat('bq1ba275xbb VOC transmission advantage (wrt BA4BA5):',bq1ba275xbb_transm,'[',bq1ba275xbb_transm_CrI,']',fill = T)
  cat('xbb15 VOC transmission advantage (wrt bq1ba275xbb):',xbb15_transm,'[',xbb15_transm_CrI,']',fill = T)
  cat('omicron VOC immune evasion (vs initial VE or protection from infection):',omicron_vereduction,'[',omicron_vereduction_CrI,']',fill = T)
  cat('ba4ba5 VOC immune evasion (vs initial VE or protection from infection):',ba4ba5_vereduction,'[',ba4ba5_vereduction_CrI,']',fill = T)
  cat('bq1ba275xbb VOC immune evasion (vs initial VE or protection from infection):',bq1ba275xbb_vereduction,'[',bq1ba275xbb_vereduction_CrI,']',fill = T)
  cat('xbb15 VOC immune evasion (vs initial VE or protection from infection):',xbb15_vereduction,'[',xbb15_vereduction_CrI,']',fill = T)
  
  cat(fill = T) 

    
  cat('Omicron exposed period (days):\t',mean(1/(gamma*omicron_gamma_factor)),fill = T) # average length of the latency period
  cat('ba4ba5 exposed period (days):\t',mean(1/(gamma*ba4ba5_gamma_factor)),fill = T) # average length of the latency period
  cat('bq1ba275xbb exposed period (days):\t',mean(1/(gamma*bq1ba275xbb_gamma_factor)),fill = T) # average length of the latency period
  cat('xbb15 exposed period (days):\t',mean(1/(gamma*xbb15_gamma_factor)),fill = T) # average length of the latency period
  cat(fill = T)
  
  cat('Omicron start:\t',format(sim_day2date(round(omicron_start))),'[',format(sim_day2date(round(omicron_start_range))),']',fill = T) # average length of the latency period
  cat('BA4BA5 start:\t',format(sim_day2date(round(ba4ba5_start))),'[',format(sim_day2date(round(ba4ba5_start_range))),']',fill = T) # average length of the latency period
  
  cat('Omicron init:\t',round(omicron_init),'[',round(omicron_init_range),']',fill = T) # average length of the latency period
  cat('BA4BA5 init:\t',round(ba4ba5_init),'[',round(ba4ba5_init_range),']',fill = T) # average length of the latency period
  
    # tbl <- table(parms_chains['VOC_omicron_start'])
  # tbl <- data.frame(date =sim_day2date(as.numeric(names(tbl))),
  #       count = (tbl))
  # print(tbl,fill = T)
  # 
  # i_sel <- 1:nrow(parms_chains)
  
  # # check VOC part 1: increased transmission
  # 100*mean(exp(parms_chains[,'log_VOC_alpha_transm'])-1);       # additional transmissibility
  # 100*quantile(exp(parms_chains[,'log_VOC_alpha_transm'])-1,c(0.025,0.975));       # additional transmissibility
  # 
  # # check VOC part 2: increased transmission
  # 100*mean(exp(parms_chains[,'log_VOC_delta_transm'])-1);       # additional transmissibility
  # 100*quantile(exp(parms_chains[,'log_VOC_delta_transm'])-1,c(0.025,0.975));       # additional transmissibility
  # 100*mean(exp(log(exp(parms_chains[,'log_VOC_alpha_transm'])*1.65))-1)
  # 100*mean(exp(log(exp(parms_chains[,'log_VOC_alpha_transm'])*1.8))-1)
  # 
  # # check VOC part 3: increased transmission
  # print(100*mean(exp(parms_chains[i_sel,'log_VOC_omicron_transm'])-1));       # additional transmissibility
  # print(100*quantile(exp(parms_chains[i_sel,'log_VOC_omicron_transm'])-1,c(0.025,0.975)));       # additional transmissibility
  # 
  # # check VOC part 4: increased transmission
  # print(100*mean(exp(parms_chains[i_sel,'log_VOC_ba4ba5_transm'])-1));       # additional transmissibility
  # print(100*quantile(exp(parms_chains[i_sel,'log_VOC_ba4ba5_transm'])-1,c(0.025,0.975),na.rm=T));       # additional transmissibility
  # 
  # print(100*mean(exp(log(exp(parms_chains[i_sel,'log_VOC_delta_transm'])*1.3))-1))
  # print(100*mean(exp(log(exp(parms_chains[i_sel,'log_VOC_delta_transm'])*1.4))-1))
  # print(100*mean(exp(log(exp(parms_chains[i_sel,'log_VOC_delta_transm'])*1.5))-1))
  # print(100*mean(exp(log(exp(parms_chains[i_sel,'log_VOC_delta_transm'])*1.6))-1))
  # print(100*mean(exp(log(exp(parms_chains[i_sel,'log_VOC_delta_transm'])*1.7))-1))
  # print(100*mean(exp(log(exp(parms_chains[i_sel,'log_VOC_delta_transm'])*1.8))-1))
  # print(100*mean(exp(log(exp(parms_chains[i_sel,'log_VOC_delta_transm'])*1.9))-1))
  # print(100*mean(exp(log(exp(parms_chains[i_sel,'log_VOC_delta_transm'])*2.0))-1))
  # print(100*mean(exp(log(exp(parms_chains[i_sel,'log_VOC_delta_transm'])*2.1))-1))
  # print(100*mean(exp(log(exp(parms_chains[i_sel,'log_VOC_delta_transm'])*2.2))-1))
  # print(100*mean(exp(log(exp(parms_chains[i_sel,'log_VOC_delta_transm'])*2.3))-1))
  # print(100*mean(exp(log(exp(parms_chains[i_sel,'log_VOC_delta_transm'])*2.4))-1))

  # omicron_init <- unlist(exp(parms_chains['log_VOC_omicron_init']))
  # print(c(mean(omicron_init),range(omicron_init)))

}

explore_invalid_param <- function(chains_param_file,n_chains=60){
  
  # load parameter file
  cat(fill = T)
  cat(chains_param_file,fill=T)
  parms_chains     = read.table(chains_param_file, sep = ",", header = T)
  dim(parms_chains)
  
  # loop over last n config, and check ll
  parms_names <- names(parms_chains)
  i <- nrow(parms_chains)
  for(i in 1:n_chains){
    invalid_param <- log_prior_model(unlist(parms_chains[i,]),parms_names)$invalid_param
    if(!is.null(invalid_param)){
      cat(c(i,invalid_param),fill = T)
    }  
  }
}

explore_param_table <- function(chains_param_file,n_chains=60){
  
  # load parameter file
  cat(fill = T)
  cat(chains_param_file,fill=T)
  parms_chains     = read.table(chains_param_file, sep = ",", header = T)
  dim(parms_chains)
  
  # loop over last n config, and check ll
  parms_names <- names(parms_chains)
  
  for(i in 1:n_chains){
    invalid_param <- log_prior_model(unlist(parms_chains[i,]),parms_names)$invalid_param
    if(!is.null(invalid_param)){
      cat(c(i,invalid_param),fill = T)
    }  
  }
}

reformat_config_file <- function(chains_param_files){
  
  for(chains_param_file in chains_param_files){
    parms      <- read.csv(chains_param_file,sep=',',header=T)  
  
    if(!any(grepl('mcmc_',names(parms)))){
      parms <- parms[nrow(parms):1,]
      parms$mcmc_core_version <- 23
    
      new_filename <- gsub('.csv','_format.csv',chains_param_file)
      write.table(parms,new_filename,sep=',',row.names = F)  
      cat(new_filename)
    } else {
      
      cat('CONFIG FILE ALREADY FORMATTED FOR v24')
    }
  }
}

# additional code, not to excecute by default ----
if(0==1){
 
 # chains_param_files = dir('data/config/',pattern='MCMCmulti_20220406.*belgium_',full.names = T,recursive = T)
  # chains_param_files = dir('output/MCMC_vaughan_v23_recap/',pattern='MCMCmulti_20220404.*e50.*belgium_',full.names = T,recursive = T)
  #chains_param_files = dir('output/MCMC_vaughan_v24_reinfection/',pattern='166945.*.csv',full.names = T,recursive = T)
  #chains_param_files = dir('output/MCMC_vaughan_v24_ba4ba5/',pattern='191740.*.csv',full.names = T,recursive = T)
  chains_param_files = dir('data/config/',pattern='MCMCmulti_20230629_d1213_e40_i150_n10_p10_crit1_workingwave63_belgium_workingend',full.names = T,recursive = T)
  
  chains_param_file <- chains_param_files
  #reformat_config_file(chains_param_files)
  
  explore_disease_history_parameters(chains_param_files[1])
  
  n_chains <- 10
  explore_disease_history_parameters(chains_param_files[1],n_chains=n_chains)
  explore_disease_history_parameters(chains_param_files[1],sel_chains=c(1:3))
  
  explore_param(chains_param_files[1])
  
  explore_invalid_param(chains_param_files[1])
  
  for(i_file in chains_param_files)
    explore_disease_history_parameters(i_file)
  
  
  ## re-use wave 40 for 41:43
  chains_param_file <- chains_param_files[1]
  parms_chains    <- read.csv(chains_param_file,sep=',',header=T)
  parms_names     <- names(parms_chains)
  parms_names_rem <- get_colnames(parms_names = parms_names,
                                    tag_list = paste0('coef_w',42:43,'_')
                                  )
                                    
  parms_chains <- parms_chains[,!parms_names %in% parms_names_rem]
  #parms_chains <- include_new_wave(parms_chains,bool_comix = TRUE) # 41
  parms_chains <- include_new_wave(parms_chains,bool_comix = TRUE) # 42
  parms_chains <- include_new_wave(parms_chains,bool_comix = TRUE) # 43
  
  write.table(parms_chains,gsub('.csv','_w41.csv',chains_param_file),sep=',',row.names = F)
  
  
  # # re-use parameter estimations
  # chains_param_files = dir('output/MCMC_vaughan_v23_recap/',pattern='MCMCmulti_20220316.*e73.*belgium_',full.names = T,recursive = T)
  # chains_param_files
  # filename_param_chain_base <- chains_param_files[3]
  # filename_param_chain_source <- chains_param_files[2]
  # 
  # nwave_Comix <- 42
  # rows_to_adopt <- c(0:1)
  # parms_to_adopt_list = c(
  #                         'log_VOC_omicron_init',
  #                          'log_VOC_omicron_transm',
  #                          'log_VOC_omicron_gamma_factor',
  #                          paste0('coef_w',36:nwave_Comix,'_'))
  # merge_tag <- '0331_sel'
  # reuse_parameter_sets(filename_param_chain_base,
  #                      filename_param_chain_source = filename_param_chain_source,
  #                      parms_to_adopt_list,
  #                      rows_to_adopt = rows_to_adopt,
  #                      merge_tag = merge_tag)
  # 
  # 
  # 
  # 
  # filename_param_chain_base  <- 'output/MCMC_vaughan_v22_age/wave4/MCMCmulti_20211118_d627_e110_i400_n10_p10_crit1_63325_belgium.csv' 
  # filename_param_chain_new   <- 'output/MCMC_vaughan_v22_age/wave4_hload/ICU_1010/MCMC_64670_final/MCMCmulti_20211123_d632_e8_i1000_n10_p10_crit3_64670_belgium.csv' 
  # parms_to_adopt <- c('delta3','phi1_add')
  # adopt_parameters(filename_param_chain_base = filename_param_chain_base,
  #                  filename_param_chain_new = filename_param_chain_new,
  #                  parms_to_adopt = parms_to_adopt)
  # 
  # 
  # # remove redundant parameters
  # chains_param_files = dir('data/config/',pattern='MCMCmulti_20211123',full.names = T,recursive = F)
  # chains_param_file <- chains_param_files[5]
  # 
  # for(chains_param_file in chains_param_files){
  #   parms      <- read.csv(chains_param_file,sep=',',header=T)  
  #   sel_remove <- grepl('imported_cases_',names(parms))
  #   sel_remove <- sel_remove | grepl('veX',names(parms))
  #   sel_remove <- sel_remove | grepl('delta4',names(parms))
  #   
  #   names(parms)[sel_remove]
  #   parms <- parms[,!sel_remove]
  #   parms <- write.table(x = parms,
  #                        file=gsub('.csv','_clean.csv',chains_param_file),sep=',',row.names = FALSE)  
  # }
  # 
  
  # ## Clean parameters and sort by LL ----
  chains_param_files = dir('output/MCMC_vaughan_v23_recap/',pattern='20220410.*.csv',full.names = T,recursive = T)
  chains_param_files = dir('output/MCMC_vaughan_v23_recap/',pattern='130137.*.csv',full.names = T,recursive = T)
  chains_param_files = dir('output/MCMC_vaughan_v24_reinfection/',pattern='191520.*.csv',full.names = T,recursive = T)
  
  chains_param_file <- chains_param_files[1]
  explore_invalid_param(chains_param_file)

  parms_chains    <- read.csv(chains_param_file,sep=',',header=T)
  #parms_chains    <- parms_chains[nrow(parms_chains)-(0:59),]
  parms_names     <- names(parms_chains)

  parms_names[grepl('omicron',parms_names)]
  
  # tmp <- parms_chains[,'log_VOC_omicron_init']
  # tmp[tmp>6.9] <- 6.9
  # parms_chains[,'log_VOC_omicron_init'] <- tmp
  # 
  # parms_chains['log_delta3_omicron']  <- -2 # 2.3
  # parms_chains[names(parms_chains)[grepl('log_mu_omicron_sev',names(parms_chains))]] <- rep(c(-9.6,-6.9,-5.7,-4.7,-3.2,-2.6,-2,-0.4),each=nrow(parms_chains))

  for(i in 1:nrow(parms_chains)){

    invalid_param <- log_prior_model(unlist(parms_chains[i,]),parms_names)$invalid_param
    if(!is.null(invalid_param)){
      cat(c(i,invalid_param,unlist(parms_chains[i,invalid_param])),fill = T)
      if(!all(grepl('coef',invalid_param))){
        cat('not all parameters are CoMix related',fill = T)
        parms_chains[i,invalid_param[!grepl('coef',invalid_param)]] <- median(parms_chains[,invalid_param[!grepl('coef',invalid_param)]])
      }
      parms_chains[i,invalid_param[grepl('coef',invalid_param)]] <- 3
    }
  }

  write.table(parms_chains,gsub('.csv','_capped.csv',chains_param_file),sep=',',row.names = F)

  # replace q-parameters by the median value
  parms_chains_median <- parms_chains[1:60,]
  col_select <- get_colnames(names(parms_chains_median),tag_list = c('log_VOC_ba4ba5_init','log_VOC_ba4ba5_transm'))
  #col_select <- get_colnames(names(parms_chains_median),tag_list = 'mu_')
  row_select <- c(41,46,39,30,2,6,9,15,21,45,47)
  row_select <- 1:5
  
  i_col <- col_select[1]
  for(i_col in col_select){
    parms_chains_median[,i_col] <- median(parms_chains_median[row_select,i_col])
  }
  parms_chains_median[,col_select]
  
  write.table(parms_chains_median,gsub('.csv','_median.csv',chains_param_file),sep=',',row.names = F)
  
  parms_chains_median[row_select,] <- parms_chains_median[row_select,]
  write.table(parms_chains_median,gsub('.csv','_select.csv',chains_param_file),sep=',',row.names = F)
  
  parms_chains_median <- parms_chains[1:60,]
  for(i_col in col_select){
    parms_chains_median[,i_col] <- median(parms_chains_median[1:3,i_col])
  }
  
  write.table(parms_chains_median,gsub('.csv','_selMedian.csv',chains_param_file),sep=',',row.names = F)
  
  explore_invalid_param(gsub('.csv','_median.csv',chains_param_file))
  
}

