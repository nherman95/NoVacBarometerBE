########################################################################### #
# This file is part of the Stochastic Compartmental Model for SARS-COV-2 
# transmission in Belgium, conceived by members of SIMID group during the 
# COVID19 pandemic.
#
# This file contains many helper function to run the stochastic model, 
# to explore the results and estimate model parameters.
#
# Copyright 2021, SIMID, University of Antwerp & Hasselt University                                        
########################################################################### #

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(zoo)) #rollmean
suppressPackageStartupMessages(library(openxlsx)) 
suppressPackageStartupMessages(library(EpiEstim)) 
suppressPackageStartupMessages(library("lubridate"))
suppressPackageStartupMessages(library(dplyr))

# requires R (≥ 4.1.0)
if("EpiLPS" %in% installed.packages()){
  suppressPackageStartupMessages(library(EpiLPS)) 
}


############################################################################ #
# PROVIDE PUBLIC HEALTH AGENCY DATA                                         ##
############################################################################ #

get_observed_incidence_data <- function(enable_data_reuse = TRUE,sel_region = NA,bool_incr_child_uptake = FALSE)
{
  
  # use (local version of) most recent SCIENSANO data (or local backup version)
  ref_data_file_name <- file.path('data/downloads',paste0('covid19_reference_data_',gsub('-','',Sys.Date()),'.csv'))
  backup_file        <- file.path('data',paste0('covid19_reference_data.csv'))
  
  # use specified region or use the default?
  if(is.na(sel_region)){
    sel_region <- get_region(1)
  }
  
  # adjust file name to region
  ref_data_file_name <- gsub('.csv',paste0('_',sel_region,'.csv'),ref_data_file_name)
  
  # if reference file is present: return data
  if(file.exists(ref_data_file_name) && enable_data_reuse){
    ref_data      <- read.table(ref_data_file_name,sep=',',header=T)
    ref_data$date <- as.Date(ref_data$date)
    return(ref_data)
  }
  
  # else download files and generate reference file
  data_dir        <- dirname(ref_data_file_name)
  hosp_ref_file   <- download_ref_file('https://epistat.sciensano.be/Data/COVID19BE_HOSP.csv',data_dir)
  cases_ref_file  <- download_ref_file('https://epistat.sciensano.be/Data/COVID19BE_CASES_AGESEX.csv',data_dir)
  tests_ref_file  <- download_ref_file('https://epistat.sciensano.be/Data/COVID19BE_tests.csv',data_dir)
  uptake_ref_file <- download_ref_file('https://epistat.sciensano.be/Data/COVID19BE_VACC.csv',data_dir)
  mort_ref_file   <- download_ref_file('https://epistat.sciensano.be/Data/COVID19BE_MORT.csv',data_dir)
  
  # if any download failed => use defaults
  if(any(is.na(c(hosp_ref_file,cases_ref_file,tests_ref_file)))){
    return(read.table(backup_file,sep=',',header=T))
  }

  # load reference hospital data
  ref_hosp_data_all     <- read.table(hosp_ref_file,sep=',',header = T)
  ref_hosp_data_all     <- select_regional_data(ref_hosp_data_all,sel_region)
  if(!'NEW_OUT' %in%  names(ref_hosp_data_all)){
    ref_hosp_data_all$NEW_OUT<-0
  }
  hosp_adm_data         <- aggregate(. ~ DATE, data= ref_hosp_data_all[,c('NEW_IN','TOTAL_IN','TOTAL_IN_ICU','DATE','NEW_OUT')],sum,na.rm=T)
  names(hosp_adm_data)  <- c('date','hospital_admissions','hospital_load','icu_load','hospital_discharges')
  hosp_adm_data$cumulative_hospital_admissions <- cumsum(hosp_adm_data$hospital_admissions)
  hosp_adm_data$region  <- sel_region
  head(hosp_adm_data)
  
  # load reference case data
  ref_case_data_all        <- read.table(cases_ref_file,sep=',',header = T)
  ref_case_data_all        <- select_regional_data(ref_case_data_all,sel_region)
  ref_case_data_table      <- data.table(ref_case_data_all[!is.na(ref_case_data_all$AGEGROUP) & !is.na(ref_case_data_all$DATE),])
  ref_case_data_age        <- data.frame(get_age_table(ref_case_data_table))
  names(ref_case_data_age) <- gsub('\\.','_',names(ref_case_data_age))
  # remove 4 latest days (starting from today)
  # or 3 before latest hospital data point (on Sunday and Monday, the info is not updated)
  #ref_case_data_age        <- ref_case_data_age[as.Date(ref_case_data_age$date) < (Sys.Date()-3),]
  ref_case_data_age        <- ref_case_data_age[as.Date(ref_case_data_age$date) < (max(as.Date(hosp_adm_data$date))-2),]
  head(ref_case_data_age)
  tail(ref_case_data_age)
  
  # load test data
  ref_test_data_all <- read.table(tests_ref_file,sep=',',header = T)
  ref_test_data_all <- select_regional_data(ref_test_data_all,sel_region)
  ref_test_data     <- aggregate(TESTS_ALL ~ DATE, data= ref_test_data_all,sum,na.rm=T)
  names(ref_test_data)  <- c('date','covid19_tests')
  
  # load mortality data
  ref_mort_data_all <- read.table(mort_ref_file,sep=',',header = T)
  ref_mort_data_all <- select_regional_data(ref_mort_data_all,sel_region)
  ref_mort_data     <- aggregate(DEATHS ~ DATE, data= ref_mort_data_all,sum,na.rm=T)
  names(ref_mort_data)  <- c('date','covid19_deaths')
  
  # load mortality data by age
  mortality_by_age <- read.table("data/COVID19BE_MORT_RAW.csv", header = T, sep = ",")
  mortality_by_age <- select_regional_data(mortality_by_age,sel_region)
  mortality_by_age$GROUP <- findInterval(mortality_by_age$AGE, c(seq(0,90,10))) # last group is open
  mortality_by_age$DATE <- as.Date(mortality_by_age$DATE,'%m/%d/%Y')

  tbl_age_dist_mort    <- table(mortality_by_age$DATE, factor(mortality_by_age$GROUP, levels = 1:10))
  age_dist_mort        <- data.frame(matrix(tbl_age_dist_mort, ncol = 10))
  names(age_dist_mort) <- paste0('covid19_deaths_age',1:ncol(age_dist_mort))
  age_dist_mort$date   <- row.names(tbl_age_dist_mort)
  ref_mort_data        <- merge(ref_mort_data,age_dist_mort,all.x = T,all.y = T)
  
  # mortality "hospital" (non nursing homes)
  mortality_by_age_hosp <- subset(mortality_by_age, PLACE != "Nursing home") # select hospital based mortality (?)
  # mortality_by_age$covid19_deaths_hospital <- 1
  # ref_mort_data_hospital                   <- aggregate(covid19_deaths_hospital ~ DATE, data= mortality_by_age,sum,na.rm=T)
  # names(ref_mort_data_hospital)[1]         <- 'date'
  # ref_mort_data_hospital$date              <- as.character(as.Date(ref_mort_data_hospital$date,'%m/%d/%Y'))
  # 
  tbl_age_dist_mort    <- table(mortality_by_age_hosp$DATE, factor(mortality_by_age_hosp$GROUP, levels = 1:10))
  age_dist_mort        <- data.frame(matrix(tbl_age_dist_mort, ncol = 10))
  names(age_dist_mort) <- paste0('covid19_deaths_hospital_age',1:ncol(age_dist_mort))
  age_dist_mort$date   <- row.names(tbl_age_dist_mort)
  age_dist_mort$covid19_deaths_hospital <- rowSums(age_dist_mort[,-11],na.rm=T)
  ref_mort_data        <- merge(ref_mort_data,age_dist_mort,all.x = T,all.y = T)
  
  # remove last 7 observations (and make sure that all dates in between are there)
  flag_last_obs                          <- as.Date(ref_mort_data$date) %in% as.Date(ref_mort_data$date[1:(nrow(age_dist_mort)-7)])
  ref_mort_data[!flag_last_obs,-1] <- NA
  
  # load vaccine uptake data
  ref_uptake_data_bel <- read.table(uptake_ref_file,sep=',',header = T)
  ref_uptake_data_all <- select_regional_data(ref_uptake_data_bel,sel_region)
  
  # make sure all possible dates with reported uptake are included
  if(max(ref_uptake_data_all$DATE)<max(ref_uptake_data_bel$DATE)){
    ref_uptake_data_extra <- ref_uptake_data_bel[!ref_uptake_data_bel$DATE %in% ref_uptake_data_all$DATE,]
    ref_uptake_data_extra$REGION <- unique(ref_uptake_data_all$REGION)
    ref_uptake_data_extra$COUNT  <- 0
    ref_uptake_data_all <- rbind(ref_uptake_data_all,
                                 ref_uptake_data_extra)
  }
  
  # # remove "E2"xtra doses of vaccine
 # ref_uptake_data_all <- ref_uptake_data_all[ref_uptake_data_all$DOSE != "E2",]
  table(ref_uptake_data_all$DOSE)

  # #temp add for child vaccine (Nicolas) !!!!
  if(bool_incr_child_uptake){
    ref_uptake_data_all <- ref_uptake_data_all[ref_uptake_data_all$AGEGROUP != "05-11",]
    addchild <- ref_uptake_data_all[ref_uptake_data_all$AGEGROUP=="12-15",]
    addchild$AGEGROUP <- "05-11"
    addchild$COUNT <- addchild$COUNT * 7/4 # 12-15y contains uptake data on 4 age bins, while we target now 5-11y (=7 age bins)
    ref_uptake_data_all <- rbind(ref_uptake_data_all,addchild)
    # addchild <- ref_uptake_data_all[ref_uptake_data_all$AGEGROUP=="16-17",]
    # addchild$AGEGROUP <- "12-15"
    # ref_uptake_data_all <- rbind(ref_uptake_data_all,addchild)
  }
  
  #Separate dedicated bivalent booster
  ref_uptake_data_all[ref_uptake_data_all$BRAND %in% c('Pfizer-BioNTech Bivalent BA1','Pfizer-BioNTech Bivalent BA4-5'),]$DOSE <- 'X'
 
  ref_uptake_data     <- aggregate(COUNT ~ DATE + DOSE, data= ref_uptake_data_all,sum,na.rm=T)
  names(ref_uptake_data) <- c('date','vaccine_dose','vaccine_uptake')

  ref_uptake_data_AC <- aggregate(vaccine_uptake ~ date, data = ref_uptake_data[ref_uptake_data$vaccine_dose %in% c('A','C'),],sum)
  ref_uptake_data_BC <- aggregate(vaccine_uptake ~ date, data = ref_uptake_data[ref_uptake_data$vaccine_dose %in% c('B','C'),],sum)
  ref_uptake_data_E  <- aggregate(vaccine_uptake ~ date, data = ref_uptake_data[ref_uptake_data$vaccine_dose %in% c('E'),],sum)
  ref_uptake_data_F  <- aggregate(vaccine_uptake ~ date, data = ref_uptake_data[ref_uptake_data$vaccine_dose %in% c('E2','E3'),],sum)
  ref_uptake_data_X  <- aggregate(vaccine_uptake ~ date, data = ref_uptake_data[ref_uptake_data$vaccine_dose %in% c('X'),],sum) #only bivalent booster
  
  ref_uptake_dose <- merge(ref_uptake_data_AC,
                           ref_uptake_data[ref_uptake_data$vaccine_dose == 'B',-2],
                           by='date',
                           all = T,
                           suffixes = c("_first", "_second"))
  ref_uptake_dose <- merge(ref_uptake_dose,
                           ref_uptake_data_BC,
                           by='date',
                           all = T)
  ref_uptake_dose <- merge(ref_uptake_dose,
                           ref_uptake_data_E,
                           by='date',
                           all = T,
                           suffixes = c("_full", "_booster"))
  ref_uptake_dose <- merge(ref_uptake_dose,
                           ref_uptake_data_F,
                           by='date',
                           all = T,
                           suffixes = c("_2ndfull", "_2ndbooster"))
  ref_uptake_dose <- merge(ref_uptake_dose,
                           ref_uptake_data_X,
                           by='date',
                           all = T,
                           suffixes = c("_extrafull", "_extrabooster"))
  # names(ref_uptake_dose)[4] <- paste0(names(ref_uptake_dose)[4],'_full')
  colSums(ref_uptake_dose[,-1],na.rm=T)
  
  # rename vaccine type
  table(ref_uptake_data_all$BRAND,useNA = 'ifany')
  ref_uptake_data_all$BRAND[ref_uptake_data_all$BRAND %in% c('Moderna Original','Pfizer-BioNTech Original','Novavax','Other','Pfizer-BioNTech Bivalent BA1')] <- 'mRNA'   # merge mRNA vaccines
  ref_uptake_data_all$BRAND[ref_uptake_data_all$BRAND %in% 'Johnson&Johnson'] <- 'JnJ'  
  ref_uptake_data_all$BRAND[ref_uptake_data_all$BRAND %in% 'AstraZeneca-Oxford'] <- 'AZ' 
 # ref_uptake_data_all$BRAND[ref_uptake_data_all$BRAND %in% 'Pfizer-BioNTech Bivalent BA1'] <- 'BIV' 
  
    
  # load vaccine uptake by dose and type
  ref_uptake_data     <- aggregate(COUNT ~ DATE + BRAND + DOSE, data= ref_uptake_data_all,sum,na.rm=T)
  names(ref_uptake_data) <- c('date','vaccine_type','vaccine_dose','vaccine_uptake')
  
  flag_AZ  <- grepl('AZ',ref_uptake_data$vaccine_type)
  flag_JnJ <- grepl('JnJ',ref_uptake_data$vaccine_type)
  flag_Vx  <- grepl('mRNA',ref_uptake_data$vaccine_type)
 # flag_BIV  <- grepl('BIV',ref_uptake_data$vaccine_type)
  
  ref_uptake_Vx <- merge(ref_uptake_data[flag_Vx & ref_uptake_data$vaccine_dose == 'A',c(1,4)],
                         ref_uptake_data[flag_Vx & ref_uptake_data$vaccine_dose == 'B',c(1,4)],
                         by='date',
                         all.x = T,
                         suffixes = c("_mRNA_first", "_mRNA_second"))
  
  ref_uptake_AZ <- merge(ref_uptake_data[flag_AZ & ref_uptake_data$vaccine_dose == 'A',c(1,4)],
                         ref_uptake_data[flag_AZ & ref_uptake_data$vaccine_dose == 'B',c(1,4)],
                         by='date',
                         all.x = T,
                         suffixes = c("_AZ_first", "_AZ_second"))
  
  ref_uptake_JnJ <- ref_uptake_data[flag_JnJ & ref_uptake_data$vaccine_dose == 'C',c(1,4)]
  names(ref_uptake_JnJ)[-1] <- paste0(names(ref_uptake_JnJ)[-1],'_JnJ_first') 
  
  ref_uptake_booster <- ref_uptake_data[ref_uptake_data$vaccine_dose == 'E',c(1,4)]
  names(ref_uptake_booster)[2] <- paste0(names(ref_uptake_booster)[2],'_mRNA_booster')
  
  ref_uptake_2ndbooster <- ref_uptake_data[ref_uptake_data$vaccine_dose == 'E2',c(1,4)]
  names(ref_uptake_2ndbooster)[2] <- paste0(names(ref_uptake_2ndbooster)[2],'_mRNA_2ndbooster')
  
  ref_uptake_extrabooster <- ref_uptake_data[ref_uptake_data$vaccine_dose == 'X',c(1,4)]
  names(ref_uptake_2ndbooster)[2] <- paste0(names(ref_uptake_2ndbooster)[2],'_mRNA_extrabooster')
  
  ref_uptake_type <- merge(ref_uptake_AZ,ref_uptake_Vx, all=T)
  ref_uptake_type <- merge(ref_uptake_type, ref_uptake_JnJ, all=T)
  ref_uptake_type <- merge(ref_uptake_type, ref_uptake_booster, all=T)
  ref_uptake_type <- merge(ref_uptake_type, ref_uptake_2ndbooster, all=T)
  ref_uptake_type <- merge(ref_uptake_type, ref_uptake_extrabooster, all=T)
  ref_uptake_type[is.na(ref_uptake_type)] <- 0
  
  # include vaccine uptake by age group ---
  ref_uptake_aggr      <- convert_uptake_age_groups(ref_uptake_data_all,sel_region)
  ref_uptake_aggr_mRNA <- convert_uptake_age_groups(ref_uptake_data_all,sel_region,sel_brand='mRNA')
  ref_uptake_aggr_AZ   <- convert_uptake_age_groups(ref_uptake_data_all,sel_region,sel_brand='AZ')
  ref_uptake_aggr_JnJ  <- convert_uptake_age_groups(ref_uptake_data_all,sel_region,sel_brand='JnJ')
  #ref_uptake_aggr_BIV  <- convert_uptake_age_groups(ref_uptake_data_all,sel_region,sel_brand='BIV')
  
    
  names(ref_uptake_aggr_mRNA) <- gsub('vaccine_','vaccine_mRNA_',names(ref_uptake_aggr_mRNA))
  names(ref_uptake_aggr_AZ)   <- gsub('vaccine_','vaccine_AZ_',names(ref_uptake_aggr_AZ))
  names(ref_uptake_aggr_JnJ)  <- gsub('vaccine_','vaccine_JnJ_',names(ref_uptake_aggr_JnJ))
#  names(ref_uptake_aggr_BIV)  <- gsub('vaccine_','vaccine_BIV_',names(ref_uptake_aggr_BIV))
  ref_uptake_aggr_JnJ         <- ref_uptake_aggr_JnJ[,grepl('_A_',names(ref_uptake_aggr_JnJ)) | grepl('date',names(ref_uptake_aggr_JnJ))]
  ref_uptake_aggr_AZ          <- ref_uptake_aggr_AZ[,!(grepl('_E_',names(ref_uptake_aggr_AZ)) | grepl('_E2_',names(ref_uptake_aggr_AZ)) | grepl('_E3_',names(ref_uptake_aggr_AZ)))]
  
  # store original age groups as well
  ref_uptake_age <- aggregate(COUNT ~ + DOSE + AGEGROUP + DATE, data=ref_uptake_data_all,sum)
  ref_uptake_age$DATE <- as.Date(ref_uptake_age$DATE)
  ref_uptake_age <- ref_uptake_age[order(ref_uptake_age$DATE),]
  
  ref_uptake_age$AGEGROUP <- gsub('\\+','-99',ref_uptake_age$AGEGROUP)
  opt_agegroup_epistat    <- sort(unique(ref_uptake_age$AGEGROUP))
  ref_uptake_age$age_factor <- as.numeric(factor(ref_uptake_age$AGEGROUP,labels = opt_agegroup_epistat))
  date_opt            <- unique(ref_uptake_age$DATE)
  AG_uptake_matrix_A <- matrix(0,nrow=length(date_opt),ncol=length(opt_agegroup_epistat))
  AG_uptake_matrix_B <- AG_uptake_matrix_A
  AG_uptake_matrix_C <- AG_uptake_matrix_A
  AG_uptake_matrix_E <- AG_uptake_matrix_A
  AG_uptake_matrix_F <- AG_uptake_matrix_A
  AG_uptake_matrix_X <- AG_uptake_matrix_A
  
  i_row <- 1
  for(i_row in 1:nrow(ref_uptake_age)){
      if(ref_uptake_age$DOSE[i_row] == 'A'){
        AG_uptake_matrix_A[ref_uptake_age$DATE[i_row] == date_opt,ref_uptake_age$age_factor[i_row]] <- ref_uptake_age$COUNT[i_row]
    } else if(ref_uptake_age$DOSE[i_row] == 'B'){
        AG_uptake_matrix_B[ref_uptake_age$DATE[i_row] == date_opt,ref_uptake_age$age_factor[i_row]] <- ref_uptake_age$COUNT[i_row]
    } else if(ref_uptake_age$DOSE[i_row] == 'C'){
      AG_uptake_matrix_C[ref_uptake_age$DATE[i_row] == date_opt,ref_uptake_age$age_factor[i_row]] <- ref_uptake_age$COUNT[i_row]
    } else if(ref_uptake_age$DOSE[i_row] == 'E'){
      AG_uptake_matrix_E[ref_uptake_age$DATE[i_row] == date_opt,ref_uptake_age$age_factor[i_row]] <- ref_uptake_age$COUNT[i_row]
    } else if(ref_uptake_age$DOSE[i_row] == 'E2' || ref_uptake_age$DOSE[i_row] == 'E3'){
      AG_uptake_matrix_F[ref_uptake_age$DATE[i_row] == date_opt,ref_uptake_age$age_factor[i_row]] <- ref_uptake_age$COUNT[i_row]
    } else if(ref_uptake_age$DOSE[i_row] == 'X'){
      AG_uptake_matrix_X[ref_uptake_age$DATE[i_row] == date_opt,ref_uptake_age$age_factor[i_row]] <- ref_uptake_age$COUNT[i_row]
    }                 
  }
  colnames(AG_uptake_matrix_A) <- paste0('epistat_uptake_A_',gsub('-','_',opt_agegroup_epistat))
  colnames(AG_uptake_matrix_B) <- paste0('epistat_uptake_B_',gsub('-','_',opt_agegroup_epistat))
  colnames(AG_uptake_matrix_C) <- paste0('epistat_uptake_C_',gsub('-','_',opt_agegroup_epistat))
  colnames(AG_uptake_matrix_E) <- paste0('epistat_uptake_E_',gsub('-','_',opt_agegroup_epistat))
  colnames(AG_uptake_matrix_F) <- paste0('epistat_uptake_F_',gsub('-','_',opt_agegroup_epistat))
  colnames(AG_uptake_matrix_X) <- paste0('epistat_uptake_X_',gsub('-','_',opt_agegroup_epistat))
  epistat_uptake_matrix <- data.frame(date = as.character(date_opt),
                                      AG_uptake_matrix_A,
                                      AG_uptake_matrix_B,
                                      AG_uptake_matrix_C,
                                      AG_uptake_matrix_E,
                                      AG_uptake_matrix_F,
                                      AG_uptake_matrix_X)
 
  # merge uptake for the total population- and by age group
  ref_uptake_dose <- merge(ref_uptake_dose,
                           ref_uptake_aggr,
                           all.x = T,
                           by='date')
  ref_uptake_dose <- merge(ref_uptake_dose,
                           ref_uptake_aggr_AZ,
                           all.x = T,
                           by='date')
  ref_uptake_dose <- merge(ref_uptake_dose,
                           ref_uptake_aggr_mRNA,
                           all.x = T,
                           by='date')
  ref_uptake_dose <- merge(ref_uptake_dose,
                           ref_uptake_aggr_JnJ,
                           all.x = T,
                           by='date')
 # ref_uptake_dose <- merge(ref_uptake_dose,
                         #  ref_uptake_aggr_BIV,
                         #  all.x = T,
                         #  by='date')
  ref_uptake_dose <- merge(ref_uptake_dose,
                           epistat_uptake_matrix,
                           all.x = T,
                           by='date')
  
  # check
  names(ref_uptake_dose)
  head(ref_uptake_dose)
  tail(ref_uptake_dose)
  ref_uptake_dose[is.na(ref_uptake_dose)] <- 0
  
  # aggregate all doses up to January 1st, 2021
  flag_prev <- ref_uptake_dose$date < as.Date('2021-01-01')
  flag_jan1 <- ref_uptake_dose$date == as.Date('2021-01-01')
  ref_uptake_dose[flag_jan1 ,-1] <- colSums(ref_uptake_dose[flag_prev| flag_jan1 ,-1])
  ref_uptake_dose[flag_prev,-1]  <- 0

  # add Rt on reported cases and hospital admissions
  ref_case_data_age$Rt_cases <- get_Rt(ref_case_data_age$cases,as.Date(ref_case_data_age$date))$Rt
  hosp_adm_data$Rt_hosp      <- get_Rt(hosp_adm_data$hospital_admissions,as.Date(hosp_adm_data$date))$Rt
  
  #tail(hosp_adm_data)
  # plot(as.Date(ref_case_data_age$date),
  #      ref_case_data_age$Rt_cases,type='l')
  # lines(as.Date(hosp_adm_data$date),
  #      hosp_adm_data$Rt_hosp,col=2)
  
  # merge hospital and case data
  ref_data <- merge(hosp_adm_data,ref_case_data_age,all=TRUE)
  ref_data <- merge(ref_data,ref_test_data,all = TRUE)
  ref_data <- merge(ref_data,ref_mort_data,all = TRUE)
  ref_data <- merge(ref_data,ref_uptake_dose,all = TRUE)
  names(ref_data)
  
  # final adjustments
  ref_data$date <- as.Date(ref_data$date)
  
  # set vaccine uptake before January 1st, 2021, to 0
  ref_data[,grepl('vaccine',names(ref_data))] <- replace(ref_data[,grepl('vaccine',names(ref_data))],
                                                         is.na(ref_data[,grepl('vaccine',names(ref_data))]),
                                                         0) 
  # add age-specific "full protection" (~ BC)
  names_vaccine_B <- names(ref_data)[grepl('vaccine_B_',names(ref_data))]
  names_vaccine_BC <- gsub('_B_','_BC_',names_vaccine_B)
  ref_data[,names_vaccine_BC] <- ref_data[,names_vaccine_B] + ref_data[,grepl('vaccine_JnJ_A_',names(ref_data))]
  
  # plot(cumsum(rowSums(ref_uptake_dose[,names_vaccine_BC])))
  # lines(cumsum(rowSums(ref_uptake_dose[,names_vaccine_B])))
  # dim(ref_uptake_dose)
  # names(ref_uptake_dose)
  #browser()
  # save as csv
  write.table(ref_data,file=ref_data_file_name,sep=',',row.names=F)
  
  # return data
  return(ref_data)
}

select_regional_data <- function(ref_data,sel_region){
  
  if(is.na(sel_region) | sel_region == 'belgium'){
    return(ref_data)
  }
  
  if(!any(names(ref_data) == 'REGION')){
    print("NO REGION SPECIFIC REFERENCE DATA")
    return(ref_data)
  }
  
  # else
  return(ref_data[!is.na(ref_data$REGION) & grepl(substr(sel_region,1,5),tolower(ref_data$REGION)),])
}

convert_uptake_age_groups <- function(ref_uptake_data_all,sel_region,sel_brand = NA){
  
  # make copy (to edit)
  ref_uptake_data_adj <- ref_uptake_data_all
  
  # get age groups based on full data set
  ref_uptake_data_adj$AGEGROUP <- gsub('\\+','-99',ref_uptake_data_adj$AGEGROUP)
  opt_agegroup_epistat         <- sort(unique(ref_uptake_data_adj$AGEGROUP))
  ref_uptake_data_adj$age_factor <- as.numeric(factor(ref_uptake_data_adj$AGEGROUP,labels = opt_agegroup_epistat))
  
  # include check on age groups, since they change frequently!
  agegroup_epistat_previous <- c("00-04","05-11","12-15","16-17","18-24","25-34",
                                 "35-44","45-54","55-64","65-74","75-84","85-99")
  if(any(opt_agegroup_epistat != agegroup_epistat_previous)){
    warning("!! EPISTAT VACCINE UPTAKE AGE GROUPS HAVE CHANGED !!")
  }
  
  # option to select a BRAND
  if(!is.na(sel_brand)){
    ref_uptake_data_adj <- ref_uptake_data_adj[ref_uptake_data_adj$BRAND == sel_brand,]
  }
  
  # reformat age groups
  unique(ref_uptake_data_adj$AGEGROUP)
  flag_age_epistat_18 <- ref_uptake_data_adj$AGEGROUP == "18-24"
  
  # get uptake 18-24y
  tmp_uptake <- ref_uptake_data_adj[flag_age_epistat_18,]
  # create list with 2/7 to add to "16-17"
  #tmp_uptake$COUNT <- round(tmp_uptake$COUNT * 2/7)
  # Update June 2022 Sciensano +1 year -> add 3/7 to 16-17
  tmp_uptake$COUNT <- round(tmp_uptake$COUNT * 3/7)
  tmp_uptake$AGEGROUP <- "16-17"
  tmp_uptake$age_factor <- 4
  
  # adjust uptake 18-24 accordingly
  ref_uptake_data_adj$COUNT[flag_age_epistat_18] <- ref_uptake_data_adj$COUNT[flag_age_epistat_18] - tmp_uptake$COUNT
  ref_uptake_data_adj <- rbind(ref_uptake_data_adj,
                               tmp_uptake)
  
  # reshift uptake 05-11y
  flag_age_epistat_11 <- ref_uptake_data_adj$AGEGROUP == "05-11"
  if(any(flag_age_epistat_11)){
    tmp_uptake <- ref_uptake_data_adj[flag_age_epistat_11,]
    # create list with 2/7 to add to "12-15"
   # tmp_uptake$COUNT <- round(tmp_uptake$COUNT * 2/7)
    # Update June 2022 Sciensano +1 year -> add 1/7 to 12-15
    tmp_uptake$COUNT <- round(tmp_uptake$COUNT * 1/7)
    tmp_uptake$AGEGROUP   <- "12-15"
    tmp_uptake$age_factor <- 3
    
    # adjust uptake 05-11 accordingly
    ref_uptake_data_adj$COUNT[flag_age_epistat_11] <- ref_uptake_data_adj$COUNT[flag_age_epistat_11] - tmp_uptake$COUNT
    ref_uptake_data_adj <- rbind(ref_uptake_data_adj,
                                 tmp_uptake)    
  }

  # include vaccine uptake by age group ---
  ref_uptake_age          <- aggregate(COUNT ~ + DOSE + AGEGROUP + age_factor + DATE, data=ref_uptake_data_adj,sum)
  ref_uptake_age$DATE     <- as.Date(ref_uptake_age$DATE)
  ref_uptake_age          <- ref_uptake_age[order(ref_uptake_age$DATE),]
  
  # load pop data, by 5y age groups, rename population column, add age group indices
  # be_pop_5y        <- read.csv('data/pop_be2020_statbel_5y.csv')
  be_pop_5y         <- get_regional_pop(sel_region,bool_10y = FALSE)
  be_pop_5y$count   <- be_pop_5y$Population.on.January.1st.2020
  be_pop_5y$age_10y  <- rep(1:10,each=2)
  be_pop_5y$age_uptake <-c(1:5,
                           rep(6:11,each=2),
                           rep(12,3))
  be_pop_5y$age_epistat <- opt_agegroup_epistat[be_pop_5y$age_uptake ]
  
  # Update June 2022 Sciensano +1 year -> add adjust coefficients #new +1 update in august 23
  be_pop_5y$age_uptake_adjustcoef <-c(rep(1,5),
                                      rep(c(7,3),6),
                           c(7,5,3))
  # calculate 5-y age-group weights per uptake age group
  # Update June 2022 Sciensano +1 year -> use adjust coefficients #new update +1 in august 23
  #be_pop_5y$age_10y == i_age
  for(i_age in unique(be_pop_5y$age_uptake)){
    flag_age <- be_pop_5y$age_uptake == i_age
    be_pop_5y$weight_age[flag_age] <- be_pop_5y$count[flag_age] * be_pop_5y$age_uptake_adjustcoef[flag_age] / sum(be_pop_5y$count[flag_age] * be_pop_5y$age_uptake_adjustcoef[flag_age])
    be_pop_5y$weight_age_original[flag_age] <- be_pop_5y$count[flag_age] / sum(be_pop_5y$count[flag_age])
  }


  # loop over original uptake data, and assign weighted first/second doses to [date,age] combinations
  date_opt            <- unique(ref_uptake_age$DATE)
  ref_uptake_matrix_A <- matrix(0,nrow=length(date_opt),ncol=20)
  ref_uptake_matrix_B <- ref_uptake_matrix_A
  ref_uptake_matrix_E <- ref_uptake_matrix_A
  ref_uptake_matrix_F <- ref_uptake_matrix_A
  ref_uptake_matrix_X <- ref_uptake_matrix_A
  i_row <- 1
  
  for(i_row in 1:nrow(ref_uptake_age)){
    uptake_age_group <- ref_uptake_age$COUNT[i_row] * be_pop_5y$weight_age * (be_pop_5y$age_uptake == ref_uptake_age$age_factor[i_row])
    if(ref_uptake_age$DOSE[i_row] == 'B'){
      ref_uptake_matrix_B[ref_uptake_age$DATE[i_row] == date_opt,] <- ref_uptake_matrix_B[ref_uptake_age$DATE[i_row] == date_opt,] + uptake_age_group
    } else if(ref_uptake_age$DOSE[i_row] == 'E'){
      ref_uptake_matrix_E[ref_uptake_age$DATE[i_row] == date_opt,] <- ref_uptake_matrix_E[ref_uptake_age$DATE[i_row] == date_opt,] + uptake_age_group
    } else if(ref_uptake_age$DOSE[i_row] %in% c('E2','E3')){
      ref_uptake_matrix_F[ref_uptake_age$DATE[i_row] == date_opt,] <- ref_uptake_matrix_F[ref_uptake_age$DATE[i_row] == date_opt,] + uptake_age_group
    } else if(ref_uptake_age$DOSE[i_row] == 'X'){
      ref_uptake_matrix_X[ref_uptake_age$DATE[i_row] == date_opt,] <- ref_uptake_matrix_X[ref_uptake_age$DATE[i_row] == date_opt,] + uptake_age_group
    }  else{
      ref_uptake_matrix_A[ref_uptake_age$DATE[i_row] == date_opt,] <- ref_uptake_matrix_A[ref_uptake_age$DATE[i_row] == date_opt,] + uptake_age_group
    }                 
  }
  
  # aggregate weighted uptake into 10y age classes
  ref_uptake_aggr <- data.frame(date=date_opt)
  age_tags        <- get_age_groups()
  
  i_age <- 1
  for(i_age in 1:length(age_tags)){
    ref_uptake_aggr[paste0('vaccine_A_ages',age_tags[i_age])] <- rowSums(ref_uptake_matrix_A[,be_pop_5y$age_10y == i_age])
    ref_uptake_aggr[paste0('vaccine_B_ages',age_tags[i_age])] <- rowSums(ref_uptake_matrix_B[,be_pop_5y$age_10y == i_age])
    ref_uptake_aggr[paste0('vaccine_E_ages',age_tags[i_age])] <- rowSums(ref_uptake_matrix_E[,be_pop_5y$age_10y == i_age])
    ref_uptake_aggr[paste0('vaccine_F_ages',age_tags[i_age])] <- rowSums(ref_uptake_matrix_F[,be_pop_5y$age_10y == i_age])
    ref_uptake_aggr[paste0('vaccine_X_ages',age_tags[i_age])] <- rowSums(ref_uptake_matrix_X[,be_pop_5y$age_10y == i_age])
    }
  ref_uptake_aggr
  ref_uptake_aggr <- ref_uptake_aggr[,order(names(ref_uptake_aggr))]
  
  # fix dates as strings
  ref_uptake_aggr$date <- as.character(ref_uptake_aggr$date)
  
  return(ref_uptake_aggr)
}


# download file from URL
# if connection is not possible: return NA
download_ref_file <- function(cases_ref_url,data_dir){
  if(!dir.exists(data_dir)){dir.create(data_dir,recursive = T)}
  case_ref_file   <- file.path(data_dir,basename(cases_ref_url))
  exit_status     <- tryCatch(download.file(cases_ref_url,case_ref_file,quiet=TRUE),
                              error = function(e){return(-1)})
  return(ifelse(exit_status==0,case_ref_file,NA))
}

get_be_pop_matrix <- function(n_row){
  
  #pop_be2020 <- pop_age(wpp_age("Belgium",2020),seq(0,90,10))$population
  pop_be2020_all <- read.table('data/pop_be2020_statbel.csv',sep=',',header=T)
  pop_be2020 <- pop_be2020_all[,3]
  
  return(matrix(rep(pop_be2020, each = n_row),
                nrow = n_row,
                ncol = length(pop_be2020)))
  
}


#region <- 'flanders'
get_regional_pop <- function(region = 'belgium',bool_10y = TRUE){
  
  pop_be2020_all <- read.table('data/pop_be2020_statbel_download.csv',sep=',',header=T)
  unique(pop_be2020_all$Region)
  if(is.na(region) | region == 'belgium'){
    pop_be2020_all <- pop_be2020_all[nchar(pop_be2020_all$Region)==0,]
  } else{
    pop_be2020_all <- pop_be2020_all[grepl(substr(region,1,5),tolower(pop_be2020_all$Region)),]
  }
  
  if(length(unique(pop_be2020_all$Region))>1){
    warning('provided "region" is not vallid, use belgium, flanders, walloon or brussels')
    return(NA)
  }
  
  # aggregate 90-99 and +100
  pop_be2020_all$Population.on.January.1st.2020[nrow(pop_be2020_all)-1] <- sum(pop_be2020_all$Population.on.January.1st.2020[nrow(pop_be2020_all)-(1:0)])
  pop_be2020_all <- pop_be2020_all[-nrow(pop_be2020_all),]
  
  # optional: aggregate 5y age groups into 10y age groups
  if(bool_10y){
    pop2020 <- colSums(matrix(c(pop_be2020_all$Population.on.January.1st.2020),ncol=10,nrow=2,byrow=F))
    return(pop2020)
  } else{
    return(pop_be2020_all)
  }

}


#OLD FUNCTION, BASED ON MERGED EXCEL FILE FROM FOD
get_latest_incidence_data_xlsx <- function(enable_data_reuse = TRUE, sel_region=NA){

  be_ref_data <- get_observed_incidence_data(enable_data_reuse=enable_data_reuse,
                                             sel_region=sel_region)

  # set gdrive on LW's computer
  gdrive_data_dir <- '~/Google Drive/nCov2019 - research/Shared COVID-19 drive with third parties (confidential information)/Data Sciensano/Hosp_surge'
  if(!dir.exists(gdrive_data_dir)){
    #return(be_ref_data)
    gdrive_data_dir <- './data/hosp_surge'
  }

  # get file names
  gdrive_data_files <- dir(gdrive_data_dir,pattern = 'Hosp_surge_overview_')
  gdrive_data_files

  # remove uncommon file names
  gdrive_data_files <- gdrive_data_files[nchar(gdrive_data_files)==38]

  # get date, and select more recent file
  file_tag_hour  <- substr(gdrive_data_files,30,31)
  file_tag_stamp <- as.Date(substr(gdrive_data_files,21,28),'%d%m%Y') + as.numeric(substr(gdrive_data_files,30,31)) / 24
  latest_data_file <- gdrive_data_files[order(file_tag_stamp,decreasing = T)[1]]

  # get file name for latest reference data
  ref_data_file_name <- file.path('data/downloads',paste0('covid19_reference_data_',gsub('-','',Sys.Date()),'.csv'))

  # add hour stamp
  if(Sys.Date() == as.Date(substr(latest_data_file,21,28),'%d%m%Y')){
    ref_data_file_name <- gsub('.csv',paste0('_',substr(latest_data_file,30,33),'.csv'),ref_data_file_name)
  } else{
    ref_data_file_name <- gsub('.csv','_0000.csv',ref_data_file_name)
  }

  # specify region?
  if(!is.na(sel_region)){
    ref_data_file_name <- gsub('.csv',paste0('_',sel_region,'.csv'),ref_data_file_name)
  }

  if(enable_data_reuse && file.exists(ref_data_file_name)){
    ref_data      <- read.table(ref_data_file_name,sep=',',header=T)
    ref_data$date <- as.Date(ref_data$date)
    return(ref_data)
  }

  # BE data: read file and check date
  db_hosp_surge <- read.xlsx(file.path(gdrive_data_dir,latest_data_file),
                             sheet = "Summary",
                             startRow = 3,
                             detectDates = T)



  # override dataset if regional data is required
  if(!is.na(sel_region) && sel_region != 'belgium'){

    db_hosp_surge_hosp <- read.xlsx(file.path(gdrive_data_dir,latest_data_file),
                               sheet = "Hospital data",
                               startRow = 1,
                               detectDates = T)
    dim(db_hosp_surge_hosp)
    names(db_hosp_surge_hosp)
    table(db_hosp_surge_hosp$Region)
    db_hosp_surge_hosp$region_en <- gsub('Brussel','brussels',db_hosp_surge_hosp$Region)
    db_hosp_surge_hosp$region_en <- gsub('Vlaanderen','flanders',db_hosp_surge_hosp$region_en)
    db_hosp_surge_hosp$region_en <- gsub('Wallonië','wallonia',db_hosp_surge_hosp$region_en)
    db_hosp_surge_hosp <- db_hosp_surge_hosp[!is.na(db_hosp_surge_hosp$region_en),]
    db_hosp_surge_hosp <- db_hosp_surge_hosp[db_hosp_surge_hosp$region_en == sel_region,]
    dim(db_hosp_surge_hosp)

    db_hosp_surge <- aggregate(. ~ Date + region_en,
                               data = db_hosp_surge_hosp[,c('Date','region_en',
                                                            'NewPatientsNotReferredHospital',
                                                            'NewPatientsOtherPathology',
                                                            'Confirmed.patients.in.hospital',
                                                            'Confirmed.patients.in.ICU',
                                                            'ICUIn',
                                                            'Mortality')],
                               FUN = sum, na.rm=T,
                               na.action = na.pass)
    dim(db_hosp_surge)
    range(db_hosp_surge$Date)
    names(db_hosp_surge)[1] <- 'Row.Labels'
    names(db_hosp_surge)[-1] <- paste0('Sum.of.',names(db_hosp_surge)[-1])
  }

  # if new hospital data available, add
  db_hosp_surge$date                 <- as.Date(db_hosp_surge$Row.Labels)
  sel_hosp_surge                     <- !db_hosp_surge$date %in% be_ref_data$date
  if(any(sel_hosp_surge)){
    ref_data_extra                     <- be_ref_data[1:sum(sel_hosp_surge),]
    ref_data_extra[]                   <- NA
    ref_data_extra$date                <- db_hosp_surge[sel_hosp_surge,'Row.Labels']
    ref_data_extra$hospital_admissions <- db_hosp_surge[sel_hosp_surge,'Sum.of.NewPatientsNotReferredHospital']
    ref_data_extra$hospital_load       <- db_hosp_surge[sel_hosp_surge,'Sum.of.Confirmed.patients.in.hospital']
    ref_data_extra$icu_load            <- db_hosp_surge[sel_hosp_surge,'Sum.of.Confirmed.patients.in.ICU']

    be_ref_data                        <- rbind(be_ref_data,ref_data_extra)
    be_ref_data$date                   <- as.Date(be_ref_data$date)
  }

  # add ICU admissions and hospital mortality
  be_ref_data <- merge(be_ref_data,
                       data.frame(date=as.Date(db_hosp_surge$Row.Labels),
                                  icu_admissions=db_hosp_surge$Sum.of.ICUIn,
                                  hospital_mortality = db_hosp_surge$Sum.of.Mortality,
                                  hospital_admissions_other = db_hosp_surge$Sum.of.NewPatientsOtherPathology
                                  ),
                       by='date',
                       all.x = T)

  # add hospital exit
  # be_ref_data$hospital_exit <- -(c(diff(be_ref_data$hospital_load),0)
  #                                    - be_ref_data$hospital_admissions)
  be_ref_data$hospital_exit <-  be_ref_data$hospital_discharges +
                                be_ref_data$hospital_mortality

  # add hospital recovery
  be_ref_data$hospital_recovery <- -(c(0,diff(be_ref_data$hospital_load))
                                     + be_ref_data$hospital_mortality
                                     - be_ref_data$hospital_admissions)

  # add 7day mean hospital recovery
  be_ref_data$hospital_recovery_7dmean <- round(rollmean(be_ref_data$hospital_recovery,k=7,align='center',fill=NA))
  be_ref_data$hospital_recovery_7dmean[be_ref_data$hospital_recovery_7dmean<0] <- 0

  # print CLI update
  print(be_ref_data[nrow(be_ref_data),1:5])

  # save updated file
  write.table(be_ref_data,file=ref_data_file_name,sep=',',row.names = F)

  # return updated dataset
  return(be_ref_data)
}

get_latest_incidence_data <- function(enable_data_reuse = TRUE, sel_region=NA){
  
  be_ref_data <- get_observed_incidence_data(enable_data_reuse=enable_data_reuse,
                                             sel_region=sel_region)
  
  # set gdrive on LW's computer
  gdrive_data_dir <- '~/Google Drive/nCov2019 - research/Shared COVID-19 drive with third parties (confidential information)/Data Sciensano/Last Update'
  if(!dir.exists(gdrive_data_dir)){
    #return(be_ref_data)
    gdrive_data_dir <- './data/hosp_surge'
  }
  
  # get file name for latest reference data
  ref_data_file_name <- file.path('data/downloads',paste0('covid19_reference_data_',gsub('-','',Sys.Date()),'.csv'))
  
  # retrieve available hosp surge data
  latest_data_file <- dir(gdrive_data_dir,pattern = 'hosp_surge_perhospital',full.names = T)

  # select most recent one
  latest_data_file <- latest_data_file[order(file.info(latest_data_file)$mtime,decreasing = T) == 1]
  
  # add hour stamp to ref file
  if(Sys.Date() == as.Date(file.info(latest_data_file)$mtime)){
    ref_data_file_name <- gsub('.csv',paste0('_',format(file.info(latest_data_file)$mtime,'%H%M'),'.csv'),ref_data_file_name)
  } else{
    ref_data_file_name <- gsub('.csv','_0000.csv',ref_data_file_name)
  }

  # specify region
  if(is.na(sel_region)){
    sel_region <- get_region(1)
  }
  ref_data_file_name <- gsub('.csv',paste0('_',sel_region,'.csv'),ref_data_file_name)
  
  if(enable_data_reuse && file.exists(ref_data_file_name)){
    ref_data      <- read.table(ref_data_file_name,sep=',',header=T)
    ref_data$date <- as.Date(ref_data$date)
    return(ref_data)
  }
  
  # read hosp surve file and set default region
  db_hosp_surge_all           <- read.csv(latest_data_file)

  if('Confirmed.patients.other.pathology' %in% names(db_hosp_surge_all)){
    db_hosp_surge_all$NewPatientsOtherPathology <- db_hosp_surge_all$Confirmed.patients.other.pathology
    db_hosp_surge_all$ICUIn <- db_hosp_surge_all$ICU.in
  }
  
  # override dataset if regional data is required  
  if(sel_region != 'belgium'){
    
    dim(db_hosp_surge_all)
    names(db_hosp_surge_all)
    table(db_hosp_surge_all$Region)
    db_hosp_surge_all$region_en <- gsub('Brussel','brussels',db_hosp_surge_all$Region)
    db_hosp_surge_all$region_en <- gsub('Vlaanderen','flanders',db_hosp_surge_all$region_en)
    db_hosp_surge_all$region_en <- gsub('Walloni\xeb','wallonia',db_hosp_surge_all$region_en)
    
    db_hosp_surge_all <- db_hosp_surge_all[!is.na(db_hosp_surge_all$region_en),]
    db_hosp_surge_all <- db_hosp_surge_all[db_hosp_surge_all$region_en == sel_region,]
    dim(db_hosp_surge_all)
  } 

    db_hosp_surge <- aggregate(. ~ Date, 
                               data = db_hosp_surge_all[,c('Date',
                                                            'NewPatientsNotReferredHospital',
                                                            'NewPatientsOtherPathology',
                                                            'Confirmed.patients.in.hospital',
                                                            'Confirmed.patients.in.ICU',
                                                            'ICUIn',
                                                            'Mortality')],
                               FUN = sum, na.rm=T,
                               na.action = na.pass)
  
    dim(db_hosp_surge)
    range(db_hosp_surge$Date)
    names(db_hosp_surge)[1] <- 'Row.Labels'
    names(db_hosp_surge)[-1] <- paste0('Sum.of.',names(db_hosp_surge)[-1])
 
  
  # if new hospital data available, add
  db_hosp_surge$date                 <- as.Date(db_hosp_surge$Row.Labels)
  sel_hosp_surge                     <- !db_hosp_surge$date %in% be_ref_data$date
  if(any(sel_hosp_surge)){
    ref_data_extra                     <- be_ref_data[1:sum(sel_hosp_surge),]
    ref_data_extra[]                   <- NA
    ref_data_extra$date                <- db_hosp_surge[sel_hosp_surge,'Row.Labels']
    ref_data_extra$hospital_admissions <- db_hosp_surge[sel_hosp_surge,'Sum.of.NewPatientsNotReferredHospital']
    ref_data_extra$hospital_load       <- db_hosp_surge[sel_hosp_surge,'Sum.of.Confirmed.patients.in.hospital']
    ref_data_extra$icu_load            <- db_hosp_surge[sel_hosp_surge,'Sum.of.Confirmed.patients.in.ICU']
    ref_data_extra$region              <- sel_region
    
    be_ref_data                        <- rbind(be_ref_data,ref_data_extra)
    be_ref_data$date                   <- as.Date(be_ref_data$date)
  }
  
  # add ICU admissions and hospital mortality
  be_ref_data <- merge(data.frame(date=as.Date(db_hosp_surge$Row.Labels),
                                  # hospital_admissions = db_hosp_surge$Sum.of.NewPatientsNotReferredHospital,
                                  # hospital_load = db_hosp_surge$Sum.of.Confirmed.patients.in.hospital,
                                  # icu_load = db_hosp_surge$Sum.of.Confirmed.patients.in.ICU,
                                  icu_admissions=db_hosp_surge$Sum.of.ICUIn,
                                  # region = db_hosp_surge$Sum.of.region_en,
                                  hospital_mortality = db_hosp_surge$Sum.of.Mortality,
                                  hospital_admissions_other = db_hosp_surge$Sum.of.NewPatientsOtherPathology
                       ),
                       be_ref_data,
                       by='date',
                       all = T)
  
  # add hospital exit
  be_ref_data$hospital_exit <-  be_ref_data$hospital_discharges +
                                  be_ref_data$hospital_mortality
  
  # add hospital recovery
  be_ref_data$hospital_recovery <- -(c(0,diff(be_ref_data$hospital_load)) 
                                     + be_ref_data$hospital_mortality
                                     - be_ref_data$hospital_admissions)
  
  # add 7day mean hospital recovery
  be_ref_data$hospital_recovery_7dmean <- round(rollmean(be_ref_data$hospital_recovery,k=7,align='center',fill=NA))
  be_ref_data$hospital_recovery_7dmean[be_ref_data$hospital_recovery_7dmean<0] <- 0
  
  # print CLI update
  print(be_ref_data[nrow(be_ref_data),c('date',
                                        'hospital_admissions',
                                        'hospital_admissions_other',
                                        'hospital_load',
                                        'icu_load',
                                        'region')])
  
  # add hospital admission age distribution
  age_dist_hosp_mat <- get_hospital_age_distr(bool_use_default = FALSE , num_days = nrow(be_ref_data))
  age_dist_hosp <- data.frame(date = be_ref_data$date[1:ncol(age_dist_hosp_mat)],
                              t(age_dist_hosp_mat))
  names(age_dist_hosp)[-1] <- paste0('prop_hospital_admission_age',1:(ncol(age_dist_hosp)-1))
  head(age_dist_hosp)
  
  be_ref_data <- merge(be_ref_data,
                       age_dist_hosp,
                       by='date',
                       all = T)
  
  # save updated file
  write.table(be_ref_data,file=ref_data_file_name,sep=',',row.names = F)
  
  # return updated dataset
  return(be_ref_data)
}

get_hospital_age_distr <- function(bool_use_default = FALSE, num_days = NA){
  
  # set gdrive on LW's computer
  gdrive_data_dir <- '~/Google Drive/nCov2019 - research/Shared COVID-19 drive with third parties (confidential information)/Data Sciensano/Last Update'
  if(!dir.exists(gdrive_data_dir)){
    #return(be_ref_data)
    gdrive_data_dir <- './data/hosp_surge'
  }
  
  # set file name and check presence
  filename_clinic <- dir(gdrive_data_dir,pattern='COVID19BE_CLINIC.csv',full.names = T)
  if(bool_use_default || length(filename_clinic)==0){
    #write.csv(get_hospital_age_distr(),file.path('./data/hosp_surge','COVID19BE_CLINIC_AGGR.csv'),row.names=F)
    return(read.csv(file.path('./data/hosp_surge','COVID19BE_CLINIC_AGGR.csv')))
  }
  
  print("Use COVID19BE_CLINIC file")
  
  # load data
  db_clinic <- read.csv(dir(gdrive_data_dir,pattern='COVID19BE_CLINIC.csv',full.names = T))
  dim(db_clinic)
  
  # # create and store placeholder for GITHUB repository
  # db_dummy <- db_clinic[as.Date(db_clinic$dt_admission) < as.Date(2020-05-01),]
  # db_dummy <- db_dummy[seq(1,nrow(db_dummy),15),]
  # write.table(db_dummy,file.path('./data/hosp_surge/COVID19BE_CLINIC_placeholder.csv'),sep=',',row.names=F)
  
  # specify key indicators
  db_clinic$AGE_GROUP     <- cut(db_clinic$age,breaks = c(seq(0,90,10),Inf),right = F)
  db_clinic$DAY           <- as.Date(db_clinic$dt_admission)
  db_clinic$N_ADMISSIONS  <- 1
  
  # aggregate
  be_clinic     <- aggregate(N_ADMISSIONS ~ DAY, data= db_clinic,sum)
  be_clinic_age <- aggregate(N_ADMISSIONS ~ DAY + AGE_GROUP, data= db_clinic,sum,drop=F)
  names(be_clinic)[2] <- "N_ADMISSIONS_TOTAL"
  be_clinic_age$N_ADMISSIONS[is.na(be_clinic_age$N_ADMISSIONS)] <- 0
  
  # merge
  be_clinic_age <- merge(be_clinic_age,
                         be_clinic)
  
  # select data after SCM start and discart last 7 days
  be_clinic_age <- be_clinic_age[be_clinic_age$DAY >= get_scm_start_date(),]
  #be_clinic_age <- be_clinic_age[be_clinic_age$DAY <= max(be_clinic_age$DAY-7),]
  range(be_clinic_age$DAY)

  # set matrix format
  be_clinic_age        <- be_clinic_age[order(be_clinic_age$DAY),]
  be_clinic_age_days   <- unique(be_clinic_age$DAY)
  be_clinic_age_groups <- sort(unique(be_clinic_age$AGE_GROUP))
  
  tail(be_clinic_age_days)
  age_dist_hosp_mat_day <- matrix(NA,nrow=10,ncol=length(be_clinic_age_days))
  for(i_age in 1:length((be_clinic_age_groups))){
    age_dist_hosp_mat_day[i_age,] <- be_clinic_age$N_ADMISSIONS[be_clinic_age$AGE_GROUP == be_clinic_age_groups[i_age]]
  }
  
  # use 7-day averages
  age_dist_hosp_mat_day <- t(apply(age_dist_hosp_mat_day,1,rollmean,k=7,align='center',fill=NA))
  
  # remove last 2 weeks
  dim(age_dist_hosp_mat_day)
  age_dist_hosp_mat_day    <- age_dist_hosp_mat_day[,1:(ncol(age_dist_hosp_mat_day)-14)]
  
  # get proportion
  age_dist_hosp_mat_day  <- t(t(age_dist_hosp_mat_day) / rowSums(t(age_dist_hosp_mat_day)))
  
  # add 7-day average to the end
  age_dist_hosp_mat_day <- cbind(age_dist_hosp_mat_day,
                                 rowMeans(age_dist_hosp_mat_day[,ncol(age_dist_hosp_mat_day)-(0:6)]))
  
  # adjust range
  dim(age_dist_hosp_mat_day)
  length(obs_total_hosp_ext)
  sel_col <- 1:ifelse(is.na(num_days),length(obs_total_hosp_ext),num_days)
  sel_col[sel_col>ncol(age_dist_hosp_mat_day)] <- ncol(age_dist_hosp_mat_day)
  age_dist_hosp_mat_day <- age_dist_hosp_mat_day[,sel_col]
  
  # return
  #write.csv(age_dist_hosp_mat_day,file.path('./data/hosp_surge','COVID19BE_CLINIC_AGGR.csv'),row.names=F)
  return(age_dist_hosp_mat_day)
}

# get_regional_hospital_age_distr <- function(sel_region = NA, date_tag = NA){
# 
#   ## NEW DATA SOURCE
#   data_age <- readRDS('data/tests_cases_mortality_hospi_vacc_per_age_and_province/tests_cases_morta_hospi_latest.Rds')
# 
#   if(!is.na(date_tag)){
#     data_age <- readRDS(paste0('data/tests_cases_mortality_hospi_vacc_per_age_and_province/tests_cases_morta_hospi_until_',date_tag,'.Rds'))
#   }
#   
#   # select region?
#   if(!(is.na(sel_region) || sel_region == 'belgium')  ){
#     
#     table(data_age$PROVINCE)
#     sel_flanders <-   data_age$PROVINCE %in% c('Antwerpen','Limburg','Oost-Vlaanderen','Vlaams-Brabant','West-Vlaanderen')
#     sel_wallonia <-   data_age$PROVINCE %in% c('Brabant Wallon','Hainaut','Liège','Luxembourg','Namur')
#     sel_brussels <-   data_age$PROVINCE %in% c('Bruxelles')
#     
#     if(sel_region == 'flanders'){
#       data_age <- data_age[sel_flanders,]
#     } else if(sel_region == 'wallonia'){
#       data_age <- data_age[sel_wallonia,]
#     } else if(sel_region == 'brussels'){
#       data_age <- data_age[sel_brussels,]
#     } else{
#       stop("REGION UNKNOWN")
#     }
#     
#   }
#   
#   be_hosp     <- aggregate(N_ADMISSIONS ~ DAY, data= data_age,sum)
#   be_hosp_age <- aggregate(N_ADMISSIONS ~ DAY + AGE_GROUP, data= data_age,sum)
#   
#   names(be_hosp)[2] <- "N_ADMISSIONS_TOTAL"
#   
#   be_hosp_age <- merge(be_hosp_age,
#                        be_hosp)
#   
#   dim(be_hosp_age)
#   range(be_hosp_age$DAY)
#   diff(range(be_hosp_age$DAY))
#   tail(be_hosp_age)
#   
#   
#   # set matrix format
#   be_hosp_age        <- be_hosp_age[order(be_hosp_age$DAY),]
#   be_hosp_age_days   <- unique(be_hosp_age$DAY)
#   be_hosp_age_groups <- sort(unique(be_hosp_age$AGE_GROUP))
#   
#   tail(be_hosp_age_days)
#   age_dist_hosp_mat_day <- matrix(NA,nrow=10,ncol=length(be_hosp_age_days))
#   for(i_age in 1:length((be_hosp_age_groups))){
#     age_dist_hosp_mat_day[i_age,] <- be_hosp_age$N_ADMISSIONS[be_hosp_age$AGE_GROUP == be_hosp_age_groups[i_age]]
#   }
#   
#   # use 7-day averages
#   age_dist_hosp_mat_day <- t(apply(age_dist_hosp_mat_day,1,rollmean,k=7,align='center',fill=NA))
#   
#   # remove last 2 weeks
#   dim(age_dist_hosp_mat_day)
#   age_dist_hosp_mat_day    <- age_dist_hosp_mat_day[,1:(ncol(age_dist_hosp_mat_day)-14)]
#   
#   # get proportion
#   age_dist_hosp_mat_day  <- t(t(age_dist_hosp_mat_day) / rowSums(t(age_dist_hosp_mat_day)))
#   
#   # adjust range
#   dim(age_dist_hosp_mat_day)
#   length(obs_total_hosp_ext)
#   sel_col <- 1:length(obs_total_hosp_ext)
#   sel_col[sel_col>ncol(age_dist_hosp_mat_day)] <- ncol(age_dist_hosp_mat_day)
#   age_dist_hosp_mat_day <- age_dist_hosp_mat_day[,sel_col]
#   
#   # return
#   return(age_dist_hosp_mat_day)
# }

# Function to efficiently generate summary tables
#data_transm <- ref_case_data_table; colname_date="DATE";colname_value='AGEGROUP';prefix='cases'
get_age_table <- function(data_transm){
  
  # overal summary
  summary_table_general         <- aggregate(formula('CASES ~ DATE'), data=data_transm, sum)
  names(summary_table_general)  <- c('date','cases')
  
  # specific summary
  summary_table                  <- dcast(data_transm, formula(paste('DATE', '~' ,'AGEGROUP')), value.var='CASES', sum)
  
  # update names
  names(summary_table)           <- c('date',paste('cases',names(summary_table)[-1],sep='_'))
  
  # merge
  summary_table <- merge(summary_table,summary_table_general)
  
  # omit NA's (default on LW's MACOS but not on VSC cluster)
  # R versions: 3.5.3 (MACOS) vs. 3.5.1 (VSC)
  summary_table <- na.omit(summary_table,cols='date')
  
  #check
  head(summary_table)
  tail(summary_table)
  
  # return
  return(summary_table)
}

#vect_cases <- incidence_run; vect_dates <- scm_dates
#vect_cases <- ref_case_data_age$cases;vect_dates <-as.Date(ref_case_data_age$date)
get_Rt <- function(vect_cases,vect_dates){
  
  # check for NA
  flag_na <- is.na(vect_cases)
  
  if(any(vect_cases<0)){
    warning('Negative number of cases! Not possible to calculate Rt')
    return(list(Rt=NA))
  }
  
  # use all data, if the last observation needs to be removed, this should be done in previous data cleaning
  cont.series = vect_cases
  cont.dates  = vect_dates
  
  # join
  dta.cont=data.frame(cont.dates,cont.series)
  names(dta.cont) <- c("dates", "I")
  
  # make sure the dates are consecutive
  dates.dummy <- data.frame(dates=seq(min(dta.cont$dates),max(dta.cont$dates),1))
  dta.cont  <- merge(dates.dummy,dta.cont,all.x = T)
  dta.cont$I[is.na(dta.cont$I)] <- 0 # set NA to "0"
  
  ## calculate effective reproduction number
  suppressMessages(
    Rt_mod <-
        estimate_R(
          dta.cont,
          method = "parametric_si",
          config = make_config(list(mean_si = 4.7, 
                                    std_si = 2.9)))
   )
  
    # add the reproduction number over the 7-day window finishing on that day.
    dta.cont$Rt <- NA
    dta.cont$Rt[Rt_mod$R$t_end] <- Rt_mod$R$`Mean(R)`
  
    # account for last observation
    dta.cont[nrow(dta.cont)+1,] <- dta.cont[nrow(dta.cont),]+1
    dta.cont[nrow(dta.cont),-1] <- NA
  
    # remove imputed dates
    if(length(vect_dates) < nrow(dta.cont)){
      dta.cont <- dta.cont[dta.cont$dates %in% vect_dates,]
    }
    
  return(dta.cont)
}

get_Rt_epilps <- function(vect_cases,vect_dates){
  
  # defensive programming
  if(!"EpiLPS" %in% installed.packages()){
    return(NA)
  }
  
  # check for NA
  flag_na <- is.na(vect_cases)
  
  if(any(vect_cases<0)){
    warning('Negative number of cases! Not possible to calculate Rt')
    return(list(Rt=NA))
  }
  
  # use all data, if the last observation needs to be removed, this should be done in previous data cleaning
  cont.series = vect_cases
  cont.dates  = vect_dates
  
  # mean 3 days and sd 2.9
  sidistr=discr_si(k=seq(0, 14), mu=3, sigma=2.9);sidistr=sidistr/sum(sidistr)
  
  # calculation Rt
  dta.rt <- epilps(incidence = cont.series, 
                     serial_interval = sidistr, verbose = FALSE)
  

  return(dta.rt)
}

recap_Rt <- function(output_dir,tag_incidence,tag_output){
  
  # get file name
  rds_filenames <- dir(output_dir,recursive = T,full.names = T)
  rds_filenames <- rds_filename_incidence[grepl(tag_incidence,rds_filename_incidence)]
  
  for(rds_filename_incidence in rds_filenames){
    # load data
    incidence <- readRDS(rds_filename_incidence)
    Rt_new    <- incidence*NA
    dim(incidence)
    for(i_sim in 1:ncol(incidence)){
      Rt_new[,i_sim] <- get_Rt(incidence[,i_sim],
                               1:nrow(incidence))$Rt
    }  
    
    # get output file name
    rds_filename_Rt <- gsub(tag_incidence,tag_output,rds_filename_incidence)
    
    # save output
    saveRDS(Rt_new,file = rds_filename_Rt)
    print(rds_filename_Rt)
  }
}

## calculate exponential growth rate
get_growth_rate <- function(f_obs, exp_window = 7){
  
  # exp_window <- 7  # days to calculate the growth
  out_r    = rep(NA,length(f_obs))
  
  for (i in exp_window:length(f_obs)){
    temp_df = as.data.frame(cbind(1:exp_window,f_obs[(i-exp_window+1):i]))
    names(temp_df)=c("day","cases")
    mylm = lm(log10(cases+1) ~ day, temp_df)
    out_r[i] <- coefficients(mylm)[[2]]
  }
  
  return(out_r)
}

get_reported_uptake_table <- function(tbl_dates,bool_proportion =  T){
  
  ref_data_be <- get_observed_incidence_data()
  
  ref_uptake_age <- ref_data_be[,grepl('_age',names(ref_data_be))]
  dim(ref_uptake_age)
  range(ref_data_be$date)
  ref_uptake_cum <- apply(ref_uptake_age,2,cumsum)
  dim(ref_uptake_cum)
  ref_uptake_tbl <- ref_uptake_cum[ref_data_be$date %in% tbl_dates,]
  ref_uptake_tbl
  dim(ref_uptake_tbl)
  
  pop_matrix  <- get_be_pop_matrix(nrow(ref_uptake_tbl))
  pop_matrix
  
  if(bool_proportion){
    ref_uptake_prop <- ref_uptake_tbl / cbind(pop_matrix,pop_matrix)
  } else{
    ref_uptake_prop <- ref_uptake_tbl
  }
  colnames(ref_uptake_prop) <- gsub('vaccine_A','dose_1',colnames(ref_uptake_prop))
  colnames(ref_uptake_prop) <- gsub('vaccine_B','dose_2',colnames(ref_uptake_prop))
  colnames(ref_uptake_prop) <- paste0(colnames(ref_uptake_prop),'_rep')
  ref_uptake_dates <- as.character(tbl_dates[tbl_dates %in% ref_data_be$date])
  return(data.frame(ref_uptake_dates,round(ref_uptake_prop,2)))
}

get_age_groups <- function(){
  age_breaks      <- seq(1,110,10)
  age_start       <- age_breaks[-length(age_breaks)]
  age_end         <- age_breaks[-1]-1
  age_tags        <- paste(age_start-1,age_end-1,sep="_")
  return(age_tags)
}

# CONVERT A CALENDAR DATE INTO A SIMULATION DAY INDEX
# default date format: YYYY-MM-DD
#note: March 1st, 2020 == "day 0"
sim_date2day <- function(x_date,x_format='%Y-%m-%d'){
  return(as.numeric(as.Date(as.character(x_date),x_format) - as.Date('2020-03-01')))
}

# CONVERT SIMULATION DAY INDEX INTO CALENDAR DATE
#note: March 1st, 2020 == "day 0"
sim_day2date <- function(x_day_index){
  return(as.Date('2020-03-01') + x_day_index)
}

# CONVERT SIMULATION TIME STEP INTO DAY INDEX
#note: the 1st step in 2020 == "step 1"
sim_step2day <- function(time_step){
  # return(floor(time_step/24))
  return(floor((time_step-1)/4))
}

# CONVERT SIMULATION TIME STEP INTO CALENDAR DATE
sim_step2date <- function(time_step){
  return(sim_day2date(sim_step2day(time_step)))
}


## Helper functions ----
##----------------- -
reduced_dates4 = c("Mar 1", "Apr 1", "May 1",
                   "Jun 1", "Jul 1", "Aug 1",
                   "Sept 1", "Oct 1", "Nov 1", 
                   "Dec 1", "Jan 1", "Feb 1",
                   "Mar 1", "Apr 1", "May 1",
                   "Jun 1", "Jul 1", "Aug 1",
                   "Sept 1")

t_col <- function(color, percent = 50, name = NULL) {
  rgb.val <- col2rgb(color)
  
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  invisible(t.col)
}

add_date_axis <- function(f_dates, side = 1, num_intervals = 12, bool_grid=FALSE){
  
  x_ticks <- pretty(as.Date(f_dates),num_intervals)
  x_labels <- format(x_ticks,'%d/%m')
  axis(side,x_ticks,x_labels)
  abline(v=x_ticks,lty=3,col='grey')
}
  
# scenario_data <- scenarioXa; col <- "blue"
add_polygon <- function(scenario_data,scenario_dates,col,alpha_val = 0.5,mean_lty=1){
  
  # get column id's
  id_mean  <- grepl('_mean',names(scenario_data))
  id_lower <- grepl('_LL',names(scenario_data))
  id_upper <- grepl('_UL',names(scenario_data))
  
  id_975  <- grepl('_975',names(scenario_data))
  id_025  <- grepl('_025',names(scenario_data))
  
  # remove NA's
  flag_na <- is.na(rowSums(scenario_data))
  if(any(flag_na)){
    scenario_data <- scenario_data[!flag_na,]
    scenario_dates <- scenario_dates[!flag_na]
  }
  
  # make sure that the dates are sorted
  scenario_data  <- scenario_data[order(scenario_dates),]
  scenario_dates <- scenario_dates[order(scenario_dates)]
  
  # add uncertainty interval
  polygon(x = c(scenario_dates,rev(scenario_dates)),
          y = c(scenario_data[,id_lower],rev(scenario_data[,id_upper])),
          col = alpha(col,alpha_val),
          border = NA)
  
  # add mean
  lines(x = scenario_dates,
        y = scenario_data[,id_mean],
        lty = mean_lty,
        col = col,
        lwd = 3)
  
  if(any(id_975) && any(id_025)){
    lines(x = scenario_dates,
          y = scenario_data[,id_975],
          lty = 2,
          col = col,
          lwd = 1)
    lines(x = scenario_dates,
          y = scenario_data[,id_025],
          lty = 2,
          col = col,
          lwd = 1)
  }
  
}

# add vertical line on given date + label on x-axis
add_vertical_line <- function(date_string,bool_text,date_tag = '',pos_factor = 0.05){
  
  plot_limits <- par("usr")
  
  v_date <- as.Date(date_string)
  abline(v=v_date,lty=2)
  #axis(1,v_date,format(v_date,'%d/%m'),
  #cex.axis=0.5,padj=-3,tck=-0.005)
  
  if(bool_text)
  {
    v_text <- ifelse(nchar(date_tag)>0,date_tag,format(v_date,'%d/%m'))
    text(x = v_date-1,
         y = mean(plot_limits[3:4])*pos_factor,
         #paste(format(v_date,'%d/%m'),date_tag),
         v_text,
         srt=90, pos=3, offset = +1.5,cex=0.6)
  }
}


# add copyright statement to the bottom right corner of your figure.
add_copyright <- function(text.cex = 0.4,bool_left = FALSE){
  plot_limits <- par("usr")
  x_value <- ifelse(bool_left,plot_limits[1],plot_limits[2])
  y_value <- plot_limits[3] + diff(plot_limits[3:4])/50
  x_adj <- ifelse(bool_left,-0.1,1)
  text(x = x_value,
       y = y_value,
       labels = ('© SIMID / UAntwerpen / UHasselt  '),
       cex = text.cex,
       adj = x_adj)
}

tocamelcase <- function(str_in){
  word_vect <- unlist(strsplit(str_in,' '))
  str_out <- ''
  for(sel_word in word_vect){
    str_out <- paste0(str_out,
                      paste0(toupper(substr(sel_word, 1, 1)), substr(sel_word, 2, nchar(sel_word))))
  }
  return(str_out)
}

totitle <- function(str_in){
  return(paste0(toupper(substr(str_in, 1, 1)), substr(str_in, 2, max(nchar(str_in)))))
}


# return model version as included in the README file 
get_scm_version <- function(filename_readme = './README.md'){
  
  readme <- readLines(filename_readme)
  str_version <- readme[grep("Version\\:", readme)]
  str_version <- unlist(strsplit(str_version,' '))
  
  return(as.double(str_version[length(str_version)]))
}

# parm_name <- i_param
get_ndays_sim <- function(parm_name,ndays_sim,extra_days=40){
  
  if(is.na(parm_name)){
    return(ndays_sim)
  }
  
  pname <- names(parm_name)
  if(!grepl('coef_w',pname)){
    return(ndays_sim)
  }
  wave_id <- gsub('.*coef_w','',pname)
  wave_id <- gsub('_.*','',wave_id)
  
  cp_wave <- ifelse(grepl('comix_coef_w',pname),
                    get_CoMix_change_day(as.numeric(wave_id)),
                    get_add_change_day(as.numeric(wave_id)))
  
  return(min(ndays_sim,cp_wave+extra_days))
}

adjust_q_param <- function(parms,parms_names,parms_names_estim,bool_setup = FALSE, 
                           opt_waves = 9:30){
  
  if(!bool_setup && !any(grepl('coef_aggr',parms_names_estim))){
    return(parms)
  }
  
  if(typeof(opt_waves) != 'list'){
    tbl_q_waves <- list(opt_waves)
  } else{
    tbl_q_waves <- opt_waves
  }
  
  length(tbl_q_waves)
  
  # option: initial setup of aggregated q-param
  if(bool_setup){
    for(i_q in 1:length(tbl_q_waves)){
      q_aggr_tag <- paste0('coef_aggr',i_q,'_age')
      sel_wave_init  <- tbl_q_waves[[i_q]][1]
      parms[paste0(q_aggr_tag,1:10)] <- parms[get_colnames(parms_names,tag_list = paste0('coef_w',sel_wave_init,'_age'))]
    }
  }
  
  # copy aggregated q-parameters to wave-specific q-parameters
  for(i_q in 1:length(tbl_q_waves)){
    q_aggr_tag <- paste0('coef_aggr',i_q,'_age')
    if(any(grepl(q_aggr_tag,parms_names_estim))){
      sel_waves <- tbl_q_waves[[i_q]]
      parms[get_colnames(parms_names,tag_list = paste0('coef_w',sel_waves,'_'))]  <- rep(parms[paste0(q_aggr_tag,1:10)],length(sel_waves))
    } 
  }
  
  return(parms)
}


## Region specific functions   ----
##---------------------- -

db_region <- c("belgium","brussels","flanders","wallonia")
get_region_id <- function(str_region){
  return(which(db_region == str_region))
}

get_region <- function(id_region){
  return(db_region[id_region])
}


## Uptake functions   ----
##---------------------- -

# note: time_step is best expressed in "days" (or at least the same unit)
# aggr_uptake<-scen_uptake_day;target_time_step<-time_step_size
#aggr_uptake <- plan_uptake_month;target_time_step<-1;dose_col<-c('A_total','B_total')
approx_uptake <- function(aggr_uptake, target_time_step, dose_col,max_date = NA){
  
  # set output dates
  x_out <- seq(min(aggr_uptake$date),max(c(aggr_uptake$date,max_date)+(1/target_time_step - 1)*target_time_step,na.rm = T),target_time_step)
  
  # get current aggregation
  aggr_time_steps <- as.numeric(c(diff(aggr_uptake$date),1))
  
  # prepare final data.frame and add dates
  target_uptake <- data.frame(matrix(0,nrow=length(x_out),ncol=(1+length(dose_col)*2)))
  names(target_uptake) <- c('date',
                            paste0('adeno_',dose_col),
                            paste0('rna_',dose_col))
  target_uptake$date <- x_out
  
  sel_column <- names(aggr_uptake)[2]
  for(sel_column in names(aggr_uptake)[-1]){
    
    aggr_uptake[,sel_column]
    aggr_uptake$date
    
    # linear extrapolate the doses over the days of the month
    approx(x = aggr_uptake$date,
           y = aggr_uptake[,sel_column] / aggr_time_steps * target_time_step,
           xout = x_out,
           yright = aggr_uptake[nrow(aggr_uptake),sel_column] / aggr_time_steps[length(aggr_time_steps)] * target_time_step,
           method='constant') -> aggr_uptake_hour
    
    target_uptake[,sel_column]   <- aggr_uptake_hour$y
  }
  
  return(target_uptake)
}
# uptake_time <- scen_out;delay_protection_rna<-delay_rna_dose2;delay_protection_adeno<-delay_adeno_dose2
uptake2protection <- function(uptake_time,
                              delay_protection_rna,
                              delay_protection_adeno){
  # format and make sure to work with unique date ids
  uptake_time$date        <- as.Date(uptake_time$date)
  if(length(unique(uptake_time$date)) != nrow(uptake_time)){
    # (re)set output dates (ERROR PRONE!)
    target_time_step        <- sum(uptake_time$date[1] == uptake_time$date)
    uptake_time$date <- seq(min(uptake_time$date),max(uptake_time$date)+(target_time_step-1)/target_time_step,1/target_time_step)
  }


  # select uptake: rna-based vaccine
  uptake_time_rna      <- uptake_time[,grepl('rna',names(uptake_time))]
  uptake_time_rna$date <- uptake_time$date + delay_protection_rna

  # select uptake: adeno-based vaccine
  uptake_time_adeno      <- uptake_time[,grepl('adeno',names(uptake_time))]
  uptake_time_adeno$date <- uptake_time$date + delay_protection_adeno

  # merge datasets
  uptake_time_delay <- data.frame(date       = c(uptake_time$date[-nrow(uptake_time)],
                                                 seq(max(uptake_time$date),max(uptake_time$date)+30,1/4)),
                                  dummy = 1)
#  uptake_time_delay <- merge(uptake_time_delay,uptake_time_adeno,by='date',all.x = T,all.y=F)
 # uptake_time_delay <- merge(uptake_time_delay,uptake_time_rna,by='date',all.x = T,all.y=F)
  # correction from R 4.3 using dplyr
   uptake_time_delay <- left_join(uptake_time_delay,uptake_time_adeno,by='date')
   uptake_time_delay <- left_join(uptake_time_delay,uptake_time_rna,by='date')  
  
  names(uptake_time_delay)
  
  uptake_time_delay$dummy <- NULL
  
  # set imputed NA's as 0
  uptake_time_delay[is.na(uptake_time_delay)] <- 0
  
  return(uptake_time_delay)
}

get_uptake_matrix <- function(vaccine_uptake, dose_tag, parms ){
  
  
  if(dose_tag == '_A_'){
    delay_protection_rna   <- parms['delay_protection_rna1']    
    delay_protection_adeno <- parms['delay_protection_adeno1']
  } else if(dose_tag == '_B_'){
    delay_protection_rna   <- parms['delay_protection_rna2']    
    delay_protection_adeno <- parms['delay_protection_adeno2']
  } else if(dose_tag == '_E_'){
    delay_protection_rna   <- parms['delay_protection_rna2']    
    delay_protection_adeno <- parms['delay_protection_adeno2']
  } else if(dose_tag == '_F_'){
    delay_protection_rna   <- parms['delay_protection_rna2']    
    delay_protection_adeno <- parms['delay_protection_adeno2']
  } else if(dose_tag == '_X_'){
    delay_protection_rna   <- parms['delay_protection_rna2']    
    delay_protection_adeno <- parms['delay_protection_adeno2']
  } else{
    stop("VACCINE UPTAKE DOSE UNKNOWN:",dose_tag)
  }

  # check of given dose_tag is present, else, use dummy and set all info to 0
  if(!any(grepl(dose_tag,names(vaccine_uptake)))){
    dose_tag <- '_A_'
    vaccine_uptake[,grepl(dose_tag,names(vaccine_uptake))] <- 0
  }
  
  # protection (after uptake)
  vaccine_protection <- uptake2protection(vaccine_uptake[,grepl('date',names(vaccine_uptake)) | grepl(dose_tag,names(vaccine_uptake))],
                                          delay_protection_rna,
                                          delay_protection_adeno)
  V_mat <- rbind(matrix(0, nrow = 305/parms['h'], ncol = ncol(vaccine_protection[,-1])),
                 as.matrix(vaccine_protection[,-1]))
  
  
  return(V_mat)
}


plot_uptake_by_age <- function(scen_out_f,
                               endpoint_reported_data = NA,
                               legend_label = 'Uptake'){
  
  age_opt     <- get_age_groups()
  ref_age_col <- paste0('vaccine_A_ages',age_opt)
  y_lim       <- c(0,1.5e6)
  x_ticks     <- pretty(range(scen_out_f$date),10)
  
  perc_breaks    <- pretty(0:100,7)
  y_ticks_perc   <- paste0(perc_breaks,'%')
  par(mar=c(5,5,4,5),mfrow=c(2,2))
  
  ref_data_be <- get_observed_incidence_data()
  
  i_age <- 5
  for(i_age in 10:3){
    
    y_ticks      <- perc_breaks * pop_data_be[i_age] / 100
    y_ticks_k    <- paste0(round(y_ticks/1e3),'k')
    y_ticks_k[1] <- 0
    y_ticks_k
    
    plot(ref_data_be$date+1,
         t(cumsum(ref_data_be[ref_age_col[i_age]])),
         ylab='Uptake',
         xlab='',
         xlim= range(scen_out_f$date),
         ylim = range(y_ticks),
         xaxt='n',
         yaxt='n',
         main=paste0(gsub('_','-',age_opt[i_age]),'y'))
    lines(scen_out_f$date+1,
          cumsum(scen_out_f[,11+i_age]),
          col=4,
          lwd=2)
    axis(1,x_ticks,format(x_ticks,'%d/%m'))
    axis(2,y_ticks,y_ticks_perc,las=2)
    axis(4,y_ticks,y_ticks_k,las=2)
    abline(v=x_ticks,lty=3,col='grey')
    abline(h=y_ticks,lty=3,col='grey')
    
    legend('topleft',
           c('Reported uptake',
             paste('SCM:', tolower(legend_label))),
           pch=c(1,NA),
           lwd=c(NA,2),
           col=c(1,4),
           cex=0.6,
           bg='white'
    )
    
    if(!is.na(endpoint_reported_data)){
      abline(v=endpoint_reported_data,lwd=2,col=3)  
      text(x = endpoint_reported_data,
           y = max(y_ticks)*0.98,
           'reported',
           pos=2,
           col=3,
           cex=0.7)
      text(x = endpoint_reported_data,
           y = max(y_ticks)*0.98,
           'scenario',
           pos=4,
           col=3,
           cex=0.7)
    }
    
    # text(x = mean(x_ticks[4:5]),
    #      y = min(y_ticks),
    #      srt=90,
    #      '65+',
    #      pos=4,
    #      col=3,
    #      cex=0.7)
    # 
    # text(x = mean(x_ticks[5:6]),
    #      y = min(y_ticks),
    #      srt=90,
    #      'risk pop',
    #      pos=4,
    #      col=3,
    #      cex=0.7)
    # 
    # text(x = mean(x_ticks[6:7]),
    #      y = min(y_ticks),
    #      srt=90,
    #      'gen pop',
    #      pos=4,
    #      col=3,
    #      cex=0.7)
  }
}

## VSC functions ----
is_vsc <- function(){
  return(nchar(system('echo $VSC_HOME',intern = T))>0)
}

get_vsc_scratch_folder <- function(){
  return(system('echo $VSC_SCRATCH',intern = T))
}

is_vsc_vaughan_cluster <- function(){
  return(system('echo $VSC_INSTITUTE_CLUSTER',intern = T) == "vaughan")
}

## WGS functions ----
approx_sequenced <- function(date,sequenced_week){
  
  out <- approx(x = date[!is.na(sequenced_week)],
                y = sequenced_week[!is.na(sequenced_week)],
                xout = seq(min(date),max(date),1))
  
  return(data.frame(x = as.character(out$x),
                    y = round(out$y)))
}

## HOSPITAL ADMISSION functions ----
# hospital admission rate: DEFAULT
# ==>> (1 - exp(-h*(1-phi0)*delta2)
#
# hospital admission rate: WITH VOC:
# ==>> (1 - exp(-h*((1-phi0) + phi0*VOC_hosp)*delta2)
#
# Odds ratio for hospital admission: ONE VOC
#     (1 - exp(-h*((1-phi0) + phi0*VOC_new)*delta2)) /
#     (1 - exp(-h*(1-phi0)*delta2))
#
# Odds ratio for hospital admission: TWO VOCs
#     (1 - exp(-h*((1-phi0) + phi0*VOC_new)*delta2)) /
#     (1 - exp(-h*((1-phi0) + phi0*VOC_previous)*delta2))
get_log_VOC_hosp <- function(phi0,delta2,h,odds_ratio,log_VOC_previous=NA){
  
  # defensive check
  if(length(odds_ratio) > 1 && length(odds_ratio) != length(phi0)){
    stop("THE PROVIDED HOSPITAL ADMISSION ODDS RATIO IS NOT VALID")
  }
  
  # define optional secondary VOC factor
  VOC_previous <- ifelse(is.na(log_VOC_previous),0,expit(log_VOC_previous))
  
  # get VOC hospital admission factor
  VOC_new <- log(1- odds_ratio * (1 - exp(-h*((1-phi0)+phi0*VOC_previous)*delta2)))/(-h*delta2*phi0) - (1-phi0)/phi0
  
  # use cutoff of 1e-15, to prevent NaN with logit
  VOC_new[VOC_new <= 1e-15] <- 1e-15
  
  # convert to logit
  log_VOC_new <- logit(VOC_new)
  
  return(log_VOC_new)
}

# help function
get_VOC_hosp <- function(phi0,delta2,h,odds_ratio,log_VOC_previous=NA){
  return(expit(get_log_VOC_hosp(phi0,delta2,h,odds_ratio,log_VOC_previous)))
}
  
# see comments with "get_log_VOC_hosp"
#log_VOC_new <- log_VOC_delta_hosp; log_VOC_previous <- log_VOC_delta_hosp
get_VOC_hosp_odds_ratio <- function(phi0,delta2,h,log_VOC_new,log_VOC_previous=NA){
  
  VOC_new       <- expit(log_VOC_new)
  
  if(any(is.na(log_VOC_previous))){
    VOC_previous <- 0
  } else {
    VOC_previous <- expit(log_VOC_previous)
  }
  
  OR <- as.numeric((1 - exp(-h*((1-phi0) + phi0*VOC_new)*delta2)) /
                     (1 - exp(-h*((1-phi0) + phi0*VOC_previous)*delta2)))
  
  return(OR)
}

