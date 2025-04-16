########################################################################### #
# This file is part of the Stochastic Compartmental Model for SARS-COV-2 
# transmission in Belgium, conceived by members of SIMID group during the 
# COVID19 pandemic.
#
# This file is used to generate vaccine uptake files with booster doses.
# Used in TechNote v20220104 and future TechNotes
#
# Copyright 2021, SIMID, University of Antwerp & Hasselt University                                        
########################################################################## # 

# clear workspace
# parse command line arguments
cl_args = commandArgs(trailingOnly=TRUE)

# clear workspace (except the command line arguments)
rm(list=ls()[ls()!='cl_args'])

library(scales)
library(data.table)
library(openxlsx)
library(foreach)
library("lubridate")

source('R/lib_stochastic_model.R')


# SETTINGS          ----
###################### # 
# This particular file is for generating an extra (2nd of 3rd) booster with different characteristics.

# set scenario tag
if(length(cl_args)>=1){ # region
cl_val <- as.numeric(unlist(strsplit(cl_args,'_')))
scen_num <- cl_val[1]
bool_bivariant_fall2023 <- as.logical(cl_val[2])
} else{
  scen_num <- 3
  bool_bivariant_fall2023 <- TRUE
}
print(scen_num)
print(bool_bivariant_fall2023)

ECDC_uptake <- read.table('data/Round_5_vaccination_uptake.csv',sep=',',header=T)
ECDC_uptake <-  ECDC_uptake[ECDC_uptake$location=="BE",]

sens_flu_pop <- numeric(10)


if (scen_num==1){
  scen_tag <- "ECDC_R5_fall23campaign_0flu"
  #sens_coef <- 0
}
if (scen_num==2){
  scen_tag <- "ECDC_R5_fall23campaign_50flu"
  sens_flu_pop[7:9] <- ECDC_uptake$value_covid_pessimistic
  sens_flu_pop[10] <- ECDC_uptake$value_covid_pessimistic[3]
#sens_coef <- 0.5
}
if (scen_num==3){
  scen_tag <- "ECDC_R5_fall23campaign_100flu"
  sens_flu_pop[7:9] <- ECDC_uptake$value_flu
  sens_flu_pop[10] <- ECDC_uptake$value_flu[3]
  #sens_coef <- 1
}
if (scen_num==4){
  scen_tag <- "ECDC_R5_fall23campaign_100flu_plus15"
  sens_flu_pop[7:9] <- ECDC_uptake$value_covid_optimistic
  sens_flu_pop[10] <- ECDC_uptake$value_covid_optimistic[3]
  #sens_coef <- 1.15
}
if (scen_num==5){
  scen_tag <- "ECDC_R5_fall23campaign_all60plus"
  sens_flu_pop[7:10] <- rep(1,4)
  #sens_coef <- 1
}

if (bool_bivariant_fall2023) {
  scen_tag <- paste0(scen_tag,"_biv")}
if (!bool_bivariant_fall2023) {
  scen_tag <- paste0(scen_tag,"_nobiv")}

#scen_tag <- "nocampaign200323"
#scen_tag <- 'onecampaign160922'
#scen_tag <- '6mcampaign160922'
#scen_tag <- '65pall280822'

#timing_booster_finished <- as.Date('2022-5-15')
timing_booster_finished <- as.Date('2024-12-31')
timing_booster_start <- as.Date('2023-10-01')

# recurrent campaign every x months
rec_campaign <- 1000

# endpoint uptake
endpoint_uptake <- timing_booster_start + 90

# output folder
output_folder <- paste0('output/uptake_',scen_tag)

# planned uptake file
#uptake_file <- 'data/OneDrive_1_07-04-2021/20200407_Input Modellers.xlsx'
#uptake_file <- 'data/OneDrive_1_28-10-2021_dummy/2021_dummy_Input Modellers.xlsx'

# endpoint to select reported data
endpoint_reported_data <- Sys.Date()
#endpoint_reported_data <- as.Date('2022-08-19')



# sensitivity on uptake
#sens_uptake_cap  <- 0.85 # default: 0.8 
# sens_uptake_rate <- 1.0
bool_threshold <- FALSE
bool_7d_avg_only <- TRUE # base future uptake only on last 7days
bool_incr_child_uptake <- FALSE

sel_region <- "belgium" # belgium, flanders, wallonia, brussels
sens_uptake_cap        <- 0#0.80
sens_uptake_cap_child  <- 0# 0.40
sens_uptake_booster    <- 0 #0.6
sens_uptake_2ndbooster    <- 0 #0.5 #0.6
sens_uptake_extrabooster_old    <- 0.25 # !!! this is absolute, not relative to previous booster
sens_uptake_extrabooster_young    <- 0.15 # !!! this is absolute, not relative to previous booster

#sens_vl_flu <- c(0.469,0.589,0.673,0.722,0.75) #percentage of flu vaccine uptake in flanders in 2019, https://www.laatjevaccineren.be/sites/default/files/2022-12/Rapport_griep65plus_MR20221012_0.pdf


# add to scen_tag 
scen_tag         <- paste0(scen_tag,
                      #     '_uadult',round(sens_uptake_cap*100),
                       #    '_uchild',round(sens_uptake_cap_child*100),
                        #   '_ubooster',round(sens_uptake_booster*100),
                         #  '_',tolower(format(timing_booster_finished,'%b%d')),
                        #   ifelse(bool_incr_child_uptake,'_INFANT',''),
                           '_',sel_region)

if(bool_threshold) { scen_tag     <- gsub('_ad','t_ad',scen_tag) }

# delay for second dose
delay_rna_dose2 <- 3*7
delay_adeno_dose2 <- 0 #12*7 # only J&J in the last months

# note: since JnJ is incorporated as AZ... it requires an arbitraty delay for  
# the virtual second JnJ dose to account for the time to build up protection 
# of the 1st dose (of e.g. 21 days).
delay_JnJ         <- 21 

# delay booster: single dose, so not applicable
delay_booster     <- 4*30 #round(365/2)

# time-step size
time_step_size <- 1/4  # 1 == per day and 1/4 == per 6 hours 

# population data
 pop_data_be  <- get_regional_pop(sel_region)
# pop_be2020_all <- read.table('data/pop_be2020_statbel_download.csv',sep=',',header=T)
# unique(pop_be2020_all$Region)
# pop_be2020_all <- pop_be2020_all[nchar(pop_be2020_all$Region)==0,]
# pop_be2020_all$Population.on.January.1st.2020[nrow(pop_be2020_all)-1] <- sum(pop_be2020_all$Population.on.January.1st.2020[nrow(pop_be2020_all)-(1:0)])
# pop_be2020_all <- pop_be2020_all[-nrow(pop_be2020_all),]
# plus65fraction=pop_be2020_all$Population.on.January.1st.2020[14]/(pop_be2020_all$Population.on.January.1st.2020[14]+pop_be2020_all$Population.on.January.1st.2020[13])



# reference data, and select uptake data
ref_data_be      <- get_observed_incidence_data(enable_data_reuse = T, #temp true
                                                sel_region=sel_region,
                                                bool_incr_child_uptake=bool_incr_child_uptake)
ref_data_be$date <- as.Date(ref_data_be$date)
db_uptake        <- ref_data_be[,grepl('vaccine_.*_._',names(ref_data_be)) | grepl('date',names(ref_data_be))]
agegroup_opt     <- unique(gsub('.*_*_ages','',names(db_uptake))[-1])
dose_col         <- paste0('A_',agegroup_opt)
dose_col_2nd     <- paste0('B_',agegroup_opt)
dose_col_booster <- paste0('E_',agegroup_opt)
dose_col_2ndbooster <- paste0('F_',agegroup_opt)
dose_col_extrabooster <- paste0('X_',agegroup_opt)   #dedicated booster bivalent
#dose_col_fallcampaign <- paste0('X_',agegroup_opt) 

# merge first dose of adeno-based vaccines
db_uptake[is.na(db_uptake)] <- 0
db_uptake[,grepl('AZ_A',names(db_uptake))] <- db_uptake[,grepl('AZ_A',names(db_uptake))] + db_uptake[,grepl('JnJ',names(db_uptake))]
db_uptake[-(1:delay_JnJ),grepl('AZ_B',names(db_uptake))] <- db_uptake[-(1:delay_JnJ),grepl('AZ_B',names(db_uptake))] + db_uptake[1:(nrow(db_uptake)-delay_JnJ),grepl('JnJ',names(db_uptake))]
db_uptake[,grepl('JnJ',names(db_uptake))] <- NULL
names(db_uptake) <- gsub('AZ','adeno',names(db_uptake))

# rename mRNA
names(db_uptake) <- gsub('mRNA','rna',names(db_uptake))
names(db_uptake)

# remove 'vaccine'
names(db_uptake) <- gsub('vaccine_','',names(db_uptake))

# select between Jan 1st and selected date
db_uptake <- db_uptake[db_uptake$date >= as.Date('2021-01-01'),]  #why is this present?
db_uptake <- db_uptake[db_uptake$date <= endpoint_uptake,]
print(range(db_uptake$date))

# change column names
names(db_uptake) <- gsub('ages','',names(db_uptake))

# split 1st and 2nd uptake
db_uptake_first   <- db_uptake[,grepl('_A_',names(db_uptake)) | grepl('date',names(db_uptake))] 
db_uptake_second  <- db_uptake[,grepl('_B_',names(db_uptake))| grepl('date',names(db_uptake))] 
db_uptake_booster <- db_uptake[,grepl('_E_',names(db_uptake))| grepl('date',names(db_uptake))] 
db_uptake_2ndbooster <- db_uptake[,(grepl('_F_',names(db_uptake)) & !grepl('adeno',names(db_uptake)))| grepl('date',names(db_uptake))] 
db_uptake_extrabooster <- db_uptake[,(grepl('_X_',names(db_uptake)) & !grepl('adeno',names(db_uptake)))| grepl('date',names(db_uptake))]
#db_uptake_fallcampaign <- db_uptake_extrabooster
#names(db_uptake_fallcampaign) <- gsub('_X_','_X2_',names(db_uptake_fallcampaign))
#db_uptake_fallcampaign[-1] <- 0

names(db_uptake_first)
names(db_uptake_second)
names(db_uptake_booster)
names(db_uptake_2ndbooster)
names(db_uptake_extrabooster)
sum(db_uptake_first[,-1])
sum(db_uptake_second[,-1])
sum(db_uptake_booster[,-1])
sum(db_uptake_2ndbooster[,-1])
sum(db_uptake_extrabooster[,-1])

sum(ref_data_be$vaccine_uptake_first)
sum(ref_data_be$vaccine_uptake_full)
sum(ref_data_be$vaccine_uptake_booster)
sum(ref_data_be$vaccine_uptake_2ndbooster)
sum(ref_data_be$db_uptake_extrabooster)

sum(db_uptake_first$rna_A_80_89) + sum(db_uptake_first$adeno_A_80_89)
sum(ref_data_be$vaccine_A_ages80_89)

sum(db_uptake_first$rna_A_90_99) + sum(db_uptake_first$adeno_A_90_99)
sum(ref_data_be$vaccine_A_ages90_99)

# PREAMBLE

# if output folder does not exit yet, create one
if(!dir.exists(output_folder)){
  dir.create(output_folder,recursive = T)
}

# UPTAKE SCENARIOS  ----
###################### # 

# parse planned uptake file

# if(nchar(uptake_file)>0){
#   uptake_file_row_start <- c(5,17,29)
#   uptake_file_nrow      <- 6
#   uptake_scenario       <- 1
#   input_uptake_month <- read.xlsx(uptake_file,sheet = "Scenarios",
#                                   rows = uptake_file_row_start[uptake_scenario] + c(1:uptake_file_nrow-1),
#                                   # cols = c(2,4:16),
#                                   cols = c(2,4:18),
#                                   detectDates=T)
#   
#   uptake_first_total    <- colSums(input_uptake_month[1:3,-1],na.rm = T)
#   uptake_second_total   <- colSums(input_uptake_month[4:5,-1],na.rm = T)
#   
#   uptake_dates <- as.Date(names(input_uptake_month)[-1])
#   uptake_dates
#   
#   # add additional point
#   uptake_dates        <- c(uptake_dates,max(uptake_dates)+31)
#   uptake_first_total  <- c(uptake_first_total,0)
#   uptake_second_total <- c(uptake_second_total,0)
#   
#   # go from 'month' to 'day'
#   plan_uptake_month    <- data.frame(date  = uptake_dates,
#                                      dose_A_total = uptake_first_total*10,
#                                      dose_B_total = uptake_second_total*10)
#   plan_uptake_day <- approx_uptake(plan_uptake_month,1,c('A_total','B_total'),max_date = endpoint_uptake)
# } 
# dim(plan_uptake_day)
# head(plan_uptake_day)
# tail(plan_uptake_day)

# # NEW: option to adapt the uptake in may/june
# if(sens_uptake_rate != 1){
# 
#   sel_date <- plan_uptake_day$date > endpoint_reported_data
#   plan_uptake_day$dose_A_total[sel_date] <- plan_uptake_day$dose_A_total[sel_date] * sens_uptake_rate
#   lines(plan_uptake_day$date,
#         cumsum(plan_uptake_day$dose_A_total),
#         lwd=2,
#         col=6,
#          type='l')
# 
#   legend('topleft',
#          c('Reported uptake: 1st or single dose',
#            'Planned uptake: 1st or single dose (month)',
#            'Planned uptake: 1st or single dose (day)',
#            scen_tag),
# 
#          col = c(1,3,4,6),
#          lwd= c(NA,2,2,2),
#          pch=c(1,NA,NA,NA),
#          cex =0.8,
#          bg='white')
#   
# }


# UPTAKE JAN-OKT: first dose  ----
###################### # 

#uptake_matrix_date <- plan_uptake_day$date
#uptake_matrix      <- data.frame(matrix(0,
#                                 ncol=ncol(db_uptake_first)-1, 
#                                 nrow = nrow(plan_uptake_day)))
#names(uptake_matrix) <- names(db_uptake_first)[-1]

# set target uptake
#target_uptake_age_fraction <- c(5/10*sens_uptake_cap_child,
#                                2/10*sens_uptake_cap_child + 8/10*sens_uptake_cap ,
#                                rep(1,8)*sens_uptake_cap)
#target_uptake_age_count    <- target_uptake_age_fraction * pop_data_be

# untill specified date: reported
#uptake_matrix[plan_uptake_day$date %in% db_uptake_first$date,] <- db_uptake_first[,grepl('_A_',names(db_uptake_first))]

# UPTAKE OKT-...  ----
###################### # 
# next (use aggregated uptake mRNA + adeno to define max)
# uptake_matrix_proj   <- uptake_matrix[1:10] + uptake_matrix[11:20]
# remaining_uptake     <- target_uptake_age_count - colSums(uptake_matrix_proj)
# pop_uptake_reported  <- colSums(uptake_matrix_proj) / pop_data_be
# 
# # if regional uptake, scale by population size
# pop_be     <- get_regional_pop(region = 'belgium')
# pop_model  <- get_regional_pop(region = sel_region)
# pop_factor   <- sum(pop_model)/sum(pop_be)
# 
# avg_7d_uptake <- colMeans(db_uptake[(-6:0)+nrow(db_uptake),-1])
# avg_7d_uptake <- avg_7d_uptake[1:10] + avg_7d_uptake[21:30]
# avg_7d_uptake[1:2] <- 1e4 * pop_factor
# 
# # reported
# i_day <- nrow(db_uptake)+1
# plan_uptake_day[i_day,]
# #for(i_day in 110:120){
# for(i_day in (nrow(db_uptake)+1):nrow(plan_uptake_day)){
# #for(i_day in 110:119){
#     
#   uptake_age             <- c(0,0,rep(0,4),0,1,1,1)
#   
#   # additional loop to shift doses if one age group has reached the max uptake
#   while(sum(uptake_matrix_proj[i_day,]) < plan_uptake_day$dose_A_total[i_day]){
#     available_doses        <- plan_uptake_day$dose_A_total[i_day] - sum(uptake_matrix_proj[i_day,])
#     remaining_uptake       <- target_uptake_age_count - colSums(uptake_matrix_proj) 
#     
#     # if reported uptake > target... make sure the maximum is used
#     remaining_uptake <- pmax(0,remaining_uptake)
#     
#     if(sum(remaining_uptake) == 0){
#       break
#     }
#     
#     # from June: General population, top down
#     if(plan_uptake_day$date[i_day] >= as.Date('2021-05-01') ||
#        sum(remaining_uptake * uptake_age) == 0){
#       
#       target_age <- min(length(remaining_uptake),which(remaining_uptake == 0)-1)
#       uptake_age <- c(rep(0,length(remaining_uptake)))
#       uptake_age[target_age:length(remaining_uptake)] <- 1
#       #print(uptake_age)
#     }
#   
#     # add default uptake ONCE, and account for this in the available doses  
#     if(!any(uptake_matrix_proj[i_day,]>0)){
#       standard_uptake_adjusted <- (avg_7d_uptake * (remaining_uptake>0))
#       available_doses          <- available_doses - sum(standard_uptake_adjusted)
#       available_doses          <- pmax(0,available_doses)
#     } else{
#       standard_uptake_adjusted[] <- 0
#     }
#   
#     if(bool_7d_avg_only){
#       plan_uptake_day$dose_A_total[i_day] <- sum(standard_uptake_adjusted)
#       available_doses <- 0
#     }
#     
#     remaining_uptake_day   <- remaining_uptake * uptake_age
#     planned_uptake         <- available_doses * (remaining_uptake_day / sum(remaining_uptake_day))
#     
#     planned_uptake         <- planned_uptake + standard_uptake_adjusted
#     planned_uptake         <- pmin(remaining_uptake,planned_uptake)  
#     
#     uptake_matrix_proj[i_day,]  <- uptake_matrix_proj[i_day,] + planned_uptake 
#   }
#   #print(uptake_matrix_proj[i_day,])
# }

# merge reported and projected uptake
#flag_projected <- plan_uptake_day$date > max(db_uptake$date)
#uptake_matrix[flag_projected,grepl('rna',names(uptake_matrix))] <- uptake_matrix_proj[flag_projected,]

# STRICT THRESHOLD ? ----
# if(bool_threshold){
#   
#   uptake_age_threshold <- pmin(sens_uptake_cap,target_uptake_age_fraction) * pop_data_be
#   
#   i_age <- 1
#   for(i_age in 1:length(uptake_age_threshold)){
#     
#     col_age <- i_age+c(0,10)
#     uptake_age_cum <- cumsum(rowSums(uptake_matrix[col_age]))
#     bool_overflow <- uptake_age_cum > uptake_age_threshold[i_age]
#     uptake_matrix[bool_overflow,col_age] <- 0
#   }
#   
# }



# REFORMAT UPTAKE   ----
###################### # 

# merge date and uptake
#scen_uptake_day    <- data.frame(date = uptake_matrix_date,
 #                                uptake_matrix)
#dim(scen_uptake_day)
#names(scen_uptake_day)
# set to 1/X time steps
#scen_out <- approx_uptake(scen_uptake_day,time_step_size,dose_col)

# write to file
#write.table(scen_out,file=file.path(output_folder,paste0('vaccine_uptake_Vx_vSCEN',scen_tag,'.csv')),sep=',',row.names=F)

# NEW UPTAKE FIRST DOSE (reported only) ----
###################### #

scen_out <- approx_uptake(db_uptake_first,time_step_size,dose_col)
scen_out_delay <- data.frame(date       = c(scen_out$date[-nrow(scen_out)],
                                               seq(max(scen_out$date),timing_booster_finished+3*time_step_size,time_step_size)),
                                dummy = 1)
scen_out <- left_join(scen_out_delay,scen_out,by='date')
names(scen_out)

scen_out$dummy <- NULL

# UPTAKE SECOND DOSE ----
###################### #

scen_out_2dose <- uptake2protection(scen_out,delay_rna_dose2,delay_adeno_dose2)
names(scen_out_2dose) <- gsub('_A_','_B_',names(scen_out_2dose))

# # untill specified date: reported
# set to 1/X time steps
reported_second <- approx_uptake(db_uptake_second,time_step_size,dose_col_2nd)
scen_out_2dose[scen_out_2dose$date %in% reported_second$date,-1] <- reported_second[,-1]

# WARNING
# In practice, some 2nd doses have been administrated before the specified timings
# in this procedure, which may cause that one dose is included as reported and scheduled.
# Therefore, limit the 2nd doses to the cumulative number of 1 doses over time.
i_col <- 22
for(i_col in 2:ncol(scen_out_2dose)){
  d_overflow <- which(cumsum(scen_out_2dose[1:nrow(scen_out),i_col]) > cumsum(scen_out[,i_col]))
  if(length(d_overflow)>0){
    print(i_col)
    scen_out[d_overflow,i_col] <- apply(cbind(scen_out[d_overflow,i_col],scen_out_2dose[d_overflow,i_col]),1,max)
    # scen_out_2dose[d_overflow,i_col] <- 0
    #scen_out_2dose[min(d_overflow),i_col] <- sum(scen_out[,i_col]) - sum(scen_out_2dose[,i_col])
  }
}

dim(scen_out)
scen_out_2dose <- left_join(scen_out,
                        scen_out_2dose,
                        by="date")
dim(scen_out_2dose)
scen_out_2dose[is.na(scen_out_2dose)] <- 0
scen_out_2dose <- scen_out_2dose[c(1,order(names(scen_out_2dose)[-1])+1)]
names(scen_out_2dose)

write.table(scen_out_2dose,file=file.path(output_folder,paste0('vaccine_uptake_2doses_Vx_vSCEN',scen_tag,'.csv')),sep=',',row.names=F)

# UPTAKE BOOSTER DOSE ----
###################### #

# define full uptake
scen_out_latest <- scen_out_2dose[,!grepl('_A_',names(scen_out_2dose))]
scen_out_latest[,grepl('rna_',names(scen_out_latest))] <- scen_out_latest[,grepl('rna_',names(scen_out_latest))] + scen_out_latest[,grepl('adeno_',names(scen_out_latest))]
scen_out_latest[,grepl('adeno_',names(scen_out_latest))] <- 0
names(scen_out_latest)

# schedule booster doses
scen_out_booster <- uptake2protection(scen_out_latest,delay_booster,delay_booster)
names(scen_out_booster) <- gsub('_B_','_E_',names(scen_out_booster))
tail(scen_out_booster)

# save copy of original planned booster uptake
scen_out_booster_db <- scen_out_booster

# # untill specified date: reported booster uptake
# set to 1/X time steps
reported_booster <- approx_uptake(db_uptake_booster,time_step_size,dose_col_booster)
scen_out_booster[scen_out_booster$date %in% reported_booster$date,-1] <- reported_booster[,-1]
scen_out_booster <- scen_out_booster[scen_out_booster$date <= max(scen_out_latest$date),]

# add booster doses until specified data
i_col <- 12
sel_booster_days <- scen_out_booster$date >= endpoint_reported_data & scen_out_booster$date < (timing_booster_finished-1)
num_time_steps <- sum(sel_booster_days)
booster_total <- 0
if(num_time_steps > 0){ # implement future booster campaign
  for(i_col in 2:ncol(scen_out_booster)){
    scen_out_booster[scen_out_booster$date > endpoint_reported_data,i_col] <- 0
    target_booster <- scen_out_latest[scen_out_latest$date<=endpoint_reported_data,i_col] * sens_uptake_booster
    
    if (grepl(paste0('_',agegroup_opt[1]),names(scen_out_booster)[i_col])){
      # target_booster[scen_out_booster$date > endpoint_reported_data] <- 0
      target_booster <- 0
    }
    remaining_booster <- pmax(0,sum(target_booster) - sum(scen_out_booster[,i_col]))
    scen_out_booster[sel_booster_days,i_col] <- remaining_booster/num_time_steps
    booster_total <- booster_total + remaining_booster
  }
} else{ # adjust historical booster campaign
  
  # # target for booster uptake = update 2nd dose
  # if(sens_uptake_booster == 1){
  #  cat("BOOSTER UPTAKE: INCREASE TO 2DOSE LEVEL")
  #   age_opt     <- get_age_groups()
  #   s_age <- age_opt[3]
  #   for(s_age in age_opt){
  #     booster_ratio <- sum(colSums(scen_out_booster[,grepl(paste0('E_',s_age),names(scen_out_booster))])) /
  #                       sum(colSums(scen_out_2dose[,grepl(paste0('B_',s_age),names(scen_out_2dose))]))
  #     # adjust
  #     scen_out_booster[,grepl(paste0('E_',s_age),names(scen_out_booster))] <- (1/booster_ratio) * scen_out_booster[,grepl(paste0('E_',s_age),names(scen_out_booster))]
  #   }
  # } else { # reduce booster uptake
  #   cat(paste0("BOOSTER UPTAKE: ADJUST TO ",sens_uptake_booster*100,"%"))
  #   scen_out_booster[,-1] <- scen_out_booster[,-1] * sens_uptake_booster
  # }

}

num_days <- length(endpoint_reported_data :(timing_booster_finished-1))
print(paste('avg booster doses / day:', round(booster_total/ num_days)))
print(paste('avg booster doses / week:', round(booster_total/ (num_days/7))))

scen_out_booster['rna_E_0_9'] <- 0
scen_out_booster[scen_out_booster$date > endpoint_reported_data,'rna_E_10_19'] <- 0

# dim(scen_out_booster)
# plot(scen_out_booster$date,
#      cumsum(scen_out_booster[,15]))
# 
# # WARNING: some scheduled doses did not take place... and have to be "rescheduled"
# i_col <- 14
# flag_projection <- !scen_out_booster$date %in% reported_booster$date
# index_latest    <- min(which(flag_projection))-(1:7)
# index_jan2021   <- scen_out_booster$date %in% seq(as.Date('2021-01-01'),as.Date('2021-02-01'),1)
# # fix for 90-99
# scen_out_booster[flag_projection,ncol(scen_out_booster)] <- scen_out_booster[flag_projection,ncol(scen_out_booster)] / 2
# for(i_col in 2:ncol(scen_out_booster)){
#   doses_removed <- sum(scen_out_booster_db[,i_col]) - sum(scen_out_booster[,i_col])
#   doses_removed
#   if(doses_removed > 0){
#     scen_out_booster[flag_projection,i_col] <- scen_out_booster[flag_projection,i_col] +  mean(scen_out_booster[index_latest,i_col])
#     #scen_out_booster[index_jan2021,i_col] <- scen_out_booster[index_jan2021,i_col] +  (doses_removed / sum(index_jan2021))
#   }
# }

# # WARNING
# # In practice, some booster doses have been administrated before the specified timing
# # in this procedure, which may cause that a booster dose is included twice as reported and as scheduled.
# # Therefore, limit the booster doses to the cumulative number of full doses over time.
# i_col <- 20
# for(i_col in 2:ncol(scen_out_booster)){
#   target_booster <- scen_out_latest[,i_col] * sens_uptake_booster
#   d_overflow <- which(cumsum(scen_out_booster[1:nrow(scen_out_latest),i_col]) > cumsum(target_booster))
#   if(length(d_overflow)>0){
#     print(i_col)
#     #scen_out_booster[d_overflow,i_col] <- apply(cbind(scen_out_latest[d_overflow,i_col],scen_out_booster[d_overflow,i_col]),1,min)
#     scen_out_booster[d_overflow,i_col] <- 0
#     scen_out_booster[min(d_overflow),i_col] <- sum(target_booster) - sum(scen_out_booster[,i_col])
#   }
# }


# drop 'adeno' columns
scen_out_booster <- scen_out_booster[!grepl('adeno_',names(scen_out_booster))]

# merge
scen_out_booster <- left_join(scen_out_2dose,
                            scen_out_booster,
                            by="date")
dim(scen_out_booster)
scen_out_booster[is.na(scen_out_booster)] <- 0

write.table(scen_out_booster,file=file.path(output_folder,paste0('vaccine_uptake_booster_vSCEN',scen_tag,'.csv')),sep=',',row.names=F)

# UPTAKE 2ndBOOSTER DOSE ----
###################### #

# define full uptake
scen_out_latest <- scen_out_booster[,grepl('_E_',names(scen_out_booster)) | grepl('date',names(scen_out_booster))]
names(scen_out_latest)

# schedule booster doses
#scen_out_2ndbooster <- uptake2protection(scen_out_latest,delay_booster,delay_booster)
#names(scen_out_booster) <- gsub('_B_','_E_',names(scen_out_booster))
#tail(scen_out_booster)

# save copy of original planned booster uptake
#scen_out_booster_db <- scen_out_booster

# # untill specified date: reported booster uptake
# set to 1/X time steps
scen_out_2ndbooster <- approx_uptake(db_uptake_2ndbooster,time_step_size,dose_col_2ndbooster,timing_booster_finished) #timing_booster_finished
scen_out_2ndbooster[-1] <- 0
reported_2ndbooster <- approx_uptake(db_uptake_2ndbooster,time_step_size,dose_col_2ndbooster) #timing_booster_finished
scen_out_2ndbooster[scen_out_2ndbooster$date %in% reported_2ndbooster$date,-1] <- reported_2ndbooster[,-1]
scen_out_2ndbooster <- scen_out_2ndbooster[scen_out_2ndbooster$date <= max(scen_out_latest$date),]


# add booster doses until specified data
i_col <- 12
sel_2ndbooster_days <- scen_out_2ndbooster$date >= timing_booster_start & scen_out_2ndbooster$date < (timing_booster_finished-1)
num_time_steps2 <- sum(sel_2ndbooster_days)
if(num_time_steps2 > 0){ # implement future booster campaign
  age_opt     <- get_age_groups()
  s_age <- age_opt[7]
  for(s_age in age_opt[7:10]){ #age_opt[7:10] age_opt[3:10]
    booster_total <- sens_uptake_2ndbooster*sum(scen_out_booster[,grepl(paste0('E_',s_age),names(scen_out_booster))]) - sum(scen_out_2ndbooster[,grepl(paste0('F_',s_age),names(scen_out_2ndbooster))])
    booster_total <- max(0,booster_total)
    scen_out_2ndbooster[sel_2ndbooster_days,grepl(paste0('F_',s_age),names(scen_out_2ndbooster))] <- booster_total/num_time_steps2
    print(paste('avg 2ndbooster doses / day:',round(sens_uptake_2ndbooster*booster_total/ num_time_steps2),' for age ',s_age))
    
      }
  
} 

num_days <- length(endpoint_reported_data :(timing_booster_finished-1))
#print(paste('avg 2ndbooster doses / day:', round(booster2nd_total/ num_days)))


scen_out_2ndbooster['rna_F_0_9'] <- 0
scen_out_2ndbooster[scen_out_2ndbooster$date > endpoint_reported_data,'rna_F_10_19'] <- 0

# dim(scen_out_booster)
# plot(scen_out_booster$date,
#      cumsum(scen_out_booster[,15]))
# 
# # WARNING: some scheduled doses did not take place... and have to be "rescheduled"
# i_col <- 14
# flag_projection <- !scen_out_booster$date %in% reported_booster$date
# index_latest    <- min(which(flag_projection))-(1:7)
# index_jan2021   <- scen_out_booster$date %in% seq(as.Date('2021-01-01'),as.Date('2021-02-01'),1)
# # fix for 90-99
# scen_out_booster[flag_projection,ncol(scen_out_booster)] <- scen_out_booster[flag_projection,ncol(scen_out_booster)] / 2
# for(i_col in 2:ncol(scen_out_booster)){
#   doses_removed <- sum(scen_out_booster_db[,i_col]) - sum(scen_out_booster[,i_col])
#   doses_removed
#   if(doses_removed > 0){
#     scen_out_booster[flag_projection,i_col] <- scen_out_booster[flag_projection,i_col] +  mean(scen_out_booster[index_latest,i_col])
#     #scen_out_booster[index_jan2021,i_col] <- scen_out_booster[index_jan2021,i_col] +  (doses_removed / sum(index_jan2021))
#   }
# }

# # WARNING
# # In practice, some booster doses have been administrated before the specified timing
# # in this procedure, which may cause that a booster dose is included twice as reported and as scheduled.
# # Therefore, limit the booster doses to the cumulative number of full doses over time.
# i_col <- 20
# for(i_col in 2:ncol(scen_out_booster)){
#   target_booster <- scen_out_latest[,i_col] * sens_uptake_booster
#   d_overflow <- which(cumsum(scen_out_booster[1:nrow(scen_out_latest),i_col]) > cumsum(target_booster))
#   if(length(d_overflow)>0){
#     print(i_col)
#     #scen_out_booster[d_overflow,i_col] <- apply(cbind(scen_out_latest[d_overflow,i_col],scen_out_booster[d_overflow,i_col]),1,min)
#     scen_out_booster[d_overflow,i_col] <- 0
#     scen_out_booster[min(d_overflow),i_col] <- sum(target_booster) - sum(scen_out_booster[,i_col])
#   }
# }


# drop 'adeno' columns
scen_out_2ndbooster <- scen_out_2ndbooster[!grepl('adeno_',names(scen_out_2ndbooster))]

# merge
scen_out_2ndbooster <- left_join(scen_out_booster,
                            scen_out_2ndbooster,
                            by="date")
dim(scen_out_2ndbooster)
scen_out_2ndbooster[is.na(scen_out_2ndbooster)] <- 0


write.table(scen_out_2ndbooster,file=file.path(output_folder,paste0('vaccine_uptake_2ndbooster_vSCEN',scen_tag,'.csv')),sep=',',row.names=F)

# UPTAKE EXTRABOOSTER DOSE (dedicated bivariant booster) ----
###################### #


# schedule booster doses
#scen_out_2ndbooster <- uptake2protection(scen_out_latest,delay_booster,delay_booster)
#names(scen_out_booster) <- gsub('_B_','_E_',names(scen_out_booster))
#tail(scen_out_booster)

# save copy of original planned booster uptake
#scen_out_booster_db <- scen_out_booster

# # untill specified date: reported booster uptake
# set to 1/X time steps
scen_out_extrabooster <- approx_uptake(db_uptake_extrabooster,time_step_size,dose_col_extrabooster,timing_booster_finished) #timing_booster_finished
scen_out_extrabooster[-1] <- 0
reported_extrabooster <- approx_uptake(db_uptake_extrabooster,time_step_size,dose_col_extrabooster) #timing_booster_finished
scen_out_extrabooster[scen_out_extrabooster$date %in% reported_extrabooster$date,-1] <- reported_extrabooster[,-1]
scen_out_extrabooster <- scen_out_extrabooster[scen_out_extrabooster$date <= max(scen_out_latest$date),]




num_days <- length(endpoint_reported_data :(timing_booster_finished-1))
#print(paste('avg 2ndbooster doses / day:', round(booster2nd_total/ num_days)))



#(bi-)annual vaccine campaign ---- 
# new_timing_booster_start <- timing_booster_start
# 
# scen_out_extrabooster <- approx_uptake(db_uptake_extrabooster,time_step_size,dose_col_extrabooster,endpoint_uptake) #timing_booster_finished
# scen_out_extrabooster[-1] <- 0
# scen_out_resset <- scen_out_extrabooster
# 
# num_campaign <- floor((endpoint_uptake-timing_booster_start)/(rec_campaign*30.5))
# print(paste("nombre campagnes :",num_campaign))
# for(year in 0:num_campaign){
#   new_timing_booster_start <- timing_booster_start %m+% months(rec_campaign*year) 
#   
# # add booster doses until specified data
# i_col <- 12
# #first campaign 70+
# sel_extrabooster_days <- scen_out_extrabooster$date >= new_timing_booster_start & scen_out_extrabooster$date < (new_timing_booster_start+30)
# num_time_steps3 <- sum(sel_extrabooster_days)
# if(num_time_steps3 > 0){ # implement future booster campaign
#   age_opt     <- get_age_groups()
#   s_age <- age_opt[7]
#   for(s_age in age_opt[8:10]){ 
#     booster_total <- sens_uptake_extrabooster_old*pop_data_be[which(age_opt==s_age)]
#     booster_total <- max(0,booster_total)
#     scen_out_extrabooster[sel_extrabooster_days,grepl(paste0('X_',s_age),names(scen_out_extrabooster))] <- booster_total/num_time_steps3
#     print(paste('***********avg extrabooster doses / day:',round(sens_uptake_extrabooster_old*booster_total/ num_time_steps3),' for age ',s_age))
#   }
# } 
# #second campaign 60+
# sel_extrabooster_days <- scen_out_extrabooster$date >= (new_timing_booster_start+30) & scen_out_extrabooster$date < (new_timing_booster_start+60)
# num_time_steps3 <- sum(sel_extrabooster_days)
# if(num_time_steps3 > 0){ # implement future booster campaign
#   age_opt     <- get_age_groups()
#   s_age <- age_opt[7]
#   for(s_age in age_opt[7:7]){ 
#     booster_total <- sens_uptake_extrabooster_old*pop_data_be[which(age_opt==s_age)]
#     booster_total <- max(0,booster_total)
#     scen_out_extrabooster[sel_extrabooster_days,grepl(paste0('X_',s_age),names(scen_out_extrabooster))] <- booster_total/num_time_steps3
#     print(paste('***********avg extrabooster doses / day:',round(sens_uptake_extrabooster_old*booster_total/ num_time_steps3),' for age ',s_age))
#   }
# } 
# #thids campaign 20-59
# sel_extrabooster_days <- scen_out_extrabooster$date >= (new_timing_booster_start+60) & scen_out_extrabooster$date < (new_timing_booster_start+90)
# num_time_steps3 <- sum(sel_extrabooster_days)
# if(num_time_steps3 > 0){ # implement future booster campaign
#   age_opt     <- get_age_groups()
#   s_age <- age_opt[7]
#   for(s_age in age_opt[3:6]){
#     booster_total <- sens_uptake_extrabooster_young*pop_data_be[which(age_opt==s_age)]
#     booster_total <- max(0,booster_total)
#     scen_out_extrabooster[sel_extrabooster_days,grepl(paste0('X_',s_age),names(scen_out_extrabooster))] <- booster_total/num_time_steps3
#     print(paste('***********avg extrabooster doses / day:',round(sens_uptake_extrabooster_young*booster_total/ num_time_steps3),' for age ',s_age))
#   }
# }
# scen_out_resset[scen_out_extrabooster$date >= new_timing_booster_start,] <- scen_out_extrabooster[scen_out_extrabooster$date >= new_timing_booster_start,]
# for(i_age in 1:10){
#   scen_out_resset[scen_out_resset$date == (new_timing_booster_start+90),paste0('rna_X_',age_opt[i_age])]  <- -sum(scen_out_resset[scen_out_resset$date < (new_timing_booster_start+90),paste0('rna_X_',age_opt[i_age])])
# }
# }
# 
# num_days <- length(endpoint_reported_data :(timing_booster_finished-1))
# #print(paste('avg 2ndbooster doses / day:', round(booster2nd_total/ num_days)))
# 
# 
# scen_out_extrabooster['rna_X_0_9'] <- 0


# drop 'adeno' columns
scen_out_extrabooster <- scen_out_extrabooster[!grepl('adeno_',names(scen_out_extrabooster))]
#scen_out_resset <- scen_out_resset[!grepl('adeno_',names(scen_out_resset))]

scen_out_extrabooster <- left_join(scen_out_2ndbooster,
                            scen_out_extrabooster,
                            by="date")

dim(scen_out_extrabooster)
scen_out_extrabooster[is.na(scen_out_extrabooster)] <- 0


write.table(scen_out_extrabooster,file=file.path(output_folder,paste0('vaccine_uptake_extrabooster_vSCEN',scen_tag,'.csv')),sep=',',row.names=F)


new_timing_booster_start <- timing_booster_start

#scen_out_fallcampaign <- approx_uptake(db_uptake_extrabooster,time_step_size,dose_col_extrabooster,endpoint_uptake) #timing_booster_finished
scen_out_fallcampaign <- scen_out_extrabooster
scen_out_resset <- scen_out_fallcampaign

# add booster doses until specified data
i_col <- 8
#campaign 70+
sel_fallcampaign_days <- scen_out_fallcampaign$date >= new_timing_booster_start & scen_out_fallcampaign$date < endpoint_uptake
num_time_steps3 <- sum(sel_fallcampaign_days)
print(num_time_steps3)
if(num_time_steps3 > 0){ # implement future booster campaign
  age_opt     <- get_age_groups()
  s_age <- age_opt[7]
  for(s_age in age_opt[8:10]){
    booster_total <- sens_flu_pop[i_col]*pop_data_be[which(age_opt==s_age)]
    booster_total <- max(0,booster_total)
    if (bool_bivariant_fall2023){
    scen_out_fallcampaign[sel_fallcampaign_days,grepl(paste0('X_',s_age),names(scen_out_fallcampaign))] <- booster_total/num_time_steps3
    }
    else{
      scen_out_fallcampaign[sel_fallcampaign_days,grepl(paste0('F_',s_age),names(scen_out_fallcampaign))] <- booster_total/num_time_steps3
    }
    print(paste('***********avg extrabooster doses / day:',round(booster_total/ num_time_steps3)/time_step_size,' for age ',s_age))
 i_col <- i_col+1
     }
}
#campaign 65+
sel_fallcampaign_days  <- scen_out_fallcampaign$date >= (new_timing_booster_start) & scen_out_extrabooster$date < endpoint_uptake
num_time_steps3 <- sum(sel_fallcampaign_days)
if(num_time_steps3 > 0){ # implement future booster campaign
  age_opt     <- get_age_groups()
  s_age <- age_opt[7]
  for(s_age in age_opt[7:7]){
    booster_total <- sens_flu_pop[7]*pop_data_be[which(age_opt==s_age)]
    booster_total <- max(0,booster_total)
    if (bool_bivariant_fall2023){
      scen_out_fallcampaign[sel_fallcampaign_days,grepl(paste0('X_',s_age),names(scen_out_fallcampaign))] <- booster_total/num_time_steps3
    }
    else{
      scen_out_fallcampaign[sel_fallcampaign_days,grepl(paste0('F_',s_age),names(scen_out_fallcampaign))] <- booster_total/num_time_steps3
    }
    print(paste('***********avg extrabooster doses / day:',round(booster_total/ num_time_steps3)/time_step_size,' for age ',s_age))
  }
}
scen_out_fallcampaign <- scen_out_fallcampaign[scen_out_fallcampaign$date < timing_booster_finished+1,]
# drop 'adeno' columns
#scen_out_fallcampaign <- scen_out_fallcampaign[!grepl('adeno_',names(scen_out_fallcampaign))]
#scen_out_resset <- scen_out_resset[!grepl('adeno_',names(scen_out_resset))]

# scen_out_fallcampaign <- left_join(scen_out_2ndbooster,
#                                    scen_out_fallcampaign,
#                                    by="date")

dim(scen_out_fallcampaign)
scen_out_fallcampaign[is.na(scen_out_fallcampaign)] <- 0

scen_out_resset[scen_out_fallcampaign$date >= new_timing_booster_start,] <- scen_out_fallcampaign[scen_out_fallcampaign$date >= new_timing_booster_start,]
for(i_age in 1:10){
  if(bool_bivariant_fall2023){
  scen_out_resset[scen_out_resset$date == (new_timing_booster_start-30),paste0('rna_X_',age_opt[i_age])]  <- -sum(scen_out_resset[scen_out_resset$date < (new_timing_booster_start-30),paste0('rna_X_',age_opt[i_age])])
  }
  else{
    scen_out_resset[scen_out_resset$date == (new_timing_booster_start-30),paste0('rna_F_',age_opt[i_age])]  <- -sum(scen_out_resset[scen_out_resset$date < (new_timing_booster_start-30),paste0('rna_F_',age_opt[i_age])])
  }
}

write.table(scen_out_fallcampaign,file=file.path(output_folder,paste0('vaccine_uptake_extrabooster_vSCEN',scen_tag,'.csv')),sep=',',row.names=F)


# explore uptake by day/week
plot(scen_out_booster$date,rowSums(scen_out_booster[,22:31]/time_step_size),main='1st doses (/day)')
plot(scen_out_booster$date,rowSums(scen_out_booster[,22:31]/time_step_size*7),main='1st doses (/week)',las=2)
plot(scen_out_booster$date,rowSums(scen_out_booster[,42:51]/time_step_size),main='booster doses (/day)')
plot(scen_out_booster$date,rowSums(scen_out_booster[,42:51])/time_step_size*7,main='booster doses (/week)',las=2)
plot(scen_out_fallcampaign$date,rowSums(scen_out_fallcampaign[,52:61])/time_step_size*7,main='booster doses (/week)',las=2)
plot(scen_out_fallcampaign$date,rowSums(scen_out_fallcampaign[,62:71])/time_step_size*7,main='booster doses (/week)',las=2)

#EXPLORE BY AGE AND TYPE ----
############################ #

# set final uptake scheme
scen_out_final <- scen_out_fallcampaign
#drop end : to remove if scenarios
#scen_out_final <- scen_out_final[scen_out_final$date < endpoint_reported_data + 60,]

age_opt     <- get_age_groups()
ref_age_col <- paste0('vaccine_A_ages',age_opt)
y_lim       <- c(0,1.5e6)
x_ticks     <- pretty(range(scen_out_final$date),10)

perc_breaks    <- pretty(0:100,7)
y_ticks_perc   <- paste0(perc_breaks,'%')

# open pdf stream
pdf(file=file.path(output_folder,paste0('vaccine_uptake_type_vSCEN',scen_tag,'.pdf')),9,12)
#png(file=file.path(output_folder,paste0('vaccine_uptake_type_vSCEN',scen_tag,'.png')),600,900)



# sum rna and adeno
scen_out_dose1_aggr <- scen_out_final[2:11] + scen_out_final[22:31]
scen_out_dose2_aggr <- scen_out_final[12:21] + scen_out_final[32:41]
scen_out_dose3_aggr <- scen_out_final[42:51]
scen_out_dose4_aggr <- scen_out_final[52:61]
scen_out_dose5_aggr <- scen_out_final[62:71]

age_opt     <- get_age_groups()
ref_age_col <- paste0('vaccine_A_ages',age_opt)
y_lim       <- c(0,1.5e6)
x_ticks     <- pretty(range(scen_out_final$date),10)

perc_breaks    <- pretty(0:100,7)
y_ticks_perc   <- paste0(perc_breaks,'%')
par(mar=c(4,4,2,4),mfrow=c(4,3))

ref_age_col_full  <- gsub('_A_','_BC_',ref_age_col)

# make selection of reference data to improve constast
ref_selection <- seq(1,nrow(ref_data_be),15)

i_age <- 3
for(i_age in 10:1){
  
  y_ticks      <- perc_breaks * pop_data_be[i_age] / 100
  y_ticks_k    <- paste0(round(y_ticks/1e3),'k')
  y_ticks_k[1] <- 0
  y_ticks_k
  
  plot(ref_data_be$date[ref_selection]+1,
       t(cumsum(ref_data_be[ref_age_col[i_age]]))[ref_selection],
       ylab='Uptake',
       xlab='',
       xlim= range(scen_out_final$date),
       ylim = range(y_ticks),
       xaxt='n',
       yaxt='n',
       #cex=0.5,
       main=paste0(gsub('_','-',age_opt[i_age]),'y'))
  
  points(ref_data_be$date[ref_selection]+1,
         t(cumsum(ref_data_be[ref_age_col_full[i_age]]))[ref_selection],
         col=1,
         #cex=0.5,
         pch=3)
  
  lines(scen_out_final$date+1,
        cumsum(scen_out_dose1_aggr[,i_age]),
        col=1,
        lty=1,
        lwd=2)
  
  lines(scen_out_final$date+1,
        cumsum(scen_out_dose2_aggr[,i_age]),
        col=1,
        lty=3,
        lwd=2)
  
  # mRNA
  points(ref_data_be$date[ref_selection]+1,
        t(cumsum(ref_data_be[gsub('_A_','_mRNA_A_',ref_age_col[i_age])]))[ref_selection],
        #cex=0.5,
        col=4)
  points(ref_data_be$date[ref_selection]+1,
         t(cumsum(ref_data_be[gsub('_A_','_mRNA_B_',ref_age_col[i_age])]))[ref_selection],
         #cex=0.5,
         col=4,
         pch=3)
  
  lines(scen_out_final$date+1,
        cumsum(scen_out_final[,paste0('rna_A_',age_opt[i_age])]),
        col=4,
        lty=1,
        lwd=2)
  lines(scen_out_final$date+1,
        cumsum(scen_out_final[,paste0('rna_B_',age_opt[i_age])]),
        col=4,
        lty=3,
        lwd=2)
  
  ## AZ
  points(ref_data_be$date[ref_selection]+1,
         t(cumsum(ref_data_be[gsub('_A_','_AZ_A_',ref_age_col[i_age])]+
                    ref_data_be[gsub('_A_','_JnJ_A_',ref_age_col[i_age])]))[ref_selection],
         #cex=0.5,
         col=5)
  points(ref_data_be$date[ref_selection]+1,
         # t(cumsum(ref_data_be[gsub('_A_','_AZ_B_',ref_age_col[i_age])])),
         t(cumsum(ref_data_be[gsub('_A_','_AZ_B_',ref_age_col[i_age])]+
                    ref_data_be[gsub('_A_','_JnJ_A_',ref_age_col[i_age])]))[ref_selection],
         #cex=0.5,
         col=5,
         pch=3)
  
  lines(scen_out_final$date+1,
        cumsum(scen_out_final[,paste0('adeno_A_',age_opt[i_age])]),
        col=5,
        lty=1,
        lwd=2)
  lines(scen_out_final$date+1,
        cumsum(scen_out_final[,paste0('adeno_B_',age_opt[i_age])]),
        col=5,
        lty=3,
        lwd=2)
  
  # booster
  points(ref_data_be$date[ref_selection]+1,
         t(cumsum(ref_data_be[gsub('_A_','_mRNA_E_',ref_age_col[i_age])]))[ref_selection],
         #cex=0.5,
         pch=6,
         col=7)
  
  lines(scen_out_final$date+1,
        cumsum(scen_out_final[,paste0('rna_E_',age_opt[i_age])]),
        col=7,
        lty=1,
        lwd=2)
  
  
  
  # 2ndbooster
  points(ref_data_be$date[ref_selection]+1,
         t(cumsum(ref_data_be[gsub('_A_','_mRNA_F_',ref_age_col[i_age])]))[ref_selection],
         #cex=0.5,
         pch=9,
         col=2)
  
  # lines(scen_out_final$date+1,
  #       cumsum(scen_out_final[,paste0('rna_F_',age_opt[i_age])]),
  #       col=2,
  #       lty=1,
  #       lwd=2)
  
  axis(1,x_ticks,format(x_ticks,'%m/%Y'),cex.axis=0.8)
  axis(2,y_ticks,y_ticks_perc,las=2)
  axis(4,y_ticks,y_ticks_k,las=2)
  abline(v=x_ticks,lty=3,col='grey')
  abline(h=y_ticks,lty=3,col='grey')
  
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
  
  
  # extrabooster
  
  points(ref_data_be$date[ref_selection]+1,
         t(cumsum(ref_data_be[gsub('_A_','_mRNA_X_',ref_age_col[i_age])]))[ref_selection],
         #cex=0.5,
         pch=2,
         col=6)
  
  if (bool_bivariant_fall2023){
  lines(scen_out_final$date+1,
          cumsum(scen_out_final[,paste0('rna_F_',age_opt[i_age])]),
          col=2,
          lty=1,
          lwd=2)
  lines(scen_out_resset$date+1,
        cumsum(scen_out_resset[,paste0('rna_X_',age_opt[i_age])]),
        col=6,
        lty=1,
        lwd=2)
  }
  else{
    lines(scen_out_resset$date+1,
          cumsum(scen_out_resset[,paste0('rna_F_',age_opt[i_age])]),
          col=2,
          lty=1,
          lwd=2)
    lines(scen_out_final$date+1,
          cumsum(scen_out_final[,paste0('rna_X_',age_opt[i_age])]),
          col=6,
          lty=1,
          lwd=2)
  }
  
  # lines(scen_out_final$date+1,
  #       cumsum(scen_out_final[,paste0('rna_X2_',age_opt[i_age])]),
  #       col=7,
  #       lty=1,
  #       lwd=2)
  
 # lines(scen_out_final$date+1,
      #  cumsum(scen_out_resset[,paste0('rna_X_',age_opt[i_age])]),
      #  col=6,
       # lty=1,
      # lwd=2)
  
  #axis(1,x_ticks,format(x_ticks,'%Y'),cex.axis=0.8)
  axis(2,y_ticks,y_ticks_perc,las=2)
  axis(4,y_ticks,y_ticks_k,las=2)
  abline(v=x_ticks,lty=3,col='grey')
  abline(h=y_ticks,lty=3,col='grey')
  
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
  
  
  
  # legend('topleft',
  #        c('Uptake: total',
  #          'Uptake: mRNA-based',
  #          'Uptake: adeno-based',
  #          'Uptake: booster',
  #          '',
  #          'Reported: 1st dose',
  #          'Reported: full scheme',
  #          'Model: 1st dose',
  #          'Model: full scheme'),
  #        col=c(1,4,5,7,NA,8,8,8,8),
  #        lty=c(1,1,1,1,NA,NA,NA,1,3),
  #        pch=c(NA,NA,NA,NA,NA,1,3,NA,NA),
  #        lwd=2,
  #        cex=0.8)

}

plot(1,bty="n",yaxt="n",xaxt="n",xlab='',ylab='',col=0)
legend('topleft',
       c('Uptake: total',
         'Uptake: mRNA-based',
         'Uptake: adeno-based',
         'Uptake: booster',
         'Uptake: 2ndbooster',
         'Uptake: bivalent booster',
         '',
         'Reported: 1st dose',
         'Reported: 2nd doses',
         'Reported: 1st booster',
         'Reported: 2nd booster',
         'Reported: bivalent booster',
         'Model: 1st dose',
         'Model: 2nd dose'),
       col=c(1,4,5,7,2,6,NA,8,8,8,8,8,8),
       lty=c(1,1,1,1,1,1,NA,NA,NA,NA,NA,1,3),
       pch=c(NA,NA,NA,NA,NA,NA,NA,1,3,6,9,2,NA,NA),
       lwd=2,
       xpd=TRUE,
       cex=1.1)


# summary statistics
names(scen_out_final)
sel_jan2022 <- scen_out_final$date >= as.Date('2022-01-01') & scen_out_final$date < as.Date('2022-02-1')
col_A <- grepl('_A_',names(scen_out_final))
col_E <- grepl('_E_',names(scen_out_final))
col_F <- grepl('_F_',names(scen_out_final))
week_avg_first   <- sum(scen_out_final[sel_jan2022,col_A]) / (sum(sel_jan2022)*time_step_size) *7
week_avg_booster <- sum(scen_out_final[sel_jan2022,col_E]) / (sum(sel_jan2022)*time_step_size) *7
week_avg_2ndbooster <- sum(scen_out_final[sel_jan2022,col_F]) / (sum(sel_jan2022)*time_step_size) *7


plot(1,bty="n",yaxt="n",xaxt="n",xlab='',ylab='',col=0)
legend('topright',
       c(#paste('Period:',paste(range(scen_out_final$date),collapse=' to ')),
         #paste0('Uptake 1st dose: ',format(week_avg_first,big.mark = ' ',scientific=F,digits=0),'/week'),
         #paste0('Uptake booster dose: ',format(week_avg_booster,big.mark = ' ',scientific=F,digits=0),'/week'),
         #'',
         #'Fall 2023 campaign',
         paste0('Start: ',timing_booster_start),
         paste0('End: ',endpoint_uptake)#,
         #paste0('Repeated each ',rec_campaign, ' months'),
         #paste0('Uptake 60+: ',sens_uptake_extrabooster_old*100,'% total population'),
         #paste0('Uptake 20-59: ',sens_uptake_extrabooster_young*100,'% total population')
       ),
       title='Fall 2023 campaign',
       col=0,
       xpd=TRUE,
       cex=1.1)


dev.off()

timing_booster_finished

# sel_region <- 'brussels'
# filename_vacc <- dir('data/uptake/uptake_region',pattern = 'vaccine_uptake_2doses_Vx_vSCENoct01',full.names = T,recursive = T)
# filename_vacc <- dir('data/uptake/uptake_region',pattern = 'vaccine_uptake_2doses_Vx_vSCENaug',full.names = T,recursive = T)
plot_historic_uptake <- function(sel_region,filename_vacc){
  
  # load vaccine uptake
  vacc_schedule      <- read.table(filename_vacc[grepl(sel_region,filename_vacc)], sep=',',header = T)
  
  # compare with current uptake
  be_ref_data <- get_latest_incidence_data(sel_region = sel_region)

  # age tags
  age_tags <- gsub('adeno_A_','',names(vacc_schedule)[2:11])
    
  # plot
  plot(be_ref_data$date,
       cumsum(be_ref_data$vaccine_uptake_first),
       xlim = range(as.Date(vacc_schedule$date)),
       col = 4,
       xaxt='n') 
  lines(as.Date(vacc_schedule$date),
        cumsum(rowSums(vacc_schedule[,grepl('_A_',names(vacc_schedule))])),
        lwd=2)
  add_date_axis(as.Date(vacc_schedule$date))  
  
  i_tag <- age_tags[5]
  par(mfrow=c(2,2))
  for(i_tag in age_tags){
    
    plot(be_ref_data$date,
         unlist(cumsum(be_ref_data[paste0('vaccine_A_ages',i_tag)])),
         xlim = range(as.Date(vacc_schedule$date)),
         col = 4,
         main= i_tag,
         xaxt='n') 
    lines(as.Date(vacc_schedule$date),
          cumsum(rowSums(vacc_schedule[,grepl(paste0('_A_',i_tag),names(vacc_schedule))])),
          lwd=2)
    add_date_axis(as.Date(vacc_schedule$date))
    
  }
  
  
}


