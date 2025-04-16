########################################################################### #
# This file is part of the Stochastic Compartmental Model for SARS-COV-2 
# transmission in Belgium, conceived by members of SIMID group during the 
# COVID19 pandemic.
#
# This file is used to generate vaccine uptake files for the Technical Notes.
#
# Copyright 2021, SIMID, University of Antwerp & Hasselt University                                        
########################################################################## # 

# clear workspace
rm(list=ls())

library(scales)
library(data.table)
library(openxlsx)
library(foreach)

source('R/lib_stochastic_model.R')


# SETTINGS          ----
###################### # 

# set scenario tag
scen_tag <- 'dec7_child'

# output folder
output_folder <- paste0('output/uptake_',scen_tag)

# planned uptake file
#uptake_file <- 'data/OneDrive_1_07-04-2021/20200407_Input Modellers.xlsx'
uptake_file <- 'data/OneDrive_1_28-10-2021_dummy/2021_dummy_Input Modellers.xlsx'
uptake_file <- 'data/OneDrive_1_28-10-2021_dummy/2021_dummy_Input Modellers_child.xlsx'

# endpoint to select reported data
endpoint_reported_data <- Sys.Date()
endpoint_reported_data <- as.Date('2021-12-06')

# endpoint uptake
endpoint_uptake <- as.Date('2022-07-01')

sel_region <- "belgium" # belgium, flanders, wallonia, brussels
sens_uptake_cap_adult  <- 0.95 #0.90
sens_uptake_cap_child  <- 0.4 

# sensitivity on uptake
#sens_uptake_cap  <- 0.85 # default: 0.8 
# sens_uptake_rate <- 1.0
#sens_uptake_adolescent <- 8/10
sens_uptake_adolescent <- 6/10 * 0.8 + 2/10 * sens_uptake_cap_child + 2/10*sens_uptake_cap_adult# 1
sens_uptake_child      <- 5/10 * sens_uptake_cap_child #  4/10

#sens_uptake_adolescent <- 6/10 * sens_uptake_cap_child + 2/10*sens_uptake_cap_adult# 1
#sens_uptake_child      <- 0/10 * sens_uptake_cap_child #  4/10

bool_incr_child_uptake <- FALSE

bool_threshold <- FALSE

#bool_7d_avg_only <- TRUE # base future uptake only on last 7days
bool_7d_avg_only <-  FALSE # base future uptake on last 7days


# target: based on 1st dose uptake 35-44y on Dec 5th, 2021
# BE: 86%
# BXL: 70%
# FL: 91%
# WL: 80%
# +5%

# # region
# sel_region <- "flanders" # belgium, flanders, wallonia, brussels
# sens_uptake_cap  <- 0.95
# 
# sel_region <- "wallonia" # belgium, flanders, wallonia, brussels
# sens_uptake_cap  <- 0.85
# 
# sel_region <- "brussels" # belgium, flanders, wallonia, brussels
# sens_uptake_cap  <- 0.75
# 

# add to scen_tag 
#scen_tag         <- paste0(scen_tag,'_uc',round(sens_uptake_cap*100),'_ur',round(sens_uptake_rate*100),'_ad',sens_uptake_adolescent*10)
#scen_tag         <- paste0(scen_tag,'_uc',round(sens_uptake_cap*100),'_ad',sens_uptake_adolescent*10,'_',sel_region)
scen_tag         <- paste0(scen_tag,'_ua',round(sens_uptake_cap_adult*100),'_uc',round(sens_uptake_cap_child*100),'_ad',sens_uptake_adolescent*10,'_child',sens_uptake_child*10,'_',ifelse(bool_incr_child_uptake,'summer_',''),sel_region)

if(bool_threshold) { scen_tag     <- gsub('_ad','t_ad',scen_tag) }

# delay for second dose
delay_rna_dose2 <- 3*7
delay_adeno_dose2 <- 0 #12*7 # only J&J in the last months

# note: since JnJ is incorporated as AZ... it requires an arbitraty delay for  
# the virtual second JnJ dose to account for the time to build up protection 
# of the 1st dose (of e.g. 21 days).
delay_JnJ         <- 21 

# delay booster: single dose, so not applicable
delay_booster     <- 6*30 #round(365/2)
if(delay_booster>12*30) { scen_tag     <- gsub(sel_region,paste0('nBoost_',sel_region),scen_tag) }

# time-step size
time_step_size <- 1/4  # 1 == per day and 1/4 == per 6 hours 

# population data
pop_data_be  <- get_regional_pop(sel_region)

# reference data, and select uptake data
ref_data_be      <- get_observed_incidence_data(enable_data_reuse = T,sel_region = sel_region,
                                                bool_incr_child_uptake = bool_incr_child_uptake)
ref_data_be$date <- as.Date(ref_data_be$date)
db_uptake        <- ref_data_be[,grepl('vaccine_.*_._',names(ref_data_be)) | grepl('date',names(ref_data_be))]
agegroup_opt     <- unique(gsub('.*_*_ages','',names(db_uptake))[-1])
dose_col         <- paste0('A_',agegroup_opt)
dose_col_2nd     <- paste0('B_',agegroup_opt)
dose_col_booster <- paste0('E_',agegroup_opt)

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
db_uptake <- db_uptake[db_uptake$date >= as.Date('2021-01-01'),]
db_uptake <- db_uptake[db_uptake$date <= endpoint_reported_data,]
# print(range(db_uptake$date))

# change column names
names(db_uptake) <- gsub('ages','',names(db_uptake))

# split 1st and 2nd uptake
db_uptake_first   <- db_uptake[,grepl('_A_',names(db_uptake)) | grepl('date',names(db_uptake))] 
db_uptake_second  <- db_uptake[,grepl('_B_',names(db_uptake))| grepl('date',names(db_uptake))] 
db_uptake_booster <- db_uptake[,grepl('_E_',names(db_uptake))| grepl('date',names(db_uptake))] 
names(db_uptake_first)
names(db_uptake_second)
names(db_uptake_booster)
sum(db_uptake_first[,-1])
sum(db_uptake_second[,-1])
sum(db_uptake_booster[,-1])

sum(ref_data_be$vaccine_uptake_first)
sum(ref_data_be$vaccine_uptake_full)
sum(ref_data_be$vaccine_uptake_booster)

sum(db_uptake_first$rna_A_80_89) + sum(db_uptake_first$adeno_A_80_89)
sum(ref_data_be$vaccine_A_ages80_89)

# PREAMBLE

# if output folder does not exit yet, create one
if(!dir.exists(output_folder)){
  dir.create(output_folder,recursive = T)
}

# UPTAKE SCENARIOS  ----
###################### # 

# parse planned uptake file

if(nchar(uptake_file)>0){
  uptake_file_row_start <- c(5,17,29)
  uptake_file_nrow      <- 6
  uptake_scenario       <- 1
  input_uptake_month <- read.xlsx(uptake_file,sheet = "Scenarios",
                                  rows = uptake_file_row_start[uptake_scenario] + c(1:uptake_file_nrow-1),
                                  # cols = c(2,4:16),
                                  cols = c(2,4:18),
                                  detectDates=T)
  
  uptake_first_total    <- colSums(input_uptake_month[1:3,-1],na.rm = T)
  uptake_second_total   <- colSums(input_uptake_month[4:5,-1],na.rm = T)
  
  uptake_dates <- as.Date(names(input_uptake_month)[-1])
  uptake_dates
  
  # add additional point
  uptake_dates        <- c(uptake_dates,max(uptake_dates)+31)
  uptake_first_total  <- c(uptake_first_total,0)
  uptake_second_total <- c(uptake_second_total,0)
  
  # go from 'month' to 'day'
  plan_uptake_month    <- data.frame(date  = uptake_dates,
                                     dose_A_total = uptake_first_total,
                                     dose_B_total = uptake_second_total)
  plan_uptake_day <- approx_uptake(plan_uptake_month,1,c('A_total','B_total'),max_date = endpoint_uptake)
  if(sens_uptake_cap_adult>=0.90){
    head(plan_uptake_day)
    flag_2022 <- plan_uptake_day$date >= as.Date('2022-01-01')
    plan_uptake_day[flag_2022 ,-1] <- plan_uptake_day[flag_2022,-1] * 30
    plan_uptake_day[!flag_2022 ,-1] <- 0
  }
} 
dim(plan_uptake_day)
head(plan_uptake_day)
tail(plan_uptake_day)

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

uptake_matrix_date <- plan_uptake_day$date
uptake_matrix      <- data.frame(matrix(0,
                                        ncol=ncol(db_uptake_first)-1, 
                                        nrow = nrow(plan_uptake_day)))
names(uptake_matrix) <- names(db_uptake_first)[-1]

# set target uptake
# target_uptake_age_fraction <- c(0.01,sens_uptake_adolescent*sens_uptake_cap_adult,rep(sens_uptake_cap_adult,3),max(0.85,sens_uptake_cap_adult),0.9,rep(0.8,3))
target_uptake_age_fraction <- c(sens_uptake_child,sens_uptake_adolescent,rep(1,8)*sens_uptake_cap_adult)
target_uptake_age_count    <- target_uptake_age_fraction * pop_data_be

# untill specified date: reported
uptake_matrix[plan_uptake_day$date %in% db_uptake_first$date,] <- db_uptake_first[,grepl('_A_',names(db_uptake_first))]

# UPTAKE OKT-...  ----
###################### # 
# next (use aggregated uptake mRNA + adeno to define max)
uptake_matrix_proj   <- uptake_matrix[1:10] + uptake_matrix[11:20]
remaining_uptake     <- target_uptake_age_count - colSums(uptake_matrix_proj)
pop_uptake_reported  <- colSums(uptake_matrix_proj) / pop_data_be

pop_uptake_reported[1:2]  <- 0

avg_7d_uptake <- colMeans(db_uptake[(-6:0)+nrow(db_uptake),-1])
avg_7d_uptake <- avg_7d_uptake[1:10] + avg_7d_uptake[21:30]

# reported
i_day <- nrow(db_uptake)+1
plan_uptake_day[i_day,]
plan_uptake_day[i_day,]
#for(i_day in 110:120){
for(i_day in (nrow(db_uptake)+1):nrow(plan_uptake_day)){
#for(i_day in 110:119){
    
  uptake_age             <- c(0,0,rep(0,4),0,1,1,1)
  
  # additional loop to shift doses if one age group has reached the max uptake
  while(sum(uptake_matrix_proj[i_day,]) < plan_uptake_day$dose_A_total[i_day]){
    available_doses        <- plan_uptake_day$dose_A_total[i_day] - sum(uptake_matrix_proj[i_day,])
    remaining_uptake       <- target_uptake_age_count - colSums(uptake_matrix_proj) 
    
    # if reported uptake > target... make sure the maximum is used
    remaining_uptake <- pmax(0,remaining_uptake)
    
    if(sum(remaining_uptake) == 0){
      break
    }
    
    # from June: General population, top down
    if(plan_uptake_day$date[i_day] >= as.Date('2021-05-01') ||
       sum(remaining_uptake * uptake_age) == 0){
      
      target_age <- min(length(remaining_uptake),which(remaining_uptake == 0)-1)
      uptake_age <- c(0,0,rep(0,length(remaining_uptake)-2))
      uptake_age[target_age:length(remaining_uptake)] <- 1
      #print(uptake_age)
    }
  
    # add default uptake ONCE, and account for this in the available doses  
    if(!any(uptake_matrix_proj[i_day,]>0)){
      standard_uptake_adjusted <- (avg_7d_uptake * (remaining_uptake>0))
      available_doses          <- available_doses - sum(standard_uptake_adjusted)
      available_doses          <- pmax(0,available_doses)
    } else{
      standard_uptake_adjusted[] <- 0
    }
  
    if(bool_7d_avg_only){
      plan_uptake_day$dose_A_total[i_day] <- sum(standard_uptake_adjusted)
      available_doses <- 0
    }
    
    remaining_uptake_day   <- remaining_uptake * uptake_age
    planned_uptake         <- available_doses * (remaining_uptake_day / sum(remaining_uptake_day))
    
    planned_uptake         <- planned_uptake + standard_uptake_adjusted
    planned_uptake         <- pmin(remaining_uptake,planned_uptake)  
    
    uptake_matrix_proj[i_day,]  <- uptake_matrix_proj[i_day,] + planned_uptake 
  }
  #print(uptake_matrix_proj[i_day,])
}

# merge reported and projected uptake
flag_projected <- plan_uptake_day$date > max(db_uptake$date)
uptake_matrix[flag_projected,grepl('rna',names(uptake_matrix))] <- uptake_matrix_proj[flag_projected,]

# STRICT THRESHOLD ? ----
if(bool_threshold){
  
  uptake_age_threshold <- pmin(sens_uptake_cap_adult,target_uptake_age_fraction) * pop_data_be
  
  i_age <- 1
  for(i_age in 1:length(uptake_age_threshold)){
    
    col_age <- i_age+c(0,10)
    uptake_age_cum <- cumsum(rowSums(uptake_matrix[col_age]))
    bool_overflow <- uptake_age_cum > uptake_age_threshold[i_age]
    uptake_matrix[bool_overflow,col_age] <- 0
  }
  
}




# REFORMAT UPTAKE   ----
###################### # 

# merge date and uptake
scen_uptake_day    <- data.frame(date = uptake_matrix_date,
                                 uptake_matrix)
dim(scen_uptake_day)
names(scen_uptake_day)
# set to 1/X time steps
scen_out <- approx_uptake(scen_uptake_day,time_step_size,dose_col,max_date=endpoint_reported_data)

# write to file
write.table(scen_out,file=file.path(output_folder,paste0('vaccine_uptake_Vx_vSCEN',scen_tag,'.csv')),sep=',',row.names=F)

# UPTAKE SECOND DOSE ----
###################### #

scen_out_2dose <- uptake2protection(scen_out,delay_rna_dose2,delay_adeno_dose2)
names(scen_out_2dose) <- gsub('_A_','_B_',names(scen_out_2dose))

# # untill specified date: reported
# set to 1/X time steps
reported_second <- approx_uptake(db_uptake_second,time_step_size,dose_col_2nd,max_date=endpoint_reported_data)
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
scen_out_2dose <- merge(scen_out,
                        scen_out_2dose,
                        all = T)
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
reported_booster <- approx_uptake(db_uptake_booster,time_step_size,dose_col_booster,max_date=endpoint_reported_data)
scen_out_booster[scen_out_booster$date %in% reported_booster$date,-1] <- reported_booster[,-1]
scen_out_booster <- scen_out_booster[scen_out_booster$date <= max(scen_out_latest$date),]

# WARNING: some scheduled doses did not take place... and have to be "rescheduled"
i_col <- 14
flag_projection <- !scen_out_booster$date %in% reported_booster$date
index_latest    <- min(which(flag_projection))-(1:7)
# fix for 90-99
scen_out_booster[flag_projection,ncol(scen_out_booster)] <- scen_out_booster[flag_projection,ncol(scen_out_booster)] / 2
for(i_col in 2:ncol(scen_out_booster)){
  doses_removed <- sum(scen_out_booster_db[,i_col]) - sum(scen_out_booster[,i_col])
  doses_removed
  if(doses_removed > 0){
    scen_out_booster[flag_projection,i_col] <- scen_out_booster[flag_projection,i_col] +  mean(scen_out_booster[index_latest,i_col])
  }
}

# WARNING
# In practice, some booster doses have been administrated before the specified timing
# in this procedure, which may cause that a booster dose is included twice as reported and as scheduled.
# Therefore, limit the booster doses to the cumulative number of full doses over time.
i_col <- 20
for(i_col in 2:ncol(scen_out_booster)){
  d_overflow <- which(cumsum(scen_out_booster[1:nrow(scen_out_latest),i_col]) > cumsum(scen_out_latest[,i_col]))
  if(length(d_overflow)>0){
    print(i_col)
    #scen_out_booster[d_overflow,i_col] <- apply(cbind(scen_out_latest[d_overflow,i_col],scen_out_booster[d_overflow,i_col]),1,min)
    scen_out_booster[d_overflow,i_col] <- 0
    scen_out_booster[min(d_overflow),i_col] <- sum(scen_out_latest[,i_col]) - sum(scen_out_booster[,i_col])
  }
}


# drop 'adeno' columns
scen_out_booster <- scen_out_booster[!grepl('adeno_',names(scen_out_booster))]

# merge
scen_out_booster <- merge(scen_out_2dose,
                          scen_out_booster,
                          all = T)
dim(scen_out_booster)
scen_out_booster[is.na(scen_out_booster)] <- 0

write.table(scen_out_booster,file=file.path(output_folder,paste0('vaccine_uptake_booster_vSCEN',scen_tag,'.csv')),sep=',',row.names=F)

#EXPLORE BY AGE AND TYPE ----
############################ #

age_opt     <- get_age_groups()
ref_age_col <- paste0('vaccine_A_ages',age_opt)
y_lim       <- c(0,1.5e6)
x_ticks     <- pretty(range(scen_out$date),10)

perc_breaks    <- pretty(0:100,7)
y_ticks_perc   <- paste0(perc_breaks,'%')

# reset reference uptake
if(bool_incr_child_uptake){
  ref_data_be      <- get_observed_incidence_data(enable_data_reuse = T,sel_region = sel_region,
                                                  bool_incr_child_uptake = FALSE)  
}

# open pdf stream
pdf(file=file.path(output_folder,paste0('vaccine_uptake_type_vSCEN',scen_tag,'.pdf')),9,6)

# set final uptake scheme
scen_out_final <- scen_out_booster

# sum rna and adeno
scen_out_dose1_aggr <- scen_out_final[2:11] + scen_out_final[22:31]
scen_out_dose2_aggr <- scen_out_final[12:21] + scen_out_final[32:41]
scen_out_dose3_aggr <- scen_out_final[42:51]

age_opt     <- get_age_groups()
ref_age_col <- paste0('vaccine_A_ages',age_opt)
y_lim       <- c(0,1.5e6)
x_ticks     <- pretty(range(scen_out_final$date),10)

perc_breaks    <- pretty(0:100,7)
y_ticks_perc   <- paste0(perc_breaks,'%')
par(mar=c(5,5,4,5),mfrow=c(1,1))

ref_age_col_full  <- gsub('_A_','_BC_',ref_age_col)

i_age <- 3
for(i_age in 10:1){
  
  y_ticks      <- perc_breaks * pop_data_be[i_age] / 100
  y_ticks_k    <- paste0(round(y_ticks/1e3),'k')
  y_ticks_k[1] <- 0
  y_ticks_k
  
  plot(ref_data_be$date+1,
       t(cumsum(ref_data_be[ref_age_col[i_age]])),
       ylab='Uptake',
       xlab='',
       xlim= range(scen_out$date),
       ylim = range(y_ticks),
       xaxt='n',
       yaxt='n',
       cex=0.5,
       main=paste0(gsub('_','-',age_opt[i_age]),'y'))
  
  points(ref_data_be$date+1,
         t(cumsum(ref_data_be[ref_age_col_full[i_age]])),
         col=1,
         cex=0.5,
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
  points(ref_data_be$date+1,
        t(cumsum(ref_data_be[gsub('_A_','_mRNA_A_',ref_age_col[i_age])])),
        cex=0.5,
        col=4)
  points(ref_data_be$date+1,
         t(cumsum(ref_data_be[gsub('_A_','_mRNA_B_',ref_age_col[i_age])])),
         cex=0.5,
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
  points(ref_data_be$date+1,
         t(cumsum(ref_data_be[gsub('_A_','_AZ_A_',ref_age_col[i_age])]+
                    ref_data_be[gsub('_A_','_JnJ_A_',ref_age_col[i_age])])),
         cex=0.5,
         col=5)
  points(ref_data_be$date+1,
         # t(cumsum(ref_data_be[gsub('_A_','_AZ_B_',ref_age_col[i_age])])),
         t(cumsum(ref_data_be[gsub('_A_','_AZ_B_',ref_age_col[i_age])]+
                    ref_data_be[gsub('_A_','_JnJ_A_',ref_age_col[i_age])])),
         cex=0.5,
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
  points(ref_data_be$date+1,
         t(cumsum(ref_data_be[gsub('_A_','_mRNA_E_',ref_age_col[i_age])])),
         cex=0.5,
         col=7)
  
  lines(scen_out_final$date+1,
        cumsum(scen_out_final[,paste0('rna_E_',age_opt[i_age])]),
        col=7,
        lty=1,
        lwd=2)
  
  axis(1,x_ticks,format(x_ticks,'%d/%m'),cex.axis=0.8)
  axis(2,y_ticks,y_ticks_perc,las=2)
  axis(4,y_ticks,y_ticks_k,las=2)
  abline(v=x_ticks,lty=3,col='grey')
  abline(h=y_ticks,lty=3,col='grey')
  
  abline_endpoint_reported_data <- ifelse(i_age <= 2 && bool_incr_child_uptake,
                                          as.Date('2021-07-01'),
                                          endpoint_reported_data)
  abline(v=abline_endpoint_reported_data,lwd=2,col=3)  
  text(x = abline_endpoint_reported_data,
       y = max(y_ticks)*0.98,
       'reported',
       pos=2,
       col=3,
       cex=0.7)
  text(x = abline_endpoint_reported_data,
       y = max(y_ticks)*0.98,
       'scenario',
       pos=4,
       col=3,
       cex=0.7)
  
  legend('topleft',
         c('Uptake: total',
           'Uptake: mRNA-based',
           'Uptake: adeno-based',
           'Uptake: booster',
           '',
           'Reported: 1st dose',
           'Reported: full scheme',
           'Model: 1st dose',
           'Model: full scheme'),
         col=c(1,4,5,7,NA,8,8,8,8),
         lty=c(1,1,1,1,NA,NA,NA,1,3),
         pch=c(NA,NA,NA,NA,NA,1,3,NA,NA),
         lwd=2,
         cex=0.8)

}

dev.off()


## ONE FIGURE ---
# open pdf stream
pdf(file=file.path(output_folder,paste0('vaccine_uptake_type_vSCEN_single',scen_tag,'.pdf')),9,6)

i_dose <- 1
for(i_dose in 2:3){

  if(i_dose == 2){
    scen_out_aggr <- scen_out_dose2_aggr
    dose_tag <- '2nd dose'
    
  }  else{
    scen_out_aggr <- scen_out_dose3_aggr
    dose_tag <- 'booster'
  }
  
  y_ticks <- perc_breaks / 100
  y_ticks_lab <- 
  
  # 2nd dose
  # scen_out_dose2_aggr
  plot(ref_data_be$date+1,
       t(cumsum(ref_data_be[ref_age_col[1]]))*0,
       ylab=paste('Uptake',dose_tag),
       xlab='',
       xlim = range(scen_out$date),
       ylim = range(y_ticks),
       xaxt='n',
       yaxt='n',
       col=0,
       cex=0.5,
       main=dose_tag)
  
  axis(1,x_ticks,format(x_ticks,'%d/%m'),cex.axis=0.8)
  axis(2,y_ticks,y_ticks_perc,las=2)
  abline(v=x_ticks,lty=3,col='grey')
  abline(h=y_ticks,lty=3,col='grey')
  
  abline_endpoint_reported_data <- ifelse(bool_incr_child_uptake,
                                          as.Date('2021-07-01'),
                                          endpoint_reported_data)
  abline(v=abline_endpoint_reported_data,lwd=2,col=3)  
  text(x = abline_endpoint_reported_data,
       y = max(y_ticks)*0.98,
       'reported',
       pos=2,
       col=3,
       cex=0.7)
  text(x = abline_endpoint_reported_data,
       y = max(y_ticks)*0.98,
       'scenario',
       pos=4,
       col=3,
       cex=0.7)
  
  legend('topleft',
         rev(paste0(gsub('_','-',age_opt),'y')),
         col=10:1,
         lty=1,
         pch=10:1,
         lwd=1,
         cex=0.8)
  
  i_age <- 3
  for(i_age in 10:1){
    
    lines(scen_out_final$date+1,
          cumsum(scen_out_aggr[,i_age]) / pop_data_be[i_age],
          col=i_age,
          lty=1,
          lwd=2)
    
    sel_points <- pretty(1:nrow(scen_out_final))
    points(scen_out_final$date[sel_points]+1,
          (cumsum(scen_out_aggr[,i_age]) / pop_data_be[i_age])[sel_points],
          col=i_age,
          pch=i_age)
  }
  
  
} # end for-loop
dev.off()

## - average coverage ----

adult_weight <- c(0,2/10,rep(1,8))
coverage_num <- colSums(scen_out_dose2_aggr)*adult_weight
pop_num      <- pop_data_be * adult_weight
print(sum(coverage_num) / sum(pop_num))
coverage_num/pop_num

adult_weight <- c(0,0,rep(1,8))
adult_weight<- 1
coverage_num <- colSums(scen_out_dose2_aggr)*adult_weight
pop_num      <- pop_data_be * adult_weight
print(sum(coverage_num) / sum(pop_num))

# explore uptake per month ----
uptake_day <- scen_out_dose1_aggr
uptake_day$date_month_lab <- format(scen_out_final$date,'%b%Y')
uptake_day$date_month <- format(scen_out_final$date,'%Y-%m')
uptake_day$uptake_sum <- rowSums(scen_out_dose1_aggr)
uptake_month <- aggregate(uptake_sum ~ date_month,data=uptake_day,sum)
#uptake_month <- uptake_month[-nrow(uptake_month),]
bplot <- barplot(uptake_month$uptake_sum,
                 col='grey',
                 ylab='Uptake 1st/single dose (month)',
                 main='Monthly uptake: 1st/single dose',
                 ylim = c(0,3.5e6))
grid(nx=NA,ny=NULL)
axis(1,bplot,format(as.Date(paste0(uptake_month$date_month,'-01')),'%b'))
print(uptake_month$uptake_sum)

print(sum(uptake_month$uptake_sum[12:15]))

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


