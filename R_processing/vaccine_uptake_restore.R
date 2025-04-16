########################################################################### #
# This file is part of the Stochastic Compartmental Model for SARS-COV-2 
# transmission in Belgium, conceived by members of SIMID group during the 
# COVID19 pandemic.
#
# This file is used to generate vaccine uptake files for RESTORE analyses.
#
# Copyright 2021, SIMID, University of Antwerp & Hasselt University                                        
########################################################################## # 

############################################################# #
# Assumptions
# ==>> all vaccines are fully protective after 4w, no waning.
# ==>> so, all uptake is assumed to be adeno-based
# ==>> so, only first dose counts for now.
#
# Output
# ==>> uptake from 1/1/2021, by age, per 1/24day
# ==>> protection from 1/1/2021, by age, per 1/24day
# ==>> columns: date, 10x RNA_A, 10x RNA_B, 10x adeno_A, 10x adeno_B
############################################################# #

# clear workspace
rm(list=ls())

library(scales)
library(data.table)
library(openxlsx)
library(foreach)

source('./R/lib_stochastic_model.R')


# SETTINGS          ----
###################### # 

# set scenario tag
scen_tag <- 'paper30'

# output folder
output_folder <- 'output/uptake_paper30'

# endpoint to select reported data
endpoint_reported_data <- as.Date('2021-06-01')

# vaccine settings
prop_rna               <- 0
prop_adeno             <- 1
delay_protection_rna   <- 21
delay_protection_adeno <- 21 #21

# uptake settings
planned_doses_day <- 60e3  # from April 2021
seq_age_groups    <- 0

# time-step size
time_step_size <- 1/24  # 1: per day and 1/24: per hour 

# population data
pop_data_be  <- get_be_pop_matrix(1)

# reference data, and select uptake data
ref_data_be   <- get_observed_incidence_data()
#ref_data_be   <- read.table('data/covid19_reference_data_20210415.csv',sep=',',header=T)
ref_data_be$date <- as.Date(ref_data_be$date)
db_uptake     <- ref_data_be[,grepl('vaccine_A',names(ref_data_be)) | grepl('date',names(ref_data_be))]
agegroup_opt  <- gsub('vaccine_A_ages','',names(db_uptake))[-1]
num_agegroup  <- length(agegroup_opt)
# dose_col      <- paste0(rep(c('A_','B_'),each=num_agegroup),rep(agegroup_opt,2))
dose_col      <- paste0('A_',agegroup_opt)

# select between Jan 1st and March 31th, 2021
db_uptake <- db_uptake[db_uptake$date >= as.Date('2021-01-01'),]
db_uptake <- db_uptake[db_uptake$date < endpoint_reported_data,]
print(range(db_uptake$date))

# change column names
names(db_uptake) <- gsub('vaccine','adeno',names(db_uptake))
names(db_uptake) <- gsub('ages','',names(db_uptake))


# PREAMBLE

# if output folder does not exit yet, create one
if(!dir.exists(output_folder)){
  dir.create(output_folder,recursive = T)
}

# ammend scenario tag with uptake settings
scen_tag <- paste0(scen_tag,'_u',planned_doses_day/1e3,'k')
scen_tag <- paste0(scen_tag,'_a',paste0(seq_age_groups,collapse=''))

# open pdf stream
pdf(file=file.path(output_folder,paste0('vaccine_uptake_vSCEN',scen_tag,'.pdf')))
print(scen_tag)

# UPTAKE SCENARIOS  ----
###################### # 

# planned uptake 
# 40k doses / day
# 60k doses / day 

uptake_dates       <- as.Date(paste0('2021-',1:9,'-01'))
uptake_num_days    <- c(diff(uptake_dates),0)
uptake_first_total <- as.numeric(uptake_num_days * planned_doses_day)
uptake_first_total

uptake_first_total[1:3] <- c(0,0,sum(db_uptake[1:max(which(db_uptake$date %in% uptake_dates)),-1]))

x_ticks <- pretty(as.Date(uptake_dates),12)
plot(db_uptake$date,
     cumsum(rowSums(db_uptake[,-1])),
     main= 'Latest Forecast based on non-confirmed estimates \n THIS CAN CHANGE ON A WEEKLY BASIS ',
     ylab='Uptake',
     xlab='',
     ylim = c(0,8e6),
     xlim=range(x_ticks),
     xaxt='n')
axis(1,x_ticks,format(x_ticks,"%d/%m"))
lines(uptake_dates[-1],
     cumsum(uptake_first_total)[-length(uptake_first_total)],
     lwd=2,
     col=4,
     type='b')

legend('topleft',
       c('Reported uptake: 1st or single dose',
         'Expected uptake: 1st or single dose'),
       col = c(1,4,4),
       lwd= c(NA,2,2),
       pch=c(1,NA,NA)
)


# UPTAKE APRIL-AUG  ----
###################### # 

# set target uptake
target_uptake_age_fraction <- c(NA,NA,rep(0.8,5),rep(0.9,3))
target_uptake_age_count    <- target_uptake_age_fraction * pop_data_be

# initiate uptake matrix
uptake_matrix_date <- data.frame(date = seq(as.Date('2021-01-01'),as.Date('2021-09-01'),1))
uptake_matrix      <- matrix(0,nrow=nrow(uptake_matrix_date),ncol=10)
dim(uptake_matrix)
colnames(uptake_matrix) <- names(db_uptake)[-1]

# add reported uptake
uptake_matrix[1:nrow(db_uptake),] <- as.matrix(db_uptake[,-1])
dim(uptake_matrix)

# afterwards, assign doses per day, by the given age sequence
target_uptake_age_count
uptake_matrix
i_day <- nrow(db_uptake)+1
i_age <- seq_age_groups[1]
for(i_day in (nrow(db_uptake)+1):nrow(uptake_matrix))
for(i_age in seq_age_groups){

  remaining_uptake              <- target_uptake_age_count - colSums(uptake_matrix)
  available_doses               <- planned_doses_day - sum(uptake_matrix[i_day,],na.rm=T)
  #print(available_doses)
  
  if(i_age == 0){ # no more vaccination
    
  } else if(i_age == 88){ # uniform by age group
    
    target_age_groups <- !is.na(remaining_uptake) & remaining_uptake > 0
    
    age_uptake <- (available_doses / sum(target_age_groups)) * target_age_groups
    uptake_matrix[i_day, ] <- uptake_matrix[i_day, ] + pmin(remaining_uptake,age_uptake)
    
  } else if(i_age == 99){ # weighted by unvaccinated individuals in the population
    age_uptake <- rep(available_doses * remaining_uptake / sum(remaining_uptake,na.rm=T))
    uptake_matrix[i_day, ] <- pmin(remaining_uptake,age_uptake)
  }else{  # sequentially by age group
    uptake_matrix[i_day, i_age] <- pmin(remaining_uptake[i_age],available_doses)
  }
}
uptake_matrix[is.na(uptake_matrix)] <- 0
range(uptake_matrix)

# quick check
uptake_matrix_cum <- apply(uptake_matrix,2,cumsum)
uptake_matrix_cum / get_be_pop_matrix(nrow(uptake_matrix_cum))

# REFORMAT UPTAKE   ----
###################### # 

# go from 'month' to 'day'
scen_uptake_day    <- data.frame(date = uptake_matrix_date,
                                   uptake_matrix)

# set to 1/24 time steps
scen_out <- approx_uptake(scen_uptake_day,time_step_size,dose_col)

# write to file
write.table(scen_out,file=file.path(output_folder,paste0('vaccine_uptake_combined_Vx_Output_vSCEN',scen_tag,'.csv')),sep=',',row.names=F)

#SET PROTECTION   ----
###################### #

scen_out_delay <- uptake2protection(scen_out,delay_protection_rna,delay_protection_adeno)

write.table(scen_out_delay,file=file.path(output_folder,paste0('vaccine_protection_Vx_Output_vSCEN',scen_tag,'.csv')),sep=',',row.names=F)

#EXPLORE           ----
###################### #

plot(scen_out$date+1,
     cumsum(rowSums(scen_out[,(2:11)+10])),type='l',ylab='uptake',
     xlim=range(scen_out_delay$date),lwd=2)
points(ref_data_be$date+1,
       cumsum(ref_data_be$vaccine_uptake_first))

lines(scen_out_delay$date+1,
       cumsum(rowSums(scen_out_delay[,(2:11)+10])),col=4,lwd=2)

legend('topleft',
       c('Reported uptake: 1st or single dose',
         'Scenario uptake: 1st or single dose',
         'Scenario protection: 1st or single dose'),
       col = c(1,4,4),
       lwd= c(NA,2,2),
       pch=c(1,NA,NA)
)

plot(scen_out$date,
     (rowSums(scen_out[,(2:11)+10])/time_step_size),
     type='l',ylab='Uptake 1st/single dose (day)',
     main='Daily uptake: 1st/single dose')
points(ref_data_be$date,
       (ref_data_be$vaccine_uptake_first))



# add uptake per month
uptake_day <- scen_out
uptake_day$date_month_lab <- format(scen_out$date,'%b')
uptake_day$date_month <- format(scen_out$date,'%m')
uptake_day$uptake_sum <- rowSums(scen_out[,(2:11)+10])
uptake_month <- aggregate(uptake_sum ~ date_month,data=uptake_day,sum)
uptake_month <- uptake_month[-nrow(uptake_month),]
bplot <- barplot(uptake_month$uptake_sum,
                 col='grey',
                 ylab='Uptake 1st/single dose (month)',
                 main='Monthly uptake: 1st/single dose',
                 ylim=c(0,2e6))
grid()
axis(1,bplot,format(as.Date(paste0('2021-',uptake_month$date_month,'-01')),'%b'))



#EXPLORE BY AGE     ----
###################### #

age_opt     <- get_age_groups()
ref_age_col <- paste0('vaccine_A_ages',age_opt)
y_lim       <- c(0,1.5e6)
x_ticks     <- pretty(range(scen_out_delay$date),10)

perc_breaks    <- pretty(0:100,7)
y_ticks_perc   <- paste0(perc_breaks,'%')
par(mar=c(5,5,4,5),mfrow=c(2,2))

i_age <- 5
for(i_age in 10:1){

  y_ticks      <- perc_breaks * pop_data_be[i_age] / 100
  y_ticks_k    <- paste0(round(y_ticks/1e3),'k')
  y_ticks_k[1] <- 0
  y_ticks_k
  
  uptake_age <- cumsum(scen_out[,11+i_age])
  plot(ref_data_be$date+1,
       t(cumsum(ref_data_be[ref_age_col[i_age]])),
       ylab='Uptake',
       xlab='',
       xlim= range(scen_out_delay$date),
       ylim = range(y_ticks),
       xaxt='n',
       yaxt='n',
       main=paste0(gsub('_','-',age_opt[i_age]),'y'))
  lines(scen_out$date+1,
        cumsum(scen_out[,11+i_age]),
        col=4,
        lwd=2)
  axis(1,x_ticks,format(x_ticks,'%d/%m'))
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
  
  
  text(x = mean(x_ticks[4:5]),
       y = min(y_ticks),
       srt=90,
       '65+',
       pos=4,
       col=3,
       cex=0.7)
  
  text(x = mean(x_ticks[5:6]),
       y = min(y_ticks),
       srt=90,
       'risk pop',
       pos=4,
       col=3,
       cex=0.7)
  
  text(x = mean(x_ticks[6:7]),
       y = min(y_ticks),
       srt=90,
       'gen pop',
       pos=4,
       col=3,
       cex=0.7)
}

dev.off()
