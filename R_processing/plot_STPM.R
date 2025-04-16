########################################################################### #
# This file is part of the Stochastic Compartmental Model for SARS-COV-2 
# transmission in Belgium, conceived by members of SIMID group during the 
# COVID19 pandemic.
#
# This file is used to update the analysis from November 2021. These results 
# have been described and presented in Technical Note v20211116.
#
# Copyright 2021, SIMID, University of Antwerp & Hasselt University                                        
########################################################################### #

rm(list=ls())

library(scales)
library(foreach)
library(parallel)
library(RColorBrewer)
source('R/lib_stochastic_model.R')
source('R/lib_social_contacts.R')
source('R/lib_combined_plot.R')



# SIMULATION DATA
output_dir <- '~/Google Drive/nCov2019 - research/STPM/'

output_files <- dir(output_dir,full.names = T,recursive = T,pattern = '.csv')
output_files <- sort(output_files,decreasing = F)

# OUTPUT FOLDER (set name, remove and re-create)
output_folder <- paste0('output/figures_update_stpm')
unlink(output_folder, recursive=TRUE)
dir.create(output_folder)
print(output_folder)

# reference data
#be_ref_data <- get_observed_incidence_data(sel_region=sel_region)
be_ref_data <- get_latest_incidence_data()

x_lim <- max(be_ref_data$date,na.rm = T) + c(-60,20)

plot(be_ref_data$date,
     be_ref_data$hospital_admissions,
     xlab='Date',
     ylab='Hospital admissions for COVID-19',
     xlim = x_lim)
add_date_axis(x_lim)

for(i in 1:length(output_files)){
  stpm_data <- read.table(output_files[i],sep=',',header = T)
  
  stpm_data$date <- as.Date(stpm_data$Date)
  
  hosp_aggr <-  data.frame(admin_mean=stpm_data$Pred.NewHosp,
                           #admin_mean=NA,
                           admin_LL = stpm_data$Min.NewHosp,
                           admin_UL = stpm_data$Max.NewHosp)
  
  add_polygon(scenario_data  = hosp_aggr,
              scenario_dates = stpm_data$date,
              col = ifelse(i == length(output_files),2,3),
              alpha_val = 0.1,
              mean_lty = 1)
}

plot(be_ref_data$date,
     be_ref_data$hospital_admissions + be_ref_data$hospital_admissions_other,
     xlab='Date',
     ylab='Hospital admissions with COVID-19',
     xlim = x_lim)
add_date_axis(x_lim)

for(i in 1:length(output_files)){
  stpm_data <- read.table(output_files[i],sep=',',header = T)
  
  stpm_data$date <- as.Date(stpm_data$Date)
  
  hosp_aggr <-  data.frame(admin_mean=stpm_data$Pred.NewHosp,
                           #admin_mean=NA,
                           admin_LL = stpm_data$Min.NewHosp,
                           admin_UL = stpm_data$Max.NewHosp)
  
  add_polygon(scenario_data  = hosp_aggr,
              scenario_dates = stpm_data$date,
              col = ifelse(i == length(output_files),2,3),
              alpha_val = 0.1,
              mean_lty = 1)
}

 
## ICU
plot(be_ref_data$date,
     be_ref_data$icu_load,
     xlim = x_lim,
     xlab='Date',
     ylab='ICU beds for COVID-19',
     xaxt='n')
add_date_axis(x_lim)

for(i in 1:length(output_files)){
  stpm_data <- read.table(output_files[i],sep=',',header = T)
  
  stpm_data$date <- as.Date(stpm_data$Date)
  
  hosp_aggr <-  data.frame(admin_mean=stpm_data$Pred.BedICU,
                           #admin_mean=NA,
                           admin_LL = stpm_data$Min.BedICU,
                           admin_UL = stpm_data$Max.BedICU)
  
  add_polygon(scenario_data  = hosp_aggr,
              scenario_dates = stpm_data$date,
              col = ifelse(i == length(output_files),2,3),
              alpha_val = 0.1,
              mean_lty = 1)
}
