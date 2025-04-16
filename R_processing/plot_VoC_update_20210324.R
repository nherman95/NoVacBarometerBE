########################################################################### #
# This file is part of the Stochastic Compartmental Model for SARS-COV-2 
# transmission in Belgium, conceived by members of SIMID group during the 
# COVID19 pandemic.
#
# This file is used to update the analysis from March 2021. These results have
# been described and presented in Technical Note v20210325.
#
# Copyright 2021, SIMID, University of Antwerp & Hasselt University                                        
########################################################################### #

# Note: this file contains a local copy of the plot functions to make sure that
# the figures do not change when the project code evolves over time. This is 
# not the most elegant way, but it is robust.

rm(list=ls())

library(scales)
library(foreach)
library(parallel)
source('R/lib_stochastic_model.R')

# SETTINGS ----

# SIMULATION DATA
#output_dir <- 'data/sm9_results/sm9_uptake2_d440_u100_ldwn100_n100/'
#output_dir <- 'data/sm9_results/sm9_uptake2_d440_u100_ldwn0_n100_oct2'
output_dirs <- c('data/results_technical_notes/sm9_results/sm9_GCS_d440_u100_ldwn100_n100/',
                 'data/results_technical_notes/sm9_results/sm9_GCS_d440_u100_ldwn0_n100',
                 'data/results_technical_notes/sm9_results/sm9_GCS_d440_u100_ldwn30_n100')

# OUTPUT FOLDER
output_folder <- 'output/figures_update_v20210325'
unlink(output_folder, recursive=TRUE)
dir.create(output_folder)

# set SCM start date
scm_start_date   <- as.Date('2020-03-01')
scm_scen_date    <- as.Date('2021-03-22')

# reference data
# be_ref_data <- get_observed_incidence_data()
be_ref_data <- get_latest_incidence_data()

# scenario options
# scen_col  <- "orange3" 
scen_cols  <- c("salmon2","darkblue","darkgrey")
scen_tags  <- paste0('scenario ',c('A','B','C'))

pdf_names <- paste0('technical note v20210325 hosp adm scenario',c('',' zoom'),'.pdf')

# PLOT ----
sel_file <- pdf_names[2]; i_scenario <- 1
# for(i_scenario in 1:length(output_dirs))
for(i_scenario in 1:2)
for(sel_file in pdf_names){
  
  # load data
  incidence_sum <- read.csv(file.path(output_dirs[i_scenario],'hosp_admin_green2.csv'),sep=',',header=T)
  scenario_file <- gsub('scenario',scen_tags[i_scenario],sel_file)
  
  pdf(file.path(output_folder,scenario_file),7,5)
  par(mar=c(4,4,2.5,0.5))
  
  # plot options
  scm_dates  <- scm_start_date +  0:(nrow(incidence_sum)-1)
  num_runs  <- ncol(incidence_sum)
  
  plot_start <- ifelse(grepl('zoom',sel_file),300,0)
  x_ticks   <- pretty(scm_start_date + c(plot_start,420),12)
  y_lim     <- c(1,ifelse(grepl('zoom',sel_file),550,1e3)) 
  
  plot(x = be_ref_data$date,
       y = be_ref_data$hospital_admissions,
       xlim = range(x_ticks),
       ylim = y_lim,
       xaxt='n',
       xlab='2020-2021',
       ylab= 'daily hospital admissions',
       main= 'Reference: technical note (v2021-03-25) at www.simid.be',
       col=0)
  axis(1,x_ticks,format(x_ticks,"%d/%m"))
  abline(v=x_ticks,col='grey',lty=3)
  grid(nx=NA,ny=NULL)
  add_copyright()
  
  # add line for callibration
  add_vertical_line(scm_scen_date,bool_text =T,date_tag= 'callibration')
  
  # add info
  legend('topleft',
         legend = c(paste0('Reported (until ',max(be_ref_data$date),')'),
                    paste0('Model projections (',scen_tags[i_scenario],')'),
                    paste0('Model average (',scen_tags[i_scenario],')')),
         col    = c('black',scen_cols[i_scenario],scen_cols[i_scenario]),
         lwd    = c(NA,1.5,4),
         pch    = c(1,NA,NA),
         bg = 'white',
         cex = 0.8)

  # add model trajectories
  for(i_run in 1:num_runs){
    lines(x   = scm_dates,
          y   = incidence_sum[,i_run],
          col = alpha(scen_cols[i_scenario],0.1))
  }
  lines(x   = scm_dates,
        y   = apply(incidence_sum,1,mean),
        lwd=3,
        col = scen_cols[i_scenario])
  
  # add reported data
  points(x = be_ref_data$date,
         y = be_ref_data$hospital_admissions)
  
  # close stream  
  dev.off()
}


## POLYGON - COMBINED


pdf(file.path(output_folder,gsub(' scenario','',pdf_names[1])),12,6)
#png(file.path(output_folder,gsub('pdf','png',gsub(' scenario','',pdf_names[1]))),width = (480*2)*1.2, height = 480*1.2)
par(mar=c(5,5,3,1))
par(mfrow=c(1,2))

sel_comparator <- 1
for(sel_comparator in c(1,3)){
  
  opt_scenario <- sort(c(2,sel_comparator))
  
  plot_start <- 370
  x_ticks   <- pretty(scm_start_date + c(plot_start,420),12)
  y_lim     <- c(1,700) 
  
  plot(x = be_ref_data$date,
       y = be_ref_data$hospital_admissions,
       xlim = range(x_ticks),
       ylim = y_lim,
       xaxt='n',
       xlab='2020-2021',
       ylab= 'Daily hospital admissions',
       main= 'Updated figure from technical note "SIMID v2021-03-25"',
       col=0)
  axis(1,x_ticks,format(x_ticks,"%d/%m"))
  abline(v=x_ticks,col='grey',lty=3)
  grid(nx=NA,ny=NULL)
  add_copyright(text.cex = 0.8)
  
  # add line for callibration
  add_vertical_line(scm_scen_date,bool_text =T,date_tag= 'callibration')
  
  # add info
  legend('topleft',
         legend = c(paste0('Reported (until ',max(be_ref_data$date),')'),
                    totitle(scen_tags[opt_scenario])),
         col    = c('black',scen_cols[opt_scenario]),
         lwd    = c(NA,rep(2,length(opt_scenario))),
         pch    = c(1,rep(NA,length(opt_scenario))),
         bg = 'white',
         cex = 0.8)
  
  for(i_scenario in rev(opt_scenario)){
    # load data
    incidence_sum <- read.csv(file.path(output_dirs[i_scenario],'hosp_admin_green2.csv'),sep=',',header=T)
    
    # plot options
    scm_dates  <- scm_start_date +  0:(nrow(incidence_sum)-1)
    num_runs  <- ncol(incidence_sum)
    
    hosp_aggr <-  data.frame(admin_mean=apply(incidence_sum,1,mean),
                             #admin_mean=NA,
                             admin_LL = apply(incidence_sum,1,min),
                             admin_UL = apply(incidence_sum,1,max))
    
    add_polygon(scenario_data  = hosp_aggr[-nrow(hosp_aggr),],
                scenario_dates = scm_dates[-length(scm_dates)],
                col = scen_cols[i_scenario],
                alpha_val = 0.5)
  }
  
  # add reported data
  points(x = be_ref_data$date,
         y = be_ref_data$hospital_admissions)

}
dev.off()

## COPY FILES TO GDRIVE
# update gdrive figures (external)
gdrive_path <- '~/Google\ Drive/nCov2019\ -\ research/Vaccine/stochastic_model/figures/shared_external/TechNote_v20210325'
if(!dir.exists(gdrive_path)){ dir.create(gdrive_path) }
output_file_names <- dir(output_folder,full.names = F)
for(sel_file in output_file_names){
  file.copy(from = file.path(output_folder,sel_file),
            to   = file.path(gdrive_path,sel_file),
            overwrite = TRUE)
  print(sel_file)
}
print('UPDATED')
print(range(be_ref_data$date))
print(pdf_names)
