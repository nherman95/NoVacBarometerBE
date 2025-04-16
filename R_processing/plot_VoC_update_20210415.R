########################################################################### #
# This file is part of the Stochastic Compartmental Model for SARS-COV-2 
# transmission in Belgium, conceived by members of SIMID group during the 
# COVID19 pandemic.
#
# This file is used to update the analysis from April 2021. These results 
# have been described and presented in Technical Note v20210415.
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
library(RColorBrewer)
source('R/lib_stochastic_model.R')

# Settings ----

# SIMULATION DATA
output_dir <- 'data/results_technical_notes/sm12_n100/'
 
output_files <- dir(output_dir,full.names = T,recursive = T)
output_files <- output_files[!grepl('pdf',output_files)]
output_files <- output_files[grepl('rds',output_files)]
output_files <- output_files[!grepl('20less',output_files)]
output_files <- output_files[!grepl('zero',output_files)]
output_files <- sort(output_files,decreasing = F)
#output_files

output_tag <- gsub('.*_incr','',output_files[2])
output_tag <- gsub('.rds','',output_tag)

# OUTPUT FOLDER (set name, remove and re-create)
output_folder <- paste0('output/figures_update_v20210415')
unlink(output_folder, recursive=TRUE)
dir.create(output_folder)
print(output_folder)

# output format
sel_output_format <- 'pdf'  # options: pdf, png, NA (=Rstudio)

# x-axis limits
x_axis_scm_day <- c(310,510)
x_axis_label <- paste(unique(format(sim_day2date(x_axis_scm_day),'Time (%Y)')),collapse='-')

# set SCM start date
scm_start_date   <- as.Date('2020-03-01')
scm_callibration_date    <- as.Date('2021-04-10')

# reference data
# be_ref_data <- get_observed_incidence_data()
be_ref_data <- get_latest_incidence_data()

col_lib <- data.frame(tags = c('xxxx','incr0_0_100','incr100_0_100','100_100_100_d426','100_100_100_d440'),
                      col  = c('black','darkred','darkgreen','darkblue','orange3'),
                      label = c(paste0('Reported (',max(be_ref_data$date),')'),
                                'No changes in behaviour',
                                'Changes on April 19-25',
                                'Changes on April 19-25 and May 1-7',
                                'Changes on April 19-25 and May 15-21'),
                      lwd = c(NA,2,2,2,2),
                      pch = c(1,NA,NA,NA,NA))

# # remove 20% less
# if(!any(grepl('20less',output_files))){
#   col_lib <- col_lib[-3,]
# } 

get_col <- function(file_name){
  flag_file <- unlist(lapply(col_lib$tags,grepl,file_name))
  return(col_lib$col[flag_file])
}

add_scen_legend <- function(bool_hosp,opt_files = NA,ref_pch = 1,ref_label = ''){
  
  legend_lib <- col_lib
  
  legend_lib$pch[1] <- ref_pch
  
  if(nchar(ref_label)>0){
    legend_lib$label[1] <- gsub(' ',paste0(c(' ',ref_label,' '),collapse ='' ),legend_lib$label[1])
  }
  
  if(!bool_hosp){
    legend_lib <- legend_lib[-1,]
  }
  
  if(!any(is.na(opt_files))){
    legend_col <- unlist(lapply(opt_files,get_col))
    legend_lib <- legend_lib[legend_lib$col %in% c(col_lib$col[1],legend_col),]
  }
  
  legend('topleft',
         legend = legend_lib$label,
         col    = legend_lib$col,
         lwd    = legend_lib$lwd,
         pch    = legend_lib$pch,
         bg = 'white',
         cex = 1)
}

# output_files_inc <- output_files_hosp; bool_hosp <- T; bool_polygon <- F
# output_files_inc <- output_files_cases; bool_hosp <- F; bool_polygon <- T
plot_incidence_time <- function(output_files_inc, bool_hosp, bool_polygon){
  
  x_ticks   <- pretty(scm_start_date + x_axis_scm_day,12)
  y_lim     <- c(1,ifelse(bool_hosp,5e2,6e4)) 
  ref_pch   <- 1
  ref_label <- ifelse(bool_hosp,'hospital admissions','cases') 
  
  plot(x = be_ref_data$date,
       y = be_ref_data$hospital_admissions,
       xlim = range(x_ticks),
       ylim = y_lim,
       xaxt='n',
       xlab=x_axis_label, 
       ylab= ifelse(bool_hosp,'Daily hospital admissions','Daily infections'),
       main= 'Updated figure from Technical Note SIMID v2021-04-15',
       col=0,
       las=2)
  axis(1,x_ticks,format(x_ticks,"%d/%m"),cex.axis=0.9)
  abline(v=x_ticks,col='grey',lty=3)
  grid(nx=NA,ny=NULL)
  add_vertical_line(scm_callibration_date,bool_text =T,date_tag= 'callibration')
  add_scen_legend(bool_hosp,output_files_inc,ref_pch = ref_pch,ref_label=ref_label)
  
  i_file <- 1
  for(i_file in 1:length(output_files_inc)){
    scen_data <- readRDS(output_files_inc[i_file])
    scen_col  <- get_col(output_files_inc[i_file])
    scm_dates  <- scm_start_date +  0:(nrow(scen_data)-1)
    
  if(bool_polygon){
      hosp_aggr <-  data.frame(admin_mean=apply(scen_data,1,mean),
                               #admin_mean=NA,
                               admin_LL = apply(scen_data,1,min),
                               admin_UL = apply(scen_data,1,max))
      
      add_polygon(scenario_data  = hosp_aggr[-nrow(hosp_aggr),],
                  scenario_dates = scm_dates[-length(scm_dates)],
                  col = scen_col,
                  alpha_val = 0.5,
                  mean_lty = 1)
    } else{
      for(i_run in 1:ncol(scen_data)){
        lines(x   = scm_dates,
              y   = scen_data[,i_run],
              col = alpha(scen_col,0.1))
      }
      lines(x   = scm_dates,
            y   = apply(scen_data,1,mean),
            lwd=3,
            col = scen_col)
    }
  }
  
  if(bool_hosp){
    points(x = be_ref_data$date,
           y = be_ref_data$hospital_admissions,
           pch=ref_pch,
           cex = 0.8)
    abline(h=100,lty=3,col='grey')
  } 
  add_copyright(text.cex = 0.8)
}

plot_load_time <- function(output_files_inc, bool_polygon){
  
  x_ticks   <- pretty(scm_start_date + x_axis_scm_day,12)
  y_lim     <- c(1,5000) 
  ref_pch   <- 2 
  ref_label <- 'hospital load'
  
  plot(x = be_ref_data$date,
       y = be_ref_data$hospital_load,
       xlim = range(x_ticks),
       ylim = y_lim,
       xaxt='n',
       xlab=x_axis_label, 
       ylab= 'Daily hospital load',
       col=0,
       las = 2)
  axis(1,x_ticks,format(x_ticks,"%d/%m"),cex.axis=0.9)
  abline(v=x_ticks,col='grey',lty=3)
  grid(nx=NA,ny=NULL)
  add_vertical_line(scm_callibration_date,bool_text =T,date_tag= 'callibration')
  add_scen_legend(TRUE,output_files_inc,ref_pch=ref_pch,ref_label=ref_label)
  
  i_file <- 1
  for(i_file in 1:length(output_files_inc)){
    scen_data <- readRDS(output_files_inc[i_file])
    scen_col  <- get_col(output_files_inc[i_file])
    scm_dates  <- scm_start_date +  0:(nrow(scen_data)-1)
    
    if(bool_polygon){
      hosp_aggr <-  data.frame(admin_mean=apply(scen_data,1,mean),
                               #admin_mean=NA,
                               admin_LL = apply(scen_data,1,min),
                               admin_UL = apply(scen_data,1,max))
      
      add_polygon(scenario_data  = hosp_aggr[-nrow(hosp_aggr),],
                  scenario_dates = scm_dates[-length(scm_dates)],
                  col = scen_col,
                  alpha_val = 0.5,
                  mean_lty = 2)
    } else{
      for(i_run in 1:ncol(scen_data)){
        lines(x   = scm_dates,
              y   = scen_data[,i_run],
              col = alpha(scen_col,0.2))
      }
      lines(x   = scm_dates,
            y   = apply(scen_data,1,mean),
            lwd=3,
            col = scen_col)
    }
  }
  
    points(x = be_ref_data$date,
           y = be_ref_data$hospital_load,
           pch = ref_pch,
           cex = 0.8)
    abline(h=1000,lty=3,col='grey')
    add_copyright(text.cex = 0.8)
}

plot_icu_load <- function(output_files_inc){
  
  x_ticks   <- pretty(scm_start_date + x_axis_scm_day,12)
  y_lim     <- c(1,1200) 
  ref_pch   <- 3 
  ref_label <- 'ICU load'
  
  plot(x = be_ref_data$date,
       y = be_ref_data$icu_load,
       xlim = range(x_ticks),
       ylim = y_lim,
       xaxt='n',
       xlab=x_axis_label,
       ylab= 'Daily ICU load',
       col=0,
       las = 2)
  axis(1,x_ticks,format(x_ticks,"%d/%m"),cex.axis=0.9)
  abline(v=x_ticks,col='grey',lty=3)
  grid(nx=NA,ny=NULL)
  add_vertical_line(scm_callibration_date,bool_text =T,date_tag= 'callibration',pos_factor = 0.1)
  add_scen_legend(TRUE,output_files_inc,ref_pch=ref_pch,ref_label=ref_label)
  abline(h=seq(100,1500,100),col='grey',lty=3)
  
  i_file <- 1
  for(i_file in 1:length(output_files_inc)){
    scen_data <- readRDS(output_files_inc[i_file])
    scen_col  <- get_col(output_files_inc[i_file])
    scm_dates  <- scm_start_date +  0:(ncol(scen_data)-1)
    
    hosp_aggr <-  data.frame(admin_mean=scen_data[2,],
                             #admin_mean=NA,
                             admin_LL = scen_data[1,],
                             admin_UL = scen_data[3,])
    
    add_polygon(scenario_data  = hosp_aggr[-nrow(hosp_aggr),],
                scenario_dates = scm_dates[-length(scm_dates)],
                col = scen_col,
                alpha_val = 0.5,
                mean_lty = 0)
  }
  
  points(x = be_ref_data$date,
         y = be_ref_data$icu_load,
         pch = ref_pch,
         cex = 0.8)
  #abline(h=c(1000),lty=3,col='grey')
  add_copyright(text.cex = 0.8)
  
}

#output_files_inc <- output_files_cases; bool_hosp <- F; bool_polygon <- T
plot_Rt <- function(output_files_inc, bool_hosp, bool_polygon){
  
  #x_ticks   <- pretty(scm_start_date + c(0,450),12)
  x_ticks   <- pretty(scm_start_date + x_axis_scm_day,12)
  y_lim     <- c(0.5,1.7) 
  ref_pch   <- 5
  ref_label <- 'cases ~ Rt'
  
  plot(x = range(x_ticks),
       y = y_lim,
       xlim = range(x_ticks),
       ylim = y_lim,
       xaxt='n',
       xlab=x_axis_label,
       ylab= 'Rt (cases)',
       main= 'Updated figure from Technical Note SIMID v2021-04-15',
       col=0)
  axis(1,x_ticks,format(x_ticks,"%d/%m")) 
  abline(v=x_ticks,col='grey',lty=3)
  grid(nx=NA,ny=NULL)
  #add_vertical_line(scm_callibration_date,bool_text =T,date_tag= 'callibration',pos_factor = 0.4)
  abline(h=1)
  
  i_file <- 2
  for(i_file in 1:length(output_files_inc)){
    scen_data <- readRDS(output_files_inc[i_file])
    scen_col  <- get_col(output_files_inc[i_file])
    scm_dates  <- scm_start_date +  0:(nrow(scen_data)-1)
    
    if(bool_polygon){
      hosp_aggr <-  data.frame(admin_mean=apply(scen_data,1,mean),
                               # admin_mean=NA,
                               admin_LL = apply(scen_data,1,min),
                               admin_UL = apply(scen_data,1,max))
      
      add_polygon(hosp_aggr[-nrow(hosp_aggr),],
                  scenario_dates = scm_dates[-length(scm_dates)],
                  col = scen_col,
                  alpha_val = 0.5,
                  mean_lty = 3)
    } else{
      for(i_run in 1:ncol(scen_data)){
        lines(x   = scm_dates,
              y   = scen_data[,i_run],
              col = alpha(scen_col,0.2))
      }
    }
  }
  
  points(x = be_ref_data$date,
         y = be_ref_data$Rt_cases,
         pch = ref_pch,
         cex = 0.8)
  
  add_scen_legend(bool_hosp,opt_files = output_files_inc,
                  ref_pch=ref_pch,ref_label=ref_label)
}

## SET FILE NAMES ----

output_files_hosp_admin <- output_files[grepl('/hosp_adm',output_files)]
output_files_hosp_load  <- output_files[grepl('hosp_load',output_files)]
output_files_icu        <- output_files[grepl('icu_load',output_files) & grepl('rds',output_files)]
output_files_cases      <- output_files[grepl('Rt_mild_sev_incr',output_files) & grepl('rds',output_files)]

## scenario A & B ----
pdf(file.path(output_folder,paste0('technical_note_v20210315_fig3_scen_AB.pdf')),9,12)
par(mar=c(4.5,4.5,2,0.5),mfrow=c(3,1))

plot_incidence_time(output_files_hosp_admin[1:2],bool_hosp = T,bool_polygon = T)
plot_load_time(output_files_hosp_load[1:2],bool_polygon = T)
plot_icu_load(output_files_icu[1:2])

dev.off()

## scenario C & D ----
pdf(file.path(output_folder,paste0('technical_note_v20210315_fig4_scen_CD.pdf')),9,12)
par(mar=c(4.5,4.5,2,0.5),mfrow=c(3,1))

plot_incidence_time(output_files_hosp_admin[4:3],bool_hosp = T,bool_polygon = T)
plot_load_time(output_files_hosp_load[4:3],bool_polygon = T)
plot_icu_load(output_files_icu[4:3])

dev.off()

# Rt ----
pdf(file.path(output_folder,paste0('technical_note_v20210315_fig8_Rt.pdf')),9,8)
par(mar=c(4.5,4.5,2,0.5),mfrow=c(2,1))
plot_Rt(output_files_cases[1:2], bool_hosp =T, bool_polygon = T)
plot_Rt(output_files_cases[4:3], bool_hosp =T, bool_polygon = T)

dev.off()


## COPY FILES TO GDRIVE ----
# update gdrive figures (external)
gdrive_path <- '~/Google\ Drive/nCov2019\ -\ research/Vaccine/stochastic_model/figures/shared_external/TechNote_v20210415'
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




