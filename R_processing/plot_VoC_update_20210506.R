########################################################################### #
# This file is part of the Stochastic Compartmental Model for SARS-COV-2 
# transmission in Belgium, conceived by members of SIMID group during the 
# COVID19 pandemic.
#
# This file is used to update the analysis from May 2021. These results 
# have been described and presented in Technical Note v2021506.
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
source('R/lib_social_contacts.R')

# Settings ----
bool_all_plots <- TRUE
bool_uptake_scen <- FALSE

# SIMULATION DATA
output_dir <- 'data/results_technical_notes/sm13_n100/'
 
output_files <- dir(output_dir,full.names = T,recursive = T)
output_files <- output_files[!grepl('pdf',output_files)]
output_files <- output_files[grepl('rds',output_files)]
output_files <- output_files[!grepl('20less',output_files)]
output_files <- output_files[!grepl('zero',output_files)]
output_files <- sort(output_files,decreasing = F)
# output_files

output_tag <- gsub('.*_incr','',output_files[2])
output_tag <- gsub('.rds','',output_tag)

# OUTPUT FOLDER (set name, remove and re-create)
output_folder <- paste0('output/figures_update_v20210506')
unlink(output_folder, recursive=TRUE)
dir.create(output_folder)
print(output_folder)

# numeric output
csv_output_file <- file.path(output_folder,'technical_note_v20210506_numeric_results.csv')
if(file.exists(csv_output_file)){file.remove(csv_output_file)}

# x-axis limits
x_axis_scm_day <- c(310,510)
x_axis_label <- paste(unique(format(sim_day2date(x_axis_scm_day),'Time (%Y)')),collapse='-')

# set SCM start date
scm_start_date   <- as.Date('2020-03-01')
scm_callibration_date    <- as.Date('2021-05-01')

bool_add_cp    <- TRUE
scm_scenario_cp <- c(as.Date('2021-04-19') + 0:6,
                    as.Date('2021-05-01') + 0:6,
                    as.Date('2021-05-15') + 0:6,
                    as.Date('2021-07-01'))

# optional
scm_scenario_cp_AB <- c(as.Date('2021-04-19') + 0:6,
                       as.Date('2021-07-01'))
#scm_scenario_cp <- scm_scenario_cp_AB

# reference data
be_ref_data <- get_latest_incidence_data()
#be_ref_data <- get_observed_incidence_data()

col_lib <- data.frame(tags = c('xxxx',
                               'cnt0_0_100',
                               'cnt100_0_100',
                               'cnt100_100_100_d426',
                               'cnt100_100_100_d440'),
                      col  = c('black','darkred','darkgreen','darkblue','orange3'),
                      label = c(paste0('Reported (',max(be_ref_data$date),')'),
                                'No changes in risk behaviour after April 18th',
                                'Changes on April 19-25 (~March2021)',
                                'Changes on April 19-25 (~March2021) and May 1-7 (~Sept2020)',
                                'Changes on April 19-25 (~March2021) and May 15-21 (~Sept2020)'),
                      lwd = c(NA,2,2,2,2),
                      pch = c(1,NA,NA,NA,NA))

if(bool_uptake_scen){
  col_lib <- data.frame(tags = c('xxxx',
                                 'uc80_ur100',
                                 'uc70_ur100',
                                 'uc60_ur100',
                                 'uc80_ur80',
                                 'uc80_ur70',
                                 'changepoint'),
                        # col  = c('black','darkred','darkgreen','darkblue','orange3','purple'),
                        col  = c('black',brewer.pal(5,"Accent"),alpha('grey',0.5)),
                        label = c(paste0('Reported (',max(be_ref_data$date),')'),
                                  'Planned uptake (80%)',
                                  'Cap uptake <70y at 70%',
                                  'Cap uptake <70y at 60%',
                                  'Slower (80%) uptake from May 1st',
                                  'Slower (70%) uptake from May 1st',
                                  'Behavioural changes possible'),
                        lwd = c(NA,2,2,2,2,2,NA),
                        pch = c(1,NA,NA,NA,NA,NA,15))
} else{
  #output_files <- output_files[grepl('uc60_ur100',output_files)]
  output_files <- output_files[grepl('uc80_ur100',output_files)]
  #output_files <- output_files[grepl('uc80_ur70',output_files)]
}

get_col <- function(file_name){
  flag_file <- unlist(lapply(col_lib$tags,grepl,file_name))
  if(!any(flag_file)) {warning(paste('NO COLOR FOR:',file_name))}
  return(col_lib$col[flag_file])
}

get_label <- function(file_name){
  flag_file <- unlist(lapply(col_lib$tags,grepl,file_name))
  if(!any(flag_file)) {warning(paste('NO LABEL FOR:',file_name))}
  return(col_lib$label[flag_file])
}

add_scen_legend <- function(opt_files = NA,ref_pch = 1,ref_label = ''){
  
  legend_lib <- col_lib
  
  legend_lib$pch[1] <- ref_pch
  
  if(nchar(ref_label)>0){
    legend_lib$label[1] <- gsub(' ',paste0(c(' ',ref_label,' '),collapse ='' ),legend_lib$label[1])
  }
  
  if(!any(is.na(opt_files))){
    legend_col <- unlist(lapply(opt_files,get_col))
    legend_lib <- legend_lib[legend_lib$col %in% c(col_lib$col[1],legend_col),]
  }
  
  if(bool_add_cp){
    legend_lib <- rbind(legend_lib,
                        col_lib[col_lib$tags == 'changepoint',])
  }
  
  legend(ifelse(min(x_axis_scm_day)==0,'topright','topleft'), # 
         legend = legend_lib$label,
         col    = legend_lib$col,
         lwd    = legend_lib$lwd,
         pch    = legend_lib$pch,
         bg = 'white',
         cex = 0.6)
}

add_change_points <- function(scm_scenario_cp){
  if(bool_add_cp){
    cp <- get_M_change_day()
    cp_date <- sim_day2date(cp)
    
    abline(v=c(rep(cp_date,each=5)+(0:4)),
           lty=1,col=alpha('grey',0.4),lwd=4)
    
    abline(v=scm_scenario_cp,
           lty=1,col=alpha('black',0.4),lwd=4)
  }
}

# output_files_inc <- output_files_hosp; bool_polygon <- T
# output_files_inc <- output_files_cases; bool_polygon <- T
plot_incidence_time <- function(output_files_inc, bool_polygon, opt_cp = NA){
  
  if(length(output_files_inc)==0){
    return(NULL)
  }
  
  y_abline <- NA
  
  # set reference
  if(any(grepl('hosp_adm',output_files_inc))){
    plot_ref_data <- data.frame(date     = be_ref_data$date,
                                reported = be_ref_data$hospital_admissions)
    ref_label    <- 'hospital admissions'
    ref_pch      <- 1
    y_lim        <- c(0,ifelse(min(x_axis_scm_day)==0,800,450))
    y_abline     <- 100
  }
  if(any(grepl('new_infect',output_files_inc))){
    plot_ref_data <- data.frame(date     = be_ref_data$date,
                                reported = be_ref_data$cases)
    ref_label    <- 'infections'
    ref_pch      <- 5 
    y_lim    <- c(0,ifelse(min(x_axis_scm_day)==0,4e4,2e4))
  }
  if(any(grepl('mortality',output_files_inc))){
    plot_ref_data <- data.frame(date     = be_ref_data$date,
                                reported = be_ref_data$covid19_deaths_hospital)
    ref_label    <- 'deaths (non NH)'
    ref_pch      <- 4
    y_lim    <- c(0,ifelse(min(x_axis_scm_day)==0,300,100))
  }
  if(any(grepl('hosp_load',output_files_inc))){
    plot_ref_data <- data.frame(date     = be_ref_data$date,
                                reported = be_ref_data$hospital_load)
    y_lim     <- c(1,ifelse(min(x_axis_scm_day)==0,7500,4500)) 
    ref_pch   <- 2 
    ref_label <- 'hospital load'
  }
  if(any(grepl('icu_load',output_files_inc))){
    plot_ref_data <- data.frame(date     = be_ref_data$date,
                                reported = be_ref_data$icu_load)
    y_lim     <- c(1,1500) 
    ref_pch   <- 3
    ref_label <- 'ICU load'
    y_abline  <- 500
  }
  if(any(grepl('Rt',output_files_inc))){
    plot_ref_data <- data.frame(date     = be_ref_data$date,
                                reported = be_ref_data$Rt_cases)
    y_lim     <- c(0.5,1.7) 
    ref_pch   <- 6
    ref_label <- 'cases ~ Rt'
    y_abline  <- 1
  }
  
  x_ticks   <- pretty(scm_start_date + x_axis_scm_day,12)
  
  plot(x = be_ref_data$date,
       y = be_ref_data$hospital_admissions,
       xlim = range(x_ticks),
       ylim = y_lim,
       xaxt='n',
       xlab=x_axis_label, 
       ylab= paste('Daily',ref_label),
       main= 'Updated figure from Technical Note SIMID v2021-05-06',
       col=0,
       las=2)
  axis(1,x_ticks,format(x_ticks,"%d/%m"),cex.axis=0.9)
  abline(v=x_ticks,col='grey',lty=3)
  grid(nx=NA,ny=NULL)
  add_vertical_line(scm_callibration_date,bool_text =T,date_tag= 'callibration')
  abline(h=y_abline,lty=2)
  if(any(is.na(opt_cp))){
    add_change_points(scm_scenario_cp)
  } else {
    add_change_points(opt_cp)
  }
  add_scen_legend(output_files_inc,ref_pch = ref_pch,ref_label=ref_label)
  
  i_file <- 2
  for(i_file in 1:length(output_files_inc)){
    scen_data <- readRDS(output_files_inc[i_file])
    scen_col  <- get_col(output_files_inc[i_file])
    scm_dates  <- scm_start_date +  0:(nrow(scen_data)-1)
    
    if(bool_polygon){
      hosp_aggr <-  data.frame(admin_mean=apply(scen_data,1,mean),
                               #admin_mean=NA,
                               admin_LL = apply(scen_data,1,quantile,0.025,na.rm=T),
                               admin_UL = apply(scen_data,1,quantile,0.975,na.rm=T))
      # fix for ICU load
      if(nrow(scen_data)==3){
        hosp_aggr <-  data.frame(admin_mean=scen_data[2,],
                                 #admin_mean=NA,
                                 admin_LL = scen_data[1,],
                                 admin_UL = scen_data[3,])
        scm_dates  <- scm_start_date +  0:(ncol(scen_data)-1)
      }
      
      add_polygon(scenario_data  = hosp_aggr[-nrow(hosp_aggr),],
                  scenario_dates = scm_dates[-length(scm_dates)],
                  col = scen_col,
                  alpha_val = 0.25,
                  mean_lty = ref_pch)
      
      
      if(!file.exists(csv_output_file)){
        write.table(t(c('scenario',
                          'clinical_endpoint',
                          as.character(scm_dates))),
                    file = csv_output_file,
                    sep=',',
                    append = T,
                    row.names=F,
                    col.names=F)    
      }
      # add (mean, LL,UL) to csv file
      write.table(cbind(get_label(output_files_inc[i_file]),
                           paste(ref_label,c(': mean',': lower limit',': upper limit')),
                           t(hosp_aggr)),
            file = csv_output_file,
            sep=',',
            append = T,
            row.names=F,
            col.names=F)      
      
    } else{
      # for(i_run in 1:ncol(scen_data)){
      #   lines(x   = scm_dates,
      #         y   = scen_data[,i_run],
      #         col = alpha(scen_col,0.1))
      # }
      if(nrow(scen_data)==3){
        scm_dates  <- scm_start_date +  0:(ncol(scen_data)-1)
        scen_data <- cbind(scen_data[2,],
                           scen_data[2,])
      }
      lines(x   = scm_dates,
            y   = apply(scen_data,1,mean),
            lwd=3,
            col = scen_col)
    }
  }
  
  # reference data
  points(x = plot_ref_data$date,
         y = plot_ref_data$reported,
         pch=ref_pch,
         cex = 0.8)
  #abline(h=100,lty=3,col='grey')
  
  add_copyright()
}

## SET FILE NAMES ----

output_files_hosp_admin <- output_files[grepl('/hosp_adm',output_files)]
output_files_hosp_load  <- output_files[grepl('hosp_load',output_files)]
output_files_icu        <- output_files[grepl('icu_load',output_files) & grepl('rds',output_files)]
output_files_rt         <- output_files[grepl('Rt_mild_sev_incr',output_files) & grepl('rds',output_files)]

## scenario A & B ----
pdf(file.path(output_folder,paste0('technical_note_v20210506_fig4_scen_AB.pdf')),9,12)
par(mar=c(4.5,4.5,2,0.5),mfrow=c(3,1))

plot_incidence_time(output_files_hosp_admin[1:2],bool_polygon = T)
plot_incidence_time(output_files_hosp_load[1:2],bool_polygon = T)
plot_incidence_time(output_files_icu[1:2],bool_polygon = T)

dev.off()

## scenario C & D ----
pdf(file.path(output_folder,paste0('technical_note_v20210506_fig5_scen_CD.pdf')),9,12)
par(mar=c(4.5,4.5,2,0.5),mfrow=c(3,1))

plot_incidence_time(output_files_hosp_admin[4:3],bool_polygon = T)
plot_incidence_time(output_files_hosp_load[4:3],bool_polygon = T)
plot_incidence_time(output_files_icu[4:3],bool_polygon = T)

dev.off()

# Rt ----
pdf(file.path(output_folder,paste0('technical_note_v20210506_fig8_Rt.pdf')),9,8)
par(mar=c(4.5,4.5,2,0.5),mfrow=c(2,1))
plot_incidence_time(output_files_rt[1:2], bool_polygon = T)
plot_incidence_time(output_files_rt[4:3],bool_polygon = T)

dev.off()


## COPY FILES TO GDRIVE ----
# update gdrive figures (external)
gdrive_path <- '~/Google\ Drive/nCov2019\ -\ research/Vaccine/stochastic_model/figures/shared_external/TechNote_v20210506'
#gdrive_path <- '~/Google\ Drive/nCov2019\ -\ research/Vaccine/stochastic_model/figures/shared_external/TechNote_latest'
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




