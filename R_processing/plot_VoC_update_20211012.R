########################################################################### #
# This file is part of the Stochastic Compartmental Model for SARS-COV-2 
# transmission in Belgium, conceived by members of SIMID group during the 
# COVID19 pandemic.
#
# This file is used to update the analysis from October 2021. These results 
# have been described and presented in Technical Note v20211012, which was 
# part of a GEMS report.
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
output_dir <- 'data/results_technical_notes/sm21_n40_bis/'

# select region
# sel_region <- 'brussels'
# sel_region <- 'belgium'

 
output_files <- dir(output_dir,full.names = T,recursive = T)
#output_files <- output_files[grepl(sel_region,output_files)]
output_files <- output_files[!grepl('_20',output_files)]
# output_files <- output_files[grepl('rds',output_files)]
# output_files <- output_files[!grepl('20less',output_files)]
# output_files <- output_files[!grepl('zero',output_files)]
# output_files <- output_files[!grepl('100_0_0',output_files)]
# output_files <- sort(output_files,decreasing = F)
# output_files

output_tag <- gsub('.*_incr','',output_files[2])
output_tag <- gsub('.rds','',output_tag)

# OUTPUT FOLDER (set name, remove and re-create)
output_folder <- paste0('output/figures_update_v20211012')
unlink(output_folder, recursive=TRUE)
dir.create(output_folder)
print(output_folder)

# numeric output
csv_output_file <- file.path(output_folder,'technical_note_v20211012_numeric_results.csv')
if(file.exists(csv_output_file)){file.remove(csv_output_file)}

# x-axis limits
x_axis_scm_day <- c(500,670)
x_axis_label <- paste(unique(format(sim_day2date(x_axis_scm_day),'Time (%Y)')),collapse='-')
default_y_axis <- TRUE

# set SCM start date
scm_start_date   <- as.Date('2020-03-01')
scm_callibration_date    <- as.Date('2021-10-08')

# number of CoMix and additional waves
scm_num_cp <- 30+6

bool_add_cp    <- TRUE
scm_scenario_cp <- c(as.Date('2021-10-01') )

# reference data
#be_ref_data <- get_observed_incidence_data(sel_region=sel_region)
be_ref_data <- get_latest_incidence_data()

col_lib <- data.frame(tags = c('xxxx',
                               'cnt0_0_0_0',
                               'cnt0_0_0_120',
                               'cnt0_0_0_130',
                               'cnt0_0_0_140',
                               'cnt0_0_0_160',
                               'short.*cnt0_0_0_0',
                               'short.*cnt0_0_0_0.*CAP',
                               'short.*cnt0_0_0_120.*CAP',
                               'final.*cnt0_0_0_0.*brussels'),
                      col  = c('black','darkred','darkblue','darkgreen','orange3','grey','darkred','darkgreen','orange3','darkred'),
                      label = c(paste0('Reported (',max(be_ref_data$date),')'),
                                'No change in behaviour and no introductions',
                                'Change on October 1st (~120% last wave)',
                                'Change on October 1st (~130% last wave)',
                                'Change on October 1st (~140% last wave)',
                                'Change on October 1st (~160% last wave)',
                                'Based on reported uptake and social contact behaviour',
                                'No additional vaccine uptake after July 1st',
                                'No additional vaccine uptake after July 1st & increased risk behaviour (120%) from September 1st',
                                'Estimated behaviour and no introductions'
                      ),
                      lwd = c(NA,2,2,2,2,2,2,2,2,2),
                      pch = c(1,NA,NA,NA,NA,NA,NA,NA,NA,NA))



# set SCM start date
get_scm_start_date <- function(){
  return(as.Date('2020-03-01'))
}

get_col <- function(file_name,col_lib){
  flag_file <- unlist(lapply(col_lib$tags,grepl,file_name))
  if(!any(flag_file)) {
    warning(paste('NO COLOR FOR:',file_name))
    return(8)
  }
  if(sum(flag_file)>0){
    flag_file <- max(which(flag_file))
  }
  return(col_lib$col[flag_file])
}

add_scen_legend <- function(col_lib, opt_files = NA,ref_pch = 1,ref_label = '',ncol = ncol, ref_date = NA){
  
  legend_lib <- col_lib
  
  legend_lib$pch[1] <- ref_pch
  
  if(!is.na(ref_date)){
    legend_lib$label[1] <- gsub('(....-..-..)',ref_date,legend_lib$label[1])
  }
  
  if(nchar(ref_label)>0){
    legend_lib$label[1] <- gsub(' ',paste0(c(' ',ref_label,' '),collapse ='' ),legend_lib$label[1])
  }
  
  if(!any(is.na(opt_files))){
    legend_col <- unlist(lapply(opt_files,get_col,col_lib))
    legend_lib <- legend_lib[legend_lib$col %in% c(col_lib$col[1],legend_col),]
  }
  
  # add changepoint if present
  if(any(col_lib$tags == 'changepoint')){
    legend_lib <- rbind(legend_lib,
                        col_lib[col_lib$tags == 'changepoint',])
  }
  
  legend(ifelse(min(x_axis_scm_day)==0,'topright','topleft'), # 
         legend = legend_lib$label,
         col    = legend_lib$col,
         lwd    = legend_lib$lwd,
         pch    = legend_lib$pch,
         bg = 'white',
         ncol = ncol,
         cex = 0.6)
}

add_change_points <- function(scm_scenario_cp){
  cp      <- get_M_change_day()
  cp_date <- sim_day2date(cp)
  
  # select
  cp_date <- cp_date[1:(scm_num_cp+1)] # add one for 2021-03-01
  
  abline(v=c(rep(cp_date,each=5)+(0:4)),
         lty=1,col=alpha('grey',0.4),lwd=4)
  
  if(!any(is.na(scm_scenario_cp))){
    abline(v=scm_scenario_cp,
           lty=1,col=alpha('black',0.4),lwd=4)      
  }
}

## POPULATION LEVEL ----
#output_tag <- '/hosp_adm'; bool_polygon <- T; opt_cp = NA
# output_files_inc <- output_files_cases; bool_polygon <- T
output_tag <- 'VOC_mild'; bool_polygon <- T; opt_cp = NA
plot_incidence_time <- function(output_files, 
                                output_tag,
                                bool_polygon,
                                col_lib,
                                x_axis_scm_day = NA,
                                scm_callibration_date = NA,
                                scm_scenario_cp = NA,
                                be_ref_data,
                                bool_legend=TRUE){
  
  output_files_inc <- output_files[grepl(gsub('xtra','',output_tag),output_files) & grepl('rds',output_files) ]
  # output_files_inc <- output_files[grepl(output_tag,output_files) & grepl('rds',output_files) ]
  
  if(length(output_files_inc)==0){
    return(NULL)
  }
  

  bool_belgium <- any(grepl('_belgium.rds',output_files))
  
  if(any(is.na(x_axis_scm_day))){
    x_axis_scm_day <- range(0,length(be_ref_data$date))
  }
  # x-axis limits
  x_ticks      <- pretty(get_scm_start_date() + x_axis_scm_day,12)
  x_ticks[x_ticks < get_scm_start_date()]      <- get_scm_start_date()
  x_axis_label <- paste0('Time (',paste(unique(format(x_ticks,'%Y')),collapse='-'),')')
  
  y_abline <- NA
  legend_ncol <- 1
  
  # set reference
  if(any(grepl('hosp_adm',output_files_inc))){
    plot_ref_data <- data.frame(date     = be_ref_data$date,
                                reported = be_ref_data$hospital_admissions + be_ref_data$hospital_admissions_other)
    ref_label    <- 'hospital admissions with COVID19'
    ref_pch      <- 1
    y_lim        <- c(0,ifelse(min(x_axis_scm_day)==0,800,300))
    y_abline     <- NA
    
    if(bool_belgium) {y_lim[2] <- 600}
  }
  if(any(grepl('hosp_adm_xtra',output_tag))){
    plot_ref_data <- data.frame(date     = be_ref_data$date,
                                reported = be_ref_data$hospital_admissions)
    ref_label    <- 'hospital admissions for COVID19'
    ref_pch      <- 5
    y_lim        <- c(0,ifelse(min(x_axis_scm_day)==0,800,300))
    y_abline     <- NA
    
    if(bool_belgium) {y_lim[2] <- 600}
  }
  if(any(grepl('new_infect',output_files_inc))){
    plot_ref_data <- data.frame(date     = be_ref_data$date,
                                reported = be_ref_data$cases)
    ref_label    <- 'infections'
    ref_pch      <- 5 
    y_lim    <- c(0,ifelse(min(x_axis_scm_day)==0,5e4,2e4))
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
    y_lim     <- c(1,700) 
    ref_pch   <- 3
    ref_label <- 'ICU load'
    y_abline  <- NA
    
    if(bool_belgium) {y_lim[2] <- 1500}
  }
  if(any(grepl('Rt',output_files_inc))){
    plot_ref_data <- data.frame(date     = be_ref_data$date,
                                reported = be_ref_data$Rt_cases)
    y_lim     <- c(0.4,1.9) 
    ref_pch   <- 6
    ref_label <- 'cases ~ Rt'
    y_abline  <- 1
  }
  if(any(grepl('VOC_mild',output_files_inc))){
    plot_ref_data <- data.frame(date     = as.Date(db_voc_raw$date),
                                reported = db_voc_raw$voc_aggr)
    y_lim     <- c(0,1.15) 
    ref_pch   <- 1
    ref_label <- 'proportion Alpha+Beta+Gamma VOC'
    y_abline  <- NA
    legend_ncol <- 2
    col_lib$label[1] <- 'Estimated (multinomial)'
  }
  
  if(any(grepl('VOC_mild_delta',output_files_inc))){
    plot_ref_data <- data.frame(date     = as.Date(db_voc_raw$date),
                                reported = db_voc_raw$voc_delta)
    ref_label <- 'proportion Delta VOC'
  }
  
  if(!default_y_axis){
    tmp_ref_data <- plot_ref_data$reported[plot_ref_data$date %in% min(x_ticks):max(x_ticks)]
    if(length(tmp_ref_data)>0 & sum(tmp_ref_data,na.rm=T)>0){
      y_lim <- range(c(tmp_ref_data*0.9,
                       tmp_ref_data*1.1),
                     na.rm=T)
      # y_abline <- NA
    }
  }
  
  plot(x = plot_ref_data$date,
       y = plot_ref_data$reported,
       xlim = range(x_ticks),
       ylim = y_lim,
       xaxt='n',
       main= 'Updated figure from Technical Note SIMID v2021-10-12',
       xlab=x_axis_label, 
       ylab= paste('Daily',ref_label),
       col=0,
       las=2)
  axis(1,x_ticks,format(x_ticks,"%d/%m"),cex.axis=0.9)
  abline(v=x_ticks,col='grey',lty=3)
  grid(nx=NA,ny=NULL)
  if(!is.na(scm_callibration_date)){
    add_vertical_line(scm_callibration_date,bool_text =T,date_tag= 'callibration',pos_factor=1.5)
  }
  abline(h=y_abline,lty=2)
  add_change_points(scm_scenario_cp)

  i_file <- 1
  for(i_file in 1:length(output_files_inc)){
    scen_data <- readRDS(output_files_inc[i_file])
    scen_col  <- get_col(output_files_inc[i_file],col_lib)
    scm_dates  <- get_scm_start_date() +  0:(nrow(scen_data)-1)
    
    if(bool_polygon){
      hosp_aggr <-  data.frame(admin_mean=apply(scen_data,1,mean),
                               #admin_mean=NA,
                               admin_LL = apply(scen_data,1,quantile,0.025,na.rm=T),
                               admin_UL = apply(scen_data,1,quantile,0.975,na.rm=T))
      # # fix for ICU load
      # if(nrow(scen_data)==3){
      #   hosp_aggr <-  data.frame(admin_mean=scen_data[2,],
      #                            #admin_mean=NA,
      #                            admin_LL = scen_data[1,],
      #                            admin_UL = scen_data[3,])
      #   scm_dates  <- get_scm_start_date() +  0:(ncol(scen_data)-1)
      # }
      
      add_polygon(scenario_data  = hosp_aggr[-nrow(hosp_aggr),],
                  scenario_dates = scm_dates[-length(scm_dates)],
                  col = scen_col,
                  alpha_val = 0.25,
                  mean_lty = ref_pch)
    } else{
      # for(i_run in 1:ncol(scen_data)){
      #   lines(x   = scm_dates,
      #         y   = scen_data[,i_run],
      #         col = alpha(scen_col,0.1))
      # }
      if(nrow(scen_data)==3){
        scm_dates  <- get_scm_start_date() +  0:(ncol(scen_data)-1)
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
  
  add_scen_legend(col_lib, 
                  output_files_inc, 
                  ref_pch = ref_pch, 
                  ref_label=ref_label,
                  ncol=legend_ncol,
                  ref_date = max(plot_ref_data$date[!is.na(plot_ref_data$reported)]))
  
  add_copyright()
}

# help function, using many variables outside the funcion scope
get_regional_figures <- function(sel_region,file_name,output_files,col_lib,bool_multi = TRUE){
  
  # pdf(file.path(output_folder,paste0('technical_note_v20210914_fig1_belgium.pdf')),9,12)
  if(bool_multi){
    pdf(file.path(output_folder,file_name),12,12)
    par(mar=c(4.5,4.5,2,0.5),mfrow=c(2,2))
  } else{
    pdf(file.path(output_folder,file_name),9,4)
    par(mar=c(4.5,4.5,2,0.5),mfrow=c(1,1))
  }
  
  output_files_x <- output_files[grepl(sel_region,output_files)]
  be_ref_data_x  <- get_latest_incidence_data(sel_region=sel_region)
  
  sel_region_title <- paste(toupper(substring(sel_region, 1,1)), substring(sel_region, 2), sep="", collapse=" ")
  if(grepl('forCOVID',file_name)){
    plot_incidence_time(output_files_x,'/hosp_adm_xtra',col_lib = col_lib ,x_axis_scm_day=x_axis_scm_day,bool_polygon = T,scm_callibration_date=scm_callibration_date,be_ref_data=be_ref_data_x)
  } else{
    plot_incidence_time(output_files_x,'/hosp_adm',col_lib = col_lib ,x_axis_scm_day=x_axis_scm_day,bool_polygon = T,scm_callibration_date=scm_callibration_date,be_ref_data=be_ref_data_x)
  }
  legend('topright',sel_region_title,bg='white')
  plot_incidence_time(output_files_x,'/hosp_load_incr',col_lib = col_lib ,x_axis_scm_day=x_axis_scm_day,bool_polygon = T,scm_callibration_date=scm_callibration_date,be_ref_data=be_ref_data_x)
  legend('topright',sel_region_title,bg='white')
  plot_incidence_time(output_files_x,'/icu_load_incr',col_lib = col_lib ,x_axis_scm_day=x_axis_scm_day,bool_polygon = T,scm_callibration_date=scm_callibration_date,be_ref_data=be_ref_data_x)
  legend('topright',sel_region_title,bg='white')
  plot_incidence_time(output_files_x,'Rt_mild',col_lib = col_lib ,x_axis_scm_day=x_axis_scm_day,bool_polygon = T,scm_callibration_date=scm_callibration_date,be_ref_data=be_ref_data_x)
  legend('topright',sel_region_title,bg='white')
  if(!bool_multi){
    plot_incidence_time(output_files_x,'/hosp_adm_xtra',col_lib = col_lib ,x_axis_scm_day=x_axis_scm_day,bool_polygon = T,scm_callibration_date=scm_callibration_date,be_ref_data=be_ref_data_x)
    legend('topright',sel_region_title,bg='white')
  }
  dev.off()
}

## SET FILE NAMES ----

# main results
output_files_be <- output_files[!grepl('CAP',output_files) & 
                                    !grepl('0_100_0_0',output_files) &
                                    !grepl('0_0_120_0',output_files)]

output_files_region <- output_files[!grepl('CAP',output_files) & 
                                    !grepl('0_100_0_0',output_files) &
                                    #!grepl('0_0_0_140',output_files) &
                                    !grepl('0_0_120_0',output_files)]

output_files_bxl <- output_files[grepl('590',output_files) &
                                   grepl('0_0_0_0',output_files)]

output_files_cap <- output_files[grepl('0_0_0_0',output_files)]

col_lib_main <- col_lib[c(1,6:2),]
col_lib_bxl <- col_lib[c(1,10),]
col_lib_cap <- col_lib[c(1,7:9),]

## FIGURE 1: BE ----
get_regional_figures('belgium','technical_note_v20211012_fig1_belgium.pdf',output_files = output_files_be,col_lib=col_lib_main)
get_regional_figures('belgium','technical_note_v20211012_xxx_belgium.pdf',output_files = output_files_be,col_lib=col_lib_main,bool_multi = FALSE)
get_regional_figures('belgium','technical_note_v20211012_forCOVID_belgium.pdf',output_files = output_files_be,col_lib=col_lib_main)

## FIGURE 2: BXL ----
get_regional_figures('brussels','technical_note_v20211012_fig2_brussels.pdf',output_files = output_files_bxl, col_lib=col_lib_bxl)
get_regional_figures('brussels','technical_note_v20211012_xxx_brussels.pdf',output_files = output_files_bxl, col_lib=col_lib_bxl,bool_multi = FALSE)
get_regional_figures('brussels','technical_note_v20211012_forCOVID_brussels.pdf',output_files = output_files_bxl, col_lib=col_lib_bxl)

## FIGURE 3: FL ----
get_regional_figures('flanders','technical_note_v20211012_fig3_flanders.pdf',output_files = output_files_region,col_lib=col_lib_main)
get_regional_figures('flanders','technical_note_v20211012_xxx_flanders.pdf',output_files = output_files_region,col_lib=col_lib_main,bool_multi = FALSE)
get_regional_figures('flanders','technical_note_v20211012_forCOVID_flanders.pdf',output_files = output_files_region,col_lib=col_lib_main)

## FIGURE 4: WAL ----
get_regional_figures('wallonia','technical_note_v20211012_fig4_wallonia.pdf',output_files = output_files_region,col_lib=col_lib_main)
get_regional_figures('wallonia','technical_note_v20211012_xxx_wallonia.pdf',output_files = output_files_region,col_lib=col_lib_main,bool_multi = FALSE)
get_regional_figures('wallonia','technical_note_v20211012_forCOVID_wallonia.pdf',output_files = output_files_region,col_lib=col_lib_main)

# COPY FILES TO GDRIVE ----
#update gdrive figures (external)
gdrive_path <- '~/Google\ Drive/nCov2019\ -\ research/Vaccine/stochastic_model/figures/shared_external/TechNote_v20211012'
#gdrive_path <- '~/Google\ Drive/nCov2019\ -\ research/Vaccine/stochastic_model/figures/shared_external/TechNote_latest'
if(!dir.exists(gdrive_path)){ dir.create(gdrive_path) }
output_file_names <- dir(output_folder,full.names = F)
for(sel_file in output_file_names[!grepl('xxx',output_file_names)]){
  file.copy(from = file.path(output_folder,sel_file),
            to   = file.path(gdrive_path,sel_file),
            overwrite = TRUE)
  print(sel_file)
}
print('UPDATED')
print(range(be_ref_data$date))


if(0==1){
  
  # additional analysis on number of admissions / region
  
  be_ref_data <- get_latest_incidence_data(sel_region = 'belgium')
  bxl_ref_data <- get_latest_incidence_data(sel_region = 'brussels')
  fld_ref_data <- get_latest_incidence_data(sel_region = 'flanders')
  wal_ref_data <- get_latest_incidence_data(sel_region = 'wallonia')
  
  k_rm <- 14
  be_main <- rollmean(be_ref_data$hospital_admissions,k=k_rm,fill = NA,align = 'center')
  be_other <- rollmean(be_ref_data$hospital_admissions_other,k=k_rm,fill = NA,align = 'center')
  be_other_rel <- be_other / be_main
  
  bxl_main <- rollmean(bxl_ref_data$hospital_admissions,k=k_rm,fill = NA,align = 'center')
  bxl_other <- rollmean(bxl_ref_data$hospital_admissions_other,k=k_rm,fill = NA,align = 'center')
  bxl_other_rel <- bxl_other / bxl_main
  
  fld_main <- rollmean(fld_ref_data$hospital_admissions,k=k_rm,fill = NA,align = 'center')
  fld_other <- rollmean(fld_ref_data$hospital_admissions_other,k=k_rm,fill = NA,align = 'center')
  fld_other_rel <- fld_other / fld_main
  
  wal_main <- rollmean(wal_ref_data$hospital_admissions,k=k_rm,fill = NA,align = 'center')
  wal_other <- rollmean(wal_ref_data$hospital_admissions_other,k=k_rm,fill = NA,align = 'center')
  wal_other_rel <- wal_other / wal_main
  
  pdf('output/hosp_other_num.pdf',9,4)
  par(mar=c(4,4,0.1,0.1))
  plot(be_ref_data$date,
       be_other,
       xlab='2020-2021',
       ylab='COVID19 hospital admissions other pathology',
       xaxt='n',
       ylim=c(0,80))  
  grid(nx=NA,ny=NULL)
  add_date_axis(be_ref_data$date)
  
  lines(bxl_ref_data$date,
        bxl_other,
       col=3,lwd=2
  )
  lines(fld_ref_data$date,
        fld_other,
        col=4,lwd=2
  )
  lines(wal_ref_data$date,
        wal_other,
        col=6,lwd=2
  )
  legend('topleft',
         c('BEL',
           'BXL',
           'FLA',
           'WAL'),
         col=c(1,3,4,6),
         lwd=2)
  dev.off()
  
  ## 
  pdf('output/hosp_other_prop.pdf',9,4)
  par(mar=c(4,4,0.1,0.1))
  plot(be_ref_data$date,
       be_other_rel,
       xlab='2020-2021',
       ylab='PROPORTION COVID19 hospital admissions other pathology',
       xaxt='n',
       ylim=c(0,1))  
  grid(nx=NA,ny=NULL)
  add_date_axis(be_ref_data$date)
  
  lines(bxl_ref_data$date,
        bxl_other_rel,
        col=3,lwd=2
  )
  lines(fld_ref_data$date,
        fld_other_rel,
        col=4,lwd=2
  )
  lines(wal_ref_data$date,
        wal_other_rel,
        col=6,lwd=2
  )
  legend('topleft',
         c('BEL',
           'BXL',
           'FLA',
           'WAL'),
         col=c(1,3,4,6),
         lwd=1)
  dev.off()
  
  
}

