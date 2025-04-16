########################################################################### #
# This file is part of the Stochastic Compartmental Model for SARS-COV-2 
# transmission in Belgium, conceived by members of SIMID group during the 
# COVID19 pandemic.
#
# This file is used to update the analysis presented at the Press Conference
# of the Prime Minister De Croo on March 22, 2021.
#
# Copyright 2021, SIMID, University of Antwerp & Hasselt University                                        
########################################################################### #

# Note: this file contains a local copy of the plot functions to make sure that
# the figures do not change when the project code evolves over time. This is 
# not the most elegant way, but it is robust.

# clear workspace
rm(list=ls())

library(scales)
library(data.table)
source('R/lib_stochastic_model.R')

# SETTINGS          ----
###################### # 

# reference data
# ref_data_be <- get_observed_incidence_data()
ref_data_be <- get_latest_incidence_data()
range(ref_data_be$date)
#tail(ref_data_be[,1:2])

# simulation data
# sim_data            <- read.table("UHasselt_predictions_v7.csv",sep=',',header=T)
sim_data              <- read.table("data/results_technical_notes/UHasselt_predictions_v7_np_full.csv",sep=',',header=T)
sim_data_callibration <- "2021-02-10"

# output folder
output_folder <- 'output/figures_update_v20210226'

# image file format
sel_output_format <- 'pdf'  # options: pdf, png, NA (=Rstudio)


# PRE-PROCESSING    ----
###################### # 

names(sim_data)[1]       <- tolower(names(sim_data)[1])
sim_data_callibration    <- as.Date(sim_data_callibration)
sim_data$date            <- as.Date(sim_data$date)
ref_data_be$date         <- as.Date(ref_data_be$date)

# if output folder does not exit yet, create one
if(!dir.exists(output_folder)){
  dir.create(output_folder,recursive = T)
}

run_tag <- 'S1' # to debug
plot_model_output <- function(sim_data,run_tag, sel_output_format=NA){

  run_tags_abc <- paste0(run_tag,c('a','b','c'),'_incidence')
  sim_exit_month <- c('','March','April','May','','March','April','May')[as.numeric(gsub('S','',run_tag))]
  sim_letter <- c('A','B','C','D','A(70)','B(70)','C(70)','D(70)')[as.numeric(gsub('S','',run_tag))]
 
  scenarioXa     <- sim_data[,grepl(run_tags_abc[1],names(sim_data))]
  scenarioXb     <- sim_data[,grepl(run_tags_abc[2],names(sim_data))]
  scenarioXc     <- sim_data[,grepl(run_tags_abc[3],names(sim_data))]
  scenario_dates <- sim_data$date
  
  if(!is.na(sel_output_format) & sel_output_format == 'pdf'){
    pdf(paste0(output_folder,"/stoch_model_hosp_20210212_",sim_letter,'_incidence.pdf'),8,4)
  }
  if(!is.na(sel_output_format) & sel_output_format == 'png'){
    png(paste0(output_folder,"/stoch_model_hosp_20210212_",sim_letter,'_incidence.png'),width = 600, height = 300)
  }
  
  # set margins
  par(mar=c(5,5,1,1))

  # set plot limits  
  y_lim    <- range(0,800,scenarioXa)
  x_lim    <- range(sim_data$date)
  x_ticks  <- pretty(x_lim,n=12)
  x_labels <- format(x_ticks,'%d/%m')
  x_lab    <- paste0('Time (',paste(unique(format(x_ticks,'%Y')),collapse='-'),')')
  
  plot(x = ref_data_be$date,
       y = ref_data_be$hospital_admissions,
       xlim = x_lim,
       ylim = y_lim,
       axes = F,
       xlab = x_lab,
       ylab = "Number of new hospitalizations",
       pch=1)
  axis(1, labels = x_labels, at = x_ticks, cex.axis = 0.8);
  axis(2, cex.axis = 0.8,las=2);
  abline(v = x_ticks,lty=3,col='grey')
  abline(h = c(1)*100,lty=3,col='grey')
  grid(nx=NA,ny=NULL)
  add_copyright()
  
  # add callibration
  add_vertical_line(sim_data_callibration,bool_text = T,date_tag='callibration')
  
  if(nchar(sim_exit_month)==0){
    legend_text <- c(paste0("Reported (until ",max(ref_data_be$date),")"),
                     paste0(sim_letter,": no changes + 70% increase VOC"),
                     paste0(sim_letter,": no changes + 50% increase VOC"),
                     paste0(sim_letter,": no changes + 30% increase VOC"))
  } else{
    legend_text <- c("Reported",
                     paste0(sim_letter,": relaxations on 1 ",sim_exit_month," + 70% increase VOC"),
                     paste0(sim_letter,": relaxations on 1 ",sim_exit_month," + 50% increase VOC"),
                     paste0(sim_letter,": relaxations on 1 ",sim_exit_month," + 30% increase VOC"))
  }
  
    legend(x = min(x_lim), y = max(y_lim)*1.05,
           legend_text,
           col = c("black","red","blue","orange"), 
           box.lwd = 0,box.col = "white",
           lty = c(NA,1,1,1),
           pch = c(1,NA,NA,NA),
           cex = 0.8, lwd = 2,
           # title = legend_title,
           bg='white')
  
  text(max(x_lim), -550, 
       "(Abrams et al. - SIMID - UHasselt & UAntwerpen)", 
       cex = 0.5,adj=1,
       xpd=T)
  
  
  add_polygon(scenarioXa, scenario_dates, col = "orange", alpha_val = 0.2)
  add_polygon(scenarioXc, scenario_dates, col = "red",    alpha_val = 0.2)
  add_polygon(scenarioXb, scenario_dates, col = "blue",   alpha_val = 0.2)
  
  if(!is.na(sel_output_format)){
    dev.off()  
  }
}



## Run functions ----

# rename output
# scenario 1-4 => S9_1, S9_2, etc
# scenario 5-8 => A, B, C, D

for(i_tag in paste0('S',c(1:8))){
  plot_model_output(sim_data = sim_data,
                    run_tag  = i_tag,
                    sel_output_format = sel_output_format)
}



## Update persconference figures ----

# pdf_persconf_filename <- paste0('persconf_stoch_model_hosp_A_',format(Sys.Date(),'%Y%m%d'),'.pdf')
pdf_persconf_filename      <- 'persconf_stoch_model_hosp_A.pdf'
pdf_persconf_filename_zoom <- 'persconf_stoch_model_hosp_A_zoom.pdf'
pdf(file.path(output_folder,pdf_persconf_filename),8,4)
plot_model_output(sim_data = sim_data,
                  run_tag  = 'S1')
dev.off()

pdf(file.path(output_folder,pdf_persconf_filename_zoom),8,4)
plot_model_output(sim_data = sim_data[sim_data$date >= as.Date('2021-01-01'),],
                  run_tag  = 'S1')
dev.off()

# update gdrive figures (external)
gdrive_path_external <- '/Users/lwillem/Google\ Drive/nCov2019\ -\ research/Vaccine/stochastic_model/figures/shared_external/PressConf_20210226'
file.copy(from = file.path(output_folder,pdf_persconf_filename),
          to   = file.path(gdrive_path_external,pdf_persconf_filename),
          overwrite = TRUE
)
file.copy(from = file.path(output_folder,pdf_persconf_filename_zoom),
          to   = file.path(gdrive_path_external,pdf_persconf_filename_zoom),
          overwrite = TRUE
)



## VACCINES ----
# planned_num_doses      <- c(0,471855, 	 857644, 	 1376805, 	 1823393, 	 3684052, 	 5432863 ) 
# planned_first_dose     <- round(planned_num_doses/2)
# planned_first_dose_cum <- cumsum(planned_first_dose)
# planned_dates          <- as.Date(paste0('2021-',1:7,'-01'))
# 
# restore_doses_day <- c(31765, 45897, 128499, 78358)/2
# restore_dates     <- as.Date(c('2021-01-01','2021-02-01','2021-03-01','2021-05-01'))

## SCIENSANO VACCINE UPTAKE
vaccin_data        <- read.table('data/downloads/COVID19BE_VACC.csv',sep=',',header = T)
vaccin_data$DATE   <- as.Date(vaccin_data$DATE)
flag_date          <- vaccin_data$DATE > as.Date("2021-01-01") & vaccin_data$DATE < as.Date("2021-02-11")
vaccin_data        <- vaccin_data[flag_date,]

# aggregate by date
vaccin_aggr <- aggregate(COUNT ~ DATE + DOSE, data = vaccin_data,sum)
vaccin_aggr <- vaccin_aggr[vaccin_aggr$DOSE == 'A',]

########## NEW UPDATE PLOT
plan_first <- cumsum(c(0,316661, 541973.7477,	579756.0323,	1125093.218,	2456540.672,	0,	2508500.921,	340841.7406))
plan_second <- cumsum(c(0,0,	316661,	178152.3548,	394401.1935,	633309.0882,	2226300,	855605.5222,	2267958.259))
plan_single <- c(0,0,	0,	0,	0,	0,	1402500,	1020000,	1020000)
plan_date <- as.Date(paste0('2021-',1:9,'-01'))

final_dates <- seq(as.Date('2021-01-01'),as.Date('2021-07-1'),1)

if(sel_output_format == 'pdf'){
  pdf(paste0(output_folder,"/stoch_model_hosp_20210212_uptake.pdf"),8,4)
}
if(sel_output_format == 'png'){
  png(paste0(output_folder,"/stoch_model_hosp_20210212_uptake.png"),width = 600, height = 300)
}

plot(final_dates,
     approx(x = plan_date,
            y = plan_first,
            xout = final_dates,
            method = "linear",
            rule = 2)$y,
     ylim=c(0,9e6),
     ylab='Vaccine uptake',
     xlab='Time',
     col = 1,
     type='l',
     lwd=4
)

lines(final_dates,
      approx(x = plan_date,
             y = plan_single + plan_second,
             xout = final_dates,
             method = "linear",
             rule = 2)$y,
      col=1,
      lty=2,
      lwd=4
)


restore_doses_day <- c(31765, 45897, 128499, 78358)
restore_dates     <- as.Date(c('2021-01-01','2021-02-01','2021-03-01','2021-05-01'))
restore_doses     <- (approx(x = restore_dates,
                             y = restore_doses_day,
                             xout = final_dates,
                             method = "constant",
                             rule = 2)$y
)

restore_first <- restore_doses / 2
restore_first[final_dates %in% vaccin_aggr$DATE] <- vaccin_aggr$COUNT

lines(final_dates,
      cumsum(restore_first),
      col='orange',
      lwd=4)

lines(final_dates+20,
      cumsum(restore_first),
       col='orange',
       lwd=4,
       lty=2)
grid(nx=NA,ny=NULL)
abline(v=plan_date,col='grey',lty=3)

legend('topleft',
       c('model: first dose',
         'model: second dose',
         'first dose',
         'full scheme'),
       col = c('orange','orange',1,1,1),
       lwd=2,
       lty=c(1,2,1,2,3))

# legend('topleft',
#        c('model: first dose',
#          'model: second dose'),
#        col = c('orange','orange',1,1,1),
#        lwd=2,
#        lty=c(1,2,1,2,3))

if(!is.na(sel_output_format)){
  dev.off()  
}


