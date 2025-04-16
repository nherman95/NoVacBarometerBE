########################################################################### #
# This file is part of the Stochastic Compartmental Model for SARS-COV-2 
# transmission in Belgium, conceived by members of SIMID group during the 
# COVID19 pandemic.
#
# This file can be used to generate output for the
# European covid-19 Scenario Hub (ECDC)
#
# This file can be sourced with command line arguments. For example:
# - terminal: Rscript R/output_EUhub.R belgium 
# - terminal: Rscript R/output_EUhub.R belgium &
# - within R: source('R/output_EUhub.R')
# - within R: system('Rscript R/output_EUhub.R belgium &')
#
# Copyright 2022, SIMID, University of Antwerp & Hasselt University                                        
########################################################################### #

args = commandArgs(trailingOnly=TRUE)

rm(list=ls()[ls()!='args'])

suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(RColorBrewer))

# load functions and data
source('R/main_vaccination.R')

################################################################ #
## SETTINGS AND PREAMBLE ----
################################################################ #

# SIMULATION DATA
output_dir <- "output/fall23"


# for Omicron
high <- ''
alter <- 0
vaccine <- '' #'ubooster90_jan01'
date <- '0301'

sel_region <- "belgium"

# set x-axis limits and date threshold
origine_date_num <- 1
x_axis_scm_day <- c(0,519)

#other possibilities
#x_axis_scm_day <- c(700,915)
x_axis_scm_day <- c(sim_date2day('2023-02-28'),sim_date2day('2024-01-10'))

# use default y axis?
default_y_axis <- TRUE

# load output file
output_files <- dir(output_dir,pattern=sel_region,full.names = T,recursive = T)

# check output files
if(length(output_files)==0){
  stop('OUTPUT SELECTION ISSUE IN "',output_dir,'" WITH "',sel_region,'"')
}


# additional figure legend library 
if(any(grepl('fall23',output_files))){
  col_lib <- data.frame(tags = c('xxxx',
                                 '0_0_0_0_100_.*_0flu_biv',
                                 '0_0_0_0_100_.*_50flu_biv',
                                 '0_0_0_0_100_.*_100flu_biv',
                                 '0_0_0_0_100_.*_100flu_plus15_biv'),
                        col  = c('black','blue','darkorchid','orange','green'),
                        label = c(paste0('Reported (',max(be_ref_data$date),')'),
                                  'Simulated contact paterns and seasonality - no vaccine campaign',
                                  'Vaccine campaign 65+ 50% of flu, XBB.1.5 bivariant booster',                                  
                                  'Vaccine campaign 65+ 100% of flu, XBB.1.5 bivariant booster',
                                  'Vaccine campaign 65+ 100% of flu + 15%, XBB.1.5 bivariant booster'),
                        lwd = c(NA,2,2,2,2),
                        pch = c(1,NA,NA,NA,NA))
  
  
  color_tags    <- get_color_multiple(output_files,col_lib = col_lib,show_warnings = FALSE)
  table(color_tags)
  color_missing <- get_color_multiple('aa',col_lib = col_lib,show_warnings = FALSE)
  output_files <- output_files[color_tags != color_missing]
  if(length(output_files) == 0){
    warning("No output files left after scenario selection")
  }
}



## MODEL PARAMETERS ----

# load model parameters and load wave and scenario info
param_file_names         <- dir(output_dir,pattern = paste0('MCMCmulti.*',sel_region),full.names = T,recursive = T)
scm_num_waves            <- vector(length = length(param_file_names))
scm_callibration_day     <- vector(length = length(param_file_names))
scm_region               <- vector(length = length(param_file_names))
scm_scenario_cp_day      <- c()

for(i_file  in 1:length(param_file_names)){
  parms                         <- read.csv(param_file_names[i_file],header=T)
  scm_num_waves[i_file]         <- identify_total_nb_waves(names(parms))
  scm_callibration_day[i_file]  <- unique(parms$ndays_calibration)
  scm_region[i_file]            <- unique(parms$region_id)
  scm_scenario_cp_day           <- c(scm_scenario_cp_day,unique(parms[grepl('cnt_adjust_day',names(parms))]))
}


#workaround to extract only some change points
#if(exists("date")){
#if(date=='0214'){  scm_scenario_cp_day <- scm_scenario_cp_day["cnt_adjust_day4"]  }     
#if(date=='0301'){   scm_scenario_cp_day <- scm_scenario_cp_day["cnt_adjust_day5"]  }     
#if(date=='0501'){   scm_scenario_cp_day <- scm_scenario_cp_day["cnt_adjust_day3"]  }  
#}



# number of CoMix and additional waves
scm_num_cp <- unique(scm_num_waves)

# set calibration and scenario dates
scm_callibration_date    <- unique(sim_day2date(scm_callibration_day))
#scm_callibration_date    <- '2022-07-23'  #temporary for correction
#scm_callibration_date    <- unique(sim_day2date(origine_date_num))
scm_scenario_cp          <- unique(sim_day2date(unlist(scm_scenario_cp_day)))

# set region from paramter file(s)
sel_region <- unique(get_region(scm_region))

# check parameter input on region, and stop if not uniform
if(length(sel_region)==0){
  stop('REGION PARAMETER IS MISSING IN THE PARAMETER FILE(S)')
} else if(length(sel_region)>1){
  stop(paste0('REGION PARAMETER IS NOT UNIQUE: "',paste0(sel_region,collapse = '", "'),'"'))
}


# sort and select output files
output_files <- sort(output_files,decreasing = T)
output_files <- output_files[grepl(sel_region,output_files)]

# get name for figures
output_dir_base <- unlist(strsplit(output_dir,'/'))
output_dir_base <- output_dir_base[length(output_dir_base)-1]

# set output
output_folder <- file.path(output_dir,paste(basename(output_dir),sel_region,sep='_'))
if(!dir.exists(output_folder)) {dir.create(output_folder)}
print(output_dir)

# load reference data
be_ref_data <- get_latest_incidence_data(sel_region = sel_region)
#be_ref_data <- get_latest_incidence_data(sel_region = sel_region,enable_data_reuse = F)
#be_ref_data <- get_observed_incidence_data(enable_data_reuse = F) # do not re-use stored data files

## PDF STREAM ----

pdf(file.path(output_folder,paste0('incidence_',basename(output_dir),
                                   '_',output_dir_base,
                                   '_',sel_region,
                                   ifelse(exists("high"),paste0('_',high),''),
                                   ifelse(exists("date"),paste0('_',date),''),
                                   ifelse(alter,paste0('_ALTER'),''),
                                   ifelse(exists("vaccine"),paste0('_',vaccine),''),
                                   #ifelse(min(x_axis_scm_day)==0,'_full',''),
                                   '.pdf')),9,6)
par(mar=c(4.5,5,1,0.5))
x_axis_vaccine = c(300,max(x_axis_scm_day))

# POPULATION LEVEL (lines and polyon) ----
if(1==1){ 
  multi_plot_incidence_time(output_files,
                            col_lib = col_lib,
                            bool_polygon_multi = T,
                            x_axis_scm_day = x_axis_scm_day,
                            be_ref_data = be_ref_data,
                            default_y_axis = default_y_axis,
                            scm_callibration_date = scm_callibration_date,
                            scm_scenario_cp = scm_scenario_cp,
                            scm_num_cp = scm_num_cp)
  
  #multi_plot_incidence_age_time(output_files = output_files,
  #                            col_lib = col_lib,
  #x_axis_scm_day = x_axis_scm_day,
  #                           x_axis_vaccine = x_axis_vaccine)
  
  
}


dev.off()


