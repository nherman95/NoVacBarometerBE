########################################################################### #
# This file is part of the Stochastic Compartmental Model for SARS-COV-2 
# transmission in Belgium, conceived by members of SIMID group during the 
# COVID19 pandemic.
#
# This file can be used to (re)create simulation output figures from one 
# or more runs. 
#
# This file can be sourced with command line arguments. For example:
# - terminal: Rscript R/plot_scm_combined.R belgium 
# - terminal: Rscript R/plot_scm_combined.R belgium &
# - within R: source('R/plot_scm_combined.R')
# - within R: system('Rscript R/plot_scm_combined.R belgium &')
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
output_dir <- "output/scm24_0707test/"
output_dir <- "output/scm24_0708/"

# temporary
sel_region <- "belgium"

# set x-axis limits
x_axis_scm_day <- c(700,915)
<<<<<<< HEAD
x_axis_scm_day <- c(700,sim_date2day('2023-01-01'))

if(!bool_combined) x_axis_scm_day <- c(740,850)

# if(grepl('0_0_100_100_0',output_dir) | grepl('114_110',output_dir) | grepl('0_0_0_0',output_dir)){
#   x_axis_scm_day <- c(550,650) 
#   x_axis_scm_day <- c(550,680) 
# }
#x_axis_scm_day <- c(0,700) # optional: full period
=======
>>>>>>> main

# use default y axis?
default_y_axis <- TRUE

# option to use command line arguments: use latest output?
if(length(args)>=2){
  sel_latest_output <- as.logical(args[2])
  if(!is.na(sel_latest_output) && sel_latest_output){
    output_folders       <- dir('output/',pattern = 'scm23',full.names = T,recursive = F)
    output_folders_mtime <- file.info(output_folders)$mtime
    output_dir           <- output_folders[order(output_folders_mtime,decreasing = T)[1]]
    print(paste("USE LATEST OUTPUT:",output_dir))
  }
}

# option to specify the region via command line arguments
if(length(args)>=1){
  sel_region   <- args[1]
} 
if(!exists('sel_region')){
  sel_region <- ''
}

# load output file
output_files <- dir(output_dir,pattern=sel_region,full.names = T,recursive = T)

# check output files
if(length(output_files)==0){
  stop('OUTPUT SELECTION ISSUE IN "',output_dir,'" WITH "',sel_region,'"')
}

# set figure legend and color libray (default = reported)
col_lib <- data.frame(tags = c('xxxx'),
                      col  = c('black'),
                      label = c(paste0('Reported (',max(be_ref_data$date),')')),
                      lwd = c(NA),
                      pch = c(1))
bool_file_sort_decreasing <- T

# add directory names
scen_dir_names <- unique(basename(dirname(dir(output_dir,recursive = T,full.names = T,pattern='.csv'))))
col_lib <- rbind(col_lib[1,],
                data.frame(tags=scen_dir_names,
                           col=topo.colors(n=length(scen_dir_names)),
                           label = scen_dir_names,
                           lwd = 2,
                           pch = NA))

# # set custom figure legend and color libray
# col_lib <- data.frame(tags = c('xxxx',
#                                "orig_915_cnt0_0_0_0_0_n.._belgium",
#                                "orig_915_cnt0_0_175_0_0_n.._belgium",
#                                "orig_915_cnt0_0_200_0_0_n.._belgium",
#                                "100_100_0_0_OPTI_vaccine_uptake_2ndbooster_vsceneu_hub_slow"),
#                       col  = c('black','darkred','orange3','blue3','darkgreen'),
#                       label = c(paste0('Reported (',max(be_ref_data$date),')'),
#                                 'Estimated dynamics, no introductions and reported vaccine uptake',
#                                 '1.75x transmission dynamics from time point A',
#                                 '2.00x transmission dynamics from time point A',
#                                 'EU Covid-19 Scenario Hub (Round 1, scenario A)'),
#                       lwd = c(NA,2,2,2,2),
#                       pch = c(1,NA,NA,NA,NA))
# # option to remove combinations that are not present in the table above
# color_indices <- get_color_index_multiple(output_files,col_lib = col_lib,show_warnings = FALSE)
# output_files  <- output_files[!is.na(color_indices)]
# col_lib <- col_lib[sort(c(1,unique(color_indices[!is.na(color_indices)]))),]

<<<<<<< HEAD
if(grepl('cnt0_aaa',output_dir)){
  
  # set figure legend and color libray
  col_lib <- data.frame(tags = c('xxxx',
                                 'cnt0_0_0_0.*ubooster60_may15',
                                 'cnt0_0_110_140_170.*ubooster60_may15',
                                 'ubooster100_may15',
                                 'x_ubooster100_apr10',
                                 'ubooster50_apr10'),
                        col  = c('black','darkred','orange3','blue3','red2','darkblue'),
                        label = c(paste0('Reported (',max(be_ref_data$date),')'),
                                  #'Estimated behaviour, no introductions and reported vaccine uptake',
                                  'Estimated dynamics, no introductions and reported vaccine uptake',
                                  'Reported vaccine uptake',
                                  'Increase booster dose uptake to 100% of 2nd dose by May 15th, 2022',
                                  'Adjusted historical booster dose uptake to 100% of 2nd dose',
                                  'Adjusted historical booster dose uptake by 50%'),
                        lwd = c(NA,2,2,2,2,2),
                        pch = c(1,NA,NA,NA,NA,NA))
  bool_file_sort_decreasing <- F
}

if(grepl('140_170',output_dir)){
  col_lib <- col_lib[-6,]
  x_axis_scm_day[1] <- 720
}
#if(grepl('0_0_0_0_0',output_dir)){
 # col_lib <- col_lib[-4,]
#}

# if(!bool_combined){
#   col_lib <- col_lib[-c(4,5),]
# } 
#color_indices <- get_color_index_multiple(output_files,col_lib = col_lib,show_warnings = FALSE)
#output_files  <- output_files[!is.na(color_indices)]

# additional figure legend library (cfr. TechNote v20211116)
if(any(grepl('0_110_110',dir(output_dir)))){

  col_lib <- data.frame(tags = c('xxxx',
                                 'cnt0_0_0_0',
                                 'cnt0_0_110_110',
                                 'cnt0_0_114_110',
                                 'cnt0_0_90_100',
                                 'cnt0_0_125_115'),
                        col  = c('black','darkred','darkgreen','darkblue','orange3','gold'),
                        label = c(paste0('Reported (',max(be_ref_data$date),')'),
                                  'Rt(infections) is continuously decreasing',
                                  'Rt(infections) remains stable until November 21th, 2021',
                                  'Rt(infections) remains stable until December 2021',
                                  'Rt(infections) decreases fast after November 21th, 2021',
                                  'Rt(infections) increases until December 2021'
                                  ),
                        lwd = c(NA,2,2,2,2,2),
                        pch = c(1,NA,NA,NA,NA,NA))
  
  # re-arrange
  col_lib <- col_lib[c(1,6,4,3,2,5),]
  
}

# additional figure legend library (cfr. TechNote v20211202)
if(any(grepl('ua80',dir(output_dir)))){
  
  col_lib <- data.frame(tags = c('xxxx',
                                 '_0_vaccine.*ua80_uc80_ad6.4_child0_belgium',
                                 '_0_vaccine.*ua80_uc80_ad8_child4_summer_belgium',
                                 '1.0_vaccine.*ua80_ua80_ad6.4_child0_belgium',
                                 '1.0_vaccine.*ua80_ua80_ad8_child4_belgium',
                                 '1.0_vaccine.*ua90_ua80_ad6.6_child0_belgium',
                                 '1.0_vaccine.*ua90_ua80_ad8.2_child4_belgium',
                                 '1.0_vaccine.*ua95_ua80_ad6.7_child0_belgium',
                                 '1.0_vaccine.*ua95_ua80_ad8.3_child4_belgium'),
                        col  = c('black','darkblue','darkgreen','darkblue','darkgreen','red2','orange3','darkred','darkgrey'),
                        label = c(paste0('Reported (',max(be_ref_data$date),')'),
                                  '(E) No increased vaccine uptake for 5-11y',
                                  '(F) With vaccine coverage of 80% for 5-11y',
                                  '(A) Additional 2022 wave, with current vaccine uptake (incl. booster)',
                                  '(B) Additional 2022 wave, + 80% vaccine coverage for 5-11y',
                                  '(C) Additional 2022 wave, with >90% vaccine coverage for +18y',
                                  '(B+C) Additional 2022 wave, with >90% vaccine coverage for +18y and 80% coverage for 5-17y',
                                  '(D) Additional 2022 wave, with 95% vaccine coverage for +18y',
                                  '(B+D) Additional 2022 wave, with 95% vaccine coverage for +18y and 80% coverage for 5-17y'
                        ),
                        lwd = c(NA,2,2,2,2,2,2,2,2),
                        pch = c(1,NA,NA,NA,NA,NA,NA,NA,NA))
 
  if(!bool_combined){
    col_lib <- col_lib[-c(7,9),]
  } else{
    col_lib <- col_lib[-c(5,6,8),]
  }
  
  
  # remove combinations that are not present in the table above
  color_indices <- get_color_index_multiple(output_files,col_lib = col_lib,show_warnings = FALSE)
  output_files  <- output_files[!is.na(color_indices)]
  col_lib <- col_lib[sort(c(1,unique(color_indices[!is.na(color_indices)]))),]
   
  if(grepl('_120_n',output_dir)){
    col_lib$label <- gsub('wave,','wave (120%),',col_lib$label)
  }
  
  if(grepl('_1.0_n',output_dir)){
    col_lib <- col_lib[-1,]
  }
}

# additional figure legend library (cfr. TechNote v20211208)
if(any(grepl('100_100_0_0_aaa',dir(output_dir)))){

  col_lib <- data.frame(tags = c('xxxx',
                                 'cnt0_0_0_0_0',
                                 'cnt100_100_0_0_0',
                                 'cnt100_100_100_0_0',
                                 'cnt100_100_100_100_0',
                                 'cnt100_100_100_100_100'),
                        col  = c('black','darkblue','darkgreen','gold','orange3','darkred'),
                        label = c(paste0('Reported (',max(be_ref_data$date),')'),
                                  '(D) No increased vaccine uptake for 5-12y',
                                  '(E) With vaccine coverage of 80% for 5-12y',
                                  '(A) Additional 2022 wave, without extra vaccine uptake for 5-12y',
                                  '(B) Additional 2022 wave, with 80% vaccine coverage for 5-12y',
                                  '(C) Additional 2022 wave, with 100% vaccine coverage for +20y'
                        ),
                        lwd = c(NA,2,2,2,2,2),
                        pch = c(1,NA,NA,NA,NA,NA))

  # TODO: fix labels!!
  col_lib$label[-1] <- col_lib$tags[-1]


  # # remove 0_0_125_115..._winter (if present)
  # # remove 350_0_125_115..._summer (if present)
  # color_tags    <- get_color_multiple(output_files,col_lib = col_lib,show_warnings = FALSE)
  # table(color_tags)
  # color_missing <- get_color_multiple('aa',col_lib = col_lib,show_warnings = FALSE)
  # output_files <- output_files[color_tags != color_missing]
}

# additional figure legend library (cfr. counterfactuals)
if(any(grepl('CAP',output_files))){
col_lib <- data.frame(tags = c('xxxx',
                               'cnt0_0_0_0',
                               'cnt0_0_0_0.*CAP',
                               'cnt0_0_0_120.*CAP',
                               'oct65_',
                               'sept65_'),
                      col  = c('black','darkred','darkgreen','orange3','darkblue','gold'),
                      label = c(paste0('Reported (',max(be_ref_data$date),')'),
                                'Based on reported uptake and social contact behaviour',
                                'No additional vaccine uptake after July 1st',
                                'No additional vaccine uptake after July 1st & increased risk behaviour (120%) from September 1st',
                                'Total uptake 65% by November 1st',
                                'Total uptake 65% by October 1st'),
                      lwd = c(NA,2,2,2,2,2),
                      pch = c(1,NA,NA,NA,NA,NA))
}

# additional figure legend library (cfr. counterfactuals)
if(any(grepl('omicron_XX',output_files))){
# # Lander (temp)
  # col_lib <- data.frame(tags = c('xxxx',
  #                                'cnt0_0_0_0_0_n',
  #                                paste0('cnt0_0_0_0.*omicronSev050.*HIGH.*',vaccine),
  #                                paste0('cnt0_0_0_0.*omicronSev025.*HIGH.*',vaccine),
  #                                paste0('cnt0_0_0_0.*omicronSev050.*LOW.*',vaccine),
  #                                paste0('cnt0_0_0_0.*omicronSev025.*LOW.*',vaccine)),
  #                       col  = c('black','darkgreen','darkred','orange3','gold','darkblue'),
  #                       label = c(paste0('Reported (',max(be_ref_data$date),')'),
  #                                 'Based on reported uptake and social contact behaviour',
  #                                 'VE HIGH and Omicron severity = 50% of Delta',
  #                                 'VE HIGH and Omicron severity = 25% of Delta',
  #                                 'VE LOW and Omicron severity = 50% of Delta',
  #                                 'VE LOW and Omicron severity = 25% of Delta'),
  #                       lwd = c(NA,2,2,2,2,2),
  #                       pch = c(1,NA,NA,NA,NA,NA))
  
  # col_lib <- data.frame(tags = c('xxxx',
  #                                'cnt0_0_0',
  #                                'cnt200_0_0',
  #                                'cnt0_80_0',
  #                                'cnt200_80_0',
  #                                'cnt0_0_4'),
  #                       col  = c('black','darkgreen','darkred','orange3','gold','darkblue'),
  #                       label = c(paste0('Reported (',max(be_ref_data$date),')'),
  #                                 'VE HIGH and Omicron severity = 25% of Delta',
  #                                 '0-19y: 200% from Jan 8th',
  #                                 '20-69y: 80% from Jan 8th',
  #                                 '0-19y: 200% and 20-69y: 80% from Jan 8th',
  #                                 'Avg Oct-Nov behaviour from Jan 8th'),
  #                       lwd = c(NA,2,2,2,2,2),
  #                       pch = c(1,NA,NA,NA,NA,NA))
  
  # col_lib <- data.frame(tags = c('xxxx',
  #                                'cnt0_0_0_0_0_n',
  #                                paste0('cnt0_0_0_0_300.*omicronSev025.*LOW*',vaccine),
  #                                paste0('cnt0_0_0_0_0.*omicronSev025.*LOW*',vaccine),
  #                                paste0('cnt0_0_0_0_300.*omicronSev025.*HIGH*',vaccine),
  #                                paste0('cnt0_0_0_0_0.*omicronSev025.*HIGH*',vaccine),
  #                                paste0('cnt0_0_0_0_200.*omicronSev025.*LOW*',vaccine),
  #                                paste0('cnt0_0_0_0_150.*omicronSev025.*LOW*',vaccine),
  #                                paste0('cnt0_0_0_0_200.*omicronSev025.*HIGH*',vaccine),
  #                                paste0('cnt0_0_0_0_150.*omicronSev025.*HIGH*',vaccine)),
  #                       col  = c('black','darkgreen','goldenrod4','darkorchid4','red4','seagreen4','goldenrod3','darkorchid2','red2','seagreen3'),
  #                       label = c(paste0('Reported (',max(be_ref_data$date),')'),
  #                                 'Based on reported uptake and social contact behaviour',
  #                                 'VE LOW & Omicron severity = 25% of Delta - 200% increase',
  #                                 'VE LOW & Omicron severity = 25% of Delta - current dynamics',
  #                                 'VE HIGH & Omicron severity = 25% of Delta - 200% increase',
  #                                 'VE HIGH & Omicron severity = 25% of Delta - current dynamics',
  #                                 'VE LOW & Omicron severity = 25% of Delta - 100% increase',
  #                                 'VE LOW & Omicron severity = 25% of Delta - 50% increase',
  #                                 'VE HIGH & Omicron severity = 25% of Delta - 100% increase',
  #                                 'VE HIGH & Omicron severity = 25% of Delta - 50% increase'),
  #                       lwd = c(NA,2,2,2,2,2,2,2,2,2),
  #                       pch = c(1,NA,NA,NA,NA,NA,NA,NA,NA,NA))
  
  
  col_lib <- data.frame(tags = c('xxxx',
                                 'cnt0_0_0_0_0_n',
                                 paste0('cnt0_0_0_0_0.*omicronSev025.*LOW*',vaccine),
                                 paste0('cnt0_0_0_0_0.*omicronSev025.*HIGH*',vaccine),
                                 paste0('cnt0_0_0_200_0.*omicronSev025.*LOW*',vaccine),
                                 paste0('cnt0_0_0_200_0.*omicronSev025.*HIGH*',vaccine),
                                 paste0('cnt0_0_0_300_0.*omicronSev025.*LOW*',vaccine),
                                 paste0('cnt0_0_0_300_0.*omicronSev025.*HIGH*',vaccine),
                                 paste0('cnt0_0_0_0_200.*omicronSev025.*LOW*',vaccine),
                                 paste0('cnt0_0_0_0_200.*omicronSev025.*HIGH*',vaccine),
                                 paste0('cnt0_0_0_0_300.*omicronSev025.*LOW*',vaccine),
                                 paste0('cnt0_0_0_0_300.*omicronSev025.*HIGH*',vaccine),
                                 paste0('cnt0_0_200_0_0.*omicronSev025.*LOW*',vaccine),
                                 paste0('cnt0_0_200_0_0.*omicronSev025.*HIGH*',vaccine),
                                 paste0('cnt0_0_300_0_0.*omicronSev025.*LOW*',vaccine),
                                 paste0('cnt0_0_300_0_0.*omicronSev025.*HIGH*',vaccine)),
                        col  = c('black','gray','blue','blue3','green','darkgreen','orange','darkorange','green','darkgreen','orange','darkorange','green','darkgreen','orange','darkorange'),
                        label = c(paste0('Reported (',max(be_ref_data$date),')'),
                                  'Based on reported uptake and social contact behaviour',
                                  'VE LOW & Omicron severity = 25% of Delta - current dynamics',
                                  'VE HIGH & Omicron severity = 25% of Delta - current dynamics',
                                  'VE LOW & Omicron severity = 25% of Delta - 100% increase on February 14',
                                  'VE HIGH & Omicron severity = 25% of Delta - 100% increase on February 14',
                                  'VE LOW & Omicron severity = 25% of Delta - 200% increase on February 14',
                                  'VE HIGH & Omicron severity = 25% of Delta - 200% increase on February 14',
                                  'VE LOW & Omicron severity = 25% of Delta - 100% increase on March 1',
                                  'VE HIGH & Omicron severity = 25% of Delta - 100% increase on March 1',
                                  'VE LOW & Omicron severity = 25% of Delta - 200% increase on March 1',
                                  'VE HIGH & Omicron severity = 25% of Delta - 200% increase on March 1',
                                  'VE LOW & Omicron severity = 25% of Delta - 100% increase on May 1',
                                  'VE HIGH & Omicron severity = 25% of Delta - 100% increase on May 1',
                                  'VE LOW & Omicron severity = 25% of Delta - 200% increase on May 1',
                                  'VE HIGH & Omicron severity = 25% of Delta - 200% increase on May 1'),
                        lwd = c(NA,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2),
                        pch = c(1,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA))
  

  
  if(high=="LOW"){
    if(date=='0214'){  col_lib <- col_lib[c(1,3,5,7),]  }     
    if(date=='0301'){  col_lib <- col_lib[c(1,3,9,11),]  }     
    if(date=='0501'){  col_lib <- col_lib[c(1,3,13,15),]  }     
  } else{
    if(date=='0214'){  col_lib <- col_lib[c(1,4,6,8),]  }     
    if(date=='0301'){  col_lib <- col_lib[c(1,4,10,12),]  }     
    if(date=='0501'){  col_lib <- col_lib[c(1,4,14,16),]  }  
  }

# Nicolas (temp)
  # col_lib <- data.frame(tags = c('xxxx',
  #                                'cnt0_0_0_0_0_n',
  #                                paste0('cnt0_0_0_0_0.*omicronSev050.*',high,'.*',vaccine),
  #                                paste0('cnt0_0_0_0_0.*omicronSev033.*',high,'.*',vaccine),
  #                                paste0('cnt0_0_0_0_0.*omicronSev025.*',high,'.*',vaccine),
  #                                paste0('cnt0_0_0_0_100.*omicronSev050.*',high,'.*',vaccine),
  #                                paste0('cnt0_0_0_0_100.*omicronSev033.*',high,'.*',vaccine),
  #                                paste0('cnt0_0_0_0_100.*omicronSev025.*',high,'.*',vaccine)),
  #                       col  = c('black','darkgreen','darkred','orange3','blue','red3','orange','dodgerblue'),
  #                       label = c(paste0('Reported (',max(be_ref_data$date),')'),
  #                                 'Based on reported uptake and social contact behaviour',
  #                                 'Omicron severity = 50% of Delta - current behaviour',
  #                                 'Omicron severity = 33% of Delta - current behaviour',
  #                                 'Omicron severity = 25% of Delta - current behaviour',
  #                                 'Omicron severity = 50% of Delta - Sept-Nov average behaviour',
  #                                 'Omicron severity = 33% of Delta - Sept-Nov average behaviour',
  #                                 'Omicron severity = 25% of Delta - Sept-Nov average behaviour'),
  #                       lwd = c(NA,2,2,2,2,2,2,2),
  #                       pch = c(1,NA,NA,NA,NA,NA,NA,NA))

   color_tags    <- get_color_multiple(output_files,col_lib = col_lib,show_warnings = FALSE)
   table(color_tags)
   color_missing <- get_color_multiple('aa',col_lib = col_lib,show_warnings = FALSE)
   output_files <- output_files[color_tags != color_missing]
  if(length(output_files) == 0){
    warning("No output files left after scenario selection")
  }
}


## MODEL PARAMETERS ----

# additional figure legend library (cfr. counterfactuals)
if(any(grepl('testopti',output_files))){
 
  col_lib <- data.frame(tags = c('xxxx',
                                 'incr0_0_0_0_0_vaccine_uptake_2ndbooster_vscennonewvac',
                                 'incr0_0_0_0_0_vaccine_uptake_2ndbooster_vscen50pc60p',
                                 'incr0_0_0_0_0_vaccine_uptake_2ndbooster_vscen50pc20p',
                                 '*cnt0_0_0_0_0.*omicronTran200.*omicronD616.*HIGH',
                                 '*cnt0_0_100_100_0.*omicronTran200.*omicronD636.*LOW'),
                        col  = c('black','darkgreen','darkred','orange3','gold','darkblue'),
                        label = c(paste0('Reported (',max(be_ref_data$date),')'),
                                  'Current behavious - no vaccine campaign',
                                  'Current behavious - vaccine campaign 60+ (Sept-Dec)',
                                  'Current behavious - vaccine campaign 20+ (Sept-Dec)',
                                  '1/1 latent period Delta (Transm x2, VE HIGH, start 8/11)',
                                  '1/4 latent period Delta (Transm x2, VE LOW, start 27/11)'),
                        lwd = c(NA,2,2,2,2,2),
                        pch = c(1,NA,NA,NA,NA,NA))
  
   color_tags    <- get_color_multiple(output_files,col_lib = col_lib,show_warnings = FALSE)
   table(color_tags)
   color_missing <- get_color_multiple('aa',col_lib = col_lib,show_warnings = FALSE)
   output_files <- output_files[color_tags != color_missing]
  
}


=======
 
>>>>>>> main
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

# number of CoMix and additional waves
scm_num_cp <- unique(scm_num_waves)

# set calibration and scenario dates
scm_callibration_date    <- unique(sim_day2date(scm_callibration_day))
scm_scenario_cp          <- unique(sim_day2date(unlist(scm_scenario_cp_day)))

# temporary disable additional change points
if(grepl('HUB_R1',output_dir)){
  scm_scenario_cp <- NULL
}

# set region from paramter file(s)
sel_region <- unique(get_region(scm_region))

# check parameter input on region, and stop if not uniform
if(length(sel_region)==0){
  stop('REGION PARAMETER IS MISSING IN THE PARAMETER FILE(S)')
} else if(length(sel_region)>1){
  stop(paste0('REGION PARAMETER IS NOT UNIQUE: "',paste0(sel_region,collapse = '", "'),'"'))
}

# sort and select output files
output_files <- sort(output_files,decreasing = bool_file_sort_decreasing)
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
#be_ref_data <- get_observed_incidence_data(enable_data_reuse = F) # do not re-use stored data files

## PDF STREAM ----
pdf(file.path(output_folder,paste0('incidence_',basename(output_dir),
                                   '_',output_dir_base,
                                   '_',sel_region,
                                   ifelse(min(x_axis_scm_day)==0,'_full',''),'.pdf')),7,4)#,9,4)
par(mar=c(4.5,5,1,1))
x_axis_vaccine = c(300,max(x_axis_scm_day))

## POPULATION LEVEL (polyons) ----
multi_plot_incidence_time(output_files,
                          col_lib = col_lib, 
                          bool_polygon_multi = T,
                          x_axis_scm_day = x_axis_scm_day,
                          be_ref_data = be_ref_data,
                          default_y_axis = default_y_axis,
                          scm_callibration_date = scm_callibration_date,
                          scm_scenario_cp = scm_scenario_cp,
                          scm_num_cp = scm_num_cp)


# # BY AGE (bars and polyon) ----
# if(grepl('_100_n',output_dir)){
# multi_plot_incidence_age_time(output_files = output_files,
#                               col_lib = col_lib,
#                               x_axis_scm_day = x_axis_scm_day,
#                               x_axis_vaccine = x_axis_vaccine)
# }


dev.off()




