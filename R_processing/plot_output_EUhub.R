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
output_dir <- "output/hubR5"

# HUB output parameters
origine_date <- "2023-06-28"
number_of_weeks <- 49 #-> to adapt




# for Omicron  -> to change
#high <- 'submitted'
alter <- 0
vaccine <- '' #'ubooster90_jan01'
date <- '0301'

sel_region <- "belgium"

# set x-axis limits and date threshold
#origine date
origine_date_num <- as.numeric(as.Date(origine_date)-get_scm_start_date())
#real start date (each period should go from Sunday to Saturday)
target_start <- origine_date_num -as.POSIXlt(origine_date)$wday-1
x_axis_scm_day <- c(target_start-0,target_start+number_of_weeks*7)


#other possibilities
#x_axis_scm_day <- c(700,915)
#x_axis_scm_day <- c(sim_date2day('2020-03-01'),sim_date2day('2023-12-31'))
#x_axis_scm_day <- c(sim_date2day('2023-03-01'),sim_date2day('2024-07-01'))

# use default y axis?
default_y_axis <- TRUE

# option to use command line arguments: use latest output?
#if(length(args)>=2){
 # sel_latest_output <- as.logical(args[2])
  #if(!is.na(sel_latest_output) && sel_latest_output){
   # output_folders       <- dir('output/',pattern = 'scm23',full.names = T,recursive = F)
    #output_folders_mtime <- file.info(output_folders)$mtime
    #output_dir           <- output_folders[order(output_folders_mtime,decreasing = T)[1]]
    #print(paste("USE LATEST OUTPUT:",output_dir))
  #}
#}

# option to specify the region via command line arguments
#if(length(args)>=1){
 # sel_region   <- args[1]
#} 
#if(!exists('sel_region')){
 # sel_region <- ''
#}

# load output file
output_files <- dir(output_dir,pattern=sel_region,full.names = T,recursive = T)

# check output files
if(length(output_files)==0){
  stop('OUTPUT SELECTION ISSUE IN "',output_dir,'" WITH "',sel_region,'"')
}


# set figure legend and color libray
if(any(grepl('hubR5',output_files))){
  
  col_lib <- data.frame(tags = c('xxxx',
                                 paste0('opti','.*0_0_0_0_100_.*_0flu_biv'),
                                 paste0('opti','.*0_0_0_0_100_.*_50flu_biv'),
                                 paste0('opti','.*0_0_0_0_100_.*_100flu_plus15_biv'),
                                 paste0('opti','.*0_0_0_0_100_.*_all60plus_biv'),
                                 paste0('pessi','.*0_0_0_0_100_.*_0flu_biv'),
                                 paste0('pessi','.*0_0_0_0_100_.*_50flu_biv'),
                                 paste0('pessi','.*0_0_0_0_100_.*_100flu_plus15_biv'),
                                 paste0('pessi','.*0_0_0_0_100_.*_all60plus_biv')),
                        col  = c('black','blue','red','orange','green','darkblue','darkred','orange4','darkgreen'),
                        label = c(paste0('Reported (',max(be_ref_data$date),')'),
                              #    'A Optimistic waning no vac',
                               #   'B Optimistic waning 50% flu',
                              #    'C Optimistic waning 100% flu + 15%',
                              #    'D Optimistic waning all 60+',
                              #    'E Pessimistic waning no vac',
                              #    'F Pessimistic waning 50% flu',
                              #    'G Pessimistic waning 100% flu + 15%',
                              #    'H Pessimistic waning all 60+'),
                        'A',
                        'B',
                        'C',
                        'D',
                        'E',
                        'F',
                        'G',
                        'H'),
                        lwd = c(NA,2,2,2,2,2,2,2,2),
                        pch = c(1,NA,NA,NA,NA,NA,NA,NA,NA))
  bool_file_sort_decreasing <- T
  color_tags    <- get_color_multiple(output_files,col_lib = col_lib,show_warnings = FALSE)
  table(color_tags)
  color_missing <- get_color_multiple('aa',col_lib = col_lib,show_warnings = FALSE)
  output_files <- output_files[color_tags != color_missing]
  
  
}

if(any(grepl('future',output_files))){

  col_lib <- data.frame(tags = c('xxxx',
                                 '.*0_0_0_0_100_.*futurepessi1.*_0flu_biv',
                                 '.*0_0_0_0_100_.*futurepessi1.*_a100flu_biv',
                                 '.*0_0_0_0_100_.*futurepessi3.*_0flu_biv',
                                 '.*0_0_0_0_100_.*futurepessi3.*_a100flu_biv',
                                 '.*0_0_0_0_100_.*hloadpessi.*_0flu_nobiv',
                                 '.*0_0_0_0_100_.*hloadpessi.*_a100flu_biv',
                                 paste0(high,'.*a0_0_0_0_100_.*_all60_biv')),
                        col  = c('black','darkblue','darkgreen','orange3','darkred','gold','red','purple'),
                        label = c(paste0('Reported (',max(be_ref_data$date),')'),
                                  'estimatation from last q',
                                  'estimatation from last q + vac',
                                  'estimatation from last 3q',
                                  'estimatation from last 3q + vac',
                                  'constant last q',
                                  'constant last q + vac',
                                  'D all65+ biv'),
                        lwd = c(NA,2,2,2,2,2,2,2),
                        pch = c(1,NA,NA,NA,NA,NA,NA,NA))
  bool_file_sort_decreasing <- T
  color_tags    <- get_color_multiple(output_files,col_lib = col_lib,show_warnings = FALSE)
  table(color_tags)
  color_missing <- get_color_multiple('aa',col_lib = col_lib,show_warnings = FALSE)
  output_files <- output_files[color_tags != color_missing]
  
  
}





if(any(grepl('aaaa',output_files))){
  col_lib <- data.frame(tags = c('xxxx',
                                 '0_0_0_0_0.*calibrated1',
                                 '0_0_0_0_100.*calibrated1',
                                 '0_0_0_0_200.*calibrated1',
                                 '0_0_0_0_100_PESSVAR'),
                        col  = c('black','darkblue','darkred','orange3','darkgreen'),
                        label = c(paste0('Reported (',max(be_ref_data$date),')'),
                                  'Constant Comix and q (wave 63 copy of wave 45 May 2022)',
                                  'Use Comix last year but constant q (wave 63)',
                                  'Use Comix and q last year',
                                  'D'),
                        lwd = c(NA,2,2,2,2),
                        pch = c(1,NA,NA,NA,NA))
  bool_file_sort_decreasing <- T
  color_tags    <- get_color_multiple(output_files,col_lib = col_lib,show_warnings = FALSE)
  table(color_tags)
  color_missing <- get_color_multiple('aa',col_lib = col_lib,show_warnings = FALSE)
  output_files <- output_files[color_tags != color_missing]
  
  
}


if(any(grepl('hub_R4',output_files))){
  col_lib <- data.frame(tags = c('xxxx',
                                 '0_0_0_0_100_OPTI',
                                 '0_0_0_0_100_OPTVAR',
                                 '0_0_0_0_100_PESSI',
                                 '0_0_0_0_100_PESSVAR'),
                        col  = c('black','darkblue','darkred','orange3','darkgreen'),
                        label = c(paste0('Reported (',max(be_ref_data$date),')'),
                                   'A',
                                   'B',
                                   'C',
                                   'D'),
                                  #'A-Optimistic waning',
                                #'B-Optimistic waning + new variant',
                                 # 'C-Pessimistic waning',
                                  #'D-Pessimistic waning + new variant'),
                        lwd = c(NA,2,2,2,2),
                        pch = c(1,NA,NA,NA,NA))
  bool_file_sort_decreasing <- T
  color_tags    <- get_color_multiple(output_files,col_lib = col_lib,show_warnings = FALSE)
  table(color_tags)
  color_missing <- get_color_multiple('aa',col_lib = col_lib,show_warnings = FALSE)
  output_files <- output_files[color_tags != color_missing]
  
  
}

# set figure legend and color libray
# if(any(grepl('bq1',output_files))){
#   
#   col_lib <- data.frame(tags = c('xxxx',
#                                  'august30aaa',
#                                  'august50aaa',
#                                  'jan23',
#                                  'sept30aaa',
#                                  'sept50aaa',
#                                  'sept70aaa'),
#                         col  = c('black','darkblue','darkgreen','orange3','darkred','gold','red'),
#                         label = c(paste0('Reported (',max(be_ref_data$date),')'),
#                                   # 'A',
#                                   # 'B',
#                                   # 'C',
#                                   # 'D',
#                                   # 'E',
#                                   # 'F'),
#                                   'Seed August 1 / -70% immune reduction',
#                                   'Seed August 1 / -50% immune reduction',
#                                   'test BQ1',
#                                   'Seed Sept 1 / -70% immune reduction',
#                                   'Seed Sept 1 / -50% immune reduction',
#                                   'Seed Sept 1 / -30% immune reduction'),
#                         lwd = c(NA,2,2,2,2,2,2),
#                         pch = c(1,NA,NA,NA,NA,NA,NA))
#   bool_file_sort_decreasing <- T
#   color_tags    <- get_color_multiple(output_files,col_lib = col_lib,show_warnings = FALSE)
#   table(color_tags)
#   color_missing <- get_color_multiple('aa',col_lib = col_lib,show_warnings = FALSE)
#   output_files <- output_files[color_tags != color_missing]
#   
#   
# }


# set figure legend and color libray
if(any(grepl('bq1aaa',output_files))){
  
  col_lib <- data.frame(tags = c('xxxx',
                                 'august30aaa',
                                 'august50aaa',
                                 'august70',
                                 'sept30aaa',
                                 'sept50aaa',
                                 'sept70aaa'),
                        col  = c('black','darkblue','darkgreen','orange3','darkred','gold','red'),
                        label = c(paste0('Reported (',max(be_ref_data$date),')'),
                                  # 'A',
                                  # 'B',
                                  # 'C',
                                  # 'D',
                                  # 'E',
                                  # 'F'),
                                  'Seed August 1 / -70% immune reduction',
                                  'Seed August 1 / -50% immune reduction',
                                  'Seed August 1 / -30% immune reduction',
                                  'Seed Sept 1 / -70% immune reduction',
                                  'Seed Sept 1 / -50% immune reduction',
                                  'Seed Sept 1 / -30% immune reduction'),
                        lwd = c(NA,2,2,2,2,2,2),
                        pch = c(1,NA,NA,NA,NA,NA,NA))
  bool_file_sort_decreasing <- T
  color_tags    <- get_color_multiple(output_files,col_lib = col_lib,show_warnings = FALSE)
  table(color_tags)
  color_missing <- get_color_multiple('aa',col_lib = col_lib,show_warnings = FALSE)
  output_files <- output_files[color_tags != color_missing]
  
  
}


# set figure legend and color libray
if(any(grepl('hubR3aa',output_files))){
  
  col_lib <- data.frame(tags = c('xxxx',
                                 'novariant.*0_0_0_0_100_.*vscenonecampaign',
                                 'opti.*0_0_0_0_100_.*vscenonecampaign',
                                 'opti.*0_0_0_0_100_.*vscen6mcampaigna',
                                 'pessi.*0_0_0_0_100_.*vscenonecampaign',
                                 'pessi.*0_0_0_0_100_.*vscen12mcampaigna',
                                 'pessi.*0_0_0_0_100_.*vscen6mcampaigna'),
                        col  = c('black','darkblue','darkgreen','orange3','darkred','gold','red'),
                        label = c(paste0('Reported (',max(be_ref_data$date),')'),
                                  # 'A',
                                  # 'B',
                                  # 'C',
                                  # 'D',
                                  # 'E',
                                  # 'F'),
                                  'No new variant',
                                  'New variant with weak reduction in immunity (ECDC hub Round 3 scenario A)',
                                  'C-vaccination campaign every 6 months',
                                  'New variant with strong reduction in immunity (ECDC hub Round 3 scenario D)',
                                  'E-vaccination campaign every 12 months (pessimistic vocs)',
                                  'F-vaccination campaign every 6 months (pessimistic vocs)'),
                        lwd = c(NA,2,2,2,2,2,2),
                        pch = c(1,NA,NA,NA,NA,NA,NA))
  bool_file_sort_decreasing <- T
  color_tags    <- get_color_multiple(output_files,col_lib = col_lib,show_warnings = FALSE)
  table(color_tags)
  color_missing <- get_color_multiple('aa',col_lib = col_lib,show_warnings = FALSE)
  output_files <- output_files[color_tags != color_missing]
  
  
}


if(any(grepl('hubR3aaa',output_files))){
  
  col_lib <- data.frame(tags = c('xxxx',
                                 'opti.*0_0_0_0_100_.*vscenonecampaigna',
                                 'opti.*0_0_0_0_100_.*vscen12mcampaigna',
                                 'opti.*0_0_0_0_100_.*vscen6mcampaigna',
                                 'pessi.*0_0_0_0_100_.*vscenonecampaign',
                                 'pessi.*0_0_0_0_100_.*vscen12mcampaign',
                                 'pessi.*0_0_0_0_100_.*vscen6mcampaign'),
                        col  = c('black','darkblue','darkred','orange3','darkgreen','gold','red'),
                        label = c(paste0('Reported (',max(be_ref_data$date),')'),
                                  # 'A',
                                  # 'B',
                                  # 'C',
                                  # 'D',
                                  # 'E',
                                  # 'F'),
                        'A-no repeated vaccination campaign',
                        'B-vaccination campaign every 12 months',
                        'C-vaccination campaign every 6 months',
                        'D-no repeated vaccination campaign (pessimistic vocs)',
                        'E-vaccination campaign every 12 months (pessimistic vocs)',
                        'F-vaccination campaign every 6 months (pessimistic vocs)'),
                        lwd = c(NA,2,2,2,2,2,2),
                        pch = c(1,NA,NA,NA,NA,NA,NA))
  bool_file_sort_decreasing <- T
  color_tags    <- get_color_multiple(output_files,col_lib = col_lib,show_warnings = FALSE)
  table(color_tags)
  color_missing <- get_color_multiple('aa',col_lib = col_lib,show_warnings = FALSE)
  output_files <- output_files[color_tags != color_missing]

  
}

# set figure legend and color libray
if(any(grepl('augustaa',output_files))){
  
  high <- 'alldelta'
  col_lib <- data.frame(tags = c('xxxx',
                                 'aug22no.*cnt0_0_0_0_0.*vscennocampaign280822',
                                 'aug22delta.*cnt100_0_0_0_0.*vscennocampaign280822',
                                 'aug22delta.*cnt100_0_0_0_0.*vscen65pall280822',
                                 'aug22delta.*cnt100_0_0_0_0.*vscen50pall280822',
                                 'aug22delta.*cnt100_0_0_0_0.*vscen18pall280822'),
                        col  = c('black','darkblue','darkred','darkorange','darkorchid','darkgreen'),#,'green3','green'),
                        label = c(paste0('Reported (',max(be_ref_data$date),')'),
                                  'Current behaviour - no vaccine campaign',
                                  'Increased transmission - no vaccine campaign',
                                  'Increased transmission - vaccine campaign 65+ dedicated booster (100% 1st booster, Sept)',
                                  'Increased transmission - vaccine campaign 50+ dedicated booster (100% 1st booster, Sept-Oct)',
                                  'Increased transmission - vaccine campaign 18+ dedicated booster (100% 1st booster, Sept-Dec)'),
                        lwd = c(NA,2,2,2,2,2),
                        pch = c(1,NA,NA,NA,NA,NA))
  
  high <- '50delta' 
  col_lib <- data.frame(tags = c('xxxx',
                                 'aug22no.*cnt0_0_0_0_0.*vscennocampaign280822',
                                 'aug22delta.*cnt100_0_0_0_0.*vscennocampaign280822',
                                 'aug22delta.*cnt100_0_0_0_0.*vscen65p50pc280822',
                                 'aug22delta.*cnt100_0_0_0_0.*vscen50p50pc280822',
                                 'aug22delta.*cnt100_0_0_0_0.*vscen18p50pc280822'),
                        col  = c('black','darkblue','darkred','darkorange','darkorchid','darkgreen'),#,'green3','green'),
                        label = c(paste0('Reported (',max(be_ref_data$date),')'),
                                  'Current behaviour - no vaccine campaign',
                                  'Increased transmission - no vaccine campaign',
                                  'Increased transmission - vaccine campaign 65+ dedicated booster (50% 1st booster, Sept)',
                                  'Increased transmission - vaccine campaign 50+ dedicated booster (50% 1st booster, Sept-Oct)',
                                  'Increased transmission - vaccine campaign 18+ dedicated booster (50% 1st booster, Sept-Dec)'),
                        lwd = c(NA,2,2,2,2,2),
                        pch = c(1,NA,NA,NA,NA,NA))
  
  high <- 'allinit'
  col_lib <- data.frame(tags = c('xxxx',
                                 'aug22no.*cnt0_0_0_0_0.*vscennocampaign280822',
                                 'aug22no.*cnt100_0_0_0_0.*vscennocampaign280822',
                                 'aug22no.*cnt100_0_0_0_0.*vscen65pall280822',
                                 'aug22no.*cnt100_0_0_0_0.*vscen50pall280822',
                                 'aug22no.*cnt100_0_0_0_0.*vscen18pall280822'),
                        col  = c('black','darkblue','darkred','darkorange','darkorchid','darkgreen'),#,'green3','green'),
                        label = c(paste0('Reported (',max(be_ref_data$date),')'),
                                  'Current behaviour - no vaccine campaign',
                                  'Increased transmission - no vaccine campaign',
                                  'Increased transmission - vaccine campaign 65+ initial booster (100% 1st booster, Sept)',
                                  'Increased transmission - vaccine campaign 50+ initial booster (100% 1st booster, Sept-Oct)',
                                  'Increased transmission - vaccine campaign 18+ initial booster (100% 1st booster, Sept-Dec)'),
                        lwd = c(NA,2,2,2,2,2),
                        pch = c(1,NA,NA,NA,NA,NA))

  high <- '50init'
  col_lib <- data.frame(tags = c('xxxx',
                                 'aug22no.*cnt0_0_0_0_0.*vscennocampaign280822',
                                 'aug22no.*cnt100_0_0_0_0.*vscennocampaign280822',
                                 'aug22no.*cnt100_0_0_0_0.*vscen65p50pc280822',
                                 'aug22no.*cnt100_0_0_0_0.*vscen50p50pc280822',
                                 'aug22no.*cnt100_0_0_0_0.*vscen18p50pc280822'),
                        col  = c('black','darkblue','darkred','darkorange','darkorchid','darkgreen'),#,'green3','green'),
                        label = c(paste0('Reported (',max(be_ref_data$date),')'),
                                  'Current behaviour - no vaccine campaign',
                                  'Increased transmission - no vaccine campaign',
                                  'Increased transmission - vaccine campaign 65+ initial booster (50% 1st booster, Sept)',
                                  'Increased transmission - vaccine campaign 50+ initial booster (50% 1st booster, Sept-Oct)',
                                  'Increased transmission - vaccine campaign 18+ initial booster (50% 1st booster, Sept-Dec)'),
                        lwd = c(NA,2,2,2,2,2),
                        pch = c(1,NA,NA,NA,NA,NA))
  
  #col_lib <- col_lib[c(1,2,5,6),] 
  #col_lib <- col_lib[c(-5),] 

  
  color_tags    <- get_color_multiple(output_files,col_lib = col_lib,show_warnings = FALSE)
  table(color_tags)
  color_missing <- get_color_multiple('aa',col_lib = col_lib,show_warnings = FALSE)
  output_files <- output_files[color_tags != color_missing]
  
  # sort(output_files)
  # color_indices <- get_color_index_multiple(output_files,col_lib = col_lib,show_warnings = FALSE)
  # table(color_indices)
  # output_files  <- output_files[!is.na(color_indices)]
  # sel_col_indices <- unique(color_indices,na.last=T)
  # sel_col_indices <- sel_col_indices[!is.na(sel_col_indices)]
  # col_lib         <- col_lib[sort(c(1,sel_col_indices)),]
  
}



# additional figure legend library 
if(any(grepl('omicronaa',output_files))){
  col_lib <- data.frame(tags = c('xxxx',
                                 paste0('copie',vaccine),
                                 paste0('1100',vaccine),
                                 'cnt0_0_0_0_0_n',
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
                        col  = c('black','blue','darkorange','blue3','green','darkgreen','orange','darkorange','green','darkgreen','orange','darkorange','green','darkgreen','orange','darkorange'),
                        label = c(paste0('Reported (',max(be_ref_data$date),')'),
                                  'Base',
                                  'New',
                                  'Waning 50% 6 months vaccine and infection',
                                  'D-2022-03-04',
                                  'E-2022-03-04',
                                  'F-2022-03-04',                                  
                                  'G-2022-03-04',
                                  'H-2022-03-04',
                                  'I-2022-03-04',                                  
                                  'J-2022-03-04',
                                  'K-2022-03-04',
                                  'L-2022-03-04',                                  
                                  'M-2022-03-04',
                                  'N-2022-03-04',
                                  'O-2022-03-04'),
                        lwd = c(NA,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2),
                        pch = c(1,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA))
  

  
  # if(high=="LOW"){
  #   if(date=='0214'){  col_lib <- col_lib[c(1,3,5,7),]  }     
  #   if(date=='0301'){  col_lib <- col_lib[c(1,3,9,11),]  }     
  #   if(date=='0501'){  col_lib <- col_lib[c(1,3,13,15),]  }     
  # } else{
  #   if(date=='0214'){  col_lib <- col_lib[c(1,4,6,8),]  }     
  #   if(date=='0301'){  col_lib <- col_lib[c(1,4,10,12),]  }     
  #   if(date=='0501'){  col_lib <- col_lib[c(1,4,14,16),]  }  
  # }

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

#################################################################@@@@@@@@@@
##csv output for EU hub ----
#################################################################@@@@@@@@@@
if(1==1){  #to avoid this part

#output_files_inc <- output_files[grepl('new_infect',output_files) | grepl('new_sympt',output_files) | grepl('hosp_adforcovid_incr',output_files) | grepl('icu_adm_incr',output_files) | grepl('scen_mortality',output_files)]
#output_files_inc <- output_files[grepl('new_infect',output_files)]
#output_files_inc <- output_files[grepl('hosp_adforcovid_incr',output_files)]
#output_files_inc <- output_files[grepl('new_infect',output_files)| grepl('hosp_adforcovid_incr',output_files) | grepl('scen_mortality',output_files)]

#output_files_inc <- output_files[grepl('scen_cases',output_files) | grepl('inc_I_hosp_for_covid',output_files) | grepl('inc_I_icu',output_files) | grepl('mortality_age',output_files)]
output_files_inc <- output_files[grepl('scen_cases',output_files) | grepl('inc_I_hosp_for_covid',output_files) | grepl('mortality_age',output_files) | grepl('scen_numvac',output_files)]
#output_files_inc <- output_files[grepl('scen_cases',output_files)]
#output_files_inc <- output_files[grepl('mortality_age',output_files)]
#output_files_inc <- output_files[grepl('inc_I_hosp_for_covid',output_files)]
#output_files_inc <- output_files[grepl('scen_numvac',output_files)]


threshold_date <- as.Date(get_scm_start_date()+target_start+7*(0:number_of_weeks))

 i_file <- 1
 foreach(i_file = 1:length(output_files_inc),
         .combine='rbind') %do%{
           
           target_name <- "unknown"
if(grepl('scen_cases',output_files_inc[i_file])){
  target_name <- "inc infection"
}  else if(grepl('mortality_age',output_files_inc[i_file])){
  target_name <- "inc death"
}  else if(grepl('inc_I_hosp_for_covid',output_files_inc[i_file])){ 
  target_name <- "inc hosp"
} else if(grepl('scen_numvac',output_files_inc[i_file])){ #???? to check
  target_name <- "inc dose"
}   
     
           scen_data <- readRDS(output_files_inc[i_file])
           scen_col  <- get_col(output_files_inc[i_file],col_lib)
           scen_label<- col_lib$label[col_lib$col == scen_col]
           scm_dates  <- get_scm_start_date() +  0:(ncol(scen_data)-1)
           sel_dates <- scm_dates %in% threshold_date
           
           if(target_name == "inc hosp"){
             scen_data_aggre <- apply(scen_data[,,],2,rowSums)
             scen_data_similar <- scen_data[,,]
           } else {
             scen_data_aggre <- apply(scen_data[,,-1],2,rowSums)
             scen_data_similar <- scen_data[,,-1]
           }
           
           scen_data_aggre[is.na(scen_data_aggre)] <- 0
           scen_data_similar[is.na(scen_data_similar)] <- 0
           
foreach(sample=c(1:dim(scen_data_aggre)[1]),.combine='rbind') %do%{
  foreach(age_id = c("child","adult","older"),.combine='rbind') %do%{
    
        #   foreach(quant = c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99),
         #          .combine='rbind') %do%{
    
    if(age_id=="child"){
      scen_data_age <- scen_data_similar[,,1] + 0.8*scen_data_similar[,,2]
    }
    if(age_id=="adult"){
      scen_data_age <- 0.2*scen_data_similar[,,2] + apply(scen_data_similar[,,3:6],2,rowSums)
    }
    if(age_id=="older"){
      scen_data_age <- apply(scen_data_similar[,,7:10],2,rowSums)
    }

  sample_data <- scen_data_age[sample,]
  sample_data_cum <- cumsum(sample_data)
  sample_data_cum <- sample_data_cum[sel_dates]
  sample_data_tbl <- round(diff(sample_data_cum))
  #scen_data_aggre_mean <- apply(scen_data_aggre,2,quantile,quant)
  #scen_data_aggre_mean_cum <- cumsum(scen_data_aggre_mean)
  #scen_data_aggre_mean_cum <- scen_data_aggre_mean_cum[sel_dates]
  #scen_data_aggre_mean_tbl <- round(diff(scen_data_aggre_mean_cum))

  scen_dates_tbl     <- format(threshold_date[-1], "%Y-%m-%d") 
  
  tbl_partial <- data.frame(
              cbind(location = "BE",
                scenario_id = scen_label,
                #origin_date = origine_date,
                #target_end_date = scen_dates_tbl,
                horizon = paste0(1:number_of_weeks," wk"),
                target_variable = target_name,
                age_id = age_id,
               #type = "quantile",
                #quantile = quant,
               sample = sample,
                value = sample_data_tbl))
   return(tbl_partial)
                   } ->  tbl_out
           return(tbl_out)  
         } ->  tbl_out2
 return(tbl_out2)  
 } -> tbl_summary
 
 write.table(tbl_summary,
             file = file.path(output_folder,paste0(origine_date,'-SIMID-SCM.csv')),
            sep=',',
            row.names = F)
}





