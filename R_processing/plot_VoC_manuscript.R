########################################################################### #
# This file is part of the Stochastic Compartmental Model for SARS-COV-2 
# transmission in Belgium, conceived by members of SIMID group during the 
# COVID19 pandemic.
#
# This file is used to plot the figures for the VoC manuscript
#
# This file can be sourced with command line arguments. For example:
# - terminal: Rscript R/plot_scm_combined.R belgium 
# - terminal: Rscript R/plot_scm_combined.R belgium &
# - within R: source('R/plot_scm_combined.R')
# - within R: system('Rscript R/plot_scm_combined.R belgium &')
#
# Copyright 2022, SIMID, University of Antwerp & Hasselt University                                        
########################################################################### #

if(0==1){

  system('Rscript R_processing/plot_VoC_manuscript.R orig &')
  system('Rscript R_processing/plot_VoC_manuscript.R delta &')
  system('Rscript R_processing/plot_VoC_manuscript.R omicron &')
  # system('Rscript R_processing/plot_VoC_manuscript.R alpha_900 &')
  # system('Rscript R_processing/plot_VoC_manuscript.R rna2_900 &')
  system('Rscript R_processing/plot_VoC_manuscript.R sep20 &')
  system('Rscript R_processing/plot_VoC_manuscript.R sep21 &')
  #system('Rscript R_processing/plot_VoC_manuscript.R oct21 &') 
  #system('Rscript R_processing/plot_VoC_manuscript.R calibr &') 
  system('Rscript R_processing/plot_VoC_manuscript.R full &')
  
}

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
output_dir <- "output/scm24_0613/"
#output_dir <- "output/scm24_0506/MCMC_140327_final_recap/omicron_760_cnt0_0_0_0_0_n20_belgium/"

prefix_subdir <- 'sep21'
#prefix_subdir <- 'sep20'
prefix_subdir <- 'omicron'
prefix_subdir <- 'delta'
prefix_subdir <- 'full'
prefix_subdir <- 'orig'

# option to specify subdirectory via command line arguments
if(length(args)>=1){
  prefix_subdir   <- args[1]
} 
print(prefix_subdir)

# temporary
sel_region <- "belgium"

bool_combined <- TRUE # option for technote v20211209


# set x-axis limits
x_axis_scm_day <- c(550,720)
# x_axis_scm_day <- c(306,720)
plot_cumulative_reference <- '2021-07-01' #'2021-09-01'
# x_axis_scm_day <- c(0,900)
if(prefix_subdir == 'sep21' | prefix_subdir == 'oct21'){
  x_axis_scm_day <- c(520,650)
}

if(prefix_subdir == 'sep20'){
  x_axis_scm_day <- c(165,300)
  plot_cumulative_reference <- '2020-03-01'
}

if(prefix_subdir == 'calibr'){
  x_axis_scm_day <- c(500,900)
}

if(grepl('alpha',prefix_subdir)){
  x_axis_scm_day <- c(500,720)
}

if(grepl('full',prefix_subdir)){
  x_axis_scm_day <- c(0,710)
}

# use default y axis?
default_y_axis <- TRUE

# load output file
output_files <- dir(output_dir,pattern=sel_region,full.names = T,recursive = T)

# check output files
if(length(output_files)==0){
  stop('OUTPUT SELECTION ISSUE IN "',output_dir,'" WITH "',sel_region,'"')
}

# VoC: 
# - all_760, 
# - omicron_760, 
# - delta_760
# - waning_760
# - w40
# - w41
scenario_legend <- data.frame(orig = '',
                              omicron = 'no Omicron VOC', 
                              delta = 'no Delta VOC',
                              alpha_900 = 'no Alpha VOC',
                              rna2_900 = 'no vaccine-related waning immunity',
                              sep21 = 'September 2021 scenarios',
                              w41_760 = 're-use q-param from wave 41 onward')

# vaccination:
# - ubooster60_jun01_belgium          # reported
# - ubooster100_may01_belgium         # increased booster uptake
# - ubooster60_jun01_belgium          # reduced booster uptake
# - ubooster60_jun01_infant_belgium   # increased infant uptake

# set figure legend and color libray
col_lib <- data.frame(tags = c('xxxx',
                               'orig.*ubooster60_jun01_belgium',
                                paste0(prefix_subdir,'.*0_x0_0_0_0.*ubooster60_jun01_belgium'),
                                paste0(prefix_subdir,'.*0_0_0_0_0.*ubooster100_may01_belgium'),
                                paste0(prefix_subdir,'.*0_0_0_0_0.*ubooster60_may01_belgium'),
                                paste0(prefix_subdir,'.*0_0_0_0_0.*ubooster60_jun01_infant_belgium'),
                                paste0(prefix_subdir,'.*0_100_0_0'),
                                paste0(prefix_subdir,'.*0_130_0'),
                                paste0(prefix_subdir,'.*0_150_0'),
                                paste0(prefix_subdir,'.*0_170_0'),
                               paste0(prefix_subdir,'.*0_0_0_100')),
                      col  = c('black','darkred',
                               'darkgreen','lightblue','orange3','darkblue',
                               'darkgreen','orange3','darkblue','gold','lightblue'),
                      label = c(paste0('Reported (',max(be_ref_data$date),')'),
                                'Estimated dynamics, no introductions and reported vaccine uptake',
 
                                'Reported vaccine uptake',
                                'Increased booster dose uptake up to 2-dose levels',
                                'Reduced booster dose uptake to 60% of 2-dose levels',
                                'Increased vaccine uptake 5-11y in 2021',
                                
                                'Extended summer holiday (from time point A)',
                                '130% of August (from time point A)',
                                '150% of August (from time point A)',
                                '170% of August (from time point A)',
                                'Dynamics Sep-Oct 2020 (from time point A)'
                                ),
                      lwd = c(NA,2,2,2,2,2,2,2,2,2,2),
                      pch = c(1,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA))
bool_file_sort_decreasing <- T

if(prefix_subdir == 'orig'){
  col_lib<- col_lib[-3,]
}

if(!prefix_subdir %in% c('sep20','sep21','oct21')){
  col_lib<- col_lib[-(7:11),]
}

if(prefix_subdir %in% c('sep20','sep21','oct21')){
  col_lib$label[2] <- 'Estimated dynamics'
}

if(prefix_subdir == 'full'){
  col_lib<- col_lib[1:2,]
}


i_scen <- names(scenario_legend)[2]
for(i_scen in names(scenario_legend)){
  if(prefix_subdir == i_scen & nchar(scenario_legend[i_scen])>0){
    col_lib$label[3:6] <- paste0(col_lib$label[3:6],' [',scenario_legend[i_scen],']')
  }
}

color_indices <- get_color_index_multiple(output_files,col_lib = col_lib,show_warnings = FALSE)
table(color_indices)
output_files  <- output_files[!is.na(color_indices)]

sel_col_indices <- unique(color_indices,na.last=T)
sel_col_indices <- sel_col_indices[!is.na(sel_col_indices)]
col_lib         <- col_lib[sort(c(1,sel_col_indices)),]



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
scm_scenario_cp          <- sort(unique(sim_day2date(unlist(scm_scenario_cp_day))))

if(prefix_subdir %in% c('sep21','oct21')){
  scm_scenario_cp <- scm_scenario_cp[scm_scenario_cp>as.Date('2021-08-01')]
  #col_lib$label <- gsub('August \\(from time point A\\)','August (from time point B)',col_lib$label)
  if(prefix_subdir %in% c('oct21')){
    col_lib$label <- gsub('August \\(from time point B\\)','August (from time point C)',col_lib$label)
  } else{
    scm_scenario_cp <- scm_scenario_cp[scm_scenario_cp<as.Date('2021-09-15')]
  }
}
if(prefix_subdir %in% c('calibr','omicron','delta','alpha_900','rna2_900','orig','full')){
  scm_scenario_cp <- c()
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

## PDF STREAM ----

pdf_name <- paste0('incidence_',prefix_subdir, #basename(output_dir),
                   '_',output_dir_base,
                   '_',sel_region,
                   ifelse(exists("high"),paste0('_',high),''),
                   ifelse(exists("vaccine"),paste0('_',vaccine),''),
                   ifelse(bool_combined,'_comb',''),
                   ifelse(min(x_axis_scm_day)==0,'_full',''),'.pdf')
print(pdf_name)
pdf(file.path(output_folder,pdf_name),7,4)#,9,4)
par(mar=c(4.5,5,1,1))
x_axis_vaccine = c(300,max(x_axis_scm_day))

# ## POPULATION LEVEL (lines and polyon) ----
multi_plot_incidence_time(output_files,
                          col_lib = col_lib,
                          bool_polygon_multi = T,
                          x_axis_scm_day = x_axis_scm_day,
                          be_ref_data = be_ref_data,
                          default_y_axis = default_y_axis,
                          include_date_in_legend = TRUE,
                          scm_callibration_date = scm_callibration_date,
                          scm_scenario_cp = scm_scenario_cp,
                          scm_num_cp = scm_num_cp)


# # BY AGE (bars and polyon) ----
# multi_plot_incidence_age_time(output_files = output_files,
#                               col_lib = col_lib,
#                               x_axis_scm_day = x_axis_scm_day,
#                               x_axis_vaccine = x_axis_vaccine)
dev.off()

# # CUMULATIVE ----
# multi_plot_incidence_age_time(out
pdf(file.path(output_folder,gsub('incidence_','cumulative_',pdf_name)),7,4)#,9,4)
par(mar=c(4.5,5,1,1))
multi_plot_cumulative_time(output_files,
                          col_lib = col_lib, 
                          x_axis_scm_day = x_axis_scm_day,
                          scen_reference_date = plot_cumulative_reference,
                          scm_callibration_date = scm_callibration_date,
                          scm_scenario_cp = scm_scenario_cp,
                          scm_num_cp = scm_num_cp)
dev.off()
