########################################################################### #
# This file is part of the Stochastic Compartmental Model for SARS-COV-2 
# transmission in Belgium, conceived by members of SIMID group during the 
# COVID19 pandemic.
#
# This file is used to generate regional vaccine uptake files to different
# uptake scenarios, based on the latest reported uptake.
#
# Copyright 2021, SIMID, University of Antwerp & Hasselt University                                        
########################################################################## # 

############################################################# #
# uptake scenarios
# ==>> reported vaccine uptake (= current)
# ==>> vaccine uptake only +20 years of age
# ==>> vaccine uptake only +60 years of age
# ==>> no vaccine uptake
#
############################################################# #

# clear workspace
rm(list=ls())


library(scales)
library(data.table)
library(openxlsx)
library(foreach)

source('./R/lib_stochastic_model.R')


# SETTINGS          ----
###################### # 

# set scenario tag
scen_tag <- 'nov05_region'

# output folder
output_folder <- paste0('output/uptake_',scen_tag)

# select region
sel_region <- "belgium" # belgium, flanders, wallonia, brussels

# uptake file
uptake_file <- dir(output_folder,pattern = paste0('vaccine_uptake_2doses_Vx.*',sel_region),full.names = T)

# uptake all
uptake_db <- read.table(uptake_file,sep=',',header=T)

# no uptake
uptake_zero <- uptake_db
uptake_zero[,-1] <- uptake_zero[,-1] *0

# no uptake -60y
uptake_elderly <- uptake_db
uptake_elderly[,c(2:7,12:18,22:28,32:38)] <- uptake_elderly[,c(2:7,12:18,22:28,32:38)]*0
colSums(uptake_elderly[,-1])

# no uptake -20y
uptake_adult <- uptake_db
age_sel <- c(2:3,12:13,22:23,32:33)
uptake_adult[,age_sel] <- uptake_adult[,age_sel]*0
colSums(uptake_adult[,-1])

# write to file
write.table(uptake_db,paste0(dirname(uptake_file),'/vaccine_uptake_2doses_scenario_current_',sel_region,'.csv'),sep=',',row.names = F)
write.table(uptake_zero,paste0(dirname(uptake_file),'/vaccine_uptake_2doses_scenario_zero_',sel_region,'.csv'),sep=',',row.names = F)
write.table(uptake_elderly,paste0(dirname(uptake_file),'/vaccine_uptake_2doses_scenario_elderly_',sel_region,'.csv'),sep=',',row.names = F)
write.table(uptake_adult,paste0(dirname(uptake_file),'/vaccine_uptake_2doses_scenario_adult_',sel_region,'.csv'),sep=',',row.names = F)

## EXPLORE
range(uptake_db[,-1])
