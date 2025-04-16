########################################################################### #
# This file is part of the Stochastic Compartmental Model for SARS-COV-2 
# transmission in Belgium, conceived by members of SIMID group during the 
# COVID19 pandemic.
#
# This file is used to cross validate the CoMix data rds files.
#
# Copyright 2022, SIMID, University of Antwerp & Hasselt University                                        
########################################################################### #

# set local file names
file_name1 <- 'data/survey_belgium2020_comix_up_to_wave_39.rds'
file_name2 <- 'data/survey_belgium2020_comix_up_to_wave_40.rds'

# load data
comix_v1 <- readRDS(file_name1)
comix_v2 <- readRDS(file_name2)

# specify number of waves
num_wave_v1 <- max(comix_v1$participants$wave)
num_wave_v2 <- max(comix_v2$participants$wave)

# loop over each wave and compare participant and 
# contact data on:
# - dimensions
# - values (if both not missing)
# - missing values
i_wave <- 12
for(i_wave in 1:min(num_wave_v1,num_wave_v2)){
  
  # set as data.frame to prevent issues with factorial variables with different levels
  pdata_v1 <- data.frame(lapply(comix_v1$participants[wave == i_wave,], as.character),stringsAsFactors = F)
  pdata_v2 <- data.frame(lapply(comix_v2$participants[wave == i_wave,], as.character),stringsAsFactors = F)
  
  cdata_v1 <- data.frame(lapply(comix_v1$contacts[part_id %in% pdata_v1$part_id,], as.character),stringsAsFactors = F)
  cdata_v2 <- data.frame(lapply(comix_v2$contacts[part_id %in% pdata_v2$part_id,], as.character),stringsAsFactors = F)
  
  # check participant data
  if(any(dim(pdata_v1) != dim(pdata_v2)) || 
     sum(pdata_v1 != pdata_v2,na.rm=T) > 0 ||
     any(table(is.na(pdata_v1)) != table(is.na(pdata_v2)))){
    
    warning('partipant data of wave ',i_wave,' changed',immediate. = T)
    
    # issue with NA?
    sel_diff <- is.na(pdata_v1) != is.na(pdata_v2)
    colnames(sel_diff)[colSums(sel_diff)>0]
    table(rowSums(sel_diff)>0)
    
  }
  
  # check contact data
  if(any(dim(cdata_v1) != dim(cdata_v2)) || 
     sum(cdata_v1 != cdata_v2,na.rm=T) > 0 ||
     any(table(is.na(cdata_v1)) != table(is.na(cdata_v2)))){
    warning('contact data of wave ',i_wave,' changed',immediate. = T)
    
    # issue with NA?
    sel_diff <- is.na(cdata_v1) != is.na(cdata_v2)
    colnames(sel_diff)[colSums(sel_diff)>0]
    table(rowSums(sel_diff)>0)
    
  }
  
} # end for loop
  