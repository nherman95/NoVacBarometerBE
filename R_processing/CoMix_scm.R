########################################################################### #
# This file is part of the Stochastic Compartmental Model for SARS-COV-2 
# transmission in Belgium, conceived by members of SIMID group during the 
# COVID19 pandemic.
#
# This file is used to pre-process the original CoMix data into matrices for 
# the Stochastic Model for healthy and symptomatic behaviour.
#
# Copyright 2021, SIMID, University of Antwerp & Hasselt University                                        
########################################################################### #

# clear workspace
rm(list=ls())

# load libraries and scripts
library(socialmixr)
library(foreach)
source('/Users/lwillem/Documents/university/research/social_contacts/repo/socrates_rshiny/R/plot_social_contact_matrix.R')
source('/Users/lwillem/Documents/university/research/social_contacts/repo/socrates_rshiny/R/npsp/simage.R')
source('/Users/lwillem/Documents/university/research/social_contacts/repo/socrates_rshiny/R/npsp/splot.R')
source('/Users/lwillem/Documents/university/research/social_contacts/repo/socrates_rshiny/R/plot_mean_number_contacts.R')
source('./R/lib_stochastic_model.R')

# source('../socrates_rshiny-master/R/plot_social_contact_matrix.R')
# source('../socrates_rshiny-master/R/npsp/simage.R')
# source('../socrates_rshiny-master/R/npsp/splot.R')
# source('../socrates_rshiny-master/R/plot_mean_number_contacts.R')
# source('./R/lib_stochastic_model.R')
# 
# source('../socrates_rshiny-master/R_socialmixr/check.r')
# source('../socrates_rshiny-master/R_socialmixr/contact_matrix.r')
# source('../socrates_rshiny-master/R_socialmixr/survey.r')
# source('../socrates_rshiny-master/R_socialmixr/wpp_age.r')
# library(httr)
# library(jsonlite)
# library(XML)
# library(curl)
# library(wpp2019)
# source('../socrates_rshiny-master/R_socialmixr/get_survey.r')
# library('countrycode')

# settings ----
wave_select <- 46
sel_waves      <- c(1:wave_select,2010)
output_dir     <- 'data/comix_matrices_socialmixr'

# load data
comix_all <- readRDS(paste0('data/survey_belgium2020_comix_up_to_wave_',wave_select,'.rds'))

# age breaks
age_breaks <- seq(0,80,10)

# set rng stream
set.seed(2021)

# clean some participants by hand
comix_all$participants <- comix_all$participants[!part_id %in% c(2179916,2179917),]

comix_all$participants$country <- as.character(comix_all$participants$country)
comix_all$participants$country <- "Belgium"
comix_all$participants$country <- as.factor(comix_all$participants$country)
table(comix_all$participants$country)

# load data
survey2010_all <- readRDS('data/survey_belgium2010.rds')
survey2010_all$participants <- survey2010_all$participants[holiday == FALSE,]
dim(survey2010_all$participants)
survey2010_all$participants[,wave:=2010]
survey2010_all$contacts[is_hh_member == FALSE & 
                          cnt_home == TRUE,cnt_otherplace := 1]
survey2010_all$contacts[is_hh_member == FALSE & 
                          cnt_home == TRUE,cnt_home := 0]
names(survey2010_all$contacts)

comix_all$participants[,holiday:=FALSE]
comix_all$participants[,weekday:=dayofweek %in% 1:5]
comix_all$contacts[,    is_imputed:=FALSE]
comix_all$contacts[,    is_hh_member:=FALSE]

table(survey2010_all$participants$part_id %in% comix_all$participants$part_id)
flag <- (survey2010_all$participants$part_id %in% comix_all$participants$part_id)
survey2010_all$participants$part_id[flag]

flag <- (comix_all$participants$part_id %in% survey2010_all$participants$part_id)
comix_all$participants$part_id[flag]
comix_all$participants[flag,]

survey2010_all$participants[,part_id := paste0(part_id,'p')]
survey2010_all$contacts[,part_id := paste0(part_id,'p')]

comix_all$participants[,part_id := paste0(part_id,'c')]
comix_all$contacts[,part_id := paste0(part_id,'c')]

table(table(comix_all$contacts$cont_id))
table(table(comix_all$contacts$cont_id))

names(survey2010_all$participants) %in% names(comix_all$participants)
names(survey2010_all$contacts) %in% names(comix_all$contacts)

sel_col_part <- names(survey2010_all$participants)
sel_col_cnt  <- names(survey2010_all$contacts)

comix_all$participants <- comix_all$participants[,..sel_col_part]
comix_all$contacts     <- comix_all$contacts[,..sel_col_cnt]

comix_all$participants <- rbind(comix_all$participants,
                                survey2010_all$participants)
comix_all$contacts     <- rbind(comix_all$contacts,
                                survey2010_all$contacts)

# pre-process
if(!dir.exists(output_dir)) dir.create(output_dir,recursive = T)
format_num_digits <- 2

names(comix_all$contacts)
comix_all$contacts[is.na(cnt_otherplace),cnt_otherplace := 1]
table(comix_all$contacts$cnt_otherplace,useNA = 'ifany')
table(comix_all$contacts$cnt_work,useNA       = 'ifany')
table(comix_all$contacts$cnt_school,useNA     = 'ifany')


# remove duplicated contact locations
# order: c("Home","Work","School","Transport","Leisure","Otherplace")
comix_all$contacts[cnt_home == 1 , 
                   cnt_work := 0]
comix_all$contacts[cnt_home == 1 | cnt_work == 1, 
                   cnt_school := 0]
comix_all$contacts[cnt_home == 1 | cnt_work == 1 | cnt_school == 1,
                   cnt_transport := 0]
comix_all$contacts[cnt_home == 1 | cnt_work == 1 | cnt_school == 1 | cnt_transport == 1,
                   cnt_leisure := 0]
comix_all$contacts[cnt_home == 1 | cnt_work == 1 | cnt_school == 1 | cnt_transport == 1 | cnt_leisure == 1,
                   cnt_otherplace := 0]

# symptomatic
locations = c("Home","Work","Transport","School","Leisure","Otherplace")
location_vec_sympt = rbind(rep(1,6),
                          c(1,0.09,0.13,0.09,0.06,0.25))
opt_filter <- data.frame('cnt_home'       = 1,
                         'cnt_work'       = 1,
                         'cnt_transport'  = 1,
                         'cnt_school'     = 1,
                         'cnt_leisure'    = 1,
                         'cnt_otherplace' = 1)

# get population data
pop_matrix_n10 <- get_be_pop_matrix(10)
pop_matrix     <- get_be_pop_matrix(9)
pop_matrix     <- cbind(pop_matrix[,1:8],rowSums(pop_matrix[,9:10]))

# location specific matrices ----
i_wave <- sel_waves[1] ; i_filter <- 1
i_wave <- 2010;i_filter <- 1
for(i_wave in sel_waves){

  set.seed(2021 + i_wave)
  
  # (re)start from all data, and select subset
  comix_sel              <- comix_all 
  comix_sel$participants <- comix_sel$participants[wave == i_wave,] 
  comix_sel$contacts     <- comix_sel$contacts[part_id %in% comix_sel$participants$part_id,]
  
  for(i_filter in 1:length(opt_filter)){
  
    print(i_filter)
    
    # calculate contact matrix with contact rate per participant
    scmxr_out              <- suppressMessages(suppressWarnings(contact_matrix(comix_sel,
                                               estimated.contact.age = "sample",
                                               missing.contact.age   = "sample",
                                               symmetric        = F,
                                               weigh.dayofweek  = TRUE,
                                               weigh.age        = TRUE,
                                               weight.threshold = 3,
                                               age.limits       = age_breaks,
                                               filter           = opt_filter[i_filter],
                                               counts           = F)))
  
    # calculate contact rate per capita (using the reference year)
    if(i_wave == 2010){
      scmxr_out$matrix_per_capita <- scmxr_out$matrix / matrix(rep(scmxr_out$demography$population,9),ncol=9,byrow=T)
    } else{
      scmxr_out$matrix_per_capita <- scmxr_out$matrix / pop_matrix
    }
    #plot_cnt_matrix(scmxr_out$matrix_per_capita)
    
    # wave 1-8: impute contacts for children 0-10y
    if(is.na(scmxr_out$matrix_per_capita[1,1])){
      scmxr_out$matrix_per_capita[1,]     <- scmxr_out$matrix_per_capita[,1]
      scmxr_out$matrix_per_capita[1,1]    <- scmxr_out$matrix_per_capita[2,2]
    }
    
    # save contact matrix
    cmat_file_name <- paste0(output_dir,'/wave_',i_wave,'_',names(opt_filter[i_filter]),'.csv')
    cmat_file_name
    write.table(scmxr_out$matrix_per_capita,
                 file=cmat_file_name,
                 sep='\t',row.names=F,col.names=F)
  }
}

#Contact matrix March 2020, lockdown, asymptomatic ---------------
# 100%home, 80% Telework, schools closed, 20% Transport, 10% leisure, 10% Other places--
locations

# (a)symptomatic behaviour ----
i_health <- 1; i_wave <- 2010
for(i_health in 1:nrow(location_vec_sympt))
for(i_wave in sel_waves){
 
  cmat_sy <- matrix(0,nrow=10,ncol=10)

  for(i_filter in 1:length(opt_filter)){
    cmat_file_name <- paste0(output_dir,'/wave_',i_wave,'_',names(opt_filter[i_filter]),'.csv')
    cmat_locations <- read.table( file = cmat_file_name,sep='\t')
    
    cmat_sy <- cmat_sy + cmat_locations * location_vec_sympt[i_health,i_filter]
  }
  cmat_sy

  # asymptomatic ==>> make symmetric
  if(i_health == 1){
    ## set c_{ij} N_i and c_{ji} N_j (which should both be equal) to
    ## 0.5 * their sum; then c_{ij} is that sum / N_i
    weighted.matrix <- as.matrix(cmat_sy) * pop_matrix
    normalised.weighted.matrix <- diag(pop_matrix[1,]) %*% weighted.matrix
    weighted.matrix <- 0.5 * diag(1/pop_matrix[1,]) %*%
      (normalised.weighted.matrix + t(normalised.weighted.matrix))
    cmat_sy <- weighted.matrix / pop_matrix
  }

  # impute contacts with/for elderly +90y
  cmat_sy_n10        <- rbind(cbind(cmat_sy,0),0)
  cmat_sy_n10[,10]   <- cmat_sy_n10[,9] # per capita rate => denominator changes when split age group 9 and 10
  cmat_sy_n10[10,10] <- cmat_sy_n10[9,9] #* pop_matrix_n10[1,9] /  pop_matrix_n10[1,10]
  cmat_sy_n10[10,]   <- cmat_sy_n10[,10] #* pop_matrix_n10[1,10] / pop_matrix_n10[1,]
  
  cmat_file_name <- paste0(output_dir,'/CoMix',i_wave,'_',ifelse(i_health==1,'asy','sy'),'.csv')
  cmat_file_name
  write.table(cmat_sy_n10,
              file=cmat_file_name,
              sep='\t',row.names=F,col.names=F)
}

## inspect dates ####
i_wave <- 1
foreach(i_wave = sel_waves,
        .combine = 'rbind') %do% {
  
  # (re)start from all data, and select subset
  p_sel              <- comix_all$participants$wave == i_wave
  
  p_dates <- as.Date(paste(comix_all$participants$year[p_sel],
                           comix_all$participants$month[p_sel],
                           comix_all$participants$day[p_sel],sep='-'))
  as.character(range(p_dates))
} -> dates_all

cmat_file_name <- paste0(output_dir,'/CoMix_dates.csv')
cmat_file_name
write.table(dates_all,
            file=cmat_file_name,
            sep='\t',row.names=F,col.names=F)


# plot matrices and contact rates ----
# get dates
# comix_start_dates <- sim_day2date(get_M_change_day()[is_M_CoMix()])
comix_start_dates <- data.frame(wave = sel_waves,
                                start = as.Date(dates_all[,1]))

# set row and column names
scmxr_out  <- suppressMessages(suppressWarnings(contact_matrix(polymod,age.limits       = c(age_breaks,90))))
mij_new    <- scmxr_out$matrix
mij_prev   <- scmxr_out$matrix

i_wave <- 22
pdf(file.path(output_dir,'all_cnt_matrices_SCM_sympt.pdf'),10,10)
par(mfrow=c(2,2))
for(i_wave in sel_waves){
  
  mij_new[]  <- as.matrix(read.table(file=paste0(output_dir,'/CoMix',i_wave,'_asy.csv'),sep='\t')) * get_be_pop_matrix(10)
  mij_prev[] <- as.matrix(read.table(file=paste0(output_dir,'/CoMix',i_wave,'_sy.csv'),sep='\t')) * get_be_pop_matrix(10)
  
  plot_cnt_matrix(mij_new)
  text(-0.2,1.1,paste('ASYMPT\nWAVE',i_wave),xpd=T)
  suppressWarnings(text(-0.2,-0.15,comix_start_dates$start[comix_start_dates$wave==as.numeric(i_wave)],xpd=T))
  plot_cnt_matrix(mij_prev)  
  text(-0.2,1.1,paste('SYMPTOMATIC\nWAVE',i_wave),xpd=T)
  suppressWarnings(text(-0.2,-0.15,comix_start_dates$start[comix_start_dates$wave==as.numeric(i_wave)],xpd=T))
  
  plot_mean_number_contacts(mij_new,scale_max = 15)
  plot_mean_number_contacts(mij_prev,scale_max = 15)
  
}
dev.off()  

# averages ----

foreach(i_wave = sel_waves) %do%{
  
  set.seed(2021 + i_wave)
  
    # (re)start from all data, and select subset
    comix_sel              <- comix_all 
    comix_sel$participants <- comix_sel$participants[wave == i_wave,] 
    
    # calculate contact matrix with contact rate per participant
    scmxr_out              <- suppressMessages(suppressWarnings(contact_matrix(comix_sel,
                                                                               estimated.contact.age = "sample",
                                                                               missing.contact.age   = "sample",
                                                                               symmetric        = F,
                                                                               weigh.dayofweek  = TRUE,
                                                                               weigh.age        = TRUE,
                                                                               weight.threshold = 3,
                                                                               age.limits       = 0,
                                                                               counts           = F)))
    scmxr_out$matrix
} -> mean_cnt_wave


data.frame(wave = sel_waves,
           date = as.character(comix_start_dates$start),
           mean_contacts = unlist(mean_cnt_wave)
          ) -> contacts_out
write.table(contacts_out,file='output/contacts_waves.csv',sep=',')


if(0==1){
  
  
  comix_prev <- readRDS('data/survey_belgium2020_comix_up_to_wave_43.rds')
  comix_new <- readRDS('data/survey_belgium2020_comix_up_to_wave_44.rds')
  
  table(comix_prev$participants$wave)
  table(comix_new$participants$wave)
  
  table(comix_new$participants$country)
  
  table(table(comix_prev$participants$part_id))
  table(table(comix_new$participants$part_id))

  table(comix_new$participants$part_id %in% comix_prev$participants$part_id)
  table(comix_prev$participants$part_id %in% comix_new$participants$part_id)
    
  dim(comix_prev$contacts)
  dim(comix_prev$contacts[part_id %in% comix_prev$participants$part_id,])
  dim(comix_new$contacts[part_id %in% comix_prev$participants$part_id,])
  sum(comix_new$contacts$part_id %in% comix_prev$participants$part_id)
  
  dim(comix_prev$participants)
  dim(comix_prev$participants[part_id %in% comix_prev$participants$part_id,])
  dim(comix_new$participants[part_id %in% comix_prev$participants$part_id,])
  
  # load data
  comix_prev <- readRDS('data/survey_belgium2020_comix_up_to_wave_43.rds')
  comix_new <- readRDS('data/survey_belgium2020_comix_up_to_wave_44.rds')
  
  # count number of contacts per participant for wave 1:35
  cnt_up_to_35 <- table(comix_prev$contacts[part_id %in% comix_prev$participants$part_id,'part_id'])
  cnt_up_to_36 <- table(comix_new$contacts[part_id %in% comix_prev$participants$part_id,'part_id'])

  # number of particpants
  length(cnt_up_to_35)
  length(cnt_up_to_36)

  # look for differences  (and select first 10)
  flag <- which(cnt_up_to_35!=cnt_up_to_36)
  flag <- flag[1:10] 
  
  # explore
  rbind(up_to_35 = cnt_up_to_35[flag],
        up_to_36 = cnt_up_to_36[flag])

  comix_prev$contacts[part_id == 2052612 ,c(1:14)]
  comix_new$contacts[part_id == 2052612 ,c(1:14)]
  
  comix_prev$participants[part_id == 2052612 ,c(1:14)]
  comix_new$participants[part_id == 2052612 ,c(1:14)]

  ii <- c(1,3)
  comix_prev$participants[wave == 43 ,3:4] == comix_new$participants[wave == 43 ,3:4]
  names(comix_prev$participants)[3:4]
    
  table(comix_prev$participants[wave == 12 ,part_gender],useNA = 'ifany')
  table(comix_new$participants[wave == 12 ,part_gender],useNA = 'ifany')
  
  prev  <- unlist(as.character(comix_prev$participants[wave == 12 ,part_gender]))
  new   <- unlist(as.character(comix_new$participants[wave == 12 ,part_gender]))
  
  table(prev,useNA = 'ifany')
  table(new,useNA = 'ifany')
  
  prev[is.na(prev)] <- 'NA'
  new[is.na(new)] <- 'NA'
  
  table(prev == new)
  
  table(comix_prev$participants$part_gender[wave == 43])
  table(comix_new$participants$part_gender)
  
  xx <- compare.list(comix_prev$participants[wave==43,],
               comix_new$participants[wave==43,])
  names(comix_prev$participants)[!xx]
  
  xx <- compare.list(comix_prev$contacts[part_id %in% (comix_prev$participants[wave==43,part_id]),],
                     comix_new$contacts[part_id %in% (comix_new$participants[wave==43,part_id]),])
  names(comix_prev$contacts)[!xx]
  
  comix_prev$participants[wave==43,c("hh_type","sday_id")]
  comix_new$participants[wave==43,c("hh_type","sday_id")]
  
  comix_prev$participants[wave==43,"hh_type"] == comix_new$participants[wave==43,"hh_type"]
  
  levels(comix_prev$participants$sday_id)
  levels(comix_new$participants$sday_id)
  
  levels(comix_prev$contacts$cnt_gender)
  levels(comix_new$contacts$cnt_gender)
  
  sum(comix_prev$contacts[part_id %in% (comix_prev$participants[wave==39,part_id]),cnt_work],na.rm=T)
  sum(comix_new$contacts[part_id %in% (comix_new$participants[wave==39,part_id]),cnt_work],na.rm=T)
  
  
}

