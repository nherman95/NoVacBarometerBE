########################################################################### #
# This file is part of the Stochastic Compartmental Model for SARS-COV-2 
# transmission in Belgium, conceived by members of SIMID group during the 
# COVID19 pandemic.
#
# This file is used to load data and help functions for the stochastic model. 
# Variables are provided via the global R environment, which should be improved
# in future versions.
#
# Copyright 2022, SIMID, University of Antwerp & Hasselt University                                        
########################################################################### #

# Libraries
#---------- -
suppressPackageStartupMessages(library(LaplacesDemon))
suppressPackageStartupMessages(library(poisbinom))
suppressPackageStartupMessages(library(MMWRweek))
suppressPackageStartupMessages(library(dplyr))

# load core and help functions
source('R/lib_model_core.R')
if(bool_vac==FALSE){
if(bool_bar==TRUE){
  source('R/lib_combined_plot.R')
} else{
  source('R/lib_combined_plot_explode.R')
}
} else{
  source('R/lib_combined_plot.R')
}
source('R/lib_model_parameters.R')
source('R/lib_stochastic_model.R')
source('R/lib_social_contacts.R')
source('R/lib_projections_vaccination.R')

# start from fixed reference hospital data
be_ref_data <- get_latest_incidence_data(enable_data_reuse = FALSE)
#be_ref_data <- read.csv('data/covid19_reference_data_20220707_1319_belgium.csv')
obs_total_hosp_ext <- be_ref_data$hospital_admissions + be_ref_data$hospital_admissions_other
obs_total_hosp_ext[is.na(obs_total_hosp_ext)] <- 0

# reported admissions for COVID
obs_forcovid_hosp     <- be_ref_data$hospital_admissions 

# fix to align previous code
obs_total_hosp <- obs_total_hosp_ext[-(1:7)]

# hospital and ICU load
obs_hosp_load_ext <- be_ref_data$hospital_load
obs_icu_load_ext  <- be_ref_data$icu_load
obs_hosp_load_ext[is.na(obs_hosp_load_ext)] <- 0
obs_icu_load_ext[is.na(obs_icu_load_ext)]   <- 0

# hospital discharges
obs_hosp_discharges_ext <- be_ref_data$hospital_discharges
obs_hosp_discharges_ext[is.na(obs_hosp_discharges_ext)] <- 0

# hospital recoveries: estimated by load, admissions and mortality
obs_hosp_recoveries_ext <- be_ref_data$hospital_recovery_7dmean
obs_hosp_recoveries_ext[is.na(obs_hosp_recoveries_ext)] <- 0

# hospital exit
obs_hosp_exit_ext <- be_ref_data$hospital_exit
obs_hosp_exit_ext[is.na(obs_hosp_exit_ext)] <- 0

# mortality
obs_age_dist_mort <- be_ref_data[,grepl('covid19_deaths_age',names(be_ref_data))]

# Serological survey data  ----
#------------------------ -
# Collection rounds 1 and 2
load("data/20200508_data_serology_collection_all.Rdata")
sero_data_all <- subset(dt_iggcat_age_cat, select = c("collection round", "age_cat","igg_cat_pos","igg_cat_total","prev_pos"))
names(sero_data_all)[1] = "cround"

# Age-dependent asymptomatic proportions  ----
#--------------------------------------- -
p_vec = c(0.94,0.90,0.84,0.61,0.49,0.21,0.02,0.02,0.02,0.02)

# # Severity of symptoms upon symptomatic infection (based on data reported by CDC)  ----
# #-------------------------------------------------------------------------------- -
phi0_vec = c(0.98,0.98,0.79,0.79,0.67,0.67,0.50,0.35,0.32,0.32);
phi1_vec = c(0.99,0.99,0.85,0.85,0.76,0.76,0.73,0.69,0.74,0.74);
phi2_vec = 1 - phi1_vec;

# Age-specific hospitalization (weekly age distribution - clinical surveillance data)  ----
#------------------------------------------------------------------------------------ -
# age_dist_hosp_mat_ext <- get_regional_hospital_age_distr()
#age_dist_hosp_mat_ext <- get_hospital_age_distr(bool_use_default = TRUE)            #NEW!! (20220704)
age_dist_hosp_mat_ext <- t(be_ref_data[,grepl('prop_hospital_admission_age',names(be_ref_data))])

# hosp_age_orig <- get_regional_hospital_age_distr()
# hosp_age_exit <- get_hospital_age_distr()
# i_age <- 5
# for(i_age in 1:10){
#   plot(hosp_age_orig[i_age,],main=i_age)
#   lines(hosp_age_exit[i_age,],lwd=2,col=4)
# }
# 
# table(hosp_age_exit==0)
# plot(colSums(hosp_age_exit==0))
# plot(rowSums(hosp_age_exit==0,na.rm = T))
# plot(colSums(hosp_age_orig==0))
# plot(rowSums(hosp_age_orig==0,na.rm = T))

## Helper functions  ----
##----------------- -
expit = function(eta){exp(eta)/(1+exp(eta))}
logit = function(pp){log(pp/(1-pp))}
proposalfunction <- function(parms, sd = 0.005){
  return(rnorm(length(parms), mean = parms, sd = sd))
}

round_tot <- function(x, digits = 0){
  up <- 10 ^ digits
  x <- x * up
  y <- floor(x)
  indices <- tail(order(x-y), round(sum(x,na.rm = T)) - sum(y,na.rm = T))
  y[indices] <- y[indices] + 1
  return(y / up)
}

## Model initialization based on confirmed cases  ----
##-------------------------------------------------- -
aggr_dat <- be_ref_data[,grepl('cases_',names(be_ref_data))]

age_dist_cc <- apply(aggr_dat[1:13,], 2, sum)
total_nr_cc <- sum(age_dist_cc)
rel_freq_cc <- age_dist_cc/total_nr_cc

## Variants of Concern (VOC)  ----
##------------------------ -
# !!! two ways : WGS excel file (old) and ECDC hub file (from GISAID and TESSy)
# set VOC types
global_lib_voc <- data.table(name = c('VOC_alpha',
                                      'VOC_delta',
                                      'VOC_omicron',
                                      'VOC_ba4ba5',
                                      'VOC_bq1ba275xbb',
                                      'VOC_xbb15'),
                      start_date_ll = c("2020-11-06", "2021-04-05", "2021-11-01","2022-1-20","2022-06-1","2022-09-1"), # lower limit
                      start_date_ul = c("2020-12-26", "2021-07-14", "2021-12-06","2022-6-30","2022-10-31","2023-02-28")  # upper limit
                      )
                      
global_lib_voc$model_strain <- rep(2:1,nrow(global_lib_voc))[1:nrow(global_lib_voc)]
global_lib_voc$start_ll <- sim_date2day(global_lib_voc$start_date_ll)
global_lib_voc$start_ul <- sim_date2day(global_lib_voc$start_date_ul)

#---------------------------------
#WGS excel file (old)
# load baseline surveillance WGS data from NRC
voc_data_file <- './data/VOC_WGS/SARS-CoV-2 WGS consortium.xlsx'

# define begin/end of baseline and active surveillence columns
db_voc_columns <- read.xlsx(voc_data_file,sheet = 2,rows = 8:9)
db_voc_column_begin <- which(db_voc_columns %in% 'Invalid results')
#db_voc_column_end   <- which(db_voc_columns %in% '% VOCs')
db_voc_column_end <- c(21,37) #temporary fix for the new file


# get raw data on baseline surveillence
db_voc_raw <- read.xlsx(voc_data_file,sheet = 2,startRow = 9,cols = 1:db_voc_column_end[1])
names(db_voc_raw)
#db_voc_raw

# get row data on active surveillence

db_voc_active <- read.xlsx(voc_data_file,sheet = 2,startRow = 9,cols = db_voc_column_begin[2]:db_voc_column_end[2])
names(db_voc_active) <- paste0(names(db_voc_active),'.active')
names(db_voc_active)

# merge baseline and active surveillance data
db_voc_raw <- cbind(db_voc_raw,db_voc_active)

# modify dates
db_voc_raw$date_start   <- as.Date('2020-11-30') + (0:(nrow(db_voc_raw)-1))*7
db_voc_raw$date_stop    <- db_voc_raw$date_start + 6
db_voc_raw$date         <- db_voc_raw$date_start + 3
#fix due to file change
db_voc_raw$Year         <- db_voc_raw[,1]
db_voc_raw$totaltest <- db_voc_raw[,20]

# select valid rows (with samples and only 2021)
#db_voc_raw <- db_voc_raw[!is.na(db_voc_raw$Year) & (db_voc_raw$`Total.tests.performed`>0) & (db_voc_raw$Year>=2021),]
db_voc_raw <- db_voc_raw[!is.na(db_voc_raw$Year) & (db_voc_raw$totaltest>0) & (db_voc_raw$Year>=2021),]

# remove last period (incomplete)
db_voc_raw <- db_voc_raw[-nrow(db_voc_raw),]

db_voc_raw$voc_alpha <- as.numeric(db_voc_raw$`501Y.V1`)
db_voc_raw$voc_beta  <- as.numeric(db_voc_raw$`501Y.V2`)
db_voc_raw$voc_gamma <- as.numeric(db_voc_raw$`501Y.V3`)
db_voc_raw$voc_delta <- as.numeric(db_voc_raw$`B.1.617.2`)
db_voc_raw$voc_omicron <- as.numeric(db_voc_raw$`BA.1`) + as.numeric(db_voc_raw$`BA.2`)
db_voc_raw$voc_ba1     <- as.numeric(db_voc_raw$`BA.1`)
db_voc_raw$voc_ba2     <- as.numeric(db_voc_raw$`BA.2`)
db_voc_raw$voc_ba4_ba5 <- as.numeric(db_voc_raw$`BA.4`) + as.numeric(db_voc_raw$`BA.5`) + as.numeric(db_voc_raw$`BA.2.75`)
db_voc_raw$voc_bq1 <- as.numeric(db_voc_raw$`BQ.1*`) + as.numeric(db_voc_raw$`Recombinants`)
db_voc_raw$n_sequenced <- as.numeric(db_voc_raw$totaltest)

db_voc_raw$voc_aggr <- db_voc_raw$voc_alpha + db_voc_raw$voc_beta + db_voc_raw$voc_gamma

# # load active surveillance WGS data from NRC
# db_voc_raw$voc_omicron_active <- db_voc_raw$`B.1.1.529.active`
# db_voc_raw$n_sequenced_active <- db_voc_raw$`Total.tests.performed.active`


#---------------------------------
# get data from ECDC hub file (from GISAID and TESSy)
# https://github.com/covid19-forecast-hub-europe/covid19-scenario-hub-europe/wiki/Data-sources
voc_data_file <- './data/VOC_WGS/data.csv'
db_voc_raw2 <- read.table(voc_data_file,sep=',',header = T)
db_voc_raw2 <- db_voc_raw2[!is.na(db_voc_raw2$country_code),]
db_voc_raw2 <- db_voc_raw2[db_voc_raw2$country_code == "BE",]
db_voc_raw2 <- db_voc_raw2[-1,]

# modify dates
db_voc_raw2$date_start   <- MMWRweek2Date(as.numeric(substr(db_voc_raw2$year_week,1,4)),as.numeric(substr(db_voc_raw2$year_week,6,7)))
db_voc_raw2$date_stop    <- db_voc_raw2$date_start + 6
db_voc_raw2$date         <- db_voc_raw2$date_start + 3
#fix due to file change
db_voc_raw2$Year         <- as.numeric(substr(db_voc_raw2$year_week,1,4))
db_voc_raw2$totaltest <- db_voc_raw2$number_sequenced
db_voc_raw2$n_sequenced <- db_voc_raw2$totaltest

db_voc_raw2 <- db_voc_raw2[!is.na(db_voc_raw2$Year) & (db_voc_raw2$totaltest>0) & (db_voc_raw2$Year>=2021),]

db_voc_raw2_extract <- db_voc_raw2[db_voc_raw2$variant == 'B.1.1.7',][ -c(1:12) ]
db_voc_raw2_extract$voc_alpha <- db_voc_raw2[db_voc_raw2$variant == 'B.1.1.7',]$number_detections_variant
db_voc_raw2_extractnew <- db_voc_raw2[db_voc_raw2$variant == 'B.1.351',][ -c(1:12) ]
db_voc_raw2_extractnew$voc_beta  <- db_voc_raw2[db_voc_raw2$variant == 'B.1.351',]$number_detections_variant
db_voc_raw2_extract <- merge(db_voc_raw2_extract,db_voc_raw2_extractnew)
db_voc_raw2_extractnew <- db_voc_raw2[db_voc_raw2$variant == 'P.1',][ -c(1:12) ]
db_voc_raw2_extractnew$voc_gamma <- db_voc_raw2[db_voc_raw2$variant == 'P.1',]$number_detections_variant
db_voc_raw2_extract <- merge(db_voc_raw2_extract,db_voc_raw2_extractnew)
db_voc_raw2_extractnew <- db_voc_raw2[db_voc_raw2$variant == 'B.1.617.2',][ -c(1:12) ]
db_voc_raw2_extractnew$voc_delta <- db_voc_raw2[db_voc_raw2$variant == 'B.1.617.2',]$number_detections_variant
db_voc_raw2_extract <- merge(db_voc_raw2_extract,db_voc_raw2_extractnew)
db_voc_raw2_extractnew <- db_voc_raw2[db_voc_raw2$variant == 'BA.1',][ -c(1:12) ]
db_voc_raw2_extractnew$voc_ba1     <- db_voc_raw2[db_voc_raw2$variant == 'BA.1',]$number_detections_variant
db_voc_raw2_extract <- merge(db_voc_raw2_extract,db_voc_raw2_extractnew)
db_voc_raw2_extractnew <- db_voc_raw2[db_voc_raw2$variant == 'BA.2',][ -c(1:12) ]
db_voc_raw2_extractnew$voc_ba2     <- db_voc_raw2[db_voc_raw2$variant == 'BA.2',]$number_detections_variant
db_voc_raw2_extract <- merge(db_voc_raw2_extract,db_voc_raw2_extractnew)
db_voc_raw2_extractnew <- db_voc_raw2[db_voc_raw2$variant == 'BA.4',][ -c(1:12) ]
db_voc_raw2_extractnew$voc_ba4 <- db_voc_raw2[db_voc_raw2$variant == 'BA.4',]$number_detections_variant
db_voc_raw2_extract <- merge(db_voc_raw2_extract,db_voc_raw2_extractnew)
db_voc_raw2_extractnew <- db_voc_raw2[db_voc_raw2$variant == 'BA.5',][ -c(1:12) ]
db_voc_raw2_extractnew$voc_ba5 <- db_voc_raw2[db_voc_raw2$variant == 'BA.5',]$number_detections_variant
db_voc_raw2_extract <- merge(db_voc_raw2_extract,db_voc_raw2_extractnew)
db_voc_raw2_extractnew <- db_voc_raw2[db_voc_raw2$variant == 'BQ.1',][ -c(1:12) ]
db_voc_raw2_extractnew$voc_bq1 <- db_voc_raw2[db_voc_raw2$variant == 'BQ.1',]$number_detections_variant
db_voc_raw2_extract <- merge(db_voc_raw2_extract,db_voc_raw2_extractnew)
db_voc_raw2_extractnew <- db_voc_raw2[db_voc_raw2$variant == 'BA.2.75',][ -c(1:12) ]
db_voc_raw2_extractnew$voc_BA275 <- db_voc_raw2[db_voc_raw2$variant == 'BA.2.75',]$number_detections_variant
db_voc_raw2_extract <- merge(db_voc_raw2_extract,db_voc_raw2_extractnew)
db_voc_raw2_extractnew <- db_voc_raw2[db_voc_raw2$variant == 'XBB',][ -c(1:12) ]
db_voc_raw2_extractnew$voc_xbb <- db_voc_raw2[db_voc_raw2$variant == 'XBB',]$number_detections_variant
db_voc_raw2_extract <- merge(db_voc_raw2_extract,db_voc_raw2_extractnew)
db_voc_raw2_extractnew <- db_voc_raw2[db_voc_raw2$variant == 'XBB.1.5' ,][ -c(1:12) ]
db_voc_raw2_extractnew$voc_xbb15 <- db_voc_raw2[db_voc_raw2$variant  == 'XBB.1.5',]$number_detections_variant
db_voc_raw2_extract <- merge(db_voc_raw2_extract,db_voc_raw2_extractnew)

db_voc_raw2_extract$voc_aggr <- db_voc_raw2_extract$voc_alpha + db_voc_raw2_extract$voc_beta + db_voc_raw2_extract$voc_gamma
db_voc_raw2_extract$voc_omicron <- db_voc_raw2_extract$voc_ba1 + db_voc_raw2_extract$voc_ba2
db_voc_raw2_extract$voc_ba4ba5 <- db_voc_raw2_extract$voc_ba4 + db_voc_raw2_extract$voc_ba5
db_voc_raw2_extract$voc_bq1ba275xbb <- db_voc_raw2_extract$voc_bq1 + db_voc_raw2_extract$voc_BA275 + db_voc_raw2_extract$voc_xbb


db_voc_raw <- db_voc_raw2_extract


##############################################################################
# load data and logistic fit from Tom W. 'deprecated)
sgtf_data <- read.table('./data/VOC_WGS/wenseleers_raw_data_share_omicron_SA_ENG_SCOT_DK_BE.csv',sep=',',header = T)
sgtf_fit  <- read.table('./data/VOC_WGS/wenseleers_emmeans_sgtf_20211229.csv',sep=',',header = T)

# select Belgian data and fit
sgtf_data <- sgtf_data[sgtf_data$country == "Belgium",]
sgtf_fit  <- sgtf_fit[sgtf_fit$country == "Belgium",]

names(sgtf_data)

sgtf_data$date <- as.Date(sgtf_data$date)
sgtf_fit$date  <- as.Date(sgtf_fit$date)

# rename columns
sgtf_data$sftg_logistic_prob_omicron <- sgtf_data$prop_omicron
sgtf_data$sftg_total   <- sgtf_data$pos_tests
sgtf_data$sftg_omicron <- round(sgtf_data$omicron)
sgtf_fit$sftg_logistic_prob_omicron <- sgtf_fit$prob


# select from December 2021
sgtf_data <- sgtf_data[sgtf_data$date >= as.Date('2021-12-01'),]
sgtf_fit <- sgtf_fit[sgtf_fit$date >= as.Date('2021-12-01'),]

# plot(sgtf_fit$date,
#        sgtf_fit$sftg_logistic_prob_omicron,
#      type='l',col=3)
# points(sgtf_data$date,
#      sgtf_data$sftg_logistic_prob_omicron)
# points(db_voc_raw$date,
#        1-db_voc_raw$`%.VOCs`,
#        pch=16)


#db_voc_raw <- merge(db_voc_raw,
                     #sgtf_fit[,c('date','sftg_logistic_prob_omicron')],
                 #    sgtf_data[,c('date','sftg_logistic_prob_omicron','sftg_total','sftg_omicron')],
                  #   all = T)



