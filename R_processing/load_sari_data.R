library(data.table); library(tidyverse)

# Script to estimate the extrapolated population SARI-numbers from the sentinel surveillance
# The method is based on the numbers reported to the sentinel surveillance and the estimated market share of these sentinel hospitals

#####################################################################################################################
## Load in datasets 

#single bootstrap sample from actual data (nrow and admidates unchanged)
sari2023 <- readRDS('data/sari_test/sim_dt.Rdata')    
#single bootstrap sample from actual data (nrow and admidates remained, no NA)
sim.noNA.dt <- readRDS('data/sari_test/sim_noNA_dt.Rdata') 

#Catchment population by age group and hosp (=population number that the hospital is representative for in this age group)
ma3 <- readRDS('data/sari_test/catchmentpop_by_ag_hosp.Rdata') 

#Total population by age group 
pop_ag3 <- readRDS('data/sari_test/pop_ag3.Rdata')


######################################################################################################################

# Zero notifications are added to the data (=> in the future hospitals should report 0)
sari2023_ag <- sari2023[, .N, by=c('admidate', 'agegr', 'hosp')] #OR use sarscov2=='Confirmed SARS-CoV-2'? corona==1 (refers to all corona-virusses)
sari_zero <- CJ(unique(sari2023$hosp), unique(sari2023$admidate), unique(sari2023$agegr), 0)
names(sari_zero) <- c('hosp', 'admidate', 'agegr', 'N')
sari2023_ag <- rbind(sari2023_ag[!is.na(admidate) & !is.na(agegr)], sari_zero[!is.na(admidate) & !is.na(agegr)], fill=T)
sari2023_ag <- sari2023_ag[!duplicated(sari2023_ag[,c('admidate', 'agegr', 'hosp')])]

# Merge Market Share/CatchmentPop with the aggregated SARI-file 
# To do this, the age groups need to be adapted
sari2023_ag[, ag:=dplyr::recode(agegr, '(-1,12]'='Jonger dan 15 jaar', '(12,17]'='Jonger dan 15 jaar', 
                                '(17,24]'='Tussen 15 en 64 jaar', '(24,44]'='Tussen 15 en 64 jaar', '(44,64]'='Tussen 15 en 64 jaar', 
                                '(64,74]'='65 jaar en ouder', '(74,84]'='65 jaar en ouder', '(84,Inf]'='65 jaar en ouder')]
sari2023_ag2 <- merge(sari2023_ag, ma3, by=c('ag', 'hosp'), allow.cartesian=TRUE)

# Aggregate out the hospitals (Estimate catchment population by only age group and date)
# I typically have several samples to allow for uncertainty in the 'catchment population', not necessary in this example
col_name_all2 <- c('N', 'sample_1') 
sari2023_ag3 <- sari2023_ag2[, lapply(.SD, sum, na.rm=TRUE), by=c('agegr', 'admidate'), .SDcols = col_name_all2]


# Merge by age group population by region 
sari2023_ag3[, ag:=dplyr::recode(agegr, '(-1,12]'='Jonger dan 15 jaar', '(12,17]'='Jonger dan 15 jaar', 
                                 '(17,24]'='Tussen 15 en 64 jaar', '(24,44]'='Tussen 15 en 64 jaar', '(44,64]'='Tussen 15 en 64 jaar', 
                                 '(64,74]'='65 jaar en ouder', '(74,84]'='65 jaar en ouder', '(84,Inf]'='65 jaar en ouder')]
sari2023_ag5 <- merge(sari2023_ag3, pop_ag3, by=c('ag'))

## EXTRAPOLATION, small numbers can be eliminated (e.g. sample_1>200) 
sari2023_ag5[, c('exh.est'):=list(N*(pop/sample_1))]
sari2023_ag_test<- sari2023_ag5[, lapply(.SD, sum, na.rm=TRUE), by=c('ag', 'admidate'), .SDcols = c('N', 'exh.est') ]
sari2023_ag_hosp<- sari2023_ag_test[, lapply(.SD, sum, na.rm=TRUE), by=c('admidate'), .SDcols = c('N','exh.est')]