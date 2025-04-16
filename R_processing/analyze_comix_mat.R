rm(list = ls())
bool_vac <- TRUE
source('R/main_vaccination.R')
## LOAD AND COMBINE COMIX DATA #####
cnt_change <-  c(0,0,0,0,1) # c(0,0,1.1,1.2,1.3)  # c(0,0,0,0,0) c(1,0,0,0,1) # 1,1.3,1.5
cnt_change_date  <- c('2022-09-01',  # restart june 2022
                      '2022-09-01',  # hub seasonality
                      '2022-09-01',  # increase in contacts, cfr. last wave
                      '2021-09-01',  # increase in contacts, cfr. last wave
                      '2023-07-01')  # reuse all matrices last year  -> careful, file lib_social_contacts must be adapted
cnt_change_pt  <- sim_date2day(cnt_change_date)
sel_region="belgium"
num_chains=1
num_stochastic_real=1
data_path     <- "data/comix_matrices_socialmixr/"
chains_param_files = dir('data/config/',pattern='MCMCmulti_20230810_d1213_e2_i1000_n10_p10_crit7_hubR5foronly_pessi_belgium_foronly',full.names = T,recursive = T)
chains_param_file <- chains_param_files[grepl(sel_region,chains_param_files) ]
print(chains_param_file)
parms_chains     = read.table(chains_param_file, sep = ",", header = T)
# select the parameter sets and make a copy for each stochastic realisations
parms_chains  <- parms_chains[rep(1:num_chains,num_stochastic_real),]
i_seq             <- round(seq(1,nrow(parms_chains),length.out = 1))
parms <- unlist(parms_chains[nrow(parms_chains)-(i_seq[1]-1), ])
## make sure that the model parameters are in line with latest requirements
# as such, update colnames and add default values if needed
##-------------------------------------------------------- -
#parms <- update2latest_model_parameter_config(parms)
parms_names <- names(parms)
nCoMix <- sum(grepl('log_comix_coef',parms_names)) / 10  # number of comix related column names, divided by 10 age groups
nOther <- sum(grepl('log_add_coef',parms_names)) / 10 # number of additional matrix related column names, divided by 10 age groups

CoMix_coef_mat = matrix(exp(parms[grepl('log_comix_coef',parms_names)]), nrow = nCoMix, ncol = 10, byrow = T);
Add_coef_mat   = matrix(exp(parms[grepl('log_add_coef',parms_names)]), nrow = nOther, ncol = 10, byrow = T);

# set file path
# set number of waves to include
comix_dates <- read.table(file.path(data_path,'CoMix_dates.csv'))
n_waves     <- nrow(comix_dates)-1 # one row contains the dates for the 2010 survey

# set file names
data_path     <- "data/comix_matrices_socialmixr/"
C_CoMix_names <- c(sprintf('CoMix%i_asy',1:n_waves),
                   sprintf('CoMix%i_sy',1:n_waves),
                   "CoMix2010_asy","CoMix2010_sy",
                   "CoMix_dates")
CoMix_mat <- list()
for(s_comix in C_CoMix_names){
  CoMix_mat[[paste0('C_',s_comix)]] <- as.matrix(read.table(file.path(data_path,paste0(s_comix,'.csv'))))  
}



rap_sy <- list()
rap_asy <- list()
eig_sy <- list()
eig_asy <- list()


M_list_all <- get_susceptiblity_matrices(CoMix_matrices=CoMix_mat, 
                                         CoMix_coef_mat, 
                                         Add_coef_mat,
                                         cnt_change_pt,cnt_change)

for (i in 1:n_waves){
  rap_sy[i]=norm(M_list_all$M_estim_sy[,,i],"F")#/norm(CoMix_mat$C_CoMix2010_sy,"F")
  rap_asy[i]=norm(M_list_all$M_estim_asy[,,i],"F")#/norm(CoMix_mat$C_CoMix2010_asy,"F")
  eig_asy[i]=eigen(M_list_all$M_estim_asy[,,i])[1]
  eig_sy[i]=eigen(M_list_all$M_estim_sy[,,i])[1]
  
}