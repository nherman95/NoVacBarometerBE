########################################################################### #
# This file is part of the Stochastic Compartmental Model for SARS-COV-2 
# transmission in Belgium, conceived by members of SIMID group during the 
# COVID19 pandemic.
#
# This file contains functions to pre-process the social contact data and to
# define the social contact scenarios.
#
# Copyright 2022, SIMID, University of Antwerp & Hasselt University                                        
########################################################################### #

################################################################ #
## LOAD DATA (GLOBAL ENVIRONMENT) ----
################################################################ #

## LOAD AND COMBINE COMIX DATA #####

# set file path
data_path     <- "data/comix_matrices_socialmixr/"

# set number of waves to include
comix_dates <- read.table(file.path(data_path,'CoMix_dates.csv'))
n_waves     <- nrow(comix_dates)-1 # one row contains the dates for the 2010 survey

# set file names
data_path     <- "data/comix_matrices_socialmixr/"
C_CoMix_names <- c(sprintf('CoMix%i_asy',1:n_waves),
                   sprintf('CoMix%i_sy',1:n_waves),
                   "CoMix2010_asy","CoMix2010_sy",
                   "CoMix_dates")


# load and combine in a single list
CoMix_mat <- list()
for(s_comix in C_CoMix_names){
  CoMix_mat[[paste0('C_',s_comix)]] <- as.matrix(read.table(file.path(data_path,paste0(s_comix,'.csv'))))  
}

# comix waves
db_M_dates <- data.frame(date = c('2020-03-01','2020-03-15','2020-04-01',
                                  CoMix_mat$C_CoMix_dates[-nrow(CoMix_mat$C_CoMix_dates),1],
                                  c('2020-09-01','2020-09-15','2020-10-01','2020-10-20','2020-10-27','2020-11-03')))
db_M_dates$sim_day  <- sim_date2day(db_M_dates$date)
db_M_dates$is_comix <- db_M_dates$date %in% CoMix_mat$C_CoMix_dates[,1]
db_M_dates          <- db_M_dates[order(db_M_dates$sim_day),]
db_M_dates


################################################################ #
## HELP FUNCTIONS W.R.T SOCIAL CONTACTS ----
################################################################ #

get_M_change_day <- function(){
  return(db_M_dates$sim_day)
} 

get_M_change_date <- function(){
  return(sim_day2date(get_M_change_day()))
} 

is_M_CoMix <- function(){
  return(db_M_dates$is_comix)
} 

get_CoMix_change_date <- function(comix_waves = NA){
  comix_dates <- sim_day2date(get_M_change_day())[is_M_CoMix()]
  if(any(is.na(comix_waves))){
    comix_waves <- 1:length(comix_dates)
  }
  return(comix_dates[comix_waves])
} 

get_CoMix_change_day <- function(comix_waves = NA){
  comix_days <- get_M_change_day()[is_M_CoMix()]
  if(any(is.na(comix_waves))){
    comix_waves <- 1:length(comix_days)
  }
  return(comix_days[comix_waves])
}

get_add_change_day <- function(add_waves = NA){
  add_days <- get_M_change_day()[!is_M_CoMix()]
  if(any(is.na(add_waves))){
    add_waves <- 1:length(add_days)
  }
  return(add_days[add_waves])
} 

# get susceptibility matrices
get_susceptiblity_matrices <- function(CoMix_matrices, 
                                       CoMix_coef_mat, 
                                       Add_coef_mat,
                                       cnt_change_pt,cnt_change=0){  # scenarios
  
  M2010_asy = CoMix_matrices$C_CoMix2010_asy
  M2010_sy  = CoMix_matrices$C_CoMix2010_sy
  
  coef_April20 = Add_coef_mat[1,];
  MApril20_asy = coef_April20*CoMix_matrices$C_CoMix1_asy;
  MApril20_sy  = coef_April20*CoMix_matrices$C_CoMix1_sy;
  
  # JULY 2020 (1)
  coef_CoMix6 <- CoMix_coef_mat[6,];
  M6_asy = coef_CoMix6*CoMix_matrices$C_CoMix6_asy;
  M6_sy = coef_CoMix6*CoMix_matrices$C_CoMix6_sy;

  # JULY 2020 (2)
  coef_CoMix7 <- CoMix_coef_mat[7,];
  M7_asy = coef_CoMix7*CoMix_matrices$C_CoMix7_asy;
  M7_sy = coef_CoMix7*CoMix_matrices$C_CoMix7_sy;

  # SEPTEMBER 2020
  coef_Sept <- (Add_coef_mat[3,] + Add_coef_mat[4,] + Add_coef_mat[5,])/3;
  MSept_asy = coef_Sept*CoMix_matrices$C_CoMix8_asy;
  MSept_sy = coef_Sept*CoMix_matrices$C_CoMix8_sy;
  
  # OCTOBER 2020
  coef_Sept <- (Add_coef_mat[5,]);
  MOct20_asy = coef_Sept*CoMix_matrices$C_CoMix8_asy;
  MOct20_sy  = coef_Sept*CoMix_matrices$C_CoMix8_sy;
  
  # DECEMBER 2020
  MDec2020_asy = (CoMix_coef_mat[11,]*CoMix_matrices$C_CoMix11_asy + CoMix_coef_mat[12,]*CoMix_matrices$C_CoMix12_asy) / 2;
  MDec2020_sy =  (CoMix_coef_mat[11,]*CoMix_matrices$C_CoMix11_sy  + CoMix_coef_mat[12,]*CoMix_matrices$C_CoMix12_sy)  / 2;
  
  # combined behaviour "July 2020
  MJuly20_asy = (M6_asy + M7_asy)/2
  MJuly20_sy  = (M6_sy + M7_sy)/2
  
  # adjusted pre-pandemic: 2010 contacts with q-param of latest two waves
  coef_CoMix_x <- rowMeans(CoMix_coef_mat[nrow(CoMix_coef_mat)-(1:0),]);  
  #coef_CoMix_x <- CoMix_coef_mat[nrow(CoMix_coef_mat),];
  M2021q_asy <- coef_CoMix_x * CoMix_matrices$C_CoMix2010_asy
  M2021q_sy <- coef_CoMix_x * CoMix_matrices$C_CoMix2010_sy
  
  
  # aggregate all estimated transmission matrices (contacts * Q) into one array
  # set array indices for comix waves
  n_comix <- nrow(CoMix_coef_mat)
  n_other_march <- 2
  n_other_sept  <- 6
  n_other <- n_other_march + n_other_sept
  n_start <- 1
  n_all <- n_comix+n_other+n_start
  index_comix <- n_start + c((1:8)+n_other_march,(9:n_comix)+n_other)
  index_other <- n_start + c(1:n_other_march,(8+n_other_march)+1:n_other_sept)
  
  # initiate arrays
  M_estim_asy <- array(dim = c(10,10,n_all))
  M_estim_sy  <- array(dim = c(10,10,n_all))
  
  # include initial contact data
  M_estim_asy[,,1] = M2010_asy;
  M_estim_sy[,,1]  = M2010_sy;

  # include CoMix waves
  for(i_wave in 1:n_comix){
    coef_CoMix <- CoMix_coef_mat[i_wave,];
    M_estim_asy[,,index_comix[i_wave]] = coef_CoMix*CoMix_matrices[[paste0('C_CoMix',i_wave,'_asy')]];
    M_estim_sy[,,index_comix[i_wave]]  = coef_CoMix*CoMix_matrices[[paste0('C_CoMix',i_wave,'_sy')]];
  }
  
  # include separated CoMix matrices and q parameters - new (temp)
  M_coef  <- array(dim = c(n_all,10))
  M_brut_asy <- array(dim = c(10,10,n_all))
  M_brut_sy  <- array(dim = c(10,10,n_all))
   for(i_wave in 1:n_comix){
     M_coef[index_comix[i_wave],] <- CoMix_coef_mat[i_wave,];
     M_brut_asy[,,index_comix[i_wave]] = CoMix_matrices[[paste0('C_CoMix',i_wave,'_asy')]];
     M_brut_sy[,,index_comix[i_wave]]  = CoMix_matrices[[paste0('C_CoMix',i_wave,'_sy')]];
   }
  
  # include additional matrices for March 2020 lockdown
  for(i_other in 1:n_other_march){
    coef_other <- Add_coef_mat[i_other,];
    M_estim_asy[,,index_other[i_other]] = coef_other*CoMix_matrices[['C_CoMix1_asy']];
    M_estim_sy[,,index_other[i_other]]  = coef_other*CoMix_matrices[['C_CoMix1_sy']];
  }
  
  # include additional waves for September - November 2020 
  for(i_other in (n_other_march+1):n_other){
    coef_other <- Add_coef_mat[i_other,];
    M_estim_asy[,,index_other[i_other]] = coef_other*CoMix_matrices[['C_CoMix8_asy']];
    M_estim_sy[,,index_other[i_other]]  = coef_other*CoMix_matrices[['C_CoMix8_sy']];
  }
  
  
  # Proxy for summer 22 using summer relative variation from waves 26-29 compared to waves 31-34
 # coefSummer21_asy  <- (M_estim_asy[,,index_comix[26]] + M_estim_asy[,,index_comix[27]] + M_estim_asy[,,index_comix[28]] + M_estim_asy[,,index_comix[29]])/(M_estim_asy[,,index_comix[31]] + M_estim_asy[,,index_comix[32]] + M_estim_asy[,,index_comix[33]] + M_estim_asy[,,index_comix[34]])
#  coefSummer21_sy  <-(M_estim_sy[,,index_comix[26]] + M_estim_sy[,,index_comix[27]] + M_estim_sy[,,index_comix[28]] + M_estim_sy[,,index_comix[29]])/(M_estim_sy[,,index_comix[31]] + M_estim_sy[,,index_comix[32]] + M_estim_sy[,,index_comix[33]] + M_estim_sy[,,index_comix[34]])
#  coefSummer21_asy[which(!is.finite(coefSummer21_asy))] <- 0
  #coefSummer21_sy[which(!is.finite(coefSummer21_sy))] <- 0
  # MSummer21_asy <- M_estim_asy[,,index_comix[46]]* coefSummer21_asy
  # MSummer21_sy <- M_estim_sy[,,index_comix[46]]*coefSummer21_sy
  
  #proxy for long term future contact rates, using jan-Aug 2022 as ref (omicron)
  #!! this is slow !!!
  # MFuture_asy <- 0
  # MFuture_sy <- 0
  # if (length(cnt_change_pt) >= 5 && cnt_change[5]>0){
  #   MFuture_asy <- array(dim = c(10,10,sim_date2day("2032-10-01")))
  #   MFuture_sy  <- array(dim = c(10,10,sim_date2day("2032-10-01")))
  #   for(calendar_time in sim_date2day("2022-09-01"):sim_date2day("2023-12-31")){
  #   M_asy <- 0*M_estim_asy[,,10]
  #   M_sy <- 0*M_estim_sy[,,10]
  #   for(days in 0:29){
  #     reference_date <- sim_day2date(calendar_time)-days
  #     year(reference_date) <- 2022
  #     reference_time <- sim_date2day(reference_date)
  #     if(month(reference_date)>= 9){  #If sept-dec, use contact from June-march backwards
  #       reference_time <- sim_date2day("2022-06-30") - (reference_time - sim_date2day("2022-09-01"))
  #     } 
  #     lastyearmat <- tail(which(get_M_change_day()<=reference_time),n=1)
  #     M_asy <- M_asy + M_estim_asy[,,lastyearmat]
  #     M_sy <- M_sy + M_estim_sy[,,lastyearmat]
  #   }
  #   MFuture_asy[,,calendar_time] <- M_asy/30
  #   MFuture_sy[,,calendar_time] <- M_sy/30
  #   }
  #   for(calendar_time in sim_date2day("2024-01-01"):sim_date2day("2032-10-01")){
  #     reference_date <- sim_day2date(calendar_time)
  #     if(day(reference_date)== 29){
  #       reference_date <- reference_date-1
  #     }
  #     year(reference_date) <- 2023
  #     reference_time <- sim_date2day(reference_date)
  #     MFuture_asy[,,calendar_time] <- MFuture_asy[,,reference_time]
  #     MFuture_sy[,,calendar_time] <- MFuture_sy[,,reference_time]
  #   }
  # }
   
  # aggregate changepoints for estimated matrices
  sel_changepoints   <- 2:n_all
  M_estim_cp         <- rep(c(get_M_change_day()[sel_changepoints]),each=6) + 0:5
  M_estim_cp_start   <- get_M_change_day()[sel_changepoints]
  M_estim_cp_end     <- M_estim_cp_start + 5
  M_estim_wave       <- rep((1:(n_comix+n_other)+1),each=6)
  
  # aggregate all changepoints (including first lockdown and for scenarios)
  M_change_days      <- unique(rep(c(M_estim_cp,cnt_change_pt),each=6) + 0:5)
  
  #warning : from 2022-08-31, CoMix matrices are duplicated using past year, with new CoMix 48 = CoMix 30
  
  return(list(M2010_asy = M2010_asy, M2010_sy = M2010_sy,
              
              MApril20_asy=MApril20_asy,MApril20_sy=MApril20_sy,
              MSept_asy=MSept_asy,MSept_sy=MSept_sy,
              MOct20_asy=MSept_asy,MOct20_sy=MSept_sy,
              MJuly20_asy=MJuly20_asy,MJuly20_sy=MJuly20_sy,
              M2021q_asy=M2021q_asy,M2021q_sy=M2021q_sy,
              MDec2020_asy=MDec2020_asy,MDec2020_sy=MDec2020_sy,
              
              M_estim_asy=M_estim_asy,M_estim_sy=M_estim_sy,
              M_estim_cp=M_estim_cp,M_estim_wave=M_estim_wave,
              M_estim_cp_start=M_estim_cp_start,
              M_estim_cp_end=M_estim_cp_end,
              
              #      MSummer21_asy=MSummer21_asy,MSummer21_sy=MSummer21_sy,
              #MFuture_asy = MFuture_asy, MFuture_sy = MFuture_sy,
              M_brut_asy=M_brut_asy,M_brut_sy=M_brut_sy,M_coef=M_coef,
            
              M_change_days = M_change_days))
}

get_contact_details <- function(calendar_time,step_index,
                                M_list_all,
                                M_asy,M_sy,
                                cnt_change_pt,cnt_change,nfull,nfull_vac_reinf,nfull_vac_reinf_oldvoc){
      # if(1==0){ #temporary removed -> to be fixed
     if( !(length(cnt_change_pt) >= 5 && calendar_time >= cnt_change_pt[5] && cnt_change[5]>0) && !round(calendar_time) %in% M_list_all$M_change_days){
      # print("warning: change days may be not correct")
     return(list(M_asy=M_asy,M_sy=M_sy))
   } else
   
    ## Social contact data  ----
    ## based on CoMix waves and proxy for lockdown in March 2020 and September-November 2020
    if(any(calendar_time >= M_list_all$M_estim_cp_start &
            calendar_time <= M_list_all$M_estim_cp_end)){

      c_wave      <- M_list_all$M_estim_wave[M_list_all$M_estim_cp == round(calendar_time)]
      wCoMix_prev <- 1 - (calendar_time - min(M_list_all$M_estim_cp[M_list_all$M_estim_wave == c_wave]))/5
      wCoMix_new  <- 1 - wCoMix_prev;
      
      M_asy       <- wCoMix_prev*M_list_all$M_estim_asy[,,c_wave-1] + wCoMix_new*M_list_all$M_estim_asy[,,c_wave]
      M_sy        <- wCoMix_prev*M_list_all$M_estim_sy[,,c_wave-1] + wCoMix_new*M_list_all$M_estim_sy[,,c_wave]
    } 
    
  ## Scenarios ----
  ##---------------------------------------------------------------- -
  ## 1. Re-use previous dynamics (last wave, child)
  ## hub: re-use June 2022 !!!progressively during the month !!!
  if (length(cnt_change_pt) >= 1 && calendar_time >= cnt_change_pt[1] & calendar_time <= (cnt_change_pt[1]+30) && cnt_change[1]>0){
   #if (length(cnt_change_pt) >= 1 && calendar_time >= cnt_change_pt[1] & calendar_time <= (cnt_change_pt[1]+5) && cnt_change[1]>0){
      
    c_wave      <- max(M_list_all$M_estim_wave)
    wCoMix_prev <- 1 - (calendar_time - cnt_change_pt[1])/30
    #wCoMix_prev <- 1 - (calendar_time - cnt_change_pt[1])/5
    wCoMix_new  <- 1 - wCoMix_prev;
    #cnt_change_matrix <- M_list_all$M_estim_asy[,,c_wave]*0+1
    #cnt_change_matrix[1:2,1:2] <- cnt_change[1] # child factor
  #  M_asy       <- wCoMix_prev*M_list_all$M_estim_asy[,,c_wave] + wCoMix_new*M_list_all$M_estim_asy[,,c_wave] *cnt_change_matrix
   # M_sy        <- wCoMix_prev*M_list_all$M_estim_sy[,,c_wave] + wCoMix_new*M_list_all$M_estim_sy[,,c_wave] *cnt_change_matrix
    
    # combined behaviour summer 2021
    #MSummer21_asy = (M_list_all$M_estim_asy[,,26] + M_list_all$M_estim_asy[,,27] + M_list_all$M_estim_asy[,,28]  + M_list_all$M_estim_asy[,,29] )/4
    #  MSummer21_sy  = (M_list_all$M_estim_sy[,,26] + M_list_all$M_estim_sy[,,27] + M_list_all$M_estim_sy[,,28]  + M_list_all$M_estim_sy[,,29] )/4
    M_asy       <- wCoMix_prev*M_list_all$M_estim_asy[,,c_wave] + wCoMix_new*M_list_all$M_estim_asy[,,c_wave-1] 
    M_sy        <- wCoMix_prev*M_list_all$M_estim_sy[,,c_wave] + wCoMix_new*M_list_all$M_estim_sy[,,c_wave-1]

  } 
  
  ## 2. Re-use previous dynamics (last wave, adult)
  ### note: possible extension of change_point 1 
  ## hub: out of average summer 2021 
  if (length(cnt_change_pt) >= 2 && calendar_time >= cnt_change_pt[2] & calendar_time <= (cnt_change_pt[2]+5) && cnt_change[2]>0){

    c_wave      <- max(M_list_all$M_estim_wave)
    wCoMix_prev <- 1 - (calendar_time - cnt_change_pt[2])/5
    wCoMix_new  <- 1 - wCoMix_prev;
    
    cnt_change1 <- ifelse(cnt_change[1]>0,cnt_change[1],1)
    M_asy       <- wCoMix_prev*M_list_all$MSummer21_asy + wCoMix_new*M_list_all$M_estim_asy[,,c_wave]
    M_sy        <- wCoMix_prev*M_list_all$MSummer21_sy + wCoMix_new*M_list_all$M_estim_sy[,,c_wave]
    
  #  cnt_change_matrix <- M_list_all$M_estim_asy[,,c_wave]*0+1 # start from dummy with all factor 1
   # cnt_change_matrix[1:2,1:2] <- cnt_change1     # re-use child factor (cnt_change 1)
    #cnt_change_matrix[3:7,3:7] <- cnt_change[2]   # adult factor

   # M_asy       <- wCoMix_prev*M_list_all$M_estim_asy[,,c_wave] + wCoMix_new*M_list_all$M_estim_asy[,,c_wave] *cnt_change_matrix
    #M_sy        <- wCoMix_prev*M_list_all$M_estim_sy[,,c_wave] + wCoMix_new*M_list_all$M_estim_sy[,,c_wave] *cnt_change_matrix
  } 

  ### 3. goal: increase of latest contacts - q
  if (length(cnt_change_pt) >= 3 && calendar_time >= cnt_change_pt[3] & calendar_time <= (cnt_change_pt[3]+5) && cnt_change[3]!=0){
    
  #  c_wave      <- max(M_list_all$M_estim_wave)
   # wCoMix_prev <- 1 - (calendar_time - cnt_change_pt[3])/5
   # wCoMix_new  <- 1 - wCoMix_prev;
    
    # get average behaviour matrices
   # num_waves <- floor(cnt_change[3])
    #M_asy_new <- rowMeans(M_list_all$M_estim_asy[,,c_wave-(1:num_waves)], dims = 2)
  #  M_sy_new  <- rowMeans(M_list_all$M_estim_sy[,,c_wave-(1:num_waves)], dims = 2)
    
    # use average behaviour
   # M_asy       <- wCoMix_prev*M_list_all$M_estim_asy[,,c_wave] + wCoMix_new*M_asy_new
  #  M_sy        <- wCoMix_prev*M_list_all$M_estim_sy[,,c_wave] + wCoMix_new*M_sy_new
    
    
    c_wave      <- max(M_list_all$M_estim_wave)
    wCoMix_prev <- 1 - (calendar_time - cnt_change_pt[3])/5
    wCoMix_new  <- 1 - wCoMix_prev;
   
    M_asy       <- wCoMix_prev*M_list_all$M_estim_asy[,,c_wave] + wCoMix_new*M_list_all$M_estim_asy[,,c_wave] *cnt_change[3]
    M_sy        <- wCoMix_prev*M_list_all$M_estim_sy[,,c_wave] + wCoMix_new*M_list_all$M_estim_sy[,,c_wave] *cnt_change[3]
    
  }
  
  ### goal: increase of latest contacts - q
  if (length(cnt_change_pt) >= 4 && calendar_time >= cnt_change_pt[4] & calendar_time <= (cnt_change_pt[4]+5) && cnt_change[4]!=0){
    
    c_wave      <- max(M_list_all$M_estim_wave)
    wCoMix_prev <- 1 - (calendar_time - cnt_change_pt[4])/5
    wCoMix_new  <- 1 - wCoMix_prev;
    
    M_asy       <- wCoMix_prev*M_list_all$M_estim_asy[,,c_wave] + wCoMix_new*M_list_all$M_estim_asy[,,c_wave] *cnt_change[4]
    M_sy        <- wCoMix_prev*M_list_all$M_estim_sy[,,c_wave] + wCoMix_new*M_list_all$M_estim_sy[,,c_wave] *cnt_change[4]
    
  }

  ### goal: increase in contacts - q - another date
  #pre-pandemic !!!
  # if (length(cnt_change_pt) >= 5 && calendar_time >= cnt_change_pt[5] & calendar_time <= (cnt_change_pt[5]+5) && cnt_change[5]>0){
  #   
  #   # # #pre-pandemic
  #    #M_asy = M_list_all$M2021q_asy ;
  #    #M_sy  = M_list_all$M2021q_sy ;
  # 
  #   c_wave      <- max(M_list_all$M_estim_wave)
  #   wCoMix_prev <- 1 - (calendar_time - cnt_change_pt[5])/5
  #   wCoMix_new  <- 1 - wCoMix_prev;
  #   
  # 
  #   M_asy       <- wCoMix_prev*M_list_all$M_estim_asy[,,c_wave] + wCoMix_new*M_list_all$M2021q_asy *cnt_change[5]
  #   M_sy        <- wCoMix_prev*M_list_all$M_estim_sy[,,c_wave] + wCoMix_new*M_list_all$M2021q_sy *cnt_change[5]
  # 
  # }
  
  # reuse of all last year matrices, but with smoother changes
  #new version using same q and previous matrices
  
  if (length(cnt_change_pt) >= 5 && calendar_time >= cnt_change_pt[5] && cnt_change[5]>0){
 #   M_asy <- M_list_all$MFuture_asy[,,calendar_time]
#    M_sy <- M_list_all$MFuture_sy[,,calendar_time]
   M_asy <- 0*M_list_all$M_estim_asy[,,10]
   M_sy <- 0*M_list_all$M_estim_sy[,,10]
   c_wave      <- max(M_list_all$M_estim_wave)
   #smoothcoef <- (3*M_list_all$M_coef[c_wave,] + 2*M_list_all$M_coef[c_wave-1,] + 1*M_list_all$M_coef[c_wave-2,])/6
   #smoothcoeflastyear <- (3*M_list_all$M_coef[c_wave-18,] + 2*M_list_all$M_coef[c_wave-1-18,] + 1*M_list_all$M_coef[c_wave-2-18,])/6
   smoothcoef <- M_list_all$M_coef[c_wave,] 
   smoothcoeflastyear <- M_list_all$M_coef[c_wave-18,] 
      for(days in 0:14){
     reference_time <- calendar_time -365-days
     lastyearmat <- tail(which(get_M_change_day()<=reference_time),n=1)
     #M_asy <- M_asy + M_list_all$M_brut_asy[,,lastyearmat]*smoothcoef
     #M_sy <- M_sy + M_list_all$M_brut_sy[,,lastyearmat]*smoothcoef
     #M_asy <- M_asy + M_list_all$M_estim_asy[,,lastyearmat]
     #M_sy <- M_sy + M_list_all$M_estim_sy[,,lastyearmat]
     changelastyearcoef <- M_list_all$M_coef[lastyearmat,] / smoothcoeflastyear
     changethisyearcoef <- changelastyearcoef * smoothcoef
     M_asy <- M_asy + M_list_all$M_brut_asy[,,lastyearmat]*changethisyearcoef
     M_sy <- M_sy + M_list_all$M_brut_sy[,,lastyearmat]*changethisyearcoef
      }
   M_asy <- M_asy/15
   M_sy <- M_sy/15
   #print(paste(sim_day2date(calendar_time),sim_day2date(reference_time)))
   
        if(calendar_time <= cnt_change_pt[5]+15){
        wCoMix_prev <- 1 - (calendar_time - cnt_change_pt[5])/15
       wCoMix_new  <- 1 - wCoMix_prev;
      M_asy       <- wCoMix_prev*M_list_all$M_estim_asy[,,c_wave] + wCoMix_new*M_asy
      M_sy        <- wCoMix_prev*M_list_all$M_estim_sy[,,c_wave] + wCoMix_new*M_sy
     }
   
  }
  
  ## VACCINE AND VOC PAPER
  if(0==1){
    ## 2. Re-use dynamics from early August 2020 (infinitly)
    if (length(cnt_change_pt) >= 2 && calendar_time >= cnt_change_pt[2] && cnt_change[2]>0){
      
      c_wave      <- max(which(get_M_change_day() == sim_date2day('2020-07-30')))
      wCoMix_prev <- pmax(0,1 - (calendar_time - cnt_change_pt[2])/5)
      wCoMix_new  <- 1 - wCoMix_prev;
      
      M_asy       <- wCoMix_prev*M_list_all$M_estim_asy[,,c_wave] + wCoMix_new*M_list_all$M_estim_asy[,,c_wave] *cnt_change[2]
      M_sy        <- wCoMix_prev*M_list_all$M_estim_sy[,,c_wave] + wCoMix_new*M_list_all$M_estim_sy[,,c_wave] *cnt_change[2]
      
    } 
    
    ## 2. Re-use dynamics from early August 2021 (infinitly)
    if (length(cnt_change_pt) >= 4 && calendar_time >= cnt_change_pt[4] && cnt_change[4]!=0){
      
      c_wave      <- max(which(get_M_change_day() == sim_date2day('2021-08-03')))
      wCoMix_prev <- pmax(0,1 - (calendar_time - cnt_change_pt[4])/5)
      wCoMix_new  <- 1 - wCoMix_prev;
      
      M_asy       <- wCoMix_prev*M_list_all$M_estim_asy[,,c_wave] + wCoMix_new*M_list_all$M_estim_asy[,,c_wave] *cnt_change[4]
      M_sy        <- wCoMix_prev*M_list_all$M_estim_sy[,,c_wave] + wCoMix_new*M_list_all$M_estim_sy[,,c_wave] *cnt_change[4]
      
    }
    
    ## 2. Re-use dynamics from September 2020 (infinitly)
    if (length(cnt_change_pt) >= 5 && calendar_time >= cnt_change_pt[5] && cnt_change[5]>0){
      
      # september 2020
      M_asy_new = M_list_all$MSept_asy ;
      M_sy_new  = M_list_all$MSept_sy ;
      
      # wave specific
      c_wave      <- max(which(get_M_change_day()<cnt_change_pt[5]))
      M_asy_prev = M_list_all$M_estim_asy[,,c_wave]
      M_sy_prev  = M_list_all$M_estim_sy[,,c_wave]
      
      wCoMix_prev <- pmax(0,1 - (calendar_time - cnt_change_pt[5])/5)
      wCoMix_new  <- 1 - wCoMix_prev;
      
      M_asy       <- wCoMix_prev*M_asy_prev + wCoMix_new*M_asy_new *cnt_change[5]
      M_sy        <- wCoMix_prev*M_sy_prev + wCoMix_new*M_sy_new *cnt_change[5]
      
    }
    
  }
  
  # if(calendar_time >sim_date2day("2021-05-01") && calendar_time < sim_date2day("2022-01-01")){
  #   #print(calendar_time)
  #   numat <- tail(which(get_M_change_day()<=calendar_time),n=1)
  #   icu_beds=sum(nfull[[step_index]][c_I_icu]+nfull[[step_index]][c_Ivoc_icu]+nfull_vac_reinf[[step_index]][c_I_icu]+nfull_vac_reinf[[step_index]][c_Ivoc_icu]+nfull_vac_reinf_oldvoc[[step_index]][c_I_icu]+nfull_vac_reinf_oldvoc[[step_index]][c_Ivoc_icu])
  #   new_hospi=sum(nfull[[step_index]][c_new_hosp_total]+nfull[[step_index]][c_new_hosp_total_voc]+nfull_vac_reinf[[step_index]][c_new_hosp_total]+nfull_vac_reinf[[step_index]][c_new_hosp_total_voc]+nfull_vac_reinf_oldvoc[[step_index]][c_new_hosp_total]+nfull_vac_reinf_oldvoc[[step_index]][c_new_hosp_total_voc])
  #   if((new_hospi > 64 && new_hospi <150) || (icu_beds>300 &icu_beds<=500)) {
  #     print("orange")
  #     M_sy=0.5*M_list_all$M_estim_sy[,,numat]
  #     M_asy=0.5*M_list_all$M_estim_asy[,,numat]
  #   }
  #   else if(new_hospi>149 || icu_beds > 500 ){
  #     print("red")
  #     M_sy=0.3*M_list_all$M_estim_sy[,,numat]
  #     M_asy=0.3*M_list_all$M_estim_asy[,,numat]
  #   }
  # }
  # if(calendar_time >sim_date2day("2022-01-01") && calendar_time < sim_date2day("2023-04-01")){
  #   #print(calendar_time)
  #   numat <- tail(which(get_M_change_day()<=calendar_time),n=1)
  #   icu_beds=sum(nfull[[step_index]][c_I_icu]+nfull[[step_index]][c_Ivoc_icu]+nfull_vac_reinf[[step_index]][c_I_icu]+nfull_vac_reinf[[step_index]][c_Ivoc_icu]+nfull_vac_reinf_oldvoc[[step_index]][c_I_icu]+nfull_vac_reinf_oldvoc[[step_index]][c_Ivoc_icu])
  #   new_hospi=sum(nfull[[step_index]][c_new_hosp_total]+nfull[[step_index]][c_new_hosp_total_voc]+nfull_vac_reinf[[step_index]][c_new_hosp_total]+nfull_vac_reinf[[step_index]][c_new_hosp_total_voc]+nfull_vac_reinf_oldvoc[[step_index]][c_new_hosp_total]+nfull_vac_reinf_oldvoc[[step_index]][c_new_hosp_total_voc])
  #   if((new_hospi > 64 && new_hospi <150) || (icu_beds>300 &icu_beds<=500)) {
  #     print("orange")
  #     M_sy=0.3*M_list_all$M_estim_sy[,,numat]
  #     M_asy=0.3*M_list_all$M_estim_asy[,,numat]
  #   }
  #   else if(new_hospi>149 || icu_beds > 500 ){
  #     print("red")
  #     M_sy=0.1*M_list_all$M_estim_sy[,,numat]
  #     M_asy=0.1*M_list_all$M_estim_asy[,,numat]
  #   }
  # }
  
  
  # include safety checks
  if(length(M_asy)==0){
    stop(paste('STOP: M_asy is NULL on day', calendar_time))
  }
  if(length(M_sy)==0){
    stop(paste('STOP: M_sy is NULL on day', calendar_time))
  }
  if(all(M_sy == M_asy)){
    stop(paste('STOP: M_sy == M_asy on day', calendar_time))
  }
  
  return(list(M_asy=M_asy,M_sy=M_sy))
}



