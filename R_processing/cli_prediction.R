########################################################################### #
# This file is part of the Stochastic Compartmental Model for SARS-COV-2 
# transmission in Belgium, conceived by members of SIMID group during the 
# COVID19 pandemic.
#
# This file contains commands to run the stochastic model projection script 
# in parallel and with pre-defined parameters. Also command to run the model 
# exploration by region are present.
#
# Copyright 2021, SIMID, University of Antwerp & Hasselt University                                        
########################################################################### #

# this if-clause is to prevent that all commands are executed together!
if(0 == 1){

  ################################################################# #
  # MODEL OUTPUT EXPLORATION SCRIPTS BY REGION ----  
  ################################################################# #

  # output directory specified in the script
  system('Rscript R/plot_scm_combined.R belgium &')
  system('Rscript R/plot_scm_combined.R brussels &')
  system('Rscript R/plot_scm_combined.R flanders &')
  system('Rscript R/plot_scm_combined.R wallonia &')
  
  # latest output
  system('Rscript R/plot_scm_combined.R belgium TRUE &')
  system('Rscript R/plot_scm_combined.R brussels TRUE &')
  system('Rscript R/plot_scm_combined.R flanders TRUE &')
  system('Rscript R/plot_scm_combined.R wallonia TRUE &')
  
  # default social contact behaviour
  system('Rscript R/projections_vaccination.R belgium &')
  system('Rscript R/projections_vaccination.R brussels &')
  system('Rscript R/projections_vaccination.R flanders &')
  system('Rscript R/projections_vaccination.R wallonia &')
  
  system('Rscript R/projections_vaccination.R belgium 0_0_0_0_0 &') # default (previous)
  system('Rscript R/projections_vaccination.R belgium 0_0_1_0_0 &') # 
  system('Rscript R/projections_vaccination.R belgium 0_0_0_1.3_0 &')
  system('Rscript R/projections_vaccination.R belgium 0_0_0_1.5_0 &')
  system('Rscript R/projections_vaccination.R belgium 0_0_0_1.7_0 &')
  system('Rscript R/projections_vaccination.R belgium 0_0_0_0_1 &')
  
  system('Rscript R/projections_vaccination.R belgium 0_1.3_0_0_0 &')
  system('Rscript R/projections_vaccination.R belgium 0_1.5_0_0_0 &')
  system('Rscript R/projections_vaccination.R belgium 0_1.7_0_0_0 &')

  system('Rscript R/projections_vaccination.R belgium 99_0_0_0_0 &')
  system('Rscript R/projections_vaccination.R belgium 999_0_0_0_0 &')
  
  system('Rscript R/projections_vaccination.R belgium 0_0_0_0_0 &')
  system('Rscript R/projections_vaccination.R belgium 0_0_1.75_0_0 &')
  system('Rscript R/projections_vaccination.R belgium 0_0_2_0_0 &')
  
  ################################################################# #
  # august 22
  ################################################################# #
  
  system('Rscript R/projections_EUhub.R belgium 0_0_0_0_0 &')
  system('Rscript R/projections_EUhub.R belgium 1_0_0_0_0 &')
  
    ################################################################# #
  # July 10
  ################################################################# #
  
  system('Rscript R/projections_EUhub.R belgium 0_0_0_0_0 &')
  system('Rscript R/projections_EUhub.R belgium 0_0_1.5_0_0 &')
  system('Rscript R/projections_EUhub.R belgium 1_1_0_0_0 &')
  
  ################################################################# #
  # SEASONALITY
  ################################################################# #
  
  system('Rscript R/projections_vaccination.R belgium 0_0_0_0_0 &')
  system('Rscript R/projections_vaccination.R belgium 2_1_0_0_0 &')
  system('Rscript R/projections_vaccination.R belgium 0_0_0.8_0.5_1 &')
  system('Rscript R/projections_vaccination.R belgium 0_0_0.5_0_0 &')
  system('Rscript R/projections_vaccination.R belgium 0_0_0.8_1_0 &')
  system('Rscript R/projections_vaccination.R belgium 0_0_0.5_1_0 &')
  

  ################################################################# #
  # MODEL PROJECTION SCRIPTS FOR OMICRON ON FEBRUARY 7th ----  
  ################################################################# #
  
  system('Rscript R/projections_vaccination.R belgium 0_0_0_0_0 HIGH_25 &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_0_0_0 LOW_25 &') 
  
 # system('Rscript R/projections_vaccination.R belgium 0_0_0_0_1.5 HIGH_25 &') 
  #system('Rscript R/projections_vaccination.R belgium 0_0_0_0_1.5 LOW_25 &') 
  
  system('Rscript R/projections_vaccination.R belgium 0_0_0_0_2 HIGH_25 &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_0_0_2 LOW_25 &') 
  
  system('Rscript R/projections_vaccination.R belgium 0_0_0_0_3 HIGH_25 &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_0_0_3 LOW_25 &') 
  
  system('Rscript R/projections_vaccination.R belgium 0_0_0_2_0 HIGH_25 &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_0_2_0 LOW_25 &') 
  
  system('Rscript R/projections_vaccination.R belgium 0_0_0_3_0 HIGH_25 &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_0_3_0 LOW_25 &')
  
  system('Rscript R/projections_vaccination.R belgium 0_0_2_0_0 HIGH_25 &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_2_0_0 LOW_25 &') 
  
  system('Rscript R/projections_vaccination.R belgium 0_0_3_0_0 HIGH_25 &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_3_0_0 LOW_25 &')
  

    
  ################################################################# #
  # MODEL PROJECTION SCRIPTS FOR OMICRON ON JANUARY 5th ----  
  ################################################################# #
  
  system('Rscript R/projections_vaccination.R belgium 0_0_0_0_0 HIGH_025 &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_0_0_0 LOW_025 &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_0_0_0 HIGH_025 &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_0_0_0 LOW_25 &') 
  
  # risk behaviour scenarios
  system('Rscript R/projections_vaccination.R belgium 0_0_0 HIGH_025 &') 
  system('Rscript R/projections_vaccination.R belgium 2_0_0 HIGH_025 &') 
  system('Rscript R/projections_vaccination.R belgium 0_0.8_0 HIGH_025 &') 
  system('Rscript R/projections_vaccination.R belgium 1_0.8_0 HIGH_025 &') 
  system('Rscript R/projections_vaccination.R belgium 2_0.7_0 HIGH_025 &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_4 HIGH_025 &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_7 HIGH_025 &') 

  system('Rscript R/projections_vaccination.R belgium 0_0_3 HIGH_025 &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_3 HIGH_050 &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_3 LOW_050 &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_3 LOW_025 &') 
  
  system('Rscript R/projections_vaccination.R belgium 0_0_0 HIGH_025 &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_0 HIGH_050 &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_0 LOW_050 &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_0 LOW_025 &') 
    
  
  ################################################################# #
  # MODEL PROJECTION SCRIPTS FOR OMICRON EXPLORATION ----  
  ################################################################# #
  
  # 05 Jan report
  system('Rscript R/projections_vaccination.R belgium 0_0_0_0_0 0_0.25_0.03_HIGH &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_0_0_0 0_0.33_0.03_HIGH &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_0_0_0 0_0.50_0.03_HIGH &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_0_0_1 0_0.25_0.03_HIGH &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_0_0_1 0_0.33_0.03_HIGH &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_0_0_1 0_0.50_0.03_HIGH &')   
  

  # transmission advantage
  system('Rscript R/projections_vaccination.R belgium 0_0_1_1_0 1.50_0.5_0.03_LOW &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_1_1_0 1.75_0.5_0.03_LOW &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_1_1_0 2.00_0.5_0.03_LOW &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_1_1_0 2.25_0.5_0.03_LOW &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_1_1_0 2.50_0.5_0.03_LOW &') 

  
  #HIGH VE
  system('Rscript R/projections_vaccination.R belgium 0_0_1_1_0 3_1_0.03_HIGH &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_1_1_0 3_0.75_0.03_HIGH &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_1_1_0 3_0.5_0.03_HIGH &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_1_1_0 3_0.25_0.03_HIGH &') 
  
  #LOW VE  
  system('Rscript R/projections_vaccination.R belgium 0_0_1_1_0 2.5_1_0.03_LOW &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_1_1_0 2.5_0.75_0.03_LOW &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_1_1_0 2.5_0.5_0.03_LOW &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_1_1_0 2.5_0.25_0.03_LOW &') 
  
  #DEFAULT (no evasion/decrease VE) 
  system('Rscript R/projections_vaccination.R belgium 0_0_1_1_0 4.3_1_0_DELTA &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_1_1_0 4.3_0.75_0_DELTA &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_1_1_0 4.3_0.5_0_DELTA &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_1_1_0 4.3_0.25_0_DELTA &') 
  
  
  #NO VE
  system('Rscript R/projections_vaccination.R belgium 0_0_1_1_0 2_1_1.0 &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_1_1_0 1.5_1_1.0 &') 
  
  #HIGH VE
  system('Rscript R/projections_vaccination.R belgium 0_0_0_0_0 3_1_0.03 &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_0_0_0 3_0.75_0.03 &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_0_0_0 3_0.5_0.03 &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_0_0_0 3_0.25_0.03 &') 
  
  #LOW VE  
  system('Rscript R/projections_vaccination.R belgium 0_0_0_0_0 2.5_1_0.03 &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_0_0_0 2.5_0.75_0.03 &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_0_0_0 2.5_0.5_0.03 &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_0_0_0 2.5_0.25_0.03 &') 
  
  #DEFAULT (no evasion/decrease VE, nor holiday) 
  system('Rscript R/projections_vaccination.R belgium 0_0_0_0_0 4.3_1_0.03 &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_0_0_0 4.3_0.75_0.03 &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_0_0_0 4.3_0.5_0.03 &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_0_0_0 4.3_0.25_0.03 &') 

  system('Rscript R/projections_vaccination.R belgium 0_0_0_0_0 4_1_0 &') 
  system('Rscript R/projections_vaccination.R belgium 0_0_0_0_0 4.3_1_0 &') 

  
  
  ################################################################# #
  # MODEL PROJECTION SCRIPTS FOR VACCINE UPTAKE ----  
  ################################################################# #

  system('Rscript R/projections_vaccination.R belgium 0_0_0_0_0 &') # default (previous)
  system('Rscript R/projections_vaccination.R belgium 1_1_0_0_0 &') # default (current notation)
  system('Rscript R/projections_vaccination.R belgium 1_1_1_0_0 &') # incl December
  system('Rscript R/projections_vaccination.R belgium 1_1_1_1_0 &') # incl December and January
  system('Rscript R/projections_vaccination.R belgium 1_1_1_1_1 &') # incl 5th wave
  
  system('Rscript R/projections_vaccination.R belgium 0_0_1_0_1 &') # incl 5th wave
  system('Rscript R/projections_vaccination.R belgium 0_0_1_1_0 &') # incl 5th wave
  
  system('Rscript R/projections_vaccination.R belgium 0_0_1_0_0 &') # 
  
  
  system('Rscript R/projections_vaccination.R belgium 0_0_1_1_0 &') # incl 5th wave
  system('Rscript R/projections_vaccination.R belgium 0_0_1_0_1') # incl 5th wave
  system('Rscript R/projections_vaccination.R belgium 0_0_1_0_1.2 &') # incl 5th wave
  
  ################################################################# #
  # MODEL PROJECTION SCRIPTS BY REGION ----  
  ################################################################# #


  
  # default social contact behaviour
  system('Rscript R/projections_vaccination.R belgium 0_0_0_0 &')
  system('Rscript R/projections_vaccination.R brussels 0_0_0_0 &')
  system('Rscript R/projections_vaccination.R flanders 0_0_0_0 &')
  system('Rscript R/projections_vaccination.R wallonia 0_0_0_0 &')
  
  # with factor 120% for last social contact option
  system('Rscript R/projections_vaccination.R belgium 0_0_0_1.2 &')
  system('Rscript R/projections_vaccination.R brussels 0_0_0_1.2 &')
  system('Rscript R/projections_vaccination.R flanders 0_0_0_1.2 &')
  system('Rscript R/projections_vaccination.R wallonia 0_0_0_1.2 &')
  
  # etc...
  system('Rscript R/projections_vaccination.R belgium 0_0_0_1.3 &')
  system('Rscript R/projections_vaccination.R brussels 0_0_0_1.3 &')
  system('Rscript R/projections_vaccination.R flanders 0_0_0_1.3 &')
  system('Rscript R/projections_vaccination.R wallonia 0_0_0_1.3 &')
  
  system('Rscript R/projections_vaccination.R belgium 0_0_0_1.4 &')
  system('Rscript R/projections_vaccination.R brussels 0_0_0_1.4 &')
  system('Rscript R/projections_vaccination.R flanders 0_0_0_1.4 &')
  system('Rscript R/projections_vaccination.R wallonia 0_0_0_1.4 &')
  
  system('Rscript R/projections_vaccination.R belgium 0_0_0_1.6 &')
  system('Rscript R/projections_vaccination.R brussels 0_0_0_1.6 &')
  system('Rscript R/projections_vaccination.R flanders 0_0_0_1.6 &')
  system('Rscript R/projections_vaccination.R wallonia 0_0_0_1.6 &')
  
  system('Rscript R/projections_vaccination.R belgium 0_0_1.2_0 &')
  system('Rscript R/projections_vaccination.R brussels 0_0_1.2_0 &')
  system('Rscript R/projections_vaccination.R flanders 0_0_1.2_0 &')
  system('Rscript R/projections_vaccination.R wallonia 0_0_1.2_0 &')
  
  # second social contact behaviour parameter
  system('Rscript R/projections_vaccination.R belgium 0_1_0_0 &')
  system('Rscript R/projections_vaccination.R brussels 0_1_0_0 &')
  system('Rscript R/projections_vaccination.R flanders 0_1_0_0 &')
  system('Rscript R/projections_vaccination.R wallonia 0_1_0_0 &')
  
  # third social contact behaviour parameter
  system('Rscript R/projections_vaccination.R belgium 0_0_1_0 &')
  system('Rscript R/projections_vaccination.R brussels 0_0_1_0 &')
  system('Rscript R/projections_vaccination.R flanders 0_0_1_0 &')
  system('Rscript R/projections_vaccination.R wallonia 0_0_1_0 &')

  # model runs for Belgium (universal steps)
  system('Rscript R/projections_vaccination.R belgium 0_0_0_0 &')
  system('Rscript R/projections_vaccination.R belgium 0_0_0_1.2 &')
  system('Rscript R/projections_vaccination.R belgium 0_0_0_1.3 &')
  system('Rscript R/projections_vaccination.R belgium 0_0_0_1.4 &')
  system('Rscript R/projections_vaccination.R belgium 0_0_0_1.6')
  system('sleep 2')
  system('Rscript R/plot_scm_combined.R belgium TRUE &')

  # model runs for Belgium (specific steps for TechNote v20211116)
  # with config file MCMCmulti_20211115_d622_e50_i400_n10_p10_crit1_62269_belgium
  system('Rscript R/projections_vaccination.R belgium 0_0_0_0 &')
  system('Rscript R/projections_vaccination.R belgium 0_0_1.1_1.1 &')
  system('Rscript R/projections_vaccination.R belgium 0_0_1.14_1.1 &')
  system('Rscript R/projections_vaccination.R belgium 0_0_0.9_1 &')
  system('Rscript R/projections_vaccination.R belgium 0_0_1.25_1.15 &')
  system('Rscript R/plot_scm_combined.R belgium TRUE &')
  
  # model runs for Belgium (specific steps for TechNote v20211202)
  # with config file MCMCmulti_20211115_d622_e50_i400_n10_p10_crit1_62269_belgium
  system('Rscript R/projections_vaccination.R belgium 0_0_1.14_1.1 &')
  system('Rscript R/projections_vaccination.R belgium 2.5_0_1.14_1.1 &')
  system('Rscript R/plot_scm_combined.R belgium TRUE &')
  
  
  # model runs for Belgium (universal steps)
  system('Rscript R/projections_vaccination.R brussels 0_0_0_0 &')
  system('Rscript R/projections_vaccination.R brussels 0_0_0_1.3 &')
  system('Rscript R/projections_vaccination.R brussels 0_0_0_1.4 &')
  system('Rscript R/projections_vaccination.R brussels 0_0_0_1.6 &')
  system('Rscript R/plot_scm_combined.R brussels TRUE &')
  
}
