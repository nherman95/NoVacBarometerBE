########################################################################### #
# This file is part of the Stochastic Compartmental Model for SARS-COV-2 
# transmission in Belgium, conceived by members of SIMID group during the 
# COVID19 pandemic.
#
# This file can be used to make sure the reference data is downloaded once
# before starting parallel MCMC chains. 
#
# This file can be excecuted with a command line argument. For example:
# - terminal: Rscript R/calibration_start.R belgium
# - terminal: Rscript R/calibration_start.R flanders &
# - within R: source('R/calibration_start.R')
# - within R: system('Rscript R/calibration_start.R brussels &')
#
# Copyright 2021, SIMID, University of Antwerp & Hasselt University                                        
########################################################################### #

args = commandArgs(trailingOnly=TRUE)

rm(list=ls()[ls()!='args'])

suppressPackageStartupMessages(library('useful')) # compare.list

# load functions and data ----
source('R/main_vaccination.R')

if(length(args)==1){
  sel_region <- args[[1]]
} else{
  sel_region <- get_region(1)
}

# (down)load reference data
be_ref_data <- get_latest_incidence_data(sel_region = sel_region)

# set logfile
cat(paste('REFERENCE DATA READY FOR:', sel_region),fill = T)

