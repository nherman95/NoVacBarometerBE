########################################################################### #
# This file is part of the Stochastic Compartmental Model for SARS-COV-2 
# transmission in Belgium, conceived by members of SIMID group during the 
# COVID19 pandemic.
#
# This file is used to install all CRAN packages and dependencies 
#
# Copyright 2021, SIMID, University of Antwerp & Hasselt University                                        
########################################################################### #

# CRAN packages
install.packages('useful')        # compare.list
install.packages('EpiEstim')      # Rt calculations
install.packages('knitr')         # to convert markdown into pdf
install.packages('RColorBrewer')  # extra colors
install.packages('scales')        # to use transparant colors (cfr. credible intervals)
install.packages('data.table')    # to use the data.table format
install.packages('zoo')           # rollmean
install.packages('openxlsx')      # to read Excel files
install.packages('LaplacesDemon') # MCMC 
install.packages('poisbinom')     # poison binomial 
install.packages('socialmixr')    # social contact data analyses
install.packages('dplyr')   # for matrix manipulation
install.packages('MMWRweek')   # conversion official weeks/days


# to install the simid.rtools package from github for help functions on 
# parallel computing etc.
# For MAC OS: make sure that "XCode Command Line Tools" are installed on your system.
install.packages('devtools')
library(devtools)
devtools::install_github("lwillem/simid_rtools",force=F,quiet=T)
