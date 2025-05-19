<base target="_blank">

# Stochastic Compartmental Model for SAR-COV-2

This project is conceived by members of SIMID group during the COVID19 pandemic. 

This version was used to gerenerate results for the manuscript: **The Impact of COVID-19 Vaccination on Mortality and Stringency Measures in Belgium: A Retrospective Analysis**

All data needed to generate the scenarios of the manuscript are contained in the data folder. Those scenarios are based on the calibration contained in the file data/config/MCMCmulti_20230810_d1213_e2_i1000_n10_p10_crit7_hubR5foronly_pessi_belgium_foronly.csv.
We are not the owner of the full unaggregate data used for the calibration. Those data can be requested from Sciensano via https://epidata.sciensano.be/epistat/dashboard/#covid.

**Project directory structure**

| Folder | Content |
|---|---|
| data    | folder with all input data and configuration files. Temporary files (e.g. reported cases or hospital admissions from Sciensano) are stored in "download"" folders, which are not included in the git repository.|
| doc   | folder with some documents (limited)|
| R   | folder with main R files to run the stochastic model, estimate parameters and plot results|
| R_processing| secondary R files to pre-process vaccine uptake or for post-processing previous results |

**Main project files**

| File | Content |
|---|---|
| R/main_vaccination.R |   script to load all principal data and functions for the stochastic model |
| R/projections_vaccination.R | script to run the stochastic model with different social contact assumptions and vaccine uptake files |
| R/plot_sm_combined.R | (re)create simulation output figures and/or combine different simulations |
| R/calibration_mcmc.R | script to start parameter estimations using MCMC |
| R_processing/vaccine_uptake_scenario.R | script to generate vaccine uptake files |
| R/lib_model_core.R |       contains the main function of the stochastic transmission model and the log-likelihood function |
| R/workbench_testing.R |   script to run different versions of the stochastic transmission model and to check the latest model output with previously stored results. |

Other files in the "R" folder contain help function for data manipulation, to visualize simulation output, or for the parameter estimation with MCMC.


**Contributors** (in alphabetical order):

* Steven Abrams 
* Philippe Beutels
* Maikel Bosschaert
* Christel Faes
* Nicolas Franco
* Nicolas Herman
* Niel Hens
* James Wambua
* Lander Willem
* The SIMID COVID-19 modelling team.
