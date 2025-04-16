<base target="_blank">

# Stochastic Compartmental Model for SAR-COV-2

This project is conceived by members of SIMID group during the COVID19 pandemic to estimate and explore the effect of social distancing and vaccine uptake in Belgium. The project is hosted on [GitHub repository](https://github.com/lwillem/stochastic_model_BE.git).

Current Version: 25.0

**Scientific output** related to this model:

* [Modeling the early phase of the Belgian COVID-19 epidemic](https://doi.org/10.1016/j.epidem.2021.100449)
* [Technical Notes](https://www.simid.be/news/technical-note-covid19/)
* [RESTORE reports](https://covid-en-wetenschap.github.io/restore)

**Project directory structure**

| Folder | Content |
|---|---|
| cluster | files to run the stochastic model on the VSC clusters using SLURM and TORQUE|
| data    | folder with all input data and configuration files. Temporary files (e.g. reported cases or hospital admissions from Sciensano) are stored in "download"" folders, which are not included in the git repository.|
| doc   | folder with some documents (limited)|
| R   | folder with main R files to run the stochastic model, estimate parameters and plot results|
| R_processing| secondary R files to pre-process vaccine uptake or for post-processing previous results |
| unit_testing | folder contains reference model output, which is used in main_workbench to contrast the latest runs |

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

**Where to start?**

The files *R/workbench_testing.R* and *R/projections_wave1.R* seem most suited to start using and exploring the Stochastic Compartmental Model.

**Contributors** (in alphabetical order):

* Steven Abrams 
* Nicolas Franco
* Philippe Beutels
* Christel Faes
* Niel Hens
* James Wambua
* Lander Willem (lander.willem@uantwerpen.be)
* The SIMID COVID-19 modelling team.
