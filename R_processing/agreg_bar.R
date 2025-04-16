# Load the necessary package
library(data.table)
# parse command line arguments
cl_args = commandArgs(trailingOnly=TRUE)
# clear workspace (except the command line arguments)
rm(list=ls()[ls()!='cl_args'])
par(cex=50)
get_arg_with_default <- function(arg_index, default_value) {
  if (length(cl_args) >= arg_index) {
    arg_value <- cl_args[arg_index]
    default_type <- typeof(default_value)
    #print(default_type)
    if (default_type == "double") {
      return(as.numeric(arg_value))
    } else if (default_type == "integer") {
      return(as.integer(arg_value))
    } else if (default_type == "logical") {
      return(as.logical(arg_value))
    } else {
      return(arg_value)
    }
  } else {
    return(default_value)
  }
}
add_polygon <- function(scenario_data, scenario_dates, col, alpha_val = 0.5, mean_lty = 1) {
  # Sélectionner les colonnes des limites inférieures et supérieures de l'intervalle de confiance
  id_lower <- grep('_LL', names(scenario_data))
  id_upper <- grep('_UL', names(scenario_data))
  
  # Supprimer les lignes avec des valeurs manquantes
  complete_cases <- complete.cases(scenario_data)
  scenario_data <- scenario_data[complete_cases, ]
  scenario_dates <- scenario_dates[complete_cases]
  
  # S'assurer que les dates sont triées
  sort_order <- order(scenario_dates)
  scenario_data <- scenario_data[sort_order, ]
  scenario_dates <- scenario_dates[sort_order]
  
  # Calculer les coordonnées pour le polygon
  polygon_x <- c(scenario_dates,rev(scenario_dates))
  print(polygon_x)
  polygon_y_lower <- c(scenario_data[, id_lower], rev(scenario_data[, id_upper]))  # Les limites inférieures
  polygon_y_upper <- c(scenario_data[, id_upper], rev(scenario_data[, id_lower]))  # Les limites supérieures
  # Ajouter le polygon à la parcelle
  polygon(x = polygon_x, y = polygon_y_lower, col = alpha(col, alpha_val), border = NA)
  #polygon(x = polygon_x, y = polygon_y_upper, col = alpha(col, alpha_val), border = NA)
  
  # Ajouter la ligne moyenne
  lines(x = scenario_dates, y = scenario_data[, "admin_mean"], lty = mean_lty, col = col, lwd = 3)
}
bool_vac<- FALSE
bool_bar<-FALSE
source('R/main_vaccination.R')

# Function to process a bar_summary file
process_bar_summary <- function(file) {
  # Read all lines from the file
  data <- readLines(file)
  
  # Extract parameters values from the first line, discard the first element and convert to numeric
  first_line <- data[1]
  parameters <- as.numeric(unlist(strsplit(first_line, ","))[-1])
  
  # Extract indicators names from the second line, discard the first element
  second_line <- data[2]
  indicator_names <- unlist(strsplit(second_line, ","))[-1]
  
  # Extract median values from the third line, discard the first element and convert to numeric
  third_line <- data[3]
  median_values <- as.numeric(unlist(strsplit(third_line, ","))[-1])
  
  # Extract lower bounds of confidence intervals from the fourth line, discard the first element and convert to numeric
  fourth_line <- data[4]
  lower_bounds <- as.numeric(unlist(strsplit(fourth_line, ","))[-1])
  
  # Extract upper bounds of confidence intervals from the fifth line, discard the first element and convert to numeric
  fifth_line <- data[5]
  upper_bounds <- as.numeric(unlist(strsplit(fifth_line, ","))[-1])
  
  # Create a dataframe for parameters
  parameters_df <- data.frame(Parameter = c("or_prop", "red_prop", "new_hosp_or", "new_hosp_red", "icu_or", "icu_red"),
                              Value = parameters)
  
  # Create a dataframe for indicators
  indicators_df <- data.frame(Indicator = indicator_names,
                              Median = median_values,
                              Lower_CI = lower_bounds,
                              Upper_CI = upper_bounds)
  
  # Combine parameters and indicators dataframes into a list
  combined_data <- list(Parameters = parameters_df, Indicators = indicators_df)
  
  # Return the combined data
  return(combined_data)
}
library(RColorBrewer)
plot_data <- function(files, ylim = c(-5000, 5000), xlab = "Date", ylab = "Valeur", legend_labels = NULL,legend_position="topleft",bool_comp=TRUE,ind_key=NULL) {
  nb_colors <- length(files)
  scen_col <- brewer.pal(n = nb_colors, name = "Set1")
  max_day <- max(sapply(files, function(file) dim(readRDS(file))[1]))
  start_date <- sim_date2day('2021-01-01')
  x_values <- sim_day2date(start_date:max_day)
  if(bool_comp==FALSE){
    plot(x = x_values, y = be_ref_data[[ind_key]][(start_date+1):(max_day+1)], ylim = ylim, xlab = xlab, ylab = ylab, xaxt="n",cex.lab=1.5,cex.axis=1.75,pch=4,cex=0.25)  
  }
  else{
  plot(x = x_values, y = rep(NA, length(x_values)), ylim = ylim, xlab = xlab, ylab = ylab, xaxt="n",cex.lab=1.5,cex.axis=1.75)}
  dates <- c('2021-01-01','2021-07-01','2022-01-01','2022-07-01','2023-01-01','2023-07-01')
  axis(side = 1, at = as.Date(dates),labels = c("Jan 21","Jul 21","Jan 22","Jul 22","Jan 23","Jul 23"),cex.axis=1.75)
  #axis(side = 1, at=sim_date2day(c('2021-01-01','2021-07-01','2022-01-01','2022-07-01','2023-01-01','2023-07-01')))
  processed_data <- lapply(files, function(file) {
    tryCatch({
      scen_data <- readRDS(file)
      df <- data.frame(admin_mean = apply(scen_data, 1, median, na.rm = TRUE),
                       admin_LL = apply(scen_data, 1, quantile, 0.025, na.rm = TRUE),
                       admin_UL = apply(scen_data, 1, quantile, 0.975, na.rm = TRUE))
      # Shift des données pour correspondre à l'axe x ajusté
      df <- df[start_date:max_day, ]
      df
    }, error = function(e) {
      message("Erreur lors du chargement du fichier ", file, ": ", e$message)
      return(NULL)
    })
  })
  lower_bounds<-list()
  upper_bounds<-list()
  mean_values<-list()
  for (i in seq_along(processed_data)) {
    if (!is.null(processed_data[[i]])) {
      df <- processed_data[[i]]
      lower_bounds[[i]] <- df$admin_LL
      upper_bounds[[i]] <- df$admin_UL
      mean_values[[i]] <- df$admin_mean
    }
  }
  for (i in seq_along(processed_data)) {
    if (!is.null(processed_data[[i]])) {
      polygon(c(x_values, rev(x_values)), c(lower_bounds[[i]], rev(upper_bounds[[i]])), col = alpha(scen_col[i], 0.1), border = NA)
    }
  }
  for (i in seq_along(processed_data)) {
    if (!is.null(processed_data[[i]])) {
  lines(x_values, mean_values[[i]], type = "l", col = scen_col[i], lwd = 3)
    }
  }
  if (!is.null(legend_labels)) {
    legend(legend_position, legend = legend_labels, col = scen_col[1:length(files)], lwd = 3)
  }
}

par(cex=1.5, cex.main=1.5, cex.lab=5, cex.axis=5)
# Specify the path to the directory containing the bar_summary files
path <-  get_arg_with_default(1,"output/sim_diff/")
all_processed_data65 <- list()
all_processed_data800 <- list()
all_processed_datanobar <- list()
legend_list <- c( "no barometer","orange=0.3, red=0.1","orange=0.6, red=0.2","orange=0.6, red=0.4","orange=0.9, red=0.2","orange=0.9, red=0.4")
# Get the list of all bar_summary CSV files in subdirectories
files_bar_summary65 <- list.files(path = path, recursive = TRUE, pattern = "^bar_summary.*65.*\\.csv$", full.names = TRUE)
files_bar_summary800 <- list.files(path = path, recursive = TRUE, pattern = "^bar_summary.*800.*\\.csv$", full.names = TRUE)
files_bar_summarynobar <- list.files(path = path, recursive = TRUE, pattern = "^bar_summary.*nobarometer.*\\.csv$", full.names = TRUE)


files_hosp_cumul65 <- list.files(path = path, recursive = TRUE, pattern = "^delta_hosp_new_cumul.*65.*\\.rds$", full.names = TRUE)
files_hosp_cumul65<- c(list.files(path = path, recursive = TRUE, pattern = "^delta_hosp_new_cumul.*nobarometer.*\\.rds$", full.names = TRUE),files_hosp_cumul65)
files_hosp_load65 <- list.files(path = path, recursive = TRUE, pattern = "^delta_hosp_load.*65.*\\.rds$", full.names = TRUE)
files_hosp_load65<- c(list.files(path = path, recursive = TRUE, pattern = "^delta_hosp_load.*nobarometer.*\\.rds$", full.names = TRUE),files_hosp_load65)
files_icu_load65 <- list.files(path = path, recursive = TRUE, pattern = "^delta_icu_load.*65.*\\.rds$", full.names = TRUE)
files_icu_load65<- c(list.files(path = path, recursive = TRUE, pattern = "^delta_icu_load.*nobarometer.*\\.rds$", full.names = TRUE),files_icu_load65)
files_cumul_deaths65 <- list.files(path = path, recursive = TRUE, pattern = "^delta_cumul_deaths.*65.*\\.rds$", full.names = TRUE)
files_cumul_deaths65<- c(list.files(path = path, recursive = TRUE, pattern = "^delta_cumul_deaths.*nobarometer.*\\.rds$", full.names = TRUE),files_cumul_deaths65)
files_cumul_cases65 <- list.files(path = path, recursive = TRUE, pattern = "^delta_cumul_cases.*65.*\\.rds$", full.names = TRUE)
files_cumul_cases65<- c(list.files(path = path, recursive = TRUE, pattern = "^delta_cumul_cases.*nobarometer.*\\.rds$", full.names = TRUE),files_cumul_cases65)
files_social_coef65 <- list.files(path = path, recursive = TRUE, pattern = "^social_coef.*65.*\\.rds$", full.names = TRUE)
files_social_coef65<- c(list.files(path = path, recursive = TRUE, pattern = "^social_coef.*nobarometer.*\\.rds$", full.names = TRUE),files_social_coef65)
# Initialize a list to store processed data for all files

#files_bar_summary800 <- list.files(path = path, recursive = TRUE, pattern = "^bar_summary.*800.*\\.csv$", full.names = TRUE)
files_hosp_cumul800 <- list.files(path = path, recursive = TRUE, pattern = "^delta_hosp_new_cumul.*800.*\\.rds$", full.names = TRUE)
files_hosp_cumul800<- c(list.files(path = path, recursive = TRUE, pattern = "^delta_hosp_new_cumul.*nobarometer.*\\.rds$", full.names = TRUE),files_hosp_cumul800)
files_hosp_load800 <- list.files(path = path, recursive = TRUE, pattern = "^delta_hosp_load.*800.*\\.rds$", full.names = TRUE)
files_hosp_load800<- c(list.files(path = path, recursive = TRUE, pattern = "^delta_hosp_load.*nobarometer.*\\.rds$", full.names = TRUE),files_hosp_load800)
files_icu_load800 <- list.files(path = path, recursive = TRUE, pattern = "^delta_icu_load.*800.*\\.rds$", full.names = TRUE)
files_icu_load800<- c(list.files(path = path, recursive = TRUE, pattern = "^delta_icu_load.*nobarometer.*\\.rds$", full.names = TRUE),files_icu_load800)
files_cumul_deaths800 <- list.files(path = path, recursive = TRUE, pattern = "^delta_cumul_deaths.*800.*\\.rds$", full.names = TRUE)
files_cumul_deaths800<- c(list.files(path = path, recursive = TRUE, pattern = "^delta_cumul_deaths.*nobarometer.*\\.rds$", full.names = TRUE),files_cumul_deaths800)
files_cumul_cases800 <- list.files(path = path, recursive = TRUE, pattern = "^delta_cumul_cases.*800.*\\.rds$", full.names = TRUE)
files_cumul_cases800<- c(list.files(path = path, recursive = TRUE, pattern = "^delta_cumul_cases.*nobarometer.*\\.rds$", full.names = TRUE),files_cumul_cases800)
files_social_coef800 <- list.files(path = path, recursive = TRUE, pattern = "^social_coef.*800.*\\.rds$", full.names = TRUE)
files_social_coef800<- c(list.files(path = path, recursive = TRUE, pattern = "^social_coef.*nobarometer.*\\.rds$", full.names = TRUE),files_social_coef800)

files_hosp_cumul_vac <- list.files(path = path, recursive = TRUE, pattern = "^hosp_adm_cum.*vac_vac.*\\.rds$", full.names = TRUE)
files_hosp_new_vac <- list.files(path = path, recursive = TRUE, pattern = "^hosp_adm_incr.*vac_vac.*\\.rds$", full.names = TRUE)
files_hosp_load_vac <- list.files(path = path, recursive = TRUE, pattern = "^hosp_load.*vac_vac.*\\.rds$", full.names = TRUE)
files_icu_load_vac <- list.files(path = path, recursive = TRUE, pattern = "^icu_load.*vac_vac.*\\.rds$", full.names = TRUE)
files_cumul_deaths_vac<- list.files(path = path, recursive = TRUE, pattern = "^cumul_deaths.*vac_vac.*\\.rds$", full.names = TRUE)
files_cumul_cases_vac<- list.files(path = path, recursive = TRUE, pattern = "^cumul_cases.*vac_vac.*\\.rds$", full.names = TRUE)


# Initialize a list to store processed data for all files
pdf(file = paste0(path,"vac_newhosp.pdf"),, width = 15, height = 10)
#par(mfrow=c(1,2))
be_ref_data$hospital_admissions <- be_ref_data$hospital_admissions + be_ref_data$hospital_admissions_other
plot_data(files_hosp_new_vac,ylab = "new hospitalisations",ylim=c(0,1000),bool_comp=FALSE,ind_key="hospital_admissions")
dev.off()
be_ref_data$hospital_admissions[is.na(be_ref_data$hospital_admissions)] <-0
be_ref_data$cumulative_hospital_admissions <- cumsum(be_ref_data$hospital_admissions)
pdf(file = paste0(path,"vac_cumhosp.pdf"),, width = 15, height = 10)
#par(mfrow=c(1,2))
plot_data(files_hosp_cumul_vac,ylab = "Cumulative new hospitalisations",ylim=c(0,250000),bool_comp=FALSE,ind_key="cumulative_hospital_admissions")
#title <- "Baseline model - New Hospitalisations"
#title(title, outer = TRUE, line = -1.5, cex.main = 1.5)
dev.off()
pdf(file = paste0(path,"recap_cumhosp.pdf"),, width = 15, height = 10)
par(mfrow=c(1,2))
plot_data(files_hosp_cumul65,ylab = "absolute difference",ylim=c(-10000,500000),legend_labels = legend_list)
mtext("a) Original barometer", side = 3, line = 1, adj = 0, cex = 1.2)
plot_data(files_hosp_cumul800,ylab = "absolute difference",ylim=c(-10000,500000),legend_labels = legend_list)
mtext("b) Less-stringent barometer", side = 3, line = 1, adj = 0, cex = 1.2)
#title <- "Cumulative difference in hospital admissions"
#title(title, outer = TRUE, line = -1.5, cex.main = 1.5)
dev.off()

pdf(file = paste0(path,"vac_hospload.pdf"),, width = 15, height = 10)
#par(mfrow=c(1,2))
plot_data(files_hosp_load_vac,ylab = "Hospital load",ylim=c(0,5000),bool_comp = FALSE,ind_key="hospital_load")
#title <- "Baseline model - Hospital load"
#title(title, outer = TRUE, line = -1.5, cex.main = 1.5)
dev.off()
pdf(file = paste0(path,"recap_hospload.pdf"),, width = 15, height = 10)
par(mfrow=c(1,2))
plot_data(files_hosp_load65,ylab = "Difference",ylim=c(-20000,100000),legend_labels = legend_list,legend_position="bottomleft")
mtext("a) Original barometer", side = 3, line = 1, adj = 0, cex = 1.2)
plot_data(files_hosp_load800,ylab = "Difference",ylim=c(-20000,100000),legend_labels = legend_list,legend_position="bottomleft")
mtext("b) Less-stringent barometer", side = 3, line = 1, adj = 0, cex = 1.2)
#title <- "Difference in hospital load"
#title(title, outer = TRUE, line = -1.5, cex.main = 1.5)
dev.off()

pdf(file = paste0(path,"vac_icuload.pdf"),, width = 15, height = 10)
#par(mfrow=c(1,2))
plot_data(files_icu_load_vac,ylab = "ICU load",ylim=c(0,1200),bool_comp=FALSE,ind_key = "icu_load")
#title <- "Baseline model - ICU load"
#title(title, outer = TRUE, line = -1.5, cex.main = 1.5)
dev.off()
pdf(file = paste0(path,"recap_icuload.pdf"),, width = 15, height = 10)
par(mfrow=c(1,2))
plot_data(files_icu_load65,ylab = "Difference",ylim=c(-1000,40000),legend_labels = legend_list)
mtext("a) Original barometer", side = 3, line = 1, adj = 0, cex = 1.2)
plot_data(files_icu_load800,ylab = "Difference",ylim=c(-1000,40000),legend_labels = legend_list)
mtext("b) Less-stringent barometer", side = 3, line = 1, adj = 0, cex = 1.2)
#title <- "Difference in ICU load"
#title(title, outer = TRUE, line = -1.5, cex.main = 1.5)
dev.off()
pdf(file = paste0(path,"recap_icuload2.pdf"),, width = 15, height = 10)
par(mfrow=c(1,2))
plot_data(files_icu_load65,ylab = "Difference",ylim=c(-1000,3000),legend_labels = legend_list,legend_position="topright")
mtext("a) Original barometer", side = 3, line = 1, adj = 0, cex = 1.2)
plot_data(files_icu_load800,ylab = "Difference",ylim=c(-1000,3000),legend_labels = legend_list,legend_position="topright")
mtext("b) Less-stringent barometer", side = 3, line = 1, adj = 0, cex = 1.2)
#title <- "Difference in ICU load"
#title(title, outer = TRUE, line = -1.5, cex.main = 1.5)
dev.off()

pdf(file = paste0(path,"vac_deaths.pdf"),, width = 15, height = 10)
#par(mfrow=c(1,2))
be_ref_data$covid19_deaths[is.na(be_ref_data$covid19_deaths)] <- 0
be_ref_data$cumul_deaths <- cumsum(be_ref_data$covid19_deaths)
plot_data(files_cumul_deaths_vac,ylab = "Cumulative Deaths",ylim=c(0,50000),bool_comp = FALSE,ind_key = "cumul_deaths")
#title <- "Baseline model - Mortality"
#title(title, outer = TRUE, line = -1.5, cex.main = 1.5)
dev.off()
pdf(file = paste0(path,"recap_avdeaths.pdf"),, width = 15, height = 10)
par(mfrow=c(1,2))
plot_data(files_cumul_deaths65,ylim=c(-5000,70000),ylab="Averted deaths",legend_labels = legend_list)
mtext("a) Original barometer", side = 3, line = 1, adj = 0, cex = 1.2)
plot_data(files_cumul_deaths800,ylim=c(-5000,70000),ylab="Averted deaths",legend_labels = legend_list)
mtext("b) Less-stringent barometer", side = 3, line = 1, adj = 0, cex = 1.2)
#title <- "Cumulative difference in mortality"
#title(title, outer = TRUE, line = -1.5, cex.main = 1.5)
dev.off()

pdf(file = paste0(path,"vac_cases.pdf"),, width = 15, height = 10)
#par(mfrow=c(1,2))
be_ref_data$cumul_cases <- cumsum(be_ref_data$cases)
plot_data(files_cumul_cases_vac,ylab = "Cumulative cases",ylim=c(0,20e6),bool_comp = FALSE,ind_key = "cumul_cases")
#title <- "Baseline model - Cumulative new infections"
#title(title, outer = TRUE, line = -1.5, cex.main = 1.5)
dev.off()
pdf(file = paste0(path,"recap_cases.pdf"),, width = 15, height = 10)
par(mfrow=c(1,2))
plot_data(files_cumul_cases65,ylim=c(-7000000,11000000),ylab="Cumulative difference of new infections",legend_labels = legend_list)
mtext("a) Original barometer", side = 3, line = 1, adj = 0, cex = 1.2)
plot_data(files_cumul_cases800,ylim=c(-7000000,11000000),ylab="Cumulative difference of new infections",legend_labels = legend_list)
mtext("b) Less-stringent barometer", side = 3, line = 1, adj = 0, cex = 1.2)
#title <- "Cumulative difference in daily infections"
#title(title, outer = TRUE, line = -1.5, cex.main = 1.5)
dev.off()

pdf(file = paste0(path,"recap_social.pdf"),, width = 15, height = 10)
par(mfrow=c(1,2))
plot_data(files_social_coef65,ylim=c(0,1.2),ylab=expression(paste("s"["bar"],"(t)")),legend_labels = legend_list)
mtext("a) Original barometer", side = 3, line = 1, adj = 0, cex = 1.2)
plot_data(files_social_coef800,ylim=c(0,1.2),ylab=expression(paste("s"["bar"],"(t)")),legend_labels = legend_list)
mtext("b) Less-stringent barometer", side = 3, line = 1, adj = 0, cex = 1.2)
#title <- ""
#title(title, outer = TRUE, line = -1.5, cex.main = 1.5)
dev.off()
exp_sbar <- readRDS("~/Documents/GitHub/stochastic_model_BE/output/sim_diff/barometer_orange_0.6_65_300_red_0.2_150_500_1217_cnt0_0_0_0_100_n75_s3_belgium/social_coef0_0_0_0_100_barometer_orange_0p6_65_300_red_0p2_150_500_vaccine_uptake_extrabooster_vscenbelgium_novac_belgium.rds")
pdf(file = paste0(path,"exp_social_coef.pdf"),, width = 15, height = 10)
scen_col <- brewer.pal(n = 3, name = "Set1")
max_day <- max(dim(exp_sbar))
start_date <- sim_date2day('2021-01-01')
x_values <- sim_day2date(start_date:max_day)
ylab <-expression(paste("s"["bar"],"(t)"))
xlab <- "Date"
plot(x = x_values, y = rep(NA, length(x_values)), ylim = c(0,1.2), xlab = xlab, ylab = ylab, xaxt="n")
dates <- c('2021-01-01','2021-07-01','2022-01-01','2022-07-01','2023-01-01','2023-07-01')
axis(side = 1, at = as.Date(dates),labels = c("Jan 21","Jul 21","Jan 22","Jul 22","Jan 23","Jul 23"))
lines(x_values, exp_sbar[start_date:max_day,5], type = "l", col = scen_col[1], lwd = 1)
#lines(x_values, exp_sbar[start_date:max_day,17], type = "l", col = scen_col[2], lwd = 2)
#lines(x_values, exp_sbar[start_date:max_day,23], type = "l", col = scen_col[3], lwd = 2)
lines(x_values, exp_sbar[start_date:max_day,32], type = "l", col = scen_col[2], lwd = 1.5)
lines(x_values, exp_sbar[start_date:max_day,61], type = "l", col = scen_col[3], lwd = 1.5)
dev.off()

# Process each bar_summary file
for (file in files_bar_summary65) {
  # Process the current file
  processed_data <- process_bar_summary(file)
  
  # Store the processed data in the list
  all_processed_data65[[file]] <- processed_data
}
for (file in files_bar_summary800) {
  # Process the current file
  processed_data <- process_bar_summary(file)
  
  # Store the processed data in the list
  all_processed_data800[[file]] <- processed_data
}
for (file in files_bar_summarynobar) {
  # Process the current file
  processed_data <- process_bar_summary(file)
  
  # Store the processed data in the list
  all_processed_datanobar[[file]] <- processed_data
}
saveRDS(c(all_processed_data65,all_processed_data800,all_processed_datanobar),file=paste0(path,"recap.rds"))

# Determine the number of parameters and indicators
num_parameters <- nrow(all_processed_data65[[files_bar_summary65[1]]]$Parameters)
num_indicators <- nrow(all_processed_data65[[files_bar_summary65[1]]]$Indicators)  
# Calculate the median values of delta_death for each file
median_delta_death65 <- sapply(all_processed_data65, function(x) {
  median(x$Indicators[x$Indicators$Indicator == "delta_death", "Median"])
})
median_delta_death800 <- sapply(all_processed_data800, function(x) {
  median(x$Indicators[x$Indicators$Indicator == "delta_death", "Median"])
})
median_delta_deathnobar <- sapply(all_processed_datanobar, function(x) {
  median(x$Indicators[x$Indicators$Indicator == "delta_death", "Median"])
})

# Sort the files based on median delta_death values
sorted_files65 <- files_bar_summary65[order(median_delta_death65)]
sorted_files800 <- files_bar_summary800[order(median_delta_death800)]
sorted_filesnobar <- files_bar_summarynobar[order(median_delta_deathnobar)]

par_slice <- c(1:2)
ind_slice <-c(1,2)
# Create or open the LaTeX file for writing
latex_file <- file(paste0(path,"recapdeath.tex"), "w")

# Write the preamble of the LaTeX document with the desired paper size
cat("\\documentclass{article}\n", file = latex_file)
cat("\\usepackage{booktabs}\n", file = latex_file)
cat("\\usepackage[a4paper, margin=1in]{geometry}\n", file = latex_file)  # Set paper size to A0 with 1-inch margins
cat("\\usepackage{pdflscape}\n", file = latex_file)  # Add pdflscape package for landscape orientation
cat("\\usepackage{makecell}\n", file = latex_file)  # Add pdflscape package for landscape orientation
cat("\\usepackage{multirow}\n", file = latex_file)  # Add pdflscape package for landscape orientation
cat("\\begin{document}\n\n", file = latex_file)

# Write the landscape environment
#cat("\\begin{landscape}\n", file = latex_file)

# Write the table environment with resizebox
cat("\\begin{table}[htbp]\n", file = latex_file)
cat("    \\centering\n", file = latex_file)
#cat("    \\caption{Summary of Parameters and Indicators}\n", file = latex_file)
#cat("    \\resizebox{\\linewidth}{!}{%\n", file = latex_file)  # Resize the table to fit the page width
cat("    \\begin{tabular}{|c", paste(rep("|c", num_parameters + num_indicators), collapse = ""), "|}\n", file = latex_file)  # Add vertical lines between columns

# Write the header line with parameter names and indicator names
cat("        \\hline\n", file = latex_file)  # Add a horizontal line above the header
cat("        ", paste(gsub("_", "\\\\\\_", all_processed_data65[[files_bar_summary65[1]]]$Parameters$Parameter[par_slice]), collapse = " & "), " & ", paste(gsub("_", "\\\\\\_", all_processed_data65[[files_bar_summary65[1]]]$Indicators$Indicator[ind_slice]), collapse = " & "), " \\\\\n", file = latex_file)
cat("        \\hline\n", file = latex_file)  # Add a horizontal line below the header
cat("        \n", file = latex_file)  # Add a horizontal line below the header


cat("        \\multicolumn{4}{|l|}{\\multirow{2}{*}{\\bf Basic Barometer}}\\\\ \n", file = latex_file)  # Add a horizontal line below the header
cat("        \\multicolumn{4}{|l|}{}\\\\ \n", file = latex_file)  # Add a horizontal line below the header
cat("        \\hline\n", file = latex_file)  # Add a horizontal line below the header


# Process each bar_summary file in sorted order
for (file in sorted_files65) {
  
  # Extract parameters values and format them as "%.2g"
  parameters_values <- sprintf("%.8g", all_processed_data65[[file]]$Parameters$Value[par_slice])
  
  # Extract indicators values (median and confidence intervals) and format them as "%.2g"
  indicators_data <- apply(all_processed_data65[[file]]$Indicators[ind_slice, -1], 1, function(x) {
    paste0("\\makecell{\\rule{0pt}{2.2ex}",sprintf("%.8g", x[1]), "\\\\ {\\footnotesize$[", sprintf("%.8g", x[2]), ";", sprintf("%.8g", x[3]), "]$}}")
  })
  
  # Combine parameter values and indicator values with "&" between each column
  data_line <- paste(c(parameters_values, unlist(indicators_data)), collapse = " & ")
  data_line <- paste(data_line, "\\\\")
  
  # Write the data line to LaTeX file
  cat("        ", data_line, "\n", file = latex_file)
  
  # Write hline
  cat("        \\hline \n ", file = latex_file)
}
cat("        \\multicolumn{4}{|l|}{\\multirow{2}{*}{\\bf Less Stringent Barometer}}\\\\ \n", file = latex_file)  # Add a horizontal line below the header
cat("        \\multicolumn{4}{|l|}{}\\\\ \n", file = latex_file)  # Add a horizontal line below the header
cat("        \\hline\n", file = latex_file)  # Add a horizontal line below the header
# Process each bar_summary file in sorted order
for (file in sorted_files800) {
  
  # Extract parameters values and format them as "%.2g"
  parameters_values <- sprintf("%.8g", all_processed_data800[[file]]$Parameters$Value[par_slice])
  
  # Extract indicators values (median and confidence intervals) and format them as "%.2g"
  indicators_data <- apply(all_processed_data800[[file]]$Indicators[ind_slice, -1], 1, function(x) {
    paste0("\\makecell{\\rule{0pt}{2.2ex}",sprintf("%.8g", x[1]), "\\\\ {\\footnotesize$[", sprintf("%.8g", x[2]), ";", sprintf("%.8g", x[3]), "]$}}")
  })
  
  # Combine parameter values and indicator values with "&" between each column
  data_line <- paste(c(parameters_values, unlist(indicators_data)), collapse = " & ")
  data_line <- paste(data_line, "\\\\")
  
  # Write the data line to LaTeX file
  cat("        ", data_line, "\n", file = latex_file)
  
  # Write hline
  cat("        \\hline \n ", file = latex_file)
}
cat("        \\multicolumn{4}{|l|}{\\multirow{2}{*}{\\bf No Barometer}}\\\\ \n", file = latex_file)  # Add a horizontal line below the header
cat("        \\multicolumn{4}{|l|}{}\\\\ \n", file = latex_file)  # Add a horizontal line below the header
cat("        \\hline\n", file = latex_file)  # Add a horizontal line below the header
# Process each bar_summary file in sorted order
for (file in sorted_filesnobar) {
  
  # Extract parameters values and format them as "%.2g"
  parameters_values <- sprintf("%.8g", all_processed_datanobar[[file]]$Parameters$Value[par_slice])
  
  # Extract indicators values (median and confidence intervals) and format them as "%.2g"
  indicators_data <- apply(all_processed_datanobar[[file]]$Indicators[ind_slice, -1], 1, function(x) {
    paste0("\\makecell{\\rule{0pt}{2.2ex}",sprintf("%.8g", x[1]), "\\\\ {\\footnotesize$[", sprintf("%.8g", x[2]), ";", sprintf("%.8g", x[3]), "]$}}")
  })
  
  # Combine parameter values and indicator values with "&" between each column
  data_line <- paste(c(parameters_values, unlist(indicators_data)), collapse = " & ")
  data_line <- paste(data_line, "\\\\")
  
  # Write the data line to LaTeX file
  cat("        ", data_line, "\n", file = latex_file)
  
  # Write hline
  cat("        \\hline \n ", file = latex_file)
}
# Write the end of the tabular environment
cat("    \\end{tabular}\n", file = latex_file)
#cat("    }%\n", file = latex_file)  # End of resizebox
cat("\\end{table}\n\n", file = latex_file)

# Write the end of the landscape environment
#cat("\\end{landscape}\n", file = latex_file)

# Write the end of the document
cat("\\end{document}\n", file = latex_file)

# Close the LaTeX file
close(latex_file)
par_slice <- c(1:2)
ind_slice <-c(6,4,7,5)
# Create or open the LaTeX file for writing
latex_file <- file(paste0(path,"recaphosp.tex"), "w")

# Write the preamble of the LaTeX document with the desired paper size
cat("\\documentclass{article}\n", file = latex_file)
cat("\\usepackage{booktabs}\n", file = latex_file)
cat("\\usepackage[a4paper, margin=1in]{geometry}\n", file = latex_file)  # Set paper size to A0 with 1-inch margins
cat("\\usepackage{pdflscape}\n", file = latex_file)  # Add pdflscape package for landscape orientation
cat("\\usepackage{makecell}\n", file = latex_file)  # Add pdflscape package for landscape orientation
cat("\\usepackage{multirow}\n", file = latex_file)  # Add pdflscape package for landscape orientation
cat("\\begin{document}\n\n", file = latex_file)

# Write the landscape environment
cat("\\begin{landscape}\n", file = latex_file)

# Write the table environment with resizebox
cat("\\begin{table}[htbp]\n", file = latex_file)
cat("    \\centering\n", file = latex_file)
#cat("    \\caption{Summary of Parameters and Indicators}\n", file = latex_file)
#cat("    \\resizebox{\\linewidth}{!}{%\n", file = latex_file)  # Resize the table to fit the page width
cat("    \\begin{tabular}{|c", paste(rep("|c", num_parameters + num_indicators), collapse = ""), "|}\n", file = latex_file)  # Add vertical lines between columns

# Write the header line with parameter names and indicator names
cat("        \\hline\n", file = latex_file)  # Add a horizontal line above the header
cat("        ", paste(gsub("_", "\\\\\\_", all_processed_data65[[files_bar_summary65[1]]]$Parameters$Parameter[par_slice]), collapse = " & "), " & ", paste(gsub("_", "\\\\\\_", all_processed_data65[[files_bar_summary65[1]]]$Indicators$Indicator[ind_slice]), collapse = " & "), " \\\\\n", file = latex_file)
cat("        \\hline\n", file = latex_file)  # Add a horizontal line below the header
cat("        \n", file = latex_file)  # Add a horizontal line below the header


cat("        \\multicolumn{6}{|l|}{\\multirow{2}{*}{\\bf Basic Barometer}}\\\\ \n", file = latex_file)  # Add a horizontal line below the header
cat("        \\multicolumn{6}{|l|}{}\\\\ \n", file = latex_file)  # Add a horizontal line below the header
cat("        \\hline\n", file = latex_file)  # Add a horizontal line below the header


# Process each bar_summary file in sorted order
for (file in sorted_files65) {
  
  # Extract parameters values and format them as "%.2g"
  parameters_values <- sprintf("%.8g", all_processed_data65[[file]]$Parameters$Value[par_slice])
  
  # Extract indicators values (median and confidence intervals) and format them as "%.2g"
  indicators_data <- apply(all_processed_data65[[file]]$Indicators[ind_slice, -1], 1, function(x) {
    paste0("\\makecell{\\rule{0pt}{2.2ex}",sprintf("%.8g", x[1]), "\\\\ {\\footnotesize$[", sprintf("%.8g", x[2]), ";", sprintf("%.8g", x[3]), "]$}}")
  })
  
  # Combine parameter values and indicator values with "&" between each column
  data_line <- paste(c(parameters_values, unlist(indicators_data)), collapse = " & ")
  data_line <- paste(data_line, "\\\\")
  
  # Write the data line to LaTeX file
  cat("        ", data_line, "\n", file = latex_file)
  
  # Write hline
  cat("        \\hline \n ", file = latex_file)
}
cat("        \\multicolumn{6}{|l|}{\\multirow{2}{*}{\\bf Less Stringent Barometer}}\\\\ \n", file = latex_file)  # Add a horizontal line below the header
cat("        \\multicolumn{6}{|l|}{}\\\\ \n", file = latex_file)  # Add a horizontal line below the header
cat("        \\hline\n", file = latex_file)  # Add a horizontal line below the header
# Process each bar_summary file in sorted order
for (file in sorted_files800) {
  
  # Extract parameters values and format them as "%.2g"
  parameters_values <- sprintf("%.8g", all_processed_data800[[file]]$Parameters$Value[par_slice])
  
  # Extract indicators values (median and confidence intervals) and format them as "%.2g"
  indicators_data <- apply(all_processed_data800[[file]]$Indicators[ind_slice, -1], 1, function(x) {
    paste0("\\makecell{\\rule{0pt}{2.2ex}",sprintf("%.8g", x[1]), "\\\\ {\\footnotesize$[", sprintf("%.8g", x[2]), ";", sprintf("%.8g", x[3]), "]$}}")
  })
  
  # Combine parameter values and indicator values with "&" between each column
  data_line <- paste(c(parameters_values, unlist(indicators_data)), collapse = " & ")
  data_line <- paste(data_line, "\\\\")
  
  # Write the data line to LaTeX file
  cat("        ", data_line, "\n", file = latex_file)
  
  # Write hline
  cat("        \\hline \n ", file = latex_file)
}
cat("        \\multicolumn{6}{|l|}{\\multirow{2}{*}{\\bf No Barometer}}\\\\ \n", file = latex_file)  # Add a horizontal line below the header
cat("        \\multicolumn{6}{|l|}{}\\\\ \n", file = latex_file)  # Add a horizontal line below the header
cat("        \\hline\n", file = latex_file)  # Add a horizontal line below the header
# Process each bar_summary file in sorted order
for (file in sorted_filesnobar) {
  
  # Extract parameters values and format them as "%.2g"
  parameters_values <- sprintf("%.8g", all_processed_datanobar[[file]]$Parameters$Value[par_slice])
  
  # Extract indicators values (median and confidence intervals) and format them as "%.2g"
  indicators_data <- apply(all_processed_datanobar[[file]]$Indicators[ind_slice, -1], 1, function(x) {
    paste0("\\makecell{\\rule{0pt}{2.2ex}",sprintf("%.8g", x[1]), "\\\\ {\\footnotesize$[", sprintf("%.8g", x[2]), ";", sprintf("%.8g", x[3]), "]$}}")
  })
  
  # Combine parameter values and indicator values with "&" between each column
  data_line <- paste(c(parameters_values, unlist(indicators_data)), collapse = " & ")
  data_line <- paste(data_line, "\\\\")
  
  # Write the data line to LaTeX file
  cat("        ", data_line, "\n", file = latex_file)
  
  # Write hline
  cat("        \\hline \n ", file = latex_file)
}
# Write the end of the tabular environment
cat("    \\end{tabular}\n", file = latex_file)
#cat("    }%\n", file = latex_file)  # End of resizebox
cat("\\end{table}\n\n", file = latex_file)

# Write the end of the landscape environment
cat("\\end{landscape}\n", file = latex_file)

# Write the end of the document
cat("\\end{document}\n", file = latex_file)

# Close the LaTeX file
close(latex_file)
par_slice <- c(1:2)
ind_slice <-c(8,9,10)
# Create or open the LaTeX file for writing
latex_file <- file(paste0(path,"recapsocial.tex"), "w")

# Write the preamble of the LaTeX document with the desired paper size
cat("\\documentclass{article}\n", file = latex_file)
cat("\\usepackage{booktabs}\n", file = latex_file)
cat("\\usepackage[a4paper, margin=1in]{geometry}\n", file = latex_file)  # Set paper size to A0 with 1-inch margins
cat("\\usepackage{pdflscape}\n", file = latex_file)  # Add pdflscape package for landscape orientation
cat("\\usepackage{makecell}\n", file = latex_file)  # Add pdflscape package for landscape orientation
cat("\\usepackage{multirow}\n", file = latex_file)  # Add pdflscape package for landscape orientation
cat("\\begin{document}\n\n", file = latex_file)

# Write the landscape environment
cat("\\begin{landscape}\n", file = latex_file)

# Write the table environment with resizebox
cat("\\begin{table}[htbp]\n", file = latex_file)
cat("    \\centering\n", file = latex_file)
#cat("    \\caption{Summary of Parameters and Indicators}\n", file = latex_file)
#cat("    \\resizebox{\\linewidth}{!}{%\n", file = latex_file)  # Resize the table to fit the page width
cat("    \\begin{tabular}{|c", paste(rep("|c", num_parameters + num_indicators), collapse = ""), "|}\n", file = latex_file)  # Add vertical lines between columns

# Write the header line with parameter names and indicator names
cat("        \\hline\n", file = latex_file)  # Add a horizontal line above the header
cat("        ", paste(gsub("_", "\\\\\\_", all_processed_data65[[files_bar_summary65[1]]]$Parameters$Parameter[par_slice]), collapse = " & "), " & ", paste(gsub("_", "\\\\\\_", all_processed_data65[[files_bar_summary65[1]]]$Indicators$Indicator[ind_slice]), collapse = " & "), " \\\\\n", file = latex_file)
cat("        \\hline\n", file = latex_file)  # Add a horizontal line below the header
cat("        \n", file = latex_file)  # Add a horizontal line below the header


cat("        \\multicolumn{5}{|l|}{\\multirow{2}{*}{\\bf Basic Barometer}}\\\\ \n", file = latex_file)  # Add a horizontal line below the header
cat("        \\multicolumn{5}{|l|}{}\\\\ \n", file = latex_file)  # Add a horizontal line below the header
cat("        \\hline\n", file = latex_file)  # Add a horizontal line below the header


# Process each bar_summary file in sorted order
for (file in sorted_files65) {
  
  # Extract parameters values and format them as "%.2g"
  parameters_values <- sprintf("%.8g", all_processed_data65[[file]]$Parameters$Value[par_slice])
  
  # Extract indicators values (median and confidence intervals) and format them as "%.2g"
  indicators_data <- apply(all_processed_data65[[file]]$Indicators[ind_slice, -1], 1, function(x) {
    paste0("\\makecell{\\rule{0pt}{2.2ex}",sprintf("%.8g", x[1]), "\\\\ {\\footnotesize$[", sprintf("%.8g", x[2]), ";", sprintf("%.8g", x[3]), "]$}}")
  })
  
  # Combine parameter values and indicator values with "&" between each column
  data_line <- paste(c(parameters_values, unlist(indicators_data)), collapse = " & ")
  data_line <- paste(data_line, "\\\\")
  
  # Write the data line to LaTeX file
  cat("        ", data_line, "\n", file = latex_file)
  
  # Write hline
  cat("        \\hline \n ", file = latex_file)
}
cat("        \\multicolumn{5}{|l|}{\\multirow{2}{*}{\\bf Less Stringent Barometer}}\\\\ \n", file = latex_file)  # Add a horizontal line below the header
cat("        \\multicolumn{5}{|l|}{}\\\\ \n", file = latex_file)  # Add a horizontal line below the header
cat("        \\hline\n", file = latex_file)  # Add a horizontal line below the header
# Process each bar_summary file in sorted order
for (file in sorted_files800) {
  
  # Extract parameters values and format them as "%.2g"
  parameters_values <- sprintf("%.8g", all_processed_data800[[file]]$Parameters$Value[par_slice])
  
  # Extract indicators values (median and confidence intervals) and format them as "%.2g"
  indicators_data <- apply(all_processed_data800[[file]]$Indicators[ind_slice, -1], 1, function(x) {
    paste0("\\makecell{\\rule{0pt}{2.2ex}",sprintf("%.8g", x[1]), "\\\\ {\\footnotesize$[", sprintf("%.8g", x[2]), ";", sprintf("%.8g", x[3]), "]$}}")
  })
  
  # Combine parameter values and indicator values with "&" between each column
  data_line <- paste(c(parameters_values, unlist(indicators_data)), collapse = " & ")
  data_line <- paste(data_line, "\\\\")
  
  # Write the data line to LaTeX file
  cat("        ", data_line, "\n", file = latex_file)
  
  # Write hline
  cat("        \\hline \n ", file = latex_file)
}
cat("        \\multicolumn{5}{|l|}{\\multirow{2}{*}{\\bf No Barometer}}\\\\ \n", file = latex_file)  # Add a horizontal line below the header
cat("        \\multicolumn{5}{|l|}{}\\\\ \n", file = latex_file)  # Add a horizontal line below the header
cat("        \\hline\n", file = latex_file)  # Add a horizontal line below the header
# Process each bar_summary file in sorted order
for (file in sorted_filesnobar) {
  
  # Extract parameters values and format them as "%.2g"
  parameters_values <- sprintf("%.8g", all_processed_datanobar[[file]]$Parameters$Value[par_slice])
  
  # Extract indicators values (median and confidence intervals) and format them as "%.2g"
  indicators_data <- apply(all_processed_datanobar[[file]]$Indicators[ind_slice, -1], 1, function(x) {
    paste0("\\makecell{\\rule{0pt}{2.2ex}",sprintf("%.8g", x[1]), "\\\\ {\\footnotesize$[", sprintf("%.8g", x[2]), ";", sprintf("%.8g", x[3]), "]$}}")
  })
  
  # Combine parameter values and indicator values with "&" between each column
  data_line <- paste(c(parameters_values, unlist(indicators_data)), collapse = " & ")
  data_line <- paste(data_line, "\\\\")
  
  # Write the data line to LaTeX file
  cat("        ", data_line, "\n", file = latex_file)
  
  # Write hline
  cat("        \\hline \n ", file = latex_file)
}
# Write the end of the tabular environment
cat("    \\end{tabular}\n", file = latex_file)
#cat("    }%\n", file = latex_file)  # End of resizebox
cat("\\end{table}\n\n", file = latex_file)

# Write the end of the landscape environment
cat("\\end{landscape}\n", file = latex_file)

# Write the end of the document
cat("\\end{document}\n", file = latex_file)

# Close the LaTeX file
close(latex_file)
# all_processed_data <- list()
# # Process each bar_summary file
# for (file in files_bar_summary800) {
#   # Process the current file
#   processed_data <- process_bar_summary(file)
#   
#   # Store the processed data in the list
#   all_processed_data[[file]] <- processed_data
# }
# saveRDS(all_processed_data,file=paste0(path,"recap800.rds"))
# # Determine the number of parameters and indicators
# num_parameters <- nrow(all_processed_data[[files_bar_summary800[1]]]$Parameters)
# num_indicators <- nrow(all_processed_data[[files_bar_summary800[1]]]$Indicators)  
# # Calculate the median values of delta_death for each file
# median_delta_death <- sapply(all_processed_data, function(x) {
#   median(x$Indicators[x$Indicators$Indicator == "delta_death", "Median"])
# })
# 
# # Sort the files based on median delta_death values
# sorted_files <- files_bar_summary800[order(median_delta_death)]
# 
# # Create or open the LaTeX file for writing
# latex_file <- file(paste0(path,"recap800.tex"), "w")
# 
# # Write the preamble of the LaTeX document with the desired paper size
# cat("\\documentclass{article}\n", file = latex_file)
# cat("\\usepackage{booktabs}\n", file = latex_file)
# cat("\\usepackage[a0paper, margin=1in]{geometry}\n", file = latex_file)  # Set paper size to A0 with 1-inch margins
# cat("\\usepackage{pdflscape}\n", file = latex_file)  # Add pdflscape package for landscape orientation
# cat("\\begin{document}\n\n", file = latex_file)
# 
# # Write the landscape environment
# cat("\\begin{landscape}\n", file = latex_file)
# 
# # Write the table environment with resizebox
# cat("\\begin{table}[htbp]\n", file = latex_file)
# cat("    \\centering\n", file = latex_file)
# cat("    \\caption{Summary of Parameters and Indicators}\n", file = latex_file)
# cat("    \\resizebox{\\linewidth}{!}{%\n", file = latex_file)  # Resize the table to fit the page width
# cat("    \\begin{tabular}{|c", paste(rep("|c", num_parameters + num_indicators), collapse = ""), "|}\n", file = latex_file)  # Add vertical lines between columns
# 
# # Write the header line with parameter names and indicator names
# cat("        \\hline\n", file = latex_file)  # Add a horizontal line above the header
# cat("        ", paste(gsub("_", "\\\\\\_", all_processed_data[[files_bar_summary800[1]]]$Parameters$Parameter), collapse = " & "), " & ", paste(gsub("_", "\\\\\\_", all_processed_data[[files_bar_summary800[1]]]$Indicators$Indicator), collapse = " & "), " \\\\\n", file = latex_file)
# cat("        \\hline\n", file = latex_file)  # Add a horizontal line below the header
# 
# # Process each bar_summary file in sorted order
# for (file in sorted_files) {
#   # Extract parameters values and format them as "%.2g"
#   parameters_values <- sprintf("%.8g", all_processed_data[[file]]$Parameters$Value)
#   
#   # Extract indicators values (median and confidence intervals) and format them as "%.2g"
#   indicators_data <- apply(all_processed_data[[file]]$Indicators[, -1], 1, function(x) {
#     paste0(sprintf("%.8g", x[1]), " ([", sprintf("%.8g", x[2]), ";", sprintf("%.8g", x[3]), "])")
#   })
#   
#   # Combine parameter values and indicator values with "&" between each column
#   data_line <- paste(c(parameters_values, unlist(indicators_data)), collapse = " & ")
#   data_line <- paste(data_line, "\\\\")
#   
#   # Write the data line to LaTeX file
#   cat("        ", data_line, "\n", file = latex_file)
#   
#   # Write hline
#   cat("        \\hline\n", file = latex_file)
# }
# 
# # Write the end of the tabular environment
# cat("    \\end{tabular}\n", file = latex_file)
# cat("    }%\n", file = latex_file)  # End of resizebox
# cat("\\end{table}\n\n", file = latex_file)
# 
# # Write the end of the landscape environment
# cat("\\end{landscape}\n", file = latex_file)
# 
# # Write the end of the document
# cat("\\end{document}\n", file = latex_file)
# 
# # Close the LaTeX file
# close(latex_file)

