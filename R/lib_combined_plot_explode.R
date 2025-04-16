########################################################################### #
# This file is part of the Stochastic Compartmental Model for SARS-COV-2 
# transmission in Belgium, conceived by members of SIMID group during the 
# COVID19 pandemic.
#
# This file contains the main functions to explore stochastic model
# results and many helper functions. The functions called "multi_plot_.XXX" 
# have to be excecuted and all other functions contribute to the actual 
# plotting.
#
# Copyright 2021, SIMID, University of Antwerp & Hasselt University                                        
########################################################################### #

library(scales)
library(RColorBrewer)

################################################################ #
# HELP FUNCTIONS ----
################################################################ #

# set SCM start date
get_scm_start_date <- function(){
  return(as.Date('2020-03-01'))
}

# check if the given file name fits the given col_lib
get_col_index <- function(file_name,col_lib,show_warnings=TRUE){
  flag_file <- unlist(lapply(col_lib$tags,grepl,file_name))
  if(!any(flag_file)) {
    return(NA)
  }
  return(max(which(flag_file)))
}

# check the given color library for ONE given file name
get_col <- function(file_name,col_lib,show_warnings=TRUE){
  col_index <- get_col_index(file_name,col_lib,show_warnings)
  if(is.na(col_index)) {
    if(show_warnings) warning(paste('NO COLOR FOR:',file_name))
    return(8)
  }
  return(col_lib$col[col_index])
}

# check the given color library for multiple file name
get_color_index_multiple <- function(file_names,col_lib,show_warnings=TRUE){
  color_index <- unlist(lapply(file_names,get_col_index,col_lib=col_lib,show_warnings=show_warnings))
  return(color_index)
}

# check the given color library for ONE given file name
get_color_multiple <- function(file_names,col_lib,show_warnings=TRUE){
  color_tag <- unlist(lapply(file_names,get_col,col_lib=col_lib,show_warnings=show_warnings))
  return(color_tag)
}

# check the given color library for multiple file names
get_label_index_multiple <- function(file_names,col_lib,show_warnings=TRUE){
  color_index <- unlist(lapply(file_names,get_label_index,col_lib=col_lib,show_warnings=show_warnings))
  return(color_index)
}

# check if the given file name fits the given col_lib
get_label_index <- function(file_name,col_lib,show_warnings=TRUE){
  flag_file <- unlist(lapply(col_lib$tags,grepl,file_name))
  if(!any(flag_file)) {
    return(NA)
  }
  return(max(which(flag_file)))
}

# add plot legend with given legend library and files
add_scen_legend <- function(col_lib, 
                            opt_files = NA, 
                            ref_pch = 1,
                            ref_label = '',
                            ncol = ncol, 
                            x_axis_scm_day = NA,
                            ref_date = NA){
  
  legend_lib <- col_lib
  
  legend_lib$pch[1] <- ref_pch[1]
  
  if(nchar(ref_label)>0 && grepl('Reported \\(',legend_lib$label[1])){
    legend_lib$label[1] <- paste(substr(legend_lib$label[1],0,nchar('Reported')),
                                  ref_label,
                                  substr(legend_lib$label[1],nchar('Reported')+2,nchar(legend_lib$label[1])))
  }
  
  if(!is.na(ref_date)){
    legend_lib$label[1] <- gsub('(....-..-..)',ref_date,legend_lib$label[1])
  } else{
    legend_lib$label[1] <- gsub('\\(....-..-..\\)','',legend_lib$label[1])
  }
  
  if(!any(is.na(opt_files))){
    legend_col <- unlist(lapply(opt_files,get_col,col_lib))
    legend_lib <- legend_lib[legend_lib$col %in% c(col_lib$col[1],legend_col),]
  }
  
  if(!grepl('Reported',legend_lib$label[1])){
    legend_lib$pch[1] <- NA
  }
  # 
  # add changepoint if present
  if(any(col_lib$tags == 'changepoint')){
    legend_lib <- rbind(legend_lib,
                        col_lib[col_lib$tags == 'changepoint',])
  }

  if(!'lty' %in% names(legend_lib)){
    legend_lib$lty <- 1
  }
  
  legend(ifelse(!any(is.na(x_axis_scm_day)) && min(x_axis_scm_day)==0,'topright','topleft'), # 
         legend = legend_lib$label,
         col    = legend_lib$col,
         lwd    = legend_lib$lwd,
         pch    = legend_lib$pch,
         lty    = legend_lib$lty,
         bg = 'white',
         ncol = ncol,
         cex = 0.6)
}

# mark behavorial change points on the current plot
add_change_points <- function(scm_scenario_cp,
                              scm_num_cp = NA){
    cp      <- get_M_change_day()
    cp_date <- sim_day2date(cp)
    
    if(!is.na(scm_num_cp)){
      cp_date <- cp_date[1:(scm_num_cp+1)] # add one for 2021-03-01
    }

    abline(v=c(rep(cp_date,each=6)+(0:5)),
           lty=1,col=alpha('lightgrey',0.4),lwd=4)
    
    n_scen_cp <- length(scm_scenario_cp)
    if(!any(is.na(scm_scenario_cp))&& n_scen_cp >0){
      abline(v=rep(scm_scenario_cp,each=6)+(0:29), #+(0:5),  temp for one month change
             lty=1,col=alpha('lightgrey',0.6),lwd=3)
      n_scen_cp <- length(scm_scenario_cp)
      text(x = scm_scenario_cp+2,
           y = rep(0,n_scen_cp),
           c('A','B','C','D','E','F')[1:n_scen_cp])
    }
}

################################################################ #
## PLOT RESULTS ON THE POPULATION LEVEL ----
################################################################ #

# NOTE:
# function to obtain all output: multi_plot_incidence_time(...)
# function that does the actual job: plot_incidence_time(...)

multi_plot_incidence_time <- function(output_files, 
                                      output_tag_multi = NA,
                                      bool_polygon_multi = TRUE,
                                      col_lib,
                                      x_axis_scm_day = NA,
                                      be_ref_data = NULL,
                                      default_y_axis = TRUE,
                                      include_date_in_legend = TRUE,
                                      scm_callibration_date = NA,
                                      scm_scenario_cp = NA,
                                      scm_num_cp = NA){

  # set all categories to plot
  if(any(is.na(output_tag_multi))){
    output_tag_multi = c('/hosp_adm',
                         '/hosp_adforcovid',
                         'hosp_load_incr',
                         'icu_load_incr',
                         'Rt_mild_sev_incr',
                         'Rt_hosp_incr',
                         #'/hosp_adm_xtra',
                         'hosp_exit_incr',
                         'new_infect',
                         'new_sympt',
                         'new_infect_xtra',
                         'new_reinfect',
                         '/scen_mortality',
                         '/scen_mortality_xtra',
                         'VOC_mild_alpha',
                         'VOC_mild_delta',
                         'VOC_mild_omicron',
                         'VOC_mild_ba4ba5',
                         'VOC_mild_bq1ba275xbb',
                         'VOC_mild_xbb15',
                         '/nvac_hosp_adm',
                         '/vac_hosp_adm',
                         'new_hosp_adm_vac_d2',
                         'new_hosp_adm_vac_booster',
                         'icu_adm_incr')
  if (bool_vac==FALSE){
      output_tag_multi <- c(output_tag_multi,'delta_hosp_new','delta_hosp_load','delta_icu_load','delta_mortality','delta_mortality_xtra','delta_new_infect','delta_Rt_hosp','delta_Rt_infection')
  }
  }
  
  # with or without polygon?
  if(any(is.na(bool_polygon_multi))){
    bool_polygon_multi = c(F,T)
  }
  
  if(any(is.na(x_axis_scm_day))){
    warning('Using a generic x-axis, since the x-limits contain NA')
  }
  
  # loop over the output tags and polygon options
  output_tag <- output_tag_multi[21];bool_polygon <- bool_polygon_multi[1]
  for(output_tag in output_tag_multi)
    for(bool_polygon in bool_polygon_multi)
      plot_incidence_time(output_files,
                          output_tag     = output_tag,
                          bool_polygon   = bool_polygon, 
                          col_lib        = col_lib, 
                          x_axis_scm_day  = x_axis_scm_day,
                          be_ref_data    = be_ref_data,
                          default_y_axis = default_y_axis,
                          include_date_in_legend = include_date_in_legend,
                          scm_callibration_date = scm_callibration_date,
                          scm_scenario_cp       = scm_scenario_cp,
                          scm_num_cp            = scm_num_cp )
  
}

#output_tag <- 'VOC_mild'; bool_polygon <- T; opt_cp = NA
plot_incidence_time <- function(output_files, 
                                output_tag,
                                bool_polygon,
                                col_lib,
                                x_axis_scm_day = NA,
                                be_ref_data = NULL,
                                default_y_axis = TRUE,
                                include_date_in_legend = TRUE,
                                scm_callibration_date = NA,
                                scm_scenario_cp = NA,
                                bool_legend=TRUE,
                                scm_num_cp = NA){

  # note: output_tags with '_xtra' make use of the same model output (file names without '_xtra') but 
  #       present the results in a different way
  output_files_inc <- output_files[grepl(gsub('_xtra','',output_tag),output_files) & grepl('rds',output_files) ]
  print(output_files_inc)
  # defensive programming check, if the 
  if(length(output_files_inc)==0){
    return(NULL)
  }
  
  # get reference data
  if(is.null(be_ref_data)){
    be_ref_data <- get_latest_incidence_data()
  }
  
  bool_hide_time_axis <- FALSE 
  if(all(is.na(x_axis_scm_day))){
    x_axis_scm_day <- range(0,length(be_ref_data$date))
  } else if(any(is.na(x_axis_scm_day))){
    x_axis_scm_day <- x_axis_scm_day[!is.na(x_axis_scm_day)]
    bool_hide_time_axis <- TRUE
  }
  
  # x-axis limits
  x_ticks         <- pretty(get_scm_start_date() + x_axis_scm_day,16)
  # x_ticks[x_ticks < get_scm_start_date()]      <- get_scm_start_date()
  x_ticks_label   <- format(x_ticks,"%d/%m")
  x_ticks_label   <- format(x_ticks,"%m/%Y")
  #x_axis_label    <- paste0('Time (',paste(unique(format(x_ticks,'%Y')),collapse='-'),')')
  x_axis_label    <- paste0('Time ')
  
  if(bool_hide_time_axis){
    x_ticks_label <- as.numeric(x_ticks - x_ticks[1])
    x_axis_label <- 'Time (days)'
  }
  
  y_abline <- NA
  legend_ncol <- 1
  
  # set reference
  if(any(grepl('hosp_adm',output_files_inc))){
    plot_ref_data <- data.frame(date     = be_ref_data$date,
                                reported = be_ref_data$hospital_admissions + be_ref_data$hospital_admissions_other)
    ref_label    <- 'hospital admissions with COVID19'
    ref_pch      <- 1
    y_lim        <- c(0,ifelse(min(x_axis_scm_day)==0,12500,400))  #700
   # y_abline     <- 100
  }
  if(any(grepl('hosp_adforcovid',output_files_inc))){
    plot_ref_data <- data.frame(date     = be_ref_data$date,
                                reported = be_ref_data$hospital_admissions)
    ref_label    <- 'hospital admissions for COVID19'
    ref_pch      <- 1
    y_lim        <- c(0,ifelse(min(x_axis_scm_day)==0,700,400))  #700
    # y_abline     <- 100
  }
  # set reference
  if(any(grepl('hosp_adm_xtra',output_tag))){
    plot_ref_data <- data.frame(date     = be_ref_data$date,
                                reported = be_ref_data$hospital_admissions)
    ref_label    <- 'hospital admissions for COVID19'
    ref_pch      <- 1
    y_lim        <- c(0,ifelse(min(x_axis_scm_day)==0,800,400))
    y_abline     <- c(65,150)
  }
  if(any(grepl('new_infect',output_files_inc))){
    plot_ref_data <- data.frame(date     = be_ref_data$date,
                                reported = be_ref_data$cases)
    ref_label    <- 'infections'
    ref_pch      <- 5 
    y_lim    <- c(0,ifelse(min(x_axis_scm_day)==0,1.4e5,4e4))
  }
  if(any(grepl('new_sympt',output_files_inc))){
    plot_ref_data <- data.frame(date     = be_ref_data$date,
                                reported = be_ref_data$cases)
    ref_label    <- 'symptomatic infections'
    ref_pch      <- 5 
    y_lim    <- c(0,ifelse(min(x_axis_scm_day)==0,15e4,4e4))
   # x_ticks <- x_ticks[6:10]
  #  x_ticks_label <- x_ticks_label[6:10]
  }
  if(any(grepl('new_infect_xtra',output_tag))){
    plot_ref_data <- data.frame(date     = be_ref_data$date,
                                reported = cumsum(be_ref_data$cases))
    ref_label    <- 'infections cumulative'
    ref_pch      <- 5 
    y_lim    <- c(1e6,ifelse(min(x_axis_scm_day)==0,1e8,1e7))
    default_y_axis <- TRUE
  }
  if(any(grepl('new_reinfect',output_files_inc))){
    plot_ref_data <- data.frame(date     = be_ref_data$date,
                                reported = be_ref_data$cases
                                )
    ref_label    <- 'reinfections'
    ref_pch      <- 5 
    y_lim    <- c(0,ifelse(min(x_axis_scm_day)==0,0.5e5,4e4))
  }
  if(any(grepl('scen_mortality',output_files_inc))){
    plot_ref_data <- data.frame(date     = be_ref_data$date,
                                reported = be_ref_data$covid19_deaths_hospital)
    # remove last 4 observations
    plot_ref_data <- plot_ref_data[1:(nrow(plot_ref_data)-4),]
    
    ref_label    <- 'deaths (excl. nursing home)' 
    ref_pch      <- 1 #4
    y_lim    <- c(0,ifelse(min(x_axis_scm_day)==0,320,20))
  }
  if(any(grepl('scen_mortality_xtra',output_tag))){
    temp <- be_ref_data$covid19_deaths
    temp[is.na(temp)] <- 0
    plot_ref_data <- data.frame(date     = be_ref_data$date,
                                reported = cumsum(temp))
    # remove last 4 observations
    plot_ref_data <- plot_ref_data[1:(nrow(plot_ref_data)-4),]
    
    ref_label    <- 'deaths cumulative' 
    ref_pch      <- 5 #4
    y_lim    <- c(0,ifelse(min(x_axis_scm_day)==0,1.2e5,1e5))
    default_y_axis <- TRUE
  }
  # if(any(grepl('scen_mortality_xtra',output_tag))){
  #   plot_ref_data <- data.frame(date     = be_ref_data$date,
  #                               reported = be_ref_data$covid19_deaths)
  #   # remove last 4 observations
  #   plot_ref_data <- plot_ref_data[1:(nrow(plot_ref_data)-4),]
  #   
  #   ref_label    <- 'deaths (incl. nursing home)' 
  #   ref_pch      <- 1 #4
  #   y_lim    <- c(0,ifelse(min(x_axis_scm_day)==0,350,20))
  #}
  if(any(grepl('hosp_load',output_files_inc))){
    plot_ref_data <- data.frame(date     = be_ref_data$date,
                                reported = be_ref_data$hospital_load)
    y_lim     <- c(1,ifelse(min(x_axis_scm_day)==0,100000,3000)) 
    ref_pch   <- 1#2 
    ref_label <- 'hospital load'
  }
  if(any(grepl('icu_load',output_files_inc))){
    plot_ref_data <- data.frame(date     = be_ref_data$date,
                                reported = be_ref_data$icu_load)
    y_lim     <- c(1,35000) 
    ref_pch   <- 1#3
    ref_label <- 'ICU load'
    # y_abline  <- c(500,600,650,700,750)
    y_abline  <- c(300,500)
  }
  if(any(grepl('icu_adm',output_files_inc))){
    plot_ref_data <- data.frame(date     = be_ref_data$date,
                                reported = be_ref_data$icu_admissions)
    y_lim     <- c(1,120) 
    ref_pch   <- 1#3
    ref_label <- 'ICU admissions'
    # y_abline  <- c(500,600,650,700,750)
    #y_abline  <- c(300,500)
  }
  # set reference
  if(any(grepl('hosp_exit_incr',output_files_inc))){
    plot_ref_data <- data.frame(date     = be_ref_data$date,
                                reported = be_ref_data$hospital_exit)
    ref_label    <- 'hospital exit (recovery + mortality)'
    ref_pch      <- 4
    y_lim        <- c(0,ifelse(min(x_axis_scm_day)==0,800,650))
    y_abline     <- NA
  }
  if(any(grepl('Rt_mild',output_files_inc))){
    plot_ref_data <- data.frame(date     = be_ref_data$date,
                                reported = be_ref_data$Rt_cases)
    y_lim     <- c(0.6,2) 
    ref_pch   <- 6
    ref_label <- 'cases ~ Rt'
    y_abline  <- 1
  }
  if(any(grepl('Rt_hosp',output_files_inc))){
    plot_ref_data <- data.frame(date     = be_ref_data$date,
                                reported = be_ref_data$Rt_hosp)
    y_lim     <- c(0.6,2) 
    ref_pch   <- 6
    ref_label <- 'hospital admissions with COVID ~ Rt'
    y_abline  <- 1
  }
  if(any(grepl('VOC_mild',output_files_inc))){
    # plot_ref_data <- data.frame(date     = as.Date(db_voc_out$date),
    #                             reported = db_voc_out$voc_aggr)
    plot_ref_data <- data.frame(date     = as.Date(db_voc_raw$date),
                                reported = db_voc_raw$voc_aggr / db_voc_raw$n_sequenced)
    #plot_ref_data <- plot_ref_data[1:(nrow(plot_ref_data)-4),]
    y_lim     <- c(0,1.25) 
    ref_pch   <- 1
    ref_label <- 'proportion Alpha+Beta+Gamma VOC'
    y_abline  <- NA
    legend_ncol <- 1
    #col_lib$label[1] <- 'Reported'
  }
  if(any(grepl('VOC_mild_alpha',output_files_inc))){
    # plot_ref_data <- data.frame(date     = as.Date(db_voc_out$date),
    #                             reported = db_voc_out$voc_aggr)
    plot_ref_data <- data.frame(date     = as.Date(db_voc_raw$date),
                                reported = db_voc_raw$voc_aggr / db_voc_raw$n_sequenced)
    #plot_ref_data <- plot_ref_data[1:(nrow(plot_ref_data)-4),]
    y_lim     <- c(0,1.25) 
    ref_pch   <- 1
    ref_label <- 'proportion Alpha+Beta+Gamma VOC'
    y_abline  <- NA
    legend_ncol <- 1
    #col_lib$label[1] <- 'Reported'
  }
  
  if(any(grepl('VOC_mild_delta',output_files_inc))){
    plot_ref_data <- data.frame(date     = as.Date(db_voc_raw$date),
                                reported = db_voc_raw$voc_delta / db_voc_raw$n_sequenced)
    y_abline  <- 0.5
    ref_label <- 'proportion Delta VOC'
  }
  
  if(any(grepl('VOC_mild_omicron',output_files_inc))){
    plot_ref_data <- data.frame(date     = as.Date(c(db_voc_raw$date,
                                                     #db_voc_raw$date,
                                                     db_voc_raw$date)),
                              #  reported = c(db_voc_raw$voc_omicron / db_voc_raw$n_sequenced)
                                             reported = c((db_voc_raw$voc_omicron) / db_voc_raw$n_sequenced)
                                                          #db_voc_raw$voc_omicron_active / db_voc_raw$n_sequenced_active,
                                            # db_voc_raw$sftg_logistic_prob_omicron
                                            
                                )
    y_abline  <- 0.5
    ref_pch   <- c(rep(1,nrow(db_voc_raw)),
                   #rep(2,nrow(db_voc_raw)),
                   rep(3,nrow(db_voc_raw)))
    ref_label <- 'proportion Omicron BA1/BA2 VOC'
  }
  if(any(grepl('VOC_mild_ba4ba5',output_files_inc))){
    plot_ref_data <- data.frame(date     = as.Date(c(db_voc_raw$date)),
                                reported = c(db_voc_raw$voc_ba4ba5 / db_voc_raw$n_sequenced)
    )
    y_abline  <- 0.5
    ref_pch   <- c(rep(1,nrow(db_voc_raw)))
    ref_label <- 'proportion Omicron BA4/BA5 VOC'
  }
  if(any(grepl('VOC_mild_bq1ba275xbb',output_files_inc))){
    plot_ref_data <- data.frame(date     = as.Date(c(db_voc_raw$date)),
                                reported = c(db_voc_raw$voc_bq1ba275xbb / db_voc_raw$n_sequenced)
    )
    y_abline  <- 0.5
    ref_pch   <- c(rep(1,nrow(db_voc_raw)))
    ref_label <- 'proportion Omicron BQ.1/BA.2.75/XBB VOC'
  }
  
  if(any(grepl('VOC_mild_xbb15',output_files_inc))){
    plot_ref_data <- data.frame(date     = as.Date(c(db_voc_raw$date)),
                                reported = c(db_voc_raw$voc_xbb15 / db_voc_raw$n_sequenced)
    )
    y_abline  <- 0.5
    ref_pch   <- c(rep(1,nrow(db_voc_raw)))
    ref_label <- 'proportion Omicron XBB.1.5 VOC'
  }
  
  if(grepl('/nvac_hosp_adm',output_tag)){
    ref_label <- paste(ref_label,'(not vaccinated)')
  }
  
  if(grepl('new_hosp_adm_vac_d2',output_tag)){
    #ref_label <- paste(ref_label,'(vaccinated (2 doses))')
    ref_label <- paste(ref_label,'(vaccinated (all >=2 doses))')
  }
  
  if(grepl('new_hosp_adm_vac_booster',output_tag)){
   # ref_label <- paste(ref_label,'(vaccinated 2 doses + booster)')
    ref_label <- paste(ref_label,'(from reinfection)')
  }
  if(any(grepl('delta_icu_load',output_files_inc))){
    plot_ref_data <- data.frame(date     = be_ref_data$date,
                                reported = rep(NA,length(be_ref_data$date)))
    y_lim     <- c(-2000,2000) 
    ref_pch   <- 1#3
    ref_label <- 'averted ICU load'
    # y_abline  <- c(500,600,650,700,750)
    y_abline  <- c(0)
  }
  if(any(grepl('delta_hosp_load',output_files_inc))){
    plot_ref_data <- data.frame(date     = be_ref_data$date,
                                reported = rep(NA,length(be_ref_data$date)))
    y_lim     <- c(-2000,2000) 
    ref_pch   <- 1#3
    ref_label <- 'averted hospital load'
    # y_abline  <- c(500,600,650,700,750)
    y_abline  <- c(0)
  }
  if(any(grepl('delta_hosp_new',output_files_inc))){
    plot_ref_data <- data.frame(date     = be_ref_data$date,
                                reported = rep(NA,length(be_ref_data$date)))
    y_lim     <- c(-1000,1000) 
    ref_pch   <- 1#3
    ref_label <- 'averted hospital admissions'
    # y_abline  <- c(500,600,650,700,750)
    y_abline  <- c(0)
  }
  if(any(grepl('delta_mortality',output_files_inc))){
    plot_ref_data <- data.frame(date     = be_ref_data$date,
                                reported = rep(NA,length(be_ref_data$date)))
    y_lim     <- c(-200,200) 
    ref_pch   <- 1#3
    ref_label <- 'averted deaths'
    # y_abline  <- c(500,600,650,700,750)
    y_abline  <- c(0)
  }
  if(any(grepl('delta_mortality_xtra',output_tag))){
    ref_label <- 'averted deaths cumulative'
    y_lim     <- c(-200,50000) 
  }
  if(any(grepl('delta_new_infect',output_files_inc))){
    plot_ref_data <- data.frame(date     = be_ref_data$date,
                                reported = rep(NA,length(be_ref_data$date)))
    y_lim     <- c(-10000,10000) 
    ref_pch   <- 1#3
    ref_label <- 'averted new infections'
    # y_abline  <- c(500,600,650,700,750)
    y_abline  <- c(0)
  }
  if(any(grepl('delta_Rt_hosp',output_files_inc))){
    plot_ref_data <- data.frame(date     = be_ref_data$date,
                                reported = rep(NA,length(be_ref_data$date)))
    y_lim     <- c(-2,2) 
    ref_pch   <- 1#3
    ref_label <- 'Rt difference (hospital)'
    # y_abline  <- c(500,600,650,700,750)
    y_abline  <- c(0)
  }
  if(any(grepl('delta_Rt_infection',output_files_inc))){
    plot_ref_data <- data.frame(date     = be_ref_data$date,
                                reported = rep(NA,length(be_ref_data$date)))
    y_lim     <- c(-2,2) 
    ref_pch   <- 1#3
    ref_label <- 'Rt difference (infection)'
    # y_abline  <- c(500,600,650,700,750)
    y_abline  <- c(0)
  }
  if(!default_y_axis){
    tmp_ref_data <- plot_ref_data$reported[plot_ref_data$date %in% min(x_ticks):max(x_ticks)]
    if(length(tmp_ref_data)>0 & sum(tmp_ref_data,na.rm=T)>0){
      y_lim <- range(c(0,tmp_ref_data*0.9,
                       tmp_ref_data*1.7),
                     na.rm=T)
      if(grepl('Rt',output_tag)){
        y_lim <- range(c(tmp_ref_data*0.8,
                         tmp_ref_data*1.15),
                       na.rm=T)
      }
      if(grepl('load',output_tag) | grepl('mort',output_tag)){
        y_lim <- range(c(0,tmp_ref_data*0.9,
                         tmp_ref_data*2.1),
                       na.rm=T)
      }
    }
  }
  
  if(any(grepl('_stpm',output_files_inc))){
    ref_label <- paste(ref_label,'(STPM)')
  }
    
  plot(x = plot_ref_data$date,
       y = plot_ref_data$reported,
       # xlim = range(x_ticks),
       xlim = get_scm_start_date() + x_axis_scm_day,
       ylim = y_lim,
       xaxt='n',
       xlab=x_axis_label, 
       ylab="",
       col=0,
       las=2)
  if(y_lim[2]>=10000){
  title(ylab = paste('Daily',ref_label), mgp = c(4, 1, 0))
  } else {
    title(ylab = paste('Daily',ref_label), mgp = c(3, 1, 0))
  }
  axis(1,x_ticks,x_ticks_label,cex.axis=0.9)
  abline(v=x_ticks,col='grey',lty=3)
  grid(nx=NA,ny=NULL)
  if(!any(is.na(scm_callibration_date))){
    add_vertical_line(scm_callibration_date,bool_text =T,date_tag= 'callibration',pos_factor=1.30)
  }
  abline(h=y_abline,lty=2)
  add_change_points(scm_scenario_cp = scm_scenario_cp, scm_num_cp = scm_num_cp)
#temp removed !!!!!!!!
  i_file <- 1
  for(i_file in 1:length(output_files_inc)){
    scen_data <- readRDS(output_files_inc[i_file])
    scen_col  <- get_col(output_files_inc[i_file],col_lib)
    scm_dates  <- get_scm_start_date() +  0:(nrow(scen_data)-1)
    
  if(bool_polygon){
      hosp_aggr <-  data.frame(admin_mean=apply(scen_data,1,mean),
                               #admin_mean=NA,
                               admin_LL = apply(scen_data,1,quantile,0.025,na.rm=T),
                               admin_UL = apply(scen_data,1,quantile,0.975,na.rm=T))

      if(any(grepl('new_infect_xtra',output_tag))){
        hosp_aggr$admin_mean <- cumsum(hosp_aggr$admin_mean)
        hosp_aggr$admin_LL   <- cumsum(hosp_aggr$admin_LL)
        hosp_aggr$admin_UL   <- cumsum(hosp_aggr$admin_UL)
      }
      
      if(any(grepl('scen_mortality_xtra',output_tag))){
        hosp_aggr$admin_mean <- cumsum(hosp_aggr$admin_mean)
        hosp_aggr$admin_LL   <- cumsum(hosp_aggr$admin_LL)
        hosp_aggr$admin_UL   <- cumsum(hosp_aggr$admin_UL)
      }
      
      if(any(grepl('delta_mortality_xtra',output_tag))){
        hosp_aggr$admin_mean <- cumsum(hosp_aggr$admin_mean)
        hosp_aggr$admin_LL   <- cumsum(hosp_aggr$admin_LL)
        hosp_aggr$admin_UL   <- cumsum(hosp_aggr$admin_UL)
      }
      
      # remove all results beyond the given date ranges
      flag_na <- scm_dates > max(sim_day2date(x_axis_scm_day))
      hosp_aggr[flag_na,] <- NA
      
      #specific add to correct for covid
     # if(any(grepl('hosp_adm_xtra',output_tag))){
       #  pop_data_be  <- get_regional_pop("belgium")
       #  output_files_inc2 <- output_files[grepl('prev_I_asymp',output_files) & grepl('rds',output_files) ]
       #  scen_data_all <- readRDS(output_files_inc2[i_file])
       #  scen_data2 <-  scen_data_all[,,6] + scen_data_all[,,7] + scen_data_all[,,8] + scen_data_all[,,9] + scen_data_all[,,10] 
       # pop60 <- pop_data_be[6] + pop_data_be[7] + pop_data_be[8] + pop_data_be[9] + pop_data_be[10]
       #   hosp_aggr$admin_mean <- hosp_aggr$admin_mean - 18000*apply(scen_data2,2,mean)/pop60
       #  hosp_aggr$admin_LL   <- hosp_aggr$admin_LL - 18000*apply(scen_data2,2,quantile,0.975,na.rm=T)/pop60
       #  hosp_aggr$admin_UL   <- hosp_aggr$admin_UL - 18000*apply(scen_data2,2,quantile,0.025,na.rm=T)/pop60
       #  hosp_aggr$admin_mean[hosp_aggr$admin_mean < 0] <- 0
       #  hosp_aggr$admin_LL[hosp_aggr$admin_LL < 0] <- 0
       #  hosp_aggr$admin_UL[hosp_aggr$admin_UL < 0] <- 0
      #}
      # if(any(grepl('hosp_adm',output_tag))){
      #   print(output_files_inc[i_file])
      #   print(max(hosp_aggr$admin_mean[scm_dates > as.Date('2022-02-01')]))
      # }
      
      
      
      # if(any(grepl('Rt',output_tag))){
      #   flag_ignore <- hosp_aggr$admin_mean < 1 & (hosp_aggr$admin_UL - hosp_aggr$admin_LL) > 0.3
      #   flag_ignore[is.na(flag_ignore)] <- FALSE
      #   flag_ignore[cumsum(flag_ignore)>5] <- TRUE
      #   flag_ignore[cumsum(flag_ignore)<5] <- FALSE
      #   hosp_aggr[flag_ignore,] <- NA
      #   
      #   tail(hosp_aggr)
      # }
      
      add_polygon(scenario_data  = hosp_aggr[-nrow(hosp_aggr),],
                  scenario_dates = scm_dates[-length(scm_dates)],
                  col = scen_col,
                  alpha_val = 0.15,
                  mean_lty = ref_pch)
    } else{
      lines(x   = scm_dates,
            y   = apply(scen_data,1,mean),
            lwd=3,
            col = scen_col)
    }
  }
  
  # reference data
  points(x = plot_ref_data$date,
         y = plot_ref_data$reported,
         pch=ref_pch,
         cex = 0.8)

  # add 7d average for hospital exit data
  if(any(grepl('hosp_exit_incr',output_files_inc))){
    lines(x = plot_ref_data$date,
           y = rollmean(plot_ref_data$reported,k=7,align = 'center',fill = NA),
           pch=ref_pch,
           lwd = 2)
  }

  if(grepl('new_infect',output_tag) || grepl('new_reinfect',output_tag) || grepl('new_sympt',output_tag)){
    ref_label <- 'positive cases'
  }
  
  add_scen_legend(col_lib, 
                  output_files_inc, 
                  ref_pch = ref_pch, 
                  ref_label=ref_label,
                  ncol=legend_ncol,
                  x_axis_scm_day=x_axis_scm_day,
                  ref_date = ifelse(include_date_in_legend,as.character(max(plot_ref_data$date[!is.na(plot_ref_data$reported)])),NA))

  add_copyright()
}





################################################################ #
## PLOT RESULTS BY AGE ----
################################################################ #

# NOTE:
# function to obtain all output: multi_plot_incidence_age_time(...)
# function that does the actual job: plot_incidence_age_time(...)

multi_plot_incidence_age_time <- function(output_files,
                                          col_lib,
                                          x_axis_scm_day,
                                          x_axis_vaccine){
  
  ## BY AGE
  plot_incidence_age_time(output_files,'scenario1', col_lib = col_lib, x_axis_scm_day = x_axis_scm_day)
  plot_incidence_age_time(output_files,'mortality_age',col_lib = col_lib, x_axis_scm_day = x_axis_scm_day)
  plot_incidence_age_time(output_files,'scen_cases',col_lib = col_lib, x_axis_scm_day = x_axis_scm_day)
  
  plot_incidence_age_time(output_files,'cumul_cases',col_lib = col_lib, x_axis_scm_day = x_axis_scm_day)
  plot_incidence_age_time(output_files,'mean_susceptible_full',col_lib = col_lib, x_axis_scm_day = x_axis_scm_day)
  plot_incidence_age_time(output_files,'mean_susceptible_initpluswaning',col_lib = col_lib, x_axis_scm_day = x_axis_scm_day)
  plot_incidence_age_time(output_files,'mean_susceptible_vac_d1',col_lib = col_lib, x_axis_scm_day = x_axis_vaccine)
  plot_incidence_age_time(output_files,'mean_susceptible_vac_d2',col_lib = col_lib, x_axis_scm_day = x_axis_vaccine)
  plot_incidence_age_time(output_files,'mean_vaccinated_age',col_lib = col_lib, x_axis_scm_day = x_axis_vaccine)
  
  plot_incidence_age_time(output_files,'mean_rna1_vaccinated',col_lib = col_lib, x_axis_scm_day = x_axis_vaccine,y_lab_extra = ' (mRNA 1 dose)')
  plot_incidence_age_time(output_files,'mean_rna2_vaccinated',col_lib = col_lib, x_axis_scm_day = x_axis_vaccine,y_lab_extra = ' (mRNA 2 doses)')
  plot_incidence_age_time(output_files,'mean_adeno1_vaccinated',col_lib = col_lib, x_axis_scm_day = x_axis_vaccine,y_lab_extra = ' (AZ 1 dose)')
  plot_incidence_age_time(output_files,'mean_adeno2_vaccinated',col_lib = col_lib, x_axis_scm_day = x_axis_vaccine,y_lab_extra = ' (AZ 2 doses)')
  plot_incidence_age_time(output_files,'mean_waning_vaccinated',col_lib = col_lib, x_axis_scm_day = x_axis_vaccine,y_lab_extra = ' (booster dose)')
  
  plot_incidence_age_time(output_files,'/nvac_age_hosp_adm',col_lib = col_lib, x_axis_scm_day = x_axis_vaccine)
  plot_incidence_age_time(output_files,'/vac_age_hosp_adm',col_lib = col_lib, x_axis_scm_day = x_axis_vaccine)
  plot_incidence_age_time(output_files,'/new_hosp_adm_vac_age_d2',col_lib = col_lib, x_axis_scm_day = x_axis_vaccine)

}


#output_tag <- 'mean_rna1_vaccinated'
plot_incidence_age_time <- function(output_files,
                                    output_tag, 
                                    col_lib,
                                    x_axis_scm_day,
                                    scm_scenario_cp = NA, 
                                    opt_cp         = NA,
                                    scm_callibration_date = NA,
                                    y_lab_extra = ""){
  
  output_files_inc <- output_files[grepl(output_tag,output_files) & grepl('rds',output_files) ]
  
  if(length(output_files_inc)==0){
    return(NULL)
  }
  
  bool_hide_time_axis <- FALSE 
  if(all(is.na(x_axis_scm_day))){
    x_axis_scm_day <- range(0,length(be_ref_data$date))
  } else if(any(is.na(x_axis_scm_day))){
    x_axis_scm_day <- x_axis_scm_day[!is.na(x_axis_scm_day)]
    bool_hide_time_axis <- TRUE
  }
  
  # x-axis limits
  x_ticks         <- pretty(get_scm_start_date() + x_axis_scm_day,12)
  x_ticks[x_ticks < get_scm_start_date()]      <- get_scm_start_date()
  x_ticks_label   <- format(x_ticks,"%m%Y")
  #x_axis_label    <- paste0('Time (',paste(unique(format(x_ticks,'%Y')),collapse='-'),')')
  x_axis_label    <- paste0('Time ')
  
  if(bool_hide_time_axis){
    x_ticks_label <- as.numeric(x_ticks - x_ticks[1])
    x_axis_label <- 'Time (days)'
  }
  
  y_lim_bar  <- c(0,400) 
  y_lim_line <- range(0,600)
  y_lab      <- 'Daily hospital admissions'
  
  if(grepl('mortality',output_tag)){
    y_lim_bar <- c(0,190)
    y_lim_line <- c(0,150)
    y_lab <- 'Daily mortality'
  }
  
  if(grepl('scen_cases',output_tag)){
    y_lim_bar <- c(0,40e3)
    y_lim_line <- c(0,15e3)
    y_lab <- 'Daily cases'
  }
  if(grepl('cumul_cases',output_tag)){
    y_lim_bar <- c(0,20e6)
    y_lim_line <- c(0,1e7)
    y_lab <- 'Total cases (cumulative)'
  }
  if(grepl('susceptible_full',output_tag)){
    y_lim_bar <- c(0,11e6)
    y_lim_line <- c(0,1e7)
    y_lab <- 'Susceptible (including vaccinated)'
  }
  if(grepl('mean_susceptible_initpluswaning',output_tag)){
    y_lim_bar <- c(0,11e6)
    y_lim_line <- c(0,1e7)
    y_lab <- 'Susceptible (only initial + waning)'
  }
  if(grepl('susceptible_vac',output_tag)){
    y_lim_bar <- c(0,11e6)
    y_lim_line <- c(0,1e7)
    y_lab <- 'Susceptible (vaccinated)'
  }
  
  if(grepl('d1',output_tag)){
    y_lab <- gsub(')',' 1 dose)',y_lab)
  } 
  if(grepl('d2',output_tag)){
    y_lab <- gsub(')',' >=2 doses)',y_lab)
  }
  
  if(grepl('vaccinated',output_tag)){
    y_lim_bar <- c(0,11e6)
    y_lim_line <- c(0,4e6)
    y_lab <- 'Vaccine protection'
  }
  
  if(grepl('mean_rna',output_tag) | grepl('mean_adeno',output_tag)){
    y_lim_bar <- c(0,6e6)
    y_lim_line <- c(0,2e6)
    y_lab <- 'Vaccine protection'
  }
  
  if(grepl('/vac_age_hosp_adm',output_tag)){
    y_lim_bar  <- c(0,5e2) 
    y_lim_line <- range(0,200)
    y_lab <- paste(y_lab,'(vaccinated)')
  }
  if(grepl('/nvac_age_hosp_adm',output_tag)){
    y_lim_bar  <- c(0,5e2) 
    y_lim_line <- range(0,400)
    y_lab <- paste(y_lab,'(not vaccinated)')
  }
  
  if(grepl('/new_hosp_adm_vac_age_d2',output_tag)){
    y_lim_bar  <- c(0,5e2) 
    y_lim_line <- range(0,200)
    y_lab <- paste(y_lab,'(vaccinated 2 doses)')
  }
  
  # optional addition to y_lab
  y_lab <- paste0(y_lab,y_lab_extra)
  
  # if(grepl('load',output_tag)){
  #   y_lim_bar <- c(0,2000)
  #   y_lim_line <- c(0,600)
  #   y_lab <- 'Daily hospital load'
  # }
  
  i_file <- 1
  for(i_file in 1:length(output_files_inc)){
    
    scen_data_all <- readRDS(output_files_inc[i_file])
    scm_dates      <- get_scm_start_date() +  0:(dim(scen_data_all)[2]-1)
    scen_label    <- col_lib$label[col_lib$col %in% get_col(output_files_inc[i_file],col_lib)]
    
    # adjust date range
    x_ticks_day    <- sim_date2day(x_ticks)
    sel_sim_day    <- min(max(x_ticks_day),dim(scen_data_all)[2])
    scen_data_all  <- scen_data_all[,min(x_ticks_day):sel_sim_day,]
    scm_dates       <- scm_dates[min(x_ticks_day):sel_sim_day]
    
    # remove day id
    scen_data_all <- scen_data_all[,,-1]
    
    # aggregate age groups
    # scen_data_all[,,1] <- scen_data_all[,,1] + scen_data_all[,,2] + scen_data_all[,,3] + scen_data_all[,,4]  # 0-39
    #scen_data_all[,,1] <- scen_data_all[,,1] + scen_data_all[,,2]  # 0-19
    scen_data_all[,,3] <- scen_data_all[,,3] + scen_data_all[,,4]  # 20-39
    scen_data_all[,,5] <- scen_data_all[,,5] + scen_data_all[,,6]  # 40-59
    scen_data_all[,,7] <- scen_data_all[,,7] + scen_data_all[,,8]  # 60-79
    scen_data_all[,,9] <- scen_data_all[,,9] + scen_data_all[,,10] # 80+
    #scen_data_all <- scen_data_all[,,-c(2,4,6,8,10)]
    scen_data_all <- scen_data_all[,,-c(4,6,8,10)]
    # scen_data_all <- scen_data_all[,,-c(2,3,4,6,8,10)]
    
    # get age-specific mean and CI
    scen_data <- apply(scen_data_all[,,],3,apply,2,mean)
    scen_data_p025 <- apply(scen_data_all[,,],3,apply,2,quantile,0.025)
    scen_data_p975 <- apply(scen_data_all[,,],3,apply,2,quantile,0.975)
    
    dim(scen_data)
    
    # set age groups
    # opt_age <- c('0-19','20-39','40-59', gsub('_','-',get_age_groups())[7:10])
    opt_age <- c('0-19','20-39','40-59', '60-79','+80')
    opt_age <- c('0-9','10-19','20-39','40-59', '60-79','+80')
    opt_col <- c('greenyellow','darkgreen','gray50','coral4','goldenrod1')
    opt_col <- brewer.pal(length(opt_age),'Paired') #[-1]
    
    # opt_age <- c('0-39','40-59', '60-79','+80')
    # opt_col <- c('darkgreen','gray50','coral4','goldenrod1')
    
    bplot <- barplot(t(scen_data),
                     col=opt_col,
                     border = NA,
                     width = 1,
                     ylim = y_lim_bar,
                     xlab=x_axis_label,
                     yaxt='n',
                     main = scen_label,
                     ylab=y_lab)
    axis(1,(x_ticks_day-min(x_ticks_day))*unique(diff(bplot))[1],x_ticks_label,cex.axis=0.9)
    abline(v=(x_ticks_day-min(x_ticks_day))*unique(diff(bplot))[1],col='grey',lty=3)
    grid(nx=NA,ny=NULL)
    axis(2,pretty(y_lim_bar,15),pretty(y_lim_bar,15),las=2,cex.axis=0.7)
    legend('topleft',
           rev(opt_age),
           fill = rev(opt_col),
           cex=0.8,
           title='Age',
           bg='white'
           )
  
    # # optional: save values as csv ?
    # x_out <- cbind(as.character(sim_day2date(min(x_ticks_day):max(x_ticks_day))),scen_data)
    # colnames(x_out) <- c('date', paste0('age_',opt_age))
    # write.table(x_out,
    #             file=gsub('.rds',paste0('_',output_tag,'_age.csv'),output_files_inc[i_file]),sep=',',row.names=F)
    
    # line plot
    plot(scm_dates,
         scen_data[,1],
         col=0,
         ylim = y_lim_line,
         xlab=x_axis_label,
         xaxt='n',
         yaxt='n',
         type='l',
         lwd=2,
         las=2,
         main = scen_label,
         cex.axis=0.9,
         ylab=y_lab)
    grid(nx=NA,ny=NULL)
    axis(2,pretty(y_lim_line),format(pretty(y_lim_line),scientific=T),las=2,cex.axis=0.7)
    axis(1,x_ticks,x_ticks_label,cex.axis=0.9)
    abline(h=0,lty=3,col='grey')
    add_vertical_line(scm_callibration_date,bool_text =T,date_tag= 'callibration')

    for(i_age in 1:ncol(scen_data)){
      scenario_data <- data.frame(inc_mean = scen_data[,i_age],
                            inc_LL = scen_data_p025[,i_age],
                            inc_UL = scen_data_p975[,i_age])
      add_polygon(scenario_data = scenario_data,
                  scenario_dates = scm_dates,
                  col = opt_col[i_age],
                  alpha_val = 0.3)
      
    }
    legend('topleft',
           rev(opt_age),
           fill = rev(opt_col),
           cex=0.8,
           title='Age',
           bg='white'
    )
    
    add_copyright()
  }

}

################################################################ #
## PLOT CUMULATIVE RESULTS OVER TIME ----
################################################################ #

# NOTE:
# function to obtain all output: multi_plot_incidence_age_time(...)
# function that does the actual job: plot_incidence_age_time(...)

multi_plot_cumulative_time <- function(output_files, 
                                        output_tag_multi = NA,
                                        col_lib,
                                        x_axis_scm_day = NA,
                                        be_ref_data = NULL,
                                        scen_reference_date = '2021-09-01',
                                        scm_callibration_date = NA,
                                        scm_scenario_cp = NA,
                                        scm_num_cp = NA){
  
  # set all categories to plot
  if(any(is.na(output_tag_multi))){
    output_tag_multi = c('/hosp_adm',
                         'hosp_load_incr',
                         'icu_load_incr',
                         '/scen_mortality',
                         'new_infect',
                         'new_sympt',
                         'new_infect_xtra',
                         'new_reinfect')
  }
  
  if(any(is.na(x_axis_scm_day))){
    warning('Using a generic x-axis, since the x-limits contain NA')
  }
  
  # loop over the output tags and polygon options
  output_tag <- output_tag_multi[4]
  for(output_tag in output_tag_multi){
    plot_cumulative_time(output_files,
                         output_tag     = output_tag,
                         col_lib        = col_lib, 
                         x_axis_scm_day  = x_axis_scm_day,
                         scen_reference_date  = scen_reference_date,
                         scm_callibration_date = scm_callibration_date,
                         scm_scenario_cp       = scm_scenario_cp,
                         scm_num_cp            = scm_num_cp )
  }
}


# output_tag <- '/hosp_adm'
plot_cumulative_time <- function(output_files,
                                    output_tag, 
                                    col_lib,
                                    x_axis_scm_day,
                                    # be_ref_data = NULL,
                                    scen_reference_date = NA,
                                    scm_callibration_date = NA,
                                    scm_scenario_cp = NA, 
                                    scm_num_cp         = NA,
                                    y_lab_extra = ""){
  
  # note: output_tags with '_xtra' make use of the same model output (file names without '_xtra') but 
  #       present the results in a different way
  output_files_inc <- output_files[grepl(gsub('_xtra','',output_tag),output_files) & grepl('rds',output_files) ]
  
  # defensive programming check, if the 
  if(length(output_files_inc)==0){
    return(NULL)
  }

  bool_hide_time_axis <- FALSE 
  if(all(is.na(x_axis_scm_day))){
    x_axis_scm_day <- range(0,length(be_ref_data$date))
  } else if(any(is.na(x_axis_scm_day))){
    x_axis_scm_day <- x_axis_scm_day[!is.na(x_axis_scm_day)]
    bool_hide_time_axis <- TRUE
  }
  
  # x-axis limits
  x_ticks         <- pretty(get_scm_start_date() + x_axis_scm_day,12)
  x_ticks[x_ticks < get_scm_start_date()]      <- get_scm_start_date()
  x_ticks_label   <- format(x_ticks,"%d/%m")
  x_ticks_label   <- format(x_ticks,"%Y")
  #x_axis_label    <- paste0('Time (',paste(unique(format(x_ticks,'%Y')),collapse='-'),')')
  x_axis_label    <- paste0('Time ')
  
  if(bool_hide_time_axis){
    x_ticks_label <- as.numeric(x_ticks - x_ticks[1])
    x_axis_label <- 'Time (days)'
  }
  
  y_abline <- NA
  legend_ncol <- 1
  
  # set plot details
  if(any(grepl('hosp_adm',output_files_inc))){
    plot_ref_data <- data.frame(date     = be_ref_data$date,
                                reported = be_ref_data$hospital_admissions + be_ref_data$hospital_admissions_other)
    ref_label    <- 'hospital adm. with COVID19'
    y_lim        <- c(ifelse(min(x_axis_scm_day)==0,0,8e4),1.5e5)  #700
  }
  if(any(grepl('hosp_adforcovid',output_files_inc))){
    ref_label    <- 'hospital admissions for COVID19 (estimated)'
    y_lim        <- c(ifelse(min(x_axis_scm_day)==0,0,8e4),1000)  #700
  }
  if(any(grepl('hosp_load',output_files_inc))){
    plot_ref_data <- data.frame(date     = be_ref_data$date,
                                reported = be_ref_data$hospital_load)
    
    y_lim     <- c(ifelse(min(x_axis_scm_day)==0,0,1e6),1.6e6) 
    ref_label <- 'hospital load'
  }
  if(any(grepl('icu_load',output_files_inc))){
    plot_ref_data <- data.frame(date     = be_ref_data$date,
                                reported = be_ref_data$icu_load)
    y_lim     <- c(ifelse(min(x_axis_scm_day)==0,0,2.3e5),3.5e5) 
    ref_label <- 'ICU load'
  }
  if(any(grepl('new_infect',output_files_inc))){
    ref_label    <- 'infections'
    y_lim    <- c(ifelse(min(x_axis_scm_day)==0,0,4e6),1.4e7)
  }
  if(any(grepl('new_sympt',output_files_inc))){
    ref_label    <- 'symptomatic infections'
    y_lim    <- c(ifelse(min(x_axis_scm_day)==0,0,4e6),1.7e7)
  }
  if(any(grepl('new_reinfect',output_files_inc))){
    ref_label    <- 'reinfections'
    y_lim    <- c(ifelse(min(x_axis_scm_day)==0,0,1e3),3e5)
  }
  if(any(grepl('mortality',output_files_inc))){
    plot_ref_data <- data.frame(date = be_ref_data$date,
                                reported = be_ref_data$covid19_deaths)
    ref_label    <- 'mortality' 
    y_lim    <- c(ifelse(min(x_axis_scm_day)==0,0,2.5e4),3.5e4)
  }

  
  if(!is.na(scen_reference_date)){
    scen_reference_data <- readRDS(output_files_inc[grepl(col_lib$tags[2],output_files_inc)])    
    
    scen_reference_day <- sim_date2day(scen_reference_date)
    scen_reference_day <- ifelse(scen_reference_day < 1,1,scen_reference_day)
    scen_reference_day <- ifelse(scen_reference_day > nrow(scen_reference_data),nrow(scen_reference_data),scen_reference_day)
  
    scen_reference_data <- scen_reference_data[-(1:scen_reference_day),]
    scen_reference_data <- apply(scen_reference_data,2,cumsum)
    scen_reference_aggr <- apply(scen_reference_data,1,mean)
    
    y_lim <- range(scen_reference_aggr[x_axis_scm_day-scen_reference_day]) * c(1,1.3)
    
    if(exists('plot_ref_data')){
      plot_ref_data <- plot_ref_data[plot_ref_data$date >= scen_reference_date,]
    }
    # scen_reference_start <- scen_reference_aggr[scen_reference_day]
    # scen_reference_range <- scen_reference_aggr[length(scen_reference_aggr)] - scen_reference_start
    # scen_reference_min   <- min(scen_reference_aggr[x_axis_scm_day]*0.9-scen_reference_start)
    # 
    # y_lim <- c(scen_reference_min,(scen_reference_range*1.1))
    
  } else{
    scen_reference_day <- 1
  }
  
  plot(x = range(x_ticks),
       y = y_lim,
       xlim = range(x_ticks),
       ylim = y_lim,
       xaxt='n',
       xlab=x_axis_label, 
       ylab="",
       col=0,
       las=2)
  if(y_lim[2]>=10000){
     title(ylab = paste('Cumulative',ref_label), mgp = c(4, 1, 0))
   } else {
    title(ylab = paste('Cumulative',ref_label), mgp = c(3, 1, 0))
   }
  axis(1,x_ticks,x_ticks_label,cex.axis=0.9)
  abline(v=x_ticks,col='grey',lty=3)
  grid(nx=NA,ny=NULL)
  if(!is.na(scm_callibration_date)){
    add_vertical_line(scm_callibration_date,bool_text =T,date_tag= 'callibration',pos_factor=1.30)
  }
  abline(h=y_abline,lty=2)
  add_change_points(scm_scenario_cp = scm_scenario_cp, scm_num_cp = scm_num_cp)
  
  if(!is.na(scen_reference_date)){
    add_vertical_line(scen_reference_date,bool_text =T,date_tag= 'reference',pos_factor=1.00)
  }
  
  if(exists('plot_ref_data')){
    plot_ref_data$reported[is.na(plot_ref_data$reported)] <- 0
    plot_ref_data$reported_cum <- cumsum(plot_ref_data$reported)
    
    lines(plot_ref_data$date,
          plot_ref_data$reported_cum,
          lty=3)
    col_lib$lty <- 1
    col_lib$lwd[grepl('xxx',col_lib$tags)] <- 1
    col_lib$pch[grepl('xxx',col_lib$tags)] <- NA
    col_lib$lty[grepl('xxx',col_lib$tags)] <- 3
    col_lib$label[1] <- paste('Reported',ref_label,'since',min(plot_ref_data$date))
  } else{
    col_lib <- col_lib[!grepl('xxx',col_lib$tags ),]
  }
  
  i_file <- 4
  for(i_file in 1:length(output_files_inc)){
    scen_data <- readRDS(output_files_inc[i_file])
    scen_col  <- get_col(output_files_inc[i_file],col_lib)
    scm_dates  <- get_scm_start_date() +  0:(nrow(scen_data)-1)
    dim(scen_data)
    
    scen_data <- scen_data[scen_reference_day:nrow(scen_data),]
    scm_dates <- scm_dates[scen_reference_day:length(scm_dates)]
    dim(scen_data)
    
    scen_data <- apply(scen_data,2,cumsum)
    
    hosp_aggr <-  data.frame(admin_mean=apply(scen_data,1,mean),
                               #admin_mean=NA,
                               admin_LL = apply(scen_data,1,quantile,0.025,na.rm=T),
                               admin_UL = apply(scen_data,1,quantile,0.975,na.rm=T))
    
    add_polygon(scenario_data  = hosp_aggr[-nrow(hosp_aggr),],
                scenario_dates = scm_dates[-length(scm_dates)],
                col = scen_col,
                alpha_val = 0.15,
                mean_lty = 1)
  }
  
  
  if(grepl('new_infect',output_tag) || grepl('new_reinfect',output_tag) || grepl('new_sympt',output_tag)){
    ref_label <- 'positive cases'
  }
  
  add_scen_legend(col_lib, 
                  output_files_inc, 
                  ref_pch = NA, 
                  ref_label=ref_label,
                  ncol=legend_ncol,
                  x_axis_scm_day= NA, # alwasy on the left-hand side 
                  ref_date = NA)
  
  add_copyright()
}




