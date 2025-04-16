########################################################################### #
# This file is part of the Stochastic Compartmental Model for SARS-COV-2 
# transmission in Belgium, conceived by members of SIMID group during the 
# COVID19 pandemic.
#
# This file contains function for unit-testing to check the model logic.
#
# Copyright 2021, SIMID, University of Antwerp & Hasselt University                                        
########################################################################### #

suppressPackageStartupMessages(library('useful')) # compare.list

# set file names for reference values
prev_file_names <- c(mean       = 'unit_testing/fit_mean.rds',
                     stochastic = 'unit_testing/fit_stochastic.rds',
                     scm_out     = 'unit_testing/scm_out.rds',
                     scm_mean    = 'unit_testing/scm_mean.rds',
                     scm_prior   = 'unit_testing/scm_prior.rds',
                     extended   = 'unit_testing/fit_extended.rds')

# compare model results with stored reference values
# fit_new <- scm_out; method <- 'scm_out'
compare_output <- function(fit_new,method,prev_file_names){

  # load previous results  
  fit_prev       <- readRDS(file=prev_file_names[method])
  
  # remove column names
  for(i_elem in 1:length(fit_prev)){
    colnames(fit_prev[[i_elem]]) <- NULL
  }
  for(i_elem in 1:length(fit_new)){
    colnames(fit_new[[i_elem]]) <- NULL
  }
  
  # compare length and names
  if(length(fit_new) != length(fit_prev) || 
     names(fit_new) != names(fit_prev)){
    warning(c('!! Model output has different length or contains different names for method = ',method))
    
    warning(paste('!! NEW: ',names(fit_new)[!names(fit_new) %in% names(fit_prev)]))
    
    if(all(names(fit_prev) %in% names(fit_new))){
      fit_new <- fit_new[names(fit_prev)]
    } else{
      
    }
  }
  
  # compare length and names
  if(length(fit_new) != length(fit_prev) || 
     names(fit_new) != names(fit_prev)){
    warning(c('!! Model output has different length or contains different names for method = ',method))
  } else{
    
    # compare lists
    bool_compare <- compare.list(fit_new,fit_prev)
    
    # check for floating point issues
    digits_cutoff   <- 6
    fit_new_fp      <- lapply(fit_new,round,digits=digits_cutoff)
    fit_prev_fp     <- lapply(fit_prev,round,digits=digits_cutoff)
    bool_compare_fp <- compare.list(fit_new_fp,fit_prev_fp)

    # report outcome of comparison
    if(all(bool_compare)){
      cat(paste0("Model output for method = '",method,"' did not change."))
    } else{
      if(all(bool_compare_fp)){
        cat(paste0("Model output for method = '",method,"' did not differ more than 1e-",digits_cutoff))
      } else{
        warning(c('Model output has substantially changed!!',
                  paste(c('\n- different results:', names(fit_new)[!bool_compare_fp]),collapse=' '))
        )
      }
    }
    
  }
  
  
  # safety check of the high-level summary statistics
  fit_prev_stats <- get_reference_stats(fit_prev)
  fit_ref_stats  <- read.table(file=gsub('.rds','_summary.csv',prev_file_names[method]),sep=',',header = F)
  if(length(fit_prev_stats) == length(fit_prev_stats) && get_reference_stats(fit_prev) == fit_prev_stats){
      cat(paste0("  [Model reference OK]\n"))
    
  } else{
    warning(c("Model reference for method = '",method,"' has changed."))
  }
}


get_reference_stats <- function(fit_xxx){
  
  fit_stats  <- c(output_names = names(fit_xxx),
                  out_dim = unlist(lapply(fit_xxx,dim)),
                  out_sum = unlist(lapply(fit_xxx,sum,na.rm=T)))
  head(fit_stats)
  
  return(fit_stats)
}

# reset reference values (optional)
rrv <- function(){
  
  # store (high-level) reference statistics
  write.table(get_reference_stats(fit_mean_scm),file=gsub('.rds','_summary.csv',prev_file_names['mean']),sep=',',col.names = F)
  write.table(get_reference_stats(fit_stochastic_scm),file=gsub('.rds','_summary.csv',prev_file_names['stochastic']),sep=',',col.names = F)
  write.table(get_reference_stats(scm_out),file=gsub('.rds','_summary.csv',prev_file_names['scm_out']),sep=',',col.names = F)
  write.table(get_reference_stats(scm_mean),file=gsub('.rds','_summary.csv',prev_file_names['scm_mean']),sep=',',col.names = F)
  write.table(get_reference_stats(scm_prior),file=gsub('.rds','_summary.csv',prev_file_names['scm_prior']),sep=',',col.names = F)
  write.table(get_reference_stats(fit_extended_scm),file=gsub('.rds','_summary.csv',prev_file_names['extended']),sep=',',col.names = F)
  
  saveRDS(fit_mean_scm,file=prev_file_names['mean'])
  saveRDS(fit_stochastic_scm,file=prev_file_names['stochastic'])
  saveRDS(scm_out,file=prev_file_names['scm_out'])
  saveRDS(scm_mean,file=prev_file_names['scm_mean'])
  saveRDS(scm_prior,file=prev_file_names['scm_prior'])
  saveRDS(fit_extended_scm,file=prev_file_names['extended'])
  
  warning('!! Reference values for UNIT TESTING have been reset !!')
}



