#!/bin/bash

# Example of couples (or_prop, red_prop)
couple_values=(
  "0.6 0.2"
  "0.9 0.4"
  "0.9 0.2"
  "0.6 0.4"
  "0.3 0.1"
)

# Example of quadruplets (new_hosp_or, new_hosp_red, icu_or, icu_red)
quadruplet_values=(
  "65 150 300 500"
  "150 300 400 800"
)

# R script to execute with the adjusted path
script="./R/projections_BE.R"
output_tag="test"
rm -r output/$output_tag
# Flag to enable or disable parallel execution
enable_parallel=true

# Number of processors to use
num_processors=6

# Function to create temporary config file
create_config_file() {
  config_file="config.txt"
  rm -f $config_file
  for couple in "${couple_values[@]}"; do
    for quadruplet in "${quadruplet_values[@]}"; do
      echo "Rscript $script FALSE TRUE $output_tag $couple $quadruplet " >> $config_file
    done
  done
  echo "Rscript $script FALSE FALSE $output_tag" >> $config_file
}

# Function to execute Rscript
execute_rscript() {
  if [ "$enable_parallel" = true ]; then
    parallel --jobs $num_processors --no-notice < config.txt
  else
    while IFS= read -r line; do
      eval "$line"
    done < config.txt
  fi
}


# Execute Rscript based on the flag
Rscript $script TRUE FALSE $output_tag
create_config_file
execute_rscript

# Gather all lines into a CSV file
rm config.txt
