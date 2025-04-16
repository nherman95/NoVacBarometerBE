#!/bin/bash

cd ..

value=1
numscen=5
while [ $value -le $numscen ]
do
    echo "$value"_1
    Rscript R_processing/vaccine_uptake_ECDC_R5.R "$value"_0
    Rscript R_processing/vaccine_uptake_ECDC_R5.R "$value"_1
    ((value++))
done