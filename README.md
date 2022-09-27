# GpGp_multi_paper

The settings folder contains the csv files vecchia.csv weather.csv model.csv many_datasets.csv and corrnug.csv 

Each csv file contains several rows with different configurations.

The weather.csv is a simple example with only 1 row.

To run this setting, first an create a R_scripts/fit.R directory. Navigate to the R_scripts directory and type

Rscript fit.R weather 1

This will create fit_weather_1.RData in R_scripts/fit.R
R_scripts/fit.R

To generate a table,  Navigate to the R_scripts directory and type

Rscipt table_weather.R
