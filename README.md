# GpGp_multi_paper

GpGpm should be in the same directory as GpGp_multi_paper. Before
starting navigate to R_scripts directory of this repository,
and run

```
Rscript install.packages.R
```

The settings folder contains the csv files 
  * vecchia.csv 
  * weather.csv 
  * model.csv 
  * many_datasets.csv  
  * corrnug.csv 

Each csv file contains several rows with different configurations.

The weather.csv is a simple example with only 1 row.

To run this setting, navigate to the R_scripts directory and type

Rscript fit.R weather 1

This will create fit_weather_1.RData in GpGp_multi_paper/fits.

To generate a table, navigate to the R_scripts directory and type

Rscript table_weather.R

As another example, navigate to the R_scripts directory and type

Rscript fit.R model 1

Rscript fit.R model 2

Rscript fit.R model 3

to run the 3 configurations in the model csv file.

This will create fit_model_1.RData fit_model_2.RData and fit_model_3.RData in GpGp_multi_paper/fits. 

To generate a table, navigate to the R_scripts directory and type

Rscript table_model.R

Note that fit_model_1.RData fit_model_2.RData and fit_model_3.RData have to be in GpGp_multi_paper/fits before table_model.R can work.

