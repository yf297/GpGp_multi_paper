# GpGp_multi_paper

The source R package GpGpm_0.4.0.tar.gz should be in the same directory as
GpGp_multi_paper. Before starting navigate to R_scripts directory of this
repository, and run

```
Rscript install.packages.R
```

The models are fit using the R script R_scripts/fit.R.
This file takes two arguments:
  1. code for a settings file
  2. which row of the settings file to run.

The settings folder contains the settings csv files 
  * vecchia.csv 
  * weather.csv 
  * model.csv 
  * many_datasets.csv  
  * corrnug.csv 

Each csv file contains several rows with different configurations.
The code for each settings file is simply the filename without
the .csv extension. 

So for example, to fit the model corresponding to the configuration
in row 7 of many_datasets.csv, navigate to the R_scripts directory
and run. 

```
Rscript fit.R many_datasets 7
```

This will create fit_many_datasets_7.RData in GpGp_multi_paper/fits.

As another example, navigate to the R_scripts directory and type

```
Rscript fit.R model 1
Rscript fit.R model 2
Rscript fit.R model 3
```

to run the 3 configurations in the model csv file.

This will create fit_model_1.RData fit_model_2.RData and fit_model_3.RData in GpGp_multi_paper/fits. 

Tables from the paper can be generated using the table_<settings_code>.R
scripts. For example, to generate the corrnug table, navigate to the
R_scripts directory and type

```
Rscript table_corrnug.R
```

Note that the fit files need to be generated before the tables can be produced.

All of the fits can be run by running the fit_all bash script.
