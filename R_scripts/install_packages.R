
#install.packages("devtools", repos="https://cloud.r-project.org")
#install.packages("rootSolve", repos="https://cloud.r-project.org")
#install.packages("fields", repos="https://cloud.r-project.org")
install.packages("matrixcalc", repos="https://cloud.r-project.org")
devtools::build("../../GpGpm")
install.packages("../../GpGpm_0.4.0.tar.gz", type = "source", repos = NULL )
