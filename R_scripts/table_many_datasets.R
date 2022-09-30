
library("xtable")

settings_code <- "many_datasets"
settings <- read.csv("../settings/many_datasets.csv")

# Make a table of the final parameter settings
results <- matrix(NA, nrow(settings), 15)

for(j in 1:nrow(settings)){
    try({
        fname <- paste0("../fits/fit_", settings_code, "_", j, ".RData" )
        load(fname)
        for(k in 1:5){
            results[j,3*(k-1) + 1] <- fit[[k]]$loglik - fit[[1]]$loglik
            results[j,3*(k-1) + 2] <- fit[[k]]$iter
            results[j,3*(k-1) + 3] <- fit[[k]]$time_elapsed/60
        }
    })
}

colnames(results) <- rep(c("loglik","iter","time"),5)

results <- as.data.frame(results)
ids <- settings$data_code
ids <- gsub("bang_", "", ids )
ids <- gsub("goes01_", "", ids )
ids <- gsub("_", ",", ids )
results$Comp <- ids

digs <- c( rep( c(2,0,1), 5 ), 1 )


row_ord <- 1:22
col_ord <- c(16,4:15)
print.xtable(
    xtable(results[row_ord,col_ord], digits = c(0,digs[col_ord])),
    sanitize.text.function = function(x) x,
    include.rownames = FALSE
)
