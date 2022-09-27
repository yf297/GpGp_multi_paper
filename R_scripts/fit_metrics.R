library(GpGpm)

args <- commandArgs(trailingOnly = TRUE)
data_code <- args[1]

# load the settings
settings <- read.csv("../settings/settings.csv")

# add enough room for 8 fits
time_elapsed <- array(NA, c(nrow(settings), 8) )
iter <- array(NA, c(nrow(settings), 8) ) 
loglik <- array(NA, c(nrow(settings), 8) )
valid <- array(NA, c(nrow(settings), 8) )

# get y,locs,X from first fit of first setting to compute NNarray with 120 neighbors (does not matter which fit or setting)
a_fit <- get(load( paste0("../fits/fit_", data_code, "_", 1, ".RData")))
y <- a_fit[[1]]$y
locs <- a_fit[[1]]$locs
X <- a_fit[[1]]$X

order_fun <- get("order_completely_random")
neighbor_fun <- get("nearest_multi_any")

locsord <- order_completely_random(locs)
ord <- order_fun(locs)

yord <- y[ord]
locsord <- locs[ord,,drop=FALSE]
Xord <- as.matrix( X[ord,,drop=FALSE] )

m <- 120
if(nrow(locs) < 400){m <- nrow(locs)}
NNarray <- nearest_multi_any(locsord, m)


for(j in 1:nrow(settings)){

    # load in the predictions
    load( paste0("../fits/fit_", data_code, "_", j, ".RData") )
    
    for(k in 1:length(fit)){
        # calculate metrics
        time_elapsed[j,k] <- fit[[k]]$time_elapsed
	iter[j,k] <- fit[[k]]$iter 
        loglik[j,k] <- vecchia_profbeta_loglik(fit[[k]]$covparms, "matern_multi", yord, Xord, locsord, NNarray)$loglik
	valid[j,k] <- fit[[k]]$valid
    }
}

