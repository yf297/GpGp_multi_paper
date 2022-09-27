
library("xtable")

settings_code <- "vecchia"
settings <- read.csv("../settings/vecchia.csv")

# Make a table of the final parameter settings
results <- matrix(NA, nrow(settings), 15)

for(j in 1:nrow(settings)){
    try({
        fname <- paste0("../fits/fit_", settings_code, "_", j, ".RData" )
        load(fname)
        for(k in 1:5){
            results[j,3*(k-1) + 1] <- fit[[k]]$loglik
            results[j,3*(k-1) + 2] <- fit[[k]]$iter
            results[j,3*(k-1) + 3] <- fit[[k]]$time_elapsed/60
        }
    })
}
results[,c(1,4,7,10,13)] <- results[,c(1,4,7,10,13)] - max(results[,c(1,4,7,10,13)])

results <- as.data.frame(results)

colnames(results) <- rep(c("loglik","iter","time"),5)

settings$ord <- NA
settings$ord[ settings$order_fun == "order_completely_random" ] <- "rand-"
settings$ord[ settings$order_fun == "order_bycomponent_random" ] <- "comp-"
settings$ord[ settings$order_fun == "order_cycle_component_random" ] <- "cyc-"
settings$nei <- NA
settings$nei[ settings$neighbor_fun == "nearest_multi_any" ] <- "any"
settings$nei[ settings$neighbor_fun == "nearest_multi_balanced" ] <- "bal"
settings$nei[ settings$neighbor_fun == "nearest_multi_pref_this" ] <- "pref"

results$Setting <- paste0( settings$ord, settings$nei )
results$m <- settings$m

digs <- c( rep( c(2,0,1), 5 ), 0, 0 )


row_ord <- 1:18
col_ord <- c(16, 17, 1,2, 4,5, 7,8, 10,11, 13,14)
print.xtable(
    xtable(results[row_ord,col_ord], digits = c(0,digs[col_ord])),
    sanitize.text.function = function(x) x,
    include.rownames = FALSE
)
