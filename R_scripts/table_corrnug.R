
library("xtable")

settings_code <- "corrnug"
settings <- read.csv("../settings/corrnug.csv")

# Make a table of the final parameter settings
parms <- matrix(NA, nrow(settings), 13)

for(j in 1:nrow(settings)){
    fname <- paste0("../fits/fit_", settings_code, "_", j, ".RData" )
    load(fname)
    #if( !fit[[5]]$no_decrease ){
        parms[j,1:12] <- fit[[4]]$covparms
        parms[j,13] <- fit[[4]]$loglik
   # }
}

for( ii in list( c(1,2), c(3,4), c(5,6), c(7,8), c(9,10), c(11,12) ) ){
    parms[ii,13] <- parms[ii,13] - max(parms[ii,13])
}

colnames(parms) <- c(
    paste0(
        rep(c("$\\sigma_","$\\alpha_","$\\nu_","$\\tau_"), each = 3 ),
        rep(c("{11}$","{12}$","{22}$"), 4 )
    ),
    "loglik"
)

parms <- as.data.frame(parms)
ids <- settings$data_code
ids <- gsub("bang_", "", ids )
ids <- gsub("_", ",", ids )
parms$Comp <- ids
parms$Comp[ seq(2,12,by=2) ] <- ""

row_ord <- 1:12
col_ord <- c(14,13,1:12)
print.xtable(
    xtable(parms[row_ord,col_ord]),
    sanitize.text.function = function(x) x,
    include.rownames = FALSE
)
