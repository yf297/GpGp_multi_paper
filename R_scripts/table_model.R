
library("xtable")

settings_code <- "model"
settings <- read.csv("../settings/model.csv")

# Make a table of the final parameter settings
parms <- matrix(NA, nrow(settings)*5, 13)

cnt <- 0
for(j in 1:nrow(settings)){
    fname <- paste0("../fits/fit_", settings_code, "_", j, ".RData" )
    load(fname)
    for(k in 1:5){
        try({
        cnt <- cnt + 1
        parms[cnt,1:12] <- fit[[k]]$covparms
        parms[cnt,13] <- fit[[k]]$loglik
	print(fit[[k]]$time_elapsed)
        })
    }
}

colnames(parms) <- c(
    paste0(
        rep(c("$\\sigma_","$\\alpha_","$\\nu_","$\\tau_"), each = 3 ),
        rep(c("{11}$","{12}$","{22}$"), 4 )
    ),
    "loglik"
)

parms <- as.data.frame(parms)

row_ord <- 1:15
col_ord <- c(1:13)
print.xtable(
    xtable(parms[row_ord,col_ord]),
    sanitize.text.function = function(x) x,
    include.rownames = FALSE
)
