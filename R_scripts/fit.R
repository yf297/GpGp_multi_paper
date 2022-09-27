library("GpGpm")
source("../R/fisher_multi.R")   
source("../R/fit_multi.R")
source("../R/link_multi.R")
source("../R/check.R")


# read the arguments
args <- commandArgs(trailingOnly = TRUE)
settings_file <- args[1]
settings_row <- as.numeric( args[2] )

# load the settings
settings <- read.csv(paste0("../settings/",settings_file,".csv"))[ settings_row, ]
fname <- paste0("../fits/fit_", settings_file, "_", settings_row, ".RData")

# define settings variables
order_fun <- get( settings$order_fun )
neighbor_fun <- get( settings$neighbor_fun )
m <- settings$m
corrnug <- settings$corrnug

# load the data
data <- get(load(paste0("../data/formatted/", settings$data_code, ".RData")))
y <- as.vector(data[,1])
locs <- as.matrix(data[,2:ncol(data)])
locs[,ncol(locs)] <- as.numeric( as.factor( locs[,ncol(locs)] ) )
X <- model.matrix(lm( y ~ -1 + as.factor(locs[,ncol(locs)])))   

# some info
ncomp <- length(unique(locs[,ncol(locs)]))
neach <- ncomp*(ncomp+1)/2
d <- ncol(locs) - 1
M <- matrix(0.5, ncomp, ncomp)
diag(M) <- 1

# start marginal parms and logparms
start_parms <- get_start_parms(y, X , locs, "matern_multi")$start_parms
start_logparms <- log(start_parms)
start_logparms <- append(start_logparms, 0, 2*neach)
start_logparms <- append(start_logparms, 0, 3*neach+1)

# some covparms indices
cross_nug_inds <- c()
for(j1 in 2:ncomp){ for(j2 in 1:(j1-1)){
    cross_nug_inds <- c( cross_nug_inds, multi_matern_parm_index(ncomp, j1, j2)$nugget )
}}
inds <- matrix(FALSE, ncomp, ncomp)    
diag(inds) <- TRUE 
marginal_ran_inds <- neach + which(t(inds)[upper.tri(inds, diag = TRUE)] == TRUE)  

# some flex-a logparms indices
inds <- matrix(FALSE, ncomp, ncomp)    
inds[upper.tri(inds, diag = FALSE)] <- TRUE
inds <- t(inds)                                      
log_cross_var_inds <- which(t(inds)[upper.tri(inds, diag = TRUE)] == TRUE)   
log_cross_ran_inds <-    neach +     log_cross_var_inds
log_delta_B_ind    <-  2*neach + 1
log_cross_smo_inds <-  2*neach + 1 + log_cross_var_inds
log_delta_A_ind    <-  3*neach + 2
log_cross_nug_inds <-  3*neach + 2 + log_cross_var_inds    
log_nug_inds <- (3*neach + 3): length(start_logparms)
inds <- matrix(FALSE, ncomp, ncomp)    
diag(inds) <- TRUE
log_marginal_ran_inds <- neach + which(t(inds)[upper.tri(inds, diag = TRUE)] == TRUE)  
log_marginal_smo_inds <- 2*neach + 1 + which(t(inds)[upper.tri(inds, diag = TRUE)] == TRUE)     

# some flex-b logparms indices
log_beta_ind <-  4*neach + 3

# some pars logparms indices
pars_log_cross_nug_inds <- (neach + ncomp + 1) + log_cross_var_inds 

# fit sequence of models

# Independent

# using flex_a link, but does not matter which link since fitting independent
linkfuns <- list()
linkfuns$link <- link_flex_a
linkfuns$dlink <- dlink_flex_a

start_logparms[log_cross_var_inds] <- 0
start_logparms[log_cross_nug_inds] <- 0 

active_logparms <- rep(TRUE, length(start_logparms))
active_logparms[c(log_cross_var_inds,
	 log_cross_ran_inds,
	 log_delta_B_ind,
	 log_cross_smo_inds,
	 log_delta_A_ind,
	 log_cross_nug_inds)] <- FALSE

print("starting independent fit")
t1 <- proc.time()
fit1 <- fit_multi(
    y, locs, X,
    neighbor_fun = neighbor_fun, 
    order_fun = order_fun,
    m = m,
    start_logparms = start_logparms, 
    linkfuns = linkfuns,
    active_logparms = active_logparms
)
t2 <- proc.time()
fit1$time_elapsed <- (t2-t1)[3]
fit1$valid <- TRUE

fit <- list( fit1 )
save(fit, file = fname)

# parsimonious

# pars link
linkfuns$link <- link_pars
linkfuns$dlink <- dlink_pars

start_logparms <- fit1$logparms
a <- log(mean(fit1$covparms[marginal_ran_inds]))
start_logparms <- c(start_logparms[1:(neach)],
		    a,
		   start_logparms[log_marginal_smo_inds],
		   start_logparms[log_nug_inds])

active_logparms <- rep(TRUE, length(start_logparms))
if( !corrnug ){ active_logparms[ pars_log_cross_nug_inds ] <- FALSE }

print("starting parsimonious fit") 
t1 <- proc.time()
fit2 <- fit_multi(
    y, locs, X,
    neighbor_fun = neighbor_fun, 
    order_fun = order_fun,
    m = m,
    start_logparms = start_logparms, 
    linkfuns = linkfuns,
    active_logparms = active_logparms
)
t2 <- proc.time()
fit2$time_elapsed <- (t2-t1)[3]   
fit2$valid <- T

fit <- list( fit1, fit2 )
save(fit, file = fname)


# flex-a

# flex-a link
linkfuns$link <- link_flex_a
linkfuns$dlink <- dlink_flex_a

start_logparms <- fit1$logparms
start_logparms[log_cross_ran_inds]  <- log(finv(M)) 
start_logparms[log_delta_B_ind] <- log(0.01)
start_logparms[log_cross_smo_inds]  <-  log(finv(M))  
start_logparms[log_delta_A_ind] <- log(0.01)
if(settings$data_code == "weather"){
start_logparms[log_delta_A_ind] <- log(0.001)}


active_logparms <- rep(TRUE, length(start_logparms))
if( !corrnug ){ active_logparms[ log_cross_nug_inds ] <- FALSE }
if(ncomp==2){
	active_logparms[log_cross_ran_inds] <- FALSE
	active_logparms[log_cross_smo_inds] <- FALSE  
}

print("starting flexible-a fit")
t1 <- proc.time()
fit3 <- fit_multi(
    y, locs, X,
    neighbor_fun = neighbor_fun, 
    order_fun = order_fun,
    m = m,
    start_logparms = start_logparms, 
    linkfuns = linkfuns,
    active_logparms = active_logparms
)
t2 <- proc.time()   
fit3$time_elapsed <- (t2-t1)[3]       
fit3$valid <- T

fit <- list( fit1, fit2, fit3)
save(fit, file = fname)


# flex-e

# flex-e link
linkfuns$link <- link_flex_e
linkfuns$dlink <- dlink_flex_e

start_logparms <- fit1$logparms
start_logparms <- append(start_logparms, 0, length(start_logparms))        
start_logparms[log_cross_ran_inds]  <- log(finv(M)) 
start_logparms[log_delta_B_ind] <- log(0.01)
start_logparms[log_cross_smo_inds]  <-  log(finv(M))  
start_logparms[log_delta_A_ind] <- log(0.01) 
start_logparms[log_beta_ind] <- log(1)


active_logparms <- rep(TRUE, length(start_logparms))
if( !corrnug ){ active_logparms[ log_cross_nug_inds ] <- FALSE }
if(ncomp==2){
	active_logparms[log_cross_ran_inds] <- FALSE
	active_logparms[log_cross_smo_inds] <- FALSE  
}

print("starting flexible-e fit")
t1 <- proc.time()
fit4 <- fit_multi(
    y, locs, X,
    neighbor_fun = neighbor_fun, 
    order_fun = order_fun,
    m = m,
    start_logparms = start_logparms, 
    linkfuns = linkfuns,
    active_logparms = active_logparms
)
t2 <- proc.time()   
fit4$time_elapsed <- (t2-t1)[3]       
fit4$valid <- T

fit <- list( fit1, fit2, fit3, fit4)
save(fit, file = fname)

# unconstrained

# unconstrained link
linkfuns$link <- link_unc03
linkfuns$dlink <- d_link_unc03
linkfuns$invlink <- invlink_unc03

start_logparms <- linkfuns$invlink(fit4$covparms)     

active_logparms <- rep(TRUE, length(start_logparms))
if( !corrnug ){ active_logparms[ cross_nug_inds ] <- FALSE }

print("starting unconstrained fit")
t1 <- proc.time()         
fit5 <- fit_multi(
    y, locs, X,
    neighbor_fun = neighbor_fun, 
    order_fun = order_fun,
    m = m,
    start_logparms = start_logparms, 
    linkfuns = linkfuns,
    active_logparms = active_logparms
)
t2 <- proc.time()   
fit5$time_elapsed <- (t2-t1)[3]       

# check condition 2 in apanasovich
c2 <- check_flex_a(fit5$covparms,d)

# check condition 1 in emery
c1 <- check_flex_b(fit5$covparms)

fit5$valid <- c1&&c2

fit <- list( fit1, fit2, fit3, fit4, fit5 )
save(fit, file = fname)


