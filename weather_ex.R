library("GpGpm")
source("R/fisher_multi.R")   
source("R/fit_multi.R")
source("R/link_multi.R")
source("R/check.R")



######################## Load data

data <- get(load(paste0("data/formatted/Weather.RData")))

######################## Format data

y <- as.vector(data[,1])
locs <- as.matrix(data[,2:ncol(data)])
locs[,ncol(locs)] <- as.numeric( as.factor( locs[,ncol(locs)] ) )
X <- model.matrix(lm( y ~ -1 + as.factor(locs[,ncol(locs)])))   

######################## Define order and neighbor functions, and number of neighbors. Just use these. Other options in the settings folder.
order_fun <- get("order_completely_random")
neighbor_fun <- get("nearest_multi_any")
m = 40
corrnug = FALSE

######################## Define some parameters. Leave this as is.

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


########################  Fit a sequence of models, starting with independent, then parsimonous, and so on
########################. This is to stabilize the optimization. We use the parameters from one model as the starting parameters for the next.
########################. Below is is the fit for the independent and parsimonious models. For others, see fit.R file
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
fit1 <- fit_multi(
  y, locs, X,
  neighbor_fun = neighbor_fun, 
  order_fun = order_fun,
  m = m,
  start_logparms = start_logparms, 
  linkfuns = linkfuns,
  active_logparms = active_logparms
)

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



