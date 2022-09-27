
predictions_multi <- function(fit = NULL, locs_pred, X_pred, 
    y_obs = fit$y, locs_obs = fit$locs, X_obs = fit$X, beta = fit$betahat,    
    covparms = fit$covparms, covfun_name = fit$covfun_name, 
    m = 60, order_fun = NULL, neighbor_fun = NULL, st_scale = NULL){
    
    n_obs <- nrow(locs_obs)
    n_pred <- nrow(locs_pred)
    
    # get orderings
    ord1 <- order_fun( locs_obs )
    ord2 <- sample( 1:n_pred )

    # reorder stuff
    X_obs <- as.matrix(X_obs)
    X_pred <- as.matrix(X_pred)
    Xord_obs  <- X_obs[ord1,,drop=FALSE]
    Xord_pred <- X_pred[ord2,,drop=FALSE]
    yord_obs  <- y_obs[ord1]
    
    # put all coordinates together
    locs_all <- rbind( locs_obs[ord1,,drop=FALSE], locs_pred[ord2,,drop=FALSE] )
    inds1 <- 1:n_obs
    inds2 <- (n_obs+1):(n_obs+n_pred)
    
    NNarray_all <- neighbor_fun(locs_all, m=m)
    
    # get entries of Linv for obs locations and pred locations
    Linv_all <- vecchia_Linv(covparms,covfun_name,locs_all, NNarray_all, n_obs+1)
    
    y_withzeros <- c(yord_obs - c(Xord_obs %*% beta), rep(0,n_pred) )
    v1 <- Linv_mult(Linv_all, y_withzeros, NNarray_all )
    v1[inds1] <- 0
    Linv_all[1:n_obs,1] <- 1.0
    v2 <- -L_mult(Linv_all,v1,NNarray_all)

    condexp <- c(v2[inds2] + Xord_pred %*% beta)
    condexp[ord2] <- condexp
    return(condexp)
}





cond_sim_multi <- function(fit = NULL, locs_pred, X_pred, 
    y_obs = fit$y, locs_obs = fit$locs, X_obs = fit$X, beta = fit$betahat,    
    covparms = fit$covparms, covfun_name = fit$covfun_name, 
    m = 60, order_fun = NULL, neighbor_fun = NULL, st_scale = NULL, nsims = 1 ){

    n_obs <- nrow(locs_obs)
    n_pred <- nrow(locs_pred)
    
    # get orderings
    ord1 <- order_fun(locs_obs)
    ord2 <- sample( 1:n_pred ) 

    # reorder stuff
    X_obs <- as.matrix(X_obs)
    X_pred <- as.matrix(X_pred)
    Xord_obs  <- X_obs[ord1,,drop=FALSE]
    Xord_pred <- X_pred[ord2,,drop=FALSE]
    yord_obs  <- y_obs[ord1]

    # put all coordinates together
    locs_all <- rbind( locs_obs[ord1,], locs_pred[ord2,] )
    inds1 <- 1:n_obs
    inds2 <- (n_obs+1):(n_obs+n_pred)
    
    # get nearest neighbor array (in space only)
    NNarray_all <- neighbor_fun(locs_all,m=m )

    # get entries of Linv for obs locations and pred locations
    Linv_all <- vecchia_Linv(covparms,covfun_name,locs_all,NNarray_all)
    
    # an unconditional simulation
    condsim <- matrix(NA, n_pred, nsims)
    for(j in 1:nsims){
        z <- L_mult(Linv_all, stats::rnorm(n_obs+n_pred), NNarray_all)
    
        y_withzeros <- c(yord_obs - Xord_obs %*% beta + z[inds1], rep(0,n_pred) )
        v1 <- Linv_mult(Linv_all, y_withzeros, NNarray_all )
        v1[inds1] <- 0
        v2 <- -L_mult(Linv_all,v1,NNarray_all)

        condsim[ord2,j] <- c(v2[inds2] + Xord_pred %*% beta) - z[inds2]
    }
    return(condsim)
}
