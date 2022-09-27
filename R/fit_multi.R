#devtools::load_all("../../GpGpm")

fit_multi <- function(y, locs, X, covfun_name = "matern_multi",
    neighbor_fun, order_fun, m, max_iter = 40,
    silent = FALSE, st_scale = NULL, convtol = 1e-4, start_logparms, linkfuns, active_logparms){

    n <- length(y)

    # check that length of observation vector same as
    # number of locations
    if( nrow(locs) != n ){
        stop("length of observation vector y not equal
              to the number of locations (rows in locs)")
    }
    locs <- as.matrix(locs)

    # detect and remove missing values
    not_missing <- apply( cbind(y,locs,X), 1,
        function(x){
            if( sum(is.na(x) | is.infinite(x)) > 0 ){
                return(FALSE)
            } else { return(TRUE) }
        }
    )
    if( sum(not_missing) < n ){
        y <- y[not_missing]
        locs <- locs[not_missing,,drop=FALSE]
        X <- X[not_missing,,drop=FALSE]
        cat(paste0( n - sum(not_missing),
            " observations removed due to missingness or Inf\n"))
    }

    # redefine n
    n <- length(y)
   

    link <- linkfuns$link
    dlink <- linkfuns$dlink
    invlink <- linkfuns$invlink
    invlink_startparms <- start_logparms
	
    active <- active_logparms

    penalty <- get_penalty(y,X,locs,covfun_name) 
    pen <- penalty$pen
    dpen <- penalty$dpen
    ddpen <- penalty$ddpen

    # get an ordering and reorder everything 
    ord <- order_fun(locs)
    yord <- y[ord]
    locsord <- locs[ord,,drop=FALSE]
    Xord <- as.matrix( X[ord,,drop=FALSE] )

    # get neighbor array if not provided
    
    NNarray <- neighbor_fun(locsord, m) 
        
    likfun <- function(logparms){

	lp <- rep(NA,length(invlink_startparms))
        lp[active] <- logparms
        lp[!active] <- invlink_startparms[!active]
                
        likobj <- vecchia_profbeta_loglik_grad_info(link(lp),covfun_name,yord,Xord,locsord,NNarray)

	likobj$loglik <- -likobj$loglik - pen(link(lp))
        likobj$grad <-  as.vector(-likobj$grad%*%dlink(lp) - dpen(link(lp))%*%dlink(lp)) 
        likobj$info <- t(dlink(lp))%*%likobj$info%*%dlink(lp) - t(dlink(lp))%*%ddpen(link(lp))%*%dlink(lp)
        likobj$grad <- likobj$grad[active]
        likobj$info <- likobj$info[active,active]



        return(likobj)    
    }

    fit <- fisher_scoring_multi( likfun,invlink_startparms[active],
            link,silent=silent, convtol = convtol, max_iter = max_iter , invlink_startparms, active)
    invlink_startparms <- fit$logparms
    #covparms <- fit$covparms
    #iter <- fit$iter
    loglik <- -fit$loglik - pen(fit$covparms)

    #fit <-list()
    fit$loglik <- loglik 
    fit$covfun_name <- covfun_name
    fit$logparms <- invlink_startparms 
    #fit$covparms <- covparms
    fit$y <- y
    fit$locs <- locs
    fit$X <- X
    #fit$iter <- iter
    class(fit) <- "GpGp_fit"
    return(fit)
}
 

# nonstandard link functions
exp2 <- function(x){ exp(x) - exp(-x) }
d_exp2 <- function(x){ exp(x) + exp(-x) }
inv_exp2 <- function(x){
    y <- rep(NA, length(x))
    y[ x >= 0] <-  log( ( x[ x>=0] + sqrt(x[ x>=0]^2+4))/2 )
    y[ x < 0 ] <- -log( (-x[ x<0 ] + sqrt(x[ x<0 ]^2+4))/2 )
    return(y)
}
tan2 <- function(x){ 2/pi*atan(x) }
d_tan2 <- function(x){ 2/pi/(1+x^2) }
inv_tan2 <- function(x){ tan(pi/2*x) }

link_unc01 <- function(x){
    # figure out indices of cross covariances
    ncomp <- ncomp_from_nparms( length(x) )
    i2 <- c()
    for(j1 in 2:ncomp){
        for(j2 in 1:(j1-1)){
            i2 <- c(i2, multi_matern_parm_index(ncomp, j1, j2)$variance )
            i2 <- c(i2, multi_matern_parm_index(ncomp, j1, j2)$nugget )
        }
    }
    i1 <- (1:length(x))[-i2]

    # calculate link
    y <- rep(NA,length(x))
    y[i1] <- exp(x[i1])
    y[i2] <- exp2(x[i2])
    return(y)
}
 
link_unc02 <- function(x){
    # figure out indices of cross covariances
    ncomp <- ncomp_from_nparms( length(x) )
    y <- rep(NA,length(x))
    for(j1 in 1:ncomp){
        for(j2 in 1:j1){
            ii <- multi_matern_parm_index(ncomp, j1, j2 )
            y[ii$range] <- exp(x[ii$range])
            y[ii$smoothness] <- exp(x[ii$smoothness])
            if(j1 == j2){
                y[ii$variance] <- exp(x[ii$variance])
                y[ii$nugget] <- exp(x[ii$nugget])
            } else {
                ii1 <- multi_matern_parm_index(ncomp, j1, j1)
                ii2 <- multi_matern_parm_index(ncomp, j2, j2)
                y[ii$variance] <- exp(0.5*(x[ii1$variance]+x[ii2$variance]))*exp2(x[ii$variance])
                y[ii$nugget] <- exp(0.5*(x[ii1$nugget]+x[ii2$nugget]))*exp2(x[ii$nugget])
            }
        }
    }
    return(y)
}
 
link_unc03 <- function(x){
    # figure out indices of cross covariances
    ncomp <- ncomp_from_nparms( length(x) )
    y <- rep(NA,length(x))
    for(j1 in 1:ncomp){
        for(j2 in 1:j1){
            ii <- multi_matern_parm_index(ncomp, j1, j2 )
            y[ii$range] <- exp(x[ii$range])
            y[ii$smoothness] <- exp(x[ii$smoothness])
            if(j1 == j2){
                y[ii$variance] <- exp(x[ii$variance])
                y[ii$nugget] <- exp(x[ii$nugget])
            } else {
                ii1 <- multi_matern_parm_index(ncomp, j1, j1)
                ii2 <- multi_matern_parm_index(ncomp, j2, j2)
                y[ii$variance] <- exp(0.5*(x[ii1$variance]+x[ii2$variance]))*tan2(x[ii$variance])
                y[ii$nugget] <- exp(0.5*(x[ii1$nugget]+x[ii2$nugget]))*tan2(x[ii$nugget])
            }
        }
    }
    return(y)
}
 
d_link_unc01 <- function(x){
    # figure out indices of cross covariances
    ncomp <- ncomp_from_nparms( length(x) )
    i2 <- c()
    for(j1 in 2:ncomp){
        for(j2 in 1:(j1-1)){
            i2 <- c(i2, multi_matern_parm_index(ncomp, j1, j2)$variance )
            i2 <- c(i2, multi_matern_parm_index(ncomp, j1, j2)$nugget )
        }
    }
    i1 <- (1:length(x))[-i2]

    # calculate d_link
    d <- diag(length(x))
    for(j in i1){ d[j,j] <- exp(x[j]) }
    for(j in i2){ d[j,j] <- d_exp2(x[j]) }
    return(d)
}
 
d_link_unc02 <- function(x){
    # figure out indices of cross covariances
    ncomp <- ncomp_from_nparms( length(x) )
    d <- diag(length(x))
    for(j1 in 1:ncomp){
        for(j2 in 1:j1){
            ii <- multi_matern_parm_index(ncomp, j1, j2 )
            d[ii$range,ii$range] <- exp(x[ii$range])
            d[ii$smoothness,ii$smoothness] <- exp(x[ii$smoothness])
            if(j1 == j2){
                d[ii$variance,ii$variance] <- exp(x[ii$variance])
                d[ii$nugget,ii$nugget] <- exp(x[ii$nugget])
            } else {
                ii1 <- multi_matern_parm_index(ncomp, j1, j1)
                ii2 <- multi_matern_parm_index(ncomp, j2, j2)
                d[ii$variance,ii$variance] <-
                    exp(0.5*(x[ii1$variance]+x[ii2$variance]))*d_exp2(x[ii$variance])
                d[ii$variance,ii1$variance] <- 
                    0.5*exp(0.5*(x[ii1$variance]+x[ii2$variance]))*exp2(x[ii$variance])
                d[ii$variance,ii2$variance] <- 
                    0.5*exp(0.5*(x[ii1$variance]+x[ii2$variance]))*exp2(x[ii$variance])
                d[ii$nugget,ii$nugget] <-
                    exp(0.5*(x[ii1$nugget]+x[ii2$nugget]))*d_exp2(x[ii$nugget])
                d[ii$nugget,ii1$nugget] <- 
                    0.5*exp(0.5*(x[ii1$nugget]+x[ii2$nugget]))*exp2(x[ii$nugget])
                d[ii$nugget,ii2$nugget] <- 
                    0.5*exp(0.5*(x[ii1$nugget]+x[ii2$nugget]))*exp2(x[ii$nugget])
            }
        }
    }
    return(d)
}
 
d_link_unc03 <- function(x){
    # figure out indices of cross covariances
    ncomp <- ncomp_from_nparms( length(x) )
    d <- diag(length(x))
    for(j1 in 1:ncomp){
        for(j2 in 1:j1){
            ii <- multi_matern_parm_index(ncomp, j1, j2 )
            d[ii$range,ii$range] <- exp(x[ii$range])
            d[ii$smoothness,ii$smoothness] <- exp(x[ii$smoothness])
            if(j1 == j2){
                d[ii$variance,ii$variance] <- exp(x[ii$variance])
                d[ii$nugget,ii$nugget] <- exp(x[ii$nugget])
            } else {
                ii1 <- multi_matern_parm_index(ncomp, j1, j1)
                ii2 <- multi_matern_parm_index(ncomp, j2, j2)
                d[ii$variance,ii$variance] <-
                    exp(0.5*(x[ii1$variance]+x[ii2$variance]))*d_tan2(x[ii$variance])
                d[ii$variance,ii1$variance] <- 
                    0.5*exp(0.5*(x[ii1$variance]+x[ii2$variance]))*tan2(x[ii$variance])
                d[ii$variance,ii2$variance] <- 
                    0.5*exp(0.5*(x[ii1$variance]+x[ii2$variance]))*tan2(x[ii$variance])
                d[ii$nugget,ii$nugget] <-
                    exp(0.5*(x[ii1$nugget]+x[ii2$nugget]))*d_tan2(x[ii$nugget])
                d[ii$nugget,ii1$nugget] <- 
                    0.5*exp(0.5*(x[ii1$nugget]+x[ii2$nugget]))*tan2(x[ii$nugget])
                d[ii$nugget,ii2$nugget] <- 
                    0.5*exp(0.5*(x[ii1$nugget]+x[ii2$nugget]))*tan2(x[ii$nugget])
            }
        }
    }
    return(d)
}
 
invlink_unc01 <- function(x){
    # figure out indices of cross covariances
    ncomp <- ncomp_from_nparms( length(x) )
    i2 <- c()
    for(j1 in 2:ncomp){
        for(j2 in 1:(j1-1)){
            i2 <- c(i2, multi_matern_parm_index(ncomp, j1, j2)$variance )
            i2 <- c(i2, multi_matern_parm_index(ncomp, j1, j2)$nugget )
        }
    }
    i1 <- (1:length(x))[-i2]

    # calculate link
    y <- rep(NA,length(x))
    y[i1] <- log(x[i1])
    y[i2] <- inv_exp2(x[i2])
    return(y)
}
 
invlink_unc02 <- function(x){
    # figure out indices of cross covariances
    ncomp <- ncomp_from_nparms( length(x) )
    y <- rep(NA,length(x))
    for(j1 in 1:ncomp){
        for(j2 in 1:j1){
            ii <- multi_matern_parm_index(ncomp, j1, j2 )
            y[ii$range] <- log(x[ii$range])
            y[ii$smoothness] <- log(x[ii$smoothness])
            if(j1 == j2){
                y[ii$variance] <- log(x[ii$variance])
                y[ii$nugget] <- log(x[ii$nugget])
            } else {
                ii1 <- multi_matern_parm_index(ncomp, j1, j1)
                ii2 <- multi_matern_parm_index(ncomp, j2, j2)
                y[ii$variance] <- inv_exp2(x[ii$variance]/sqrt(x[ii1$variance]*x[ii2$variance]))
                y[ii$nugget] <- inv_exp2(x[ii$nugget]/sqrt(x[ii1$nugget]*x[ii2$nugget]))
            }
        }
    }
    return(y)
}
 
invlink_unc03 <- function(x){
    # figure out indices of cross covariances
    ncomp <- ncomp_from_nparms( length(x) )
    y <- rep(NA,length(x))
    for(j1 in 1:ncomp){
        for(j2 in 1:j1){
            ii <- multi_matern_parm_index(ncomp, j1, j2 )
            y[ii$range] <- log(x[ii$range])
            y[ii$smoothness] <- log(x[ii$smoothness])
            if(j1 == j2){
                y[ii$variance] <- log(x[ii$variance])
                y[ii$nugget] <- log(x[ii$nugget])
            } else {
                ii1 <- multi_matern_parm_index(ncomp, j1, j1)
                ii2 <- multi_matern_parm_index(ncomp, j2, j2)
                y[ii$variance] <- inv_tan2(x[ii$variance]/sqrt(x[ii1$variance]*x[ii2$variance]))
                y[ii$nugget] <- inv_tan2(x[ii$nugget]/sqrt(x[ii1$nugget]*x[ii2$nugget]))
            }
        }
    }
    return(y)
}
 

    

