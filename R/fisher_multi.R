

check_step <- function(step,silent){

    tol <- 1
    # if step size large, then make it smaller
    if (sum(step^2) > tol) {
        if(!silent) cat(sprintf("step too large, norm = %4.2f\n", sqrt(sum(step^2))))
        step <- sqrt(tol)*step/sqrt(sum(step^2))
    }
    return(step)

}

check_likelihood <- function(likobj, likobj0){

    pass <- TRUE
    allvals <- c( likobj$loglik, likobj$grad, c(likobj$info) )
    if( sum(is.na(allvals)) > 0  ||  sum( abs(allvals) == Inf ) > 0 ){
        print(likobj$loglik)
        pass <- FALSE
    }
    tol <- 0.0
    if( likobj$loglik > likobj0$loglik + tol ){
        print(likobj$loglik)
        pass <- FALSE 
    }
    return(pass)

}


fisher_scoring_multi0 <- function( likfun, start_parms, link, 
    silent = FALSE, convtol = 1e-4, max_iter = 40, invlink_startparms, active ){
    
    # function for checking wolfe conditions
    wolfe_check <- function(likobj0,likobj1,logparms,step,both){
        c1 <- 1e-4
        c2 <- 0.9
        tol <- 0.1
        ll0 <- likobj0$loglik
        gr0 <- likobj0$grad
        ll1 <- likobj1$loglik
        gr1 <- likobj1$grad
        if(!both){
            satfd <- ll1 <= ll0 + c1*crossprod(step,gr0) + tol
        } else {
            satfd <- ll1 <= ll0 + c1*crossprod(step,gr0) + tol &&
                 -crossprod(step,gr1) <= -c2*crossprod(step,gr0) + tol
        }
        return(satfd)
    }
    
    # evaluate function at initial values
    logparms <- start_parms
    likobj <- likfun(logparms)
    
    # test likelihood object    
    if( !test_likelihood_object(likobj) ){
        stop("Invalid Starting values")
        logparms <- 0.1*logparms
        likobj <- likfun(logparms)
    }
    
    # assign loglik, grad, and info
    loglik <- likobj$loglik        
    grad <- likobj$grad
    info <- as.matrix(likobj$info)
    # add a small amount of regularization
    info <- info + 0.1*min(diag(info))*diag(nrow(info))

    # print some stuff out
    lp <- rep(NA,length(invlink_startparms))
    lp[active] <- logparms
    lp[!active] <- invlink_startparms[!active]
    
    if(!silent){
        cat(paste0("Iter ",0,": \n"))
        cat("pars = ",  paste0(round(link(lp),4)), "  \n" )
        cat(paste0("loglik = ", round(-loglik,6),         "  \n"))
        cat("grad = ")
        cat(as.character(round(-grad,3)))
        cat("\n\n")
    }
    
    iter = 0
    for(j in seq_len(max_iter)){
        
        likobj0 <- likobj
        
        # if condition number of info matrix large, 
        # then gradient descent
        tol <- 1e-10
        if ( 1/condition_number(info) < tol) {
            if (!silent) cat("Cond # of info matrix > 1/tol \n")
            #info <- 1.0*max(likobj0$info)*diag(nrow(likobj0$info))
            # regularize
            ee <- eigen(info)
            ee_ratios <- ee$values/max(ee$values)
            ee_ratios[ ee_ratios < 1e-8 ] <- 1e-8
            ee$values <- max(ee$values)*ee_ratios
            info <- ee$vectors %*% diag(ee$values) %*% t(ee$vectors)
            #diag(info) <- diag(info) + tol*max(diag(info))
        }

        # calculate fisher step 
        step <- -solve(info, grad)
        
        # if step size large, then make it smaller
        if (sum(step^2) > 1) {
            if(!silent) cat("##\n")
            step <- step/sqrt(sum(step^2))
        }
        
        # take step and calculate loglik, grad, and info
        newlogparms <- logparms + step
        likobj <- likfun(newlogparms)
        
        # check for Inf, NA, or NaN
        cnt <- 1
        while (!test_likelihood_object(likobj)) {
            if (!silent) cat("inf or na or nan in likobj\n")
            step <- 0.5 * step
            newlogparms <- logparms + step
            likobj <- likfun(newlogparms)
            if (cnt == 10) { stop("could not avoid inf, na or nan\n") }
        }

        # check whether negative likelihood decreased
        # take smaller stepsize
        if( likobj$loglik > likobj0$loglik ){
            step <- 0.25*step
            newlogparms <- logparms + step
	    if(!silent){print("loglik decreased")}
            likobj <- likfun(newlogparms)
        }
            
        # check again, move along gradient
        #if( likobj$loglik > likobj0$loglik ){
        #    info0 <- diag( rep(mean(diag(info)),nrow(info)) )
        #    step <- -solve(info0,grad)
        #    newlogparms <- logparms + step
	#    if(!silent){print("moving along gradient")}       
        #    likobj <- likfun(newlogparms)
        #}
            
        # check once move, take smaller step along gradient
        if( likobj$loglik > likobj0$loglik ){
            info0 <- diag( rep(max(diag(info)),nrow(info)) )
            step <- -solve(info0,grad)
            if(!silent){print("moving along gradient")}              
            newlogparms <- logparms + step
            likobj <- likfun(newlogparms)
        }
        
        # Check the wolfe conditions
        # if not satisfied, shrink fisher step
        # after some iterations, switch to gradient
        cnt <- 1
        no_decrease <- FALSE
        both <- FALSE
        mult <- 1.0
        #if(!wolfe_check(likobj0,likobj,logparms,newlogparms-logparms,both) &&
        #        !no_decrease ){
        #    cat("@@\n")
        #    step <- -0.5*sqrt(mean(step^2))*grad/sqrt(sum(grad^2))
        #}
        #while (!wolfe_check(likobj0,likobj,logparms,newlogparms-logparms,both) &&
        #        !no_decrease ){
        #    info <- 1/mult*max(likobj$info)*diag(nrow(likobj0$info))
        #    step <- -solve(info,grad)
        #    if(!silent) cat("**\n") 
        #    if ( sqrt(sum(step^2)) < 1e-4 ){ no_decrease <- TRUE }  # maybe we should throw error here?
        #    newlogparms <- logparms + step
        #    likobj <- likfun(newlogparms)
        #    cnt <- cnt + 1
        #    mult <- mult*0.5
        #}
        stepgrad <- c(crossprod(step,grad))
        
        # redefine logparms, loglik, grad, info
        logparms <- logparms + step
        loglik <- likobj$loglik        
        grad <- likobj$grad
        info <- likobj$info

        lp[active] <- logparms
        # print some stuff out
        if(!silent){
            cat(paste0("Iter ",j,": \n"))
            cat("pars = ",  paste0(round(link(lp),4)), "  \n" )
            cat(paste0("loglik = ", round(-loglik,6),         "  \n"))
            cat("grad = ")
            cat(as.character(round(-grad,4)),"\n")
            cat("step dot grad = ",stepgrad,"\n")
            cat("\n")
        }
        
        # if gradient is small, then break and return the results        
        iter <- j
        if( abs(stepgrad) < convtol || no_decrease ){
	    break
        }
    }

    # collect information to return
    betahatinfo <- likobj        
    betahat <- as.vector(betahatinfo$betahat)
    betacov <- solve(betahatinfo$betainfo)
    sebeta <- sqrt(diag(betacov))
    tbeta <- betahat/sebeta

    ret <- list(
        covparms = link(lp), 
        logparms = lp,
        betahat = betahat, 
        sebeta = sebeta,
        betacov = betacov,
        tbeta = tbeta,
        loglik = loglik,
        no_decrease = no_decrease,
        grad = likobj$grad,
        info = likobj$info,
        conv = ( abs(stepgrad) < convtol || no_decrease ),
	iter = iter
    )
    return(ret)
}


fisher_scoring_multi <- function( likfun, start_parms, link, 
    silent = FALSE, convtol = 1e-4, max_iter = 40, start_logparms, active ){
    
    # evaluate function at initial values
    logparms <- start_logparms[active]
    likobj <- likfun(logparms)
    
    # test likelihood object    
    if( !test_likelihood_object(likobj) ){
        stop("Invalid Starting values")
        logparms <- 0.1*logparms
        likobj <- likfun(logparms)
    }
    
    # assign loglik, grad, and info
    loglik <- likobj$loglik        
    grad <- likobj$grad
    info <- as.matrix(likobj$info)
    # add a small amount of regularization
    info <- info + 0.1*min(diag(info))*diag(nrow(info))

    # print some stuff out
    lp <- rep(NA,length(start_logparms))
    lp[active] <- logparms
    lp[!active] <- start_logparms[!active]
    
    if(!silent){
        cat(paste0("Iter ",0,": \n"))
        cat("pars = ",  paste0(round(link(lp),4)), "  \n" )
        cat(paste0("loglik = ", round(-loglik,6),         "  \n"))
        cat("grad = ")
        cat(as.character(round(-grad,3)))
        cat("\n\n")
    }
    
    iter <- 0
    for(j in seq_len(max_iter)){

        iter <- j
        
        likobj0 <- likobj
        
        # if condition number of info matrix large, regularize
        tol <- 1e-10
        if ( 1/condition_number(info) < tol) {
            if (!silent) cat("Cond # of info matrix > 1/tol \n")
            # regularize
            ee <- eigen(info)
            ee_ratios <- ee$values/max(ee$values)
            ee_ratios[ ee_ratios < 1e-5 ] <- 1e-5
            ee$values <- max(ee$values)*ee_ratios
            info <- ee$vectors %*% diag(ee$values) %*% t(ee$vectors)
        }

        # calculate fisher step 
        step <- -solve(info, grad)
        step <- check_step(step,silent)

        # take step and calculate loglik, grad, and info
        newlogparms <- logparms + step
        likobj <- likfun(newlogparms)
        
        # take smaller step along fisher direction if bad step
        if( !check_likelihood(likobj, likobj0) ){
            if (!silent) cat("increase or inf or na or nan in likobj\n")
            step <- 0.25 * step
            newlogparms <- logparms + step
            likobj <- likfun(newlogparms)
        }

        # take step along gradient if bad step
        if( !check_likelihood(likobj, likobj0) ){
            info0 <- diag( rep(max(diag(info)),nrow(info)) )
            step <- -0.25*solve(info0,grad)
            step <- check_step(step,silent)
            if(!silent){print("moving along gradient")}              
            newlogparms <- logparms + step
            likobj <- likfun(newlogparms)
        }

        # take step along gradient if bad step
        if( !check_likelihood(likobj, likobj0) ){
            info0 <- diag( rep(max(diag(info)),nrow(info)) )
            step <- -0.05*solve(info0,grad)
            step <- check_step(step,silent)
            if(!silent){print("moving along gradient")}              
            newlogparms <- logparms + step
            likobj <- likfun(newlogparms)
        }

        # take step along gradient if bad step
        if( !check_likelihood(likobj, likobj0) ){
            info0 <- diag( rep(max(diag(info)),nrow(info)) )
            step <- -0.01*solve(info0,grad)
            step <- check_step(step,silent)
            if(!silent){print("moving along gradient")}              
            newlogparms <- logparms + step
            likobj <- likfun(newlogparms)
        }

        no_decrease <- FALSE
        # set no_decrease to true if we can't move
        if( !check_likelihood(likobj, likobj0) ){
            if(!silent){print("could not increase along gradient")}              
            no_decrease <- TRUE
	    # revert back to previous step
	    likobj <- likobj0
	    step <- rep(0,length(logparms))
        }

        stepgrad <- c(crossprod(step,grad))
        
        # redefine logparms, loglik, grad, info
        logparms <- logparms + step
        loglik <- likobj$loglik        
        grad <- likobj$grad
        info <- likobj$info

        lp[active] <- logparms
        # print some stuff out
        if(!silent){
            cat(paste0("Iter ",j,": \n"))
            cat("pars = ",  paste0(round(link(lp),4)), "  \n" )
            cat(paste0("loglik = ", round(-loglik,6),         "  \n"))
            cat("grad = ")
            cat(as.character(round(-grad,4)),"\n")
            cat("step dot grad = ",stepgrad,"\n")
            cat("\n")
        }
        
        # if gradient is small, then break and return the results        
        if( abs(stepgrad) < convtol || no_decrease ){
	    break
        }
    }

    # collect information to return
    betahatinfo <- likobj        
    betahat <- as.vector(betahatinfo$betahat)
    betacov <- solve(betahatinfo$betainfo)
    sebeta <- sqrt(diag(betacov))
    tbeta <- betahat/sebeta

    ret <- list(
        covparms = link(lp), 
        logparms = lp,
        betahat = betahat, 
        sebeta = sebeta,
        betacov = betacov,
        tbeta = tbeta,
        loglik = loglik,
        no_decrease = no_decrease,
        grad = likobj$grad,
        info = likobj$info,
        conv = ( abs(stepgrad) < convtol || no_decrease ),
	iter = iter
    )
    return(ret)
}

