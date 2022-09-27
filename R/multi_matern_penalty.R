
#plot_specden <- function(x){
#    par(mfrow=c(1,2))
#    fields::image.plot( Re(x) )
#    fields::image.plot( Im(x) )
#    par(mfrow=c(1,1))
#}

##############################################
##### penalties, derivatives, hessian
##############################################

pen_logdet_cross_spec <- function(covparms, effrange, nwrap = 2 ){
    # compute the cross spectral densities and their cholesky decomps
    cross_spec <- compute_cross_spec( covparms, effrange, nwrap )
    chol_spec <- cholesky_array( cross_spec )
    # return the negative sum_logdet penalty
    return( -sum_logdet( chol_spec ) )
}

dpen_logdet_cross_spec <- function(covparms, effrange, wrap = 2 ){
    # compute the cross spectral densities and their cholesky decomps
    nparms <- length(covparms)
    cross_spec <- compute_cross_spec( covparms, effrange, nwrap )
    chol_spec <- cholesky_array( cross_spec )
    # compute the derivative of the cross spectral densities
    d_cross_spec <- compute_d_cross_spec( covparms, effrange, nwrap = 2 )
    # compute Li_Sj_Lit
    LSL <- array(NA, dim(d_cross_spec))
    for(j in 1:nparms){ LSL[,,,,j] <- Li_Sj_Lit( chol_spec, d_cross_spec[,,,,j] ) }
    # compute derivative of logdets using LSL
    dpn <- rep(NA, nparms )
    for(j in 1:nparms){ dpn[j] <- d_sum_logdet( LSL[,,,,j] ) }
    return(dpn)
}

ddpen_logdet_cross_spec <- function(covparms, effrange, wrap = 2 ){
    # compute the cross spectral densities and their cholesky decomps
    nparms <- length(covparms)
    cross_spec <- compute_cross_spec( covparms, effrange, nwrap )
    chol_spec <- cholesky_array( cross_spec )
    # compute the derivative of the cross spectral densities
    d_cross_spec <- compute_d_cross_spec( covparms, effrange, nwrap = 2 )
    # compute Li_Sj_Lit
    LSL <- array(NA, dim(d_cross_spec))
    for(j in 1:nparms){ LSL[,,,,j] <- Li_Sj_Lit( chol_spec, d_cross_spec[,,,,j] ) }
    # compute "second derivative" of logdets using LSL
    # (not exactly the right expression)
    ddpn <- array(NA, c(nparms,nparms) )
    for(j1 in 1:nparms){ for(j2 in 1:nparms){
        ddpn[j1,j2] <- dd_sum_logdet( LSL[,,,,j1], LSL[,,,,j2] )
    }}
    return(ddpn)
}


approx_max_range <- function(locs){
    locs0 <- locs[ sample(1:nrow(locs),400), ]
    return( max( fields::rdist(locs0) ) )
}

##############################################
##### cross spectral density and derivatives
##############################################

compute_cross_spec <- function(covparms, effrange, nwrap = 2){

    # matrix of covparms and number of components
    parm_mat <- multi_parms_mat( covparms )
    ncomp <- nrow( parm_mat$var )

    # number of gridpoints in the integration
    # also deserves some scrutiny
    ngrid <- 100

    # initialize array of ncomp x ncomp cross spectral density matrices
    cross_spec <- array( NA, c(ngrid,ngrid,ncomp,ncomp) )

    for( j1 in 1:ncomp){ for(j2 in 1:ncomp){
        # use a nugget of zero because we want to evaluate whether
        # the underlying matern process is valid, not allow an
        # invalid matern process to be "saved" by adding a nugget
        matern_parms <- c(parm_mat$var[j1,j2], parm_mat$range[j1,j2], parm_mat$smoothness[j1,j2], 0)
        cov_arr <- wrap_cov( matern_parms, ngrid, effrange, 2 )
        cross_spec[,,j1,j2] <- fft(cov_arr)
    }}
    return(cross_spec)
}


compute_d_cross_spec <- function( covparms, effrange, nwrap = 2 ){

    # matrix of covparms and number of components
    parm_mat <- multi_parms_mat( covparms )
    ncomp <- nrow( parm_mat$var )

    # number of gridpoints in the integration
    # also deserves some scrutiny
    ngrid <- 100
    d_cross_spec <- array(0, c( ngrid, ngrid, ncomp, ncomp, length(covparms) ) )

    # take fft of derivatives of covariances with respect to each parameter
    # for every pair of components
    # this one only has to loop over every unique pair, i.e. j2 <= j1
    for( j1 in 1:ncomp){ for(j2 in 1:j1){
        matern_parms <- c(parm_mat$var[j1,j2], parm_mat$range[j1,j2], parm_mat$smoothness[j1,j2], 0)
        d_cov_arr <- d_wrap_cov( matern_parms, ngrid, effrange, 2 )
        # figure how to assign to the right index
        index <- multi_matern_parm_index( ncomp, j1, j2 )
        for( k in 1:3){
            i <- index[[k]]
            d_cross_spec[,,j1,j2,i] <- fft(d_cov_arr[,,k])
            if(j1 != j2){ d_cross_spec[,,j2,j1,i] <- Conj(d_cross_spec[,,j1,j2,i]) }
        }
    }}
    return(d_cross_spec)
}



##############################################
##### lower-level penalties
##############################################

sum_logdet <- function( cholcs ){
    ss <- 0
    dimcs <- dim(cholcs)
    ncomp <- dimcs[ length(dimcs) ]

    for(j in 1:ncomp ){
        ss <- ss + 2*sum( log( Re(cholcs[,,j,j]) ) )
    }
    return(ss)
}

d_sum_logdet <- function( LSLj ){
    return( array_trace( LSLj ) )
}

dd_sum_logdet <- function( LSLj, LSLk ){
    LSLjk <- array_mat_mult( LSLj, LSLk )
    return( array_trace( LSLjk ) )
}



##############################################
##### linear algebra functions
##############################################

array_trace <- function(A){
    # array of matrices A
    dimA <- dim(A)
    ncomp <- dimA[ length(dimA) ]
    ss <- 0
    for(k in 1:ncomp){ ss <- ss + sum(Re(A[,,k,k])) }
    return(ss)
}
    


forward_solve <- function( mat, vec ){
    n <- length(vec)
    sol <- rep(NA, n)
    sol[1] <- vec[1]/mat[1,1]

    for(j in 2:n){
        sol[j] <- ( vec[j] - sol[1:(j-1)] %*% mat[j,1:(j-1)] )/mat[j,j]
    }
    return(c(sol))
}


array_mat_mult <- function( A, B ){

    # arrays of matrices A and B
    # get the dimensions
    dimA <- dim(A)
    ncomp <- dimA[ length(dimA) ]
    
    # initialize array of products
    P <- array(0, dimA )

    # do the multiplication
    for(j1 in 1:ncomp){ for(j2 in 1:ncomp){
        temp <- P[,,j1,j2]
        for(j3 in 1:ncomp){
            temp <- temp + A[,,j1,j3]*B[,,j3,j2]
        }
        P[,,j1,j2] <- temp
    }}
    return(P)
}


forward_solve_array <- function( L, B ){

    # array of lower triangular matrices L, array of right hand side matrices B

    # get the dimensions
    dimL <- dim(L)
    ncomp <- dimL[ length(dimL) ]
    
    # array of solutions
    sol <- array(0, dimL )

    # do each solve separately. sol and B should always have last index set to k
    for(k in 1:ncomp){
        
        sol[,,1,k] <- B[,,1,k]/L[,,1,1]

        for(j in 2:ncomp){
            temp <- B[,,j,k]
            for( el in 1:(j-1) ){
                temp <- temp - sol[,,el,k]*L[,,j,el]
            }
            sol[,,j,k] <- temp/L[,,j,j]
        }
    }

    return(sol)
}


cholesky_array <- function( cs ){

    # get the dimensions
    dimcs <- dim(cs)
    ncomp <- dimcs[ length(dimcs) ]
    
    # initialize cholesky array
    cholcs <- array(0, dimcs) 

    # (1,1) entry
    cholcs[,,1,1] <- sqrt( cs[,,1,1] )
    if( ncomp == 1 ){ return( cholcs ) }

    # (2,1) entry
    cholcs[,,2,1] <- cs[,,2,1]/cholcs[,,1,1]
    # (2,2) entry
    cholcs[,,2,2] <- sqrt( cs[,,2,2] - cholcs[,,2,1]*Conj( cholcs[,,2,1] ) )
    if( ncomp == 2 ){ return( cholcs ) }

    # loop over rows
    for(j in 3:ncomp){

        # do a forward solve to get the first j-1 columns of row j # a^j = L^(j-1) ell^j
        # get entry (j,1)
        cholcs[,,j,1] <- cs[,,1,j]/cholcs[,,1,1]

        # entries (j,2) through (j,j-1)
        for( k in 2:(j-1) ){

            temp <- cs[,,j,k]

            # subtract off inner product, doing each term in a loop
            for( el in 1:(k-1) ){
                temp <- temp - cholcs[,,j,el]*cholcs[,,k,el]
            }

            # divide by the diagonal
            cholcs[,,j,k] <- temp/cholcs[,,k,k]
        }

        # get entry (j,j)
        temp <- cs[,,j,j]
        for( k in 1:(j-1)){
            temp <- temp - cholcs[,,j,k]*Conj( cholcs[,,j,k] )
        }
        cholcs[,,j,j] <- sqrt( temp )
    }
    
    return(cholcs)
}


Li_Sj_Lit <- function( cholcs, d_cross_specden_j ){

    # get the dimensions
    dimcs <- dim(cholcs)
    ncomp <- dimcs[ length(dimcs) ]
    
    A <- forward_solve_array( cholcs, d_cross_specden_j )
    A <- aperm( A, c(1,2,4,3) )
    A <- forward_solve_array( cholcs, A )
    A <- aperm( A, c(1,2,4,3) )
    return(A)

}








##############################################
##### parameter index operations
##############################################

multi_parms_mat <- function( covparms ){
    # length(covparms) is 3*ncomp*(ncomp+1)/2 + ncomp
    for(j in 1:10){
        if( length(covparms) == 3*j*(j+1)/2 + j ){
            ncomp <- j
        }
    }
    vmat <- matrix(NA, ncomp, ncomp)
    vst <- 0
    rmat <- matrix(NA, ncomp, ncomp)
    rst <- ncomp*(ncomp+1)/2
    smat <- matrix(NA, ncomp, ncomp)
    sst <- 2*rst
    cnt <- 0
    for(j1 in 1:ncomp){
        for(j2 in 1:j1){
            k <- (j1-1)*j1/2 + j2         
            vmat[j1,j2] <- covparms[k]
            vmat[j2,j1] <- covparms[k]
            rmat[j1,j2] <- covparms[rst + k]
            rmat[j2,j1] <- covparms[rst + k]
            smat[j1,j2] <- covparms[sst + k]
            smat[j2,j1] <- covparms[sst + k]
        }
    }
    return( list( variance=vmat, range=rmat, smoothness=smat ) )
}


multi_matern_parm_index <- function( ncomp, j1, j2 ){

    v1 <- 0
    r1 <- ncomp*(ncomp+1)/2
    s1 <- 2*r1

    k <- (j1-1)*j1/2 + j2
    index <- list()
    index$variance   <- v1 + k
    index$range      <- r1 + k
    index$smoothness <- s1 + k
    
    return(index)
}


##############################################
##### wrapping functions
##############################################

wrap_cov <- function( matern_parms, ngrid, effrange, nwrap = 1 ){
    xgrid <- as.matrix( expand.grid( (0:(ngrid-1))*effrange/ngrid, (0:(ngrid-1))*effrange/ngrid) )
    cov_vec <- array(0, nrow(xgrid) )
    for( j1 in (-nwrap):(nwrap-1) ){
        for( j2 in (-nwrap):(nwrap-1) ){
            xx <- xgrid
            xx[,1] <- xgrid[,1] + j1*effrange
            xx[,2] <- xgrid[,2] + j2*effrange
            cov_vec <- cov_vec + c(matern_isotropic2( matern_parms, matrix(0,1,2), xx ))
        }
    }
    cov_arr <- array( cov_vec, c(ngrid,ngrid) )
    return(cov_arr)
}

    
d_wrap_cov <- function( matern_parms, ngrid, effrange, nwrap = 1 ){
    nparms <- length(matern_parms)
    xgrid <- as.matrix( expand.grid( (0:(ngrid-1))*effrange/ngrid, (0:(ngrid-1))*effrange/ngrid) )
    d_cov_vec <- array(0, c(nrow(xgrid),nparms))
    for( j1 in (-nwrap):(nwrap-1) ){
        for( j2 in (-nwrap):(nwrap-1) ){
            xx <- xgrid
            xx[,1] <- xgrid[,1] + j1*effrange
            xx[,2] <- xgrid[,2] + j2*effrange
            d_cov_vec <- d_cov_vec + c(d_matern_isotropic2( matern_parms, matrix(0,1,2), xx ))
        }
    }
    d_cov_arr <- array( NA, c(ngrid,ngrid,nparms) )
    for(j in 1:nparms){
        d_cov_arr[ , , j ] <- array( d_cov_vec[,j], c(ngrid,ngrid) )
    }
    return(d_cov_arr)
}



