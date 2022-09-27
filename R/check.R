library(matrixcalc)

makeSymm <- function(m) {
   m[upper.tri(m)] <- t(m)[upper.tri(m)]
   return(m)
}

is.cnd <- function(M){
	p <- nrow(M)
	new <- matrix(0,p,p)
	for(i in 1:p){
	    for(j in 1:p){
		new[i,j] <- M[i,p] + M[p,j] - M[i,j] - M[p,p]
	    }
	}
	return(is.positive.semi.definite(new))
}


check_flex_a <- function(covparms,d){

   ncomp <- ncomp_from_nparms(length(covparms))
   neach <- ncomp*(ncomp+1)/2
   ran_inds <- (neach + 1):(2*neach)
   Z <- matrix(NA,ncomp,ncomp)	
   Z[upper.tri(Z, diag=T)] <- covparms[ran_inds]
   ran <- t(Z)
   ran <- makeSymm(ran)
   rann2 <- ran^(-2) 

  b2 <-  is.cnd(rann2)

  return(b2)
}


check_flex_b <- function(covparms){

   ncomp <- ncomp_from_nparms(length(covparms))
   neach <- ncomp*(ncomp+1)/2
   smo_inds <- (2*neach + 1):(3*neach)  	
   Z <- matrix(NA,ncomp,ncomp)	
   Z[upper.tri(Z, diag=T)] <- covparms[smo_inds]
   smo <- t(Z)
   smo <- makeSymm(smo)

   b0 <- is.cnd(smo) 

  return(b0)
}
