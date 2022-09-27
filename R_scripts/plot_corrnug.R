
library("GpGp")

load("../fits/fit_corrnug_1.RData")
c1 <- fit[[4]]$covparms

load("../fits/fit_corrnug_2.RData")
c2 <- fit[[4]]$covparms

i0 <- c(1,4,7,10)
c111 <- c1[i0+0]
c112 <- c1[i0+1]
c122 <- c1[i0+2]

c211 <- c2[i0+0]
c212 <- c2[i0+1]
c222 <- c2[i0+2]

x <- seq(-1.5,1.5,by=0.01)
j0 <- which( x==0 )
locs <- matrix(x, length(x), 1)

pdf("../figures/corrnug.pdf", width=9,height=3)
par(mfrow=c(1,3))

plot( 0, type="n", xlim = range(x), ylim = c(0,8), xlab = "", ylab="" )
mtext("Marginal As", line = 1)
parms <- c111
nug <- parms[4]
parms[4] <- 0
covvec <- matern_isotropic( covparms = parms, locs = locs )[,j0]
lines(x,covvec, type="l")
points(x[j0], covvec[j0] + nug)

parms <- c211
nug <- parms[4]
parms[4] <- 0
covvec <- matern_isotropic( covparms = parms, locs = locs )[,j0]
lines(x,covvec, type="l", col = "magenta")
points(x[j0], covvec[j0] + nug, col = "magenta")

plot( 0, type="n", xlim = range(x), ylim = c(0,4.5), xlab = "", ylab="" )
mtext("Cross As,Fe", line = 1)
parms <- c112
nug <- parms[4]
parms[4] <- 0
covvec <- matern_isotropic( covparms = parms, locs = locs )[,j0]
lines(x,covvec, type="l")
points(x[j0], covvec[j0] + nug)

parms <- c212
nug <- parms[4]
parms[4] <- 0
covvec <- matern_isotropic( covparms = parms, locs = locs )[,j0]
lines(x,covvec, type="l", col = "magenta")
points(x[j0], covvec[j0] + nug, col = "magenta")

plot( 0, type="n", xlim = range(x), ylim = c(0,4.5), xlab = "", ylab="" )
mtext("Marginal Fe", line = 1)
parms <- c122
nug <- parms[4]
parms[4] <- 0
covvec <- matern_isotropic( covparms = parms, locs = locs )[,j0]
lines(x,covvec, type="l")
points(x[j0], covvec[j0] + nug)

parms <- c222
nug <- parms[4]
parms[4] <- 0
covvec <- matern_isotropic( covparms = parms, locs = locs )[,j0]
lines(x,covvec, type="l", col = "magenta")
points(x[j0], covvec[j0] + nug, col = "magenta")
dev.off()



pdf("../figures/corrnug22.pdf", width=7,height=7)
par(mfrow=c(2,2),mar=c(2,4,4,1))

plot( 0, type="n", xlim = range(x), ylim = c(0,8), xlab = "", ylab="" )
mtext("Marginal As", line = 1)
parms <- c111
nug <- parms[4]
parms[4] <- 0
covvec <- matern_isotropic( covparms = parms, locs = locs )[,j0]
lines(x,covvec, type="l",lwd=2)
points(x[j0], covvec[j0] + nug, cex = 2)


parms <- c211
nug <- parms[4]
parms[4] <- 0
covvec <- matern_isotropic( covparms = parms, locs = locs )[,j0]
lines(x,covvec, type="l", col = "magenta",lwd=2)
points(x[j0], covvec[j0] + nug, col = "magenta", cex = 2)

plot(0,type="n",axes=FALSE,ann=FALSE)

plot( 0, type="n", xlim = range(x), ylim = c(0,4.5), xlab = "", ylab="" )
mtext("Cross As,Fe", line = 1)
parms <- c112
nug <- parms[4]
parms[4] <- 0
covvec <- matern_isotropic( covparms = parms, locs = locs )[,j0]
lines(x,covvec, type="l",lwd=2)
points(x[j0], covvec[j0] + nug, cex = 2)

parms <- c212
nug <- parms[4]
parms[4] <- 0
covvec <- matern_isotropic( covparms = parms, locs = locs )[,j0]
lines(x,covvec, type="l", col = "magenta",lwd=2)
points(x[j0], covvec[j0] + nug, col = "magenta", cex = 2)

plot( 0, type="n", xlim = range(x), ylim = c(0,4.5), xlab = "", ylab="" )
mtext("Marginal Fe", line = 1)
parms <- c122
nug <- parms[4]
parms[4] <- 0
covvec <- matern_isotropic( covparms = parms, locs = locs )[,j0]
lines(x,covvec, type="l",lwd=2)
points(x[j0], covvec[j0] + nug, cex = 2)

parms <- c222
nug <- parms[4]
parms[4] <- 0
covvec <- matern_isotropic( covparms = parms, locs = locs )[,j0]
lines(x,covvec, type="l", col = "magenta",lwd=2)
points(x[j0], covvec[j0] + nug, col = "magenta", cex = 2)
dev.off()




