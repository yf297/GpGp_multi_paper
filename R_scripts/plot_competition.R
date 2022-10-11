
library("viridis")

load("../data/formatted/3a1.RData")
dat <- as.data.frame(a1)

data_range <- range( dat$response )

pdf("../figures/plot_competition.pdf", width=6,height=2.5)
par(mfrow=c(1,2),mar=c(1,1,1,1),oma=c(0,0,0,2))
ii <- dat$comp == 1
fields::quilt.plot( dat$x[ii], dat$y[ii], dat$response[ii], axes = FALSE, col = viridis(64), zlim = data_range )
ii <- dat$comp == 2
fields::quilt.plot( dat$x[ii], dat$y[ii], dat$response[ii], axes = FALSE, col = viridis(64), zlim = data_range )
dev.off()

