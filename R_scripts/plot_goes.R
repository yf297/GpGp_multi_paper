
library(fields)
library(viridis)

load("../data/formatted/goes01_1_2_3_4.RData")
dat <- as.data.frame( goes01_1_2_3_4 )

titles <- c("Band 1","Band 6", "Band 7", "Band 9")

pdf("../figures/plot_goes.pdf",width=8,height=2.0)
par(mfrow=c(1,4), family = "serif", mar=c(1,2,3,1), oma = c(0,0,0,2) )
for(j in 1:4){
    ii <- dat$j == j
    mat <- matrix(NA, 44, 54)
    mat[ cbind( dat$row[ii] - 37, dat$col[ii] - 27 ) ] <- dat$y[ii]
    image.plot( mat, axes = FALSE, col = viridis(64) )
    mtext(titles[j], side=3, line = 1)
}
dev.off()
    
    
