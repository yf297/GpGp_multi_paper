
library(fields)
library(viridis)
library(maps)
library(mapdata)

load("../data/formatted/weather.RData")
dat <- as.data.frame(weather)
dat$long <- dat$X
dat$lat <- dat$Y
dat$component <- dat[,5]
dat$response <- dat[,1]


titles <- c("","")

pdf("../figures/plot_weather.pdf",width=8,height=2.0)
par(mfrow=c(1,2), family = "serif", mar=c(1,2,3,1), oma = c(0,9,0,11) )
for(j in 1:2){
    ii <- dat$component == j
    quilt.plot( dat$long[ii], dat$lat[ii], dat$response[ii],
        axes = FALSE, col = viridis(64) )
    box()
    mtext(titles[j], side=3, line = 1)
    #map("worldHires", add = TRUE )
}
dev.off()
    
    
