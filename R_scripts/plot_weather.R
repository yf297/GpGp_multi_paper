
library(fields)
library(viridis)
library(maps)
library(mapdata)

data("weather",package="RandomFields")
weather <- as.data.frame(weather)

#load("../data/formatted/weather.RData")
#dat <- as.data.frame(weather)
#dat$long <- dat$X
#dat$lat <- dat$Y
#dat$component <- dat[,5]
#dat$response <- dat[,1]


titles <- c("Pressure (Pa)","Temperature (C)")

pdf("../figures/plot_weather.pdf",width=9,height=2.50)
par(mfrow=c(1,2), family = "serif", mar=c(1,2,2,3), oma = c(0,7,0,9) )
for(j in 1:2){
    quilt.plot( weather$lon, weather$lat, weather[,j],
        axes = FALSE, col = viridis(64) )
    box()
    mtext(titles[j], side=3, line = 0.5, cex = 1.25)
    map("worldHires", add = TRUE )
}
dev.off()
    
    
