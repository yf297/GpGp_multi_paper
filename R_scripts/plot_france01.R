
library(fields)
library(viridis)
library(maps)
library(mapdata)

load("../data/formatted/france01.RData")
dat <- as.data.frame(france01)


# read in original data to figure out how to transform coordinates back
load("../data/processed/france_composite_analysis.RData")

# subset to these two 
comp1 <- "n_tot_31_1"
comp2 <- "zn_tot_hf"

france$date <- as.Date( france$date, format = "%d/%m/%Y" )
i0 <- france$date >= "2006-01-01" & france$date <= "2008-12-31"
france <- france[i0,]

i1 <- france$component == comp1
france1 <- france[i1,]
i2 <- france$component == comp2
france2 <- france[i2,]

# locs
lonlat <- france1[, c("long","lat")]
meanlat <- mean( lonlat$lat )

dat$long <- dat[,2]/cos( 2*pi*meanlat/360 ) + mean(lonlat[,1])
dat$lat <- dat[,3] + mean(lonlat[,2] )

# define response and component
dat$logvalue <- dat[,1]
dat$component <- dat[,4]


titles <- c("log N (g/kg)","log Zn (mg/kg)")

pdf("../figures/plot_france01.pdf",width=8.0,height=2.0)
par(mfrow=c(1,2), family = "serif", mar=c(1,2,3,1), oma = c(0,9,0,11) )
for(j in 1:2){
    ii <- dat$component == j
    quilt.plot( dat$long[ii], dat$lat[ii], dat$logvalue[ii],
        axes = FALSE, col = viridis(64) )
    box()
    mtext(titles[j], side=3, line = 1)
    map("worldHires", add = TRUE )
}
dev.off()
    
    
