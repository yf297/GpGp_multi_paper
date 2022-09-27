
library(fields)
library(viridis)
library(maps)
library(mapdata)

load("../data/formatted/bang_As_Fe_Mn_P.RData")
dat <- bang_As_Fe_Mn_P


titles <- c("log As (ug/L)","log Fe (mg/L)", "log Mn (mg/L)", "log P (mg/L)")

pdf("../figures/plot_bang.pdf",width=8,height=2.0)
par(mfrow=c(1,4), family = "serif", mar=c(1,2,3,1), oma = c(0,0,0,2) )
for(j in 1:4){
    ii <- dat$component == j
    quilt.plot( dat$long[ii], dat$lat[ii], dat$logvalue[ii],
        axes = FALSE, col = viridis(64) )
    box()
    mtext(titles[j], side=3, line = 1)
    map("worldHires", add = TRUE )
}
dev.off()
    
    
