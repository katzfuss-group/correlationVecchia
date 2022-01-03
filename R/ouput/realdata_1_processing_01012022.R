####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview:
###
####################################################################################

options("rgdal_show_exportToProj4_warnings"="none")

rm(list = ls())

#########################################################################

library(rgdal)
library(raster)
library(sp)
library(ncdf4)
library(maptools)
library(correlationVecchia)

#########################################################################

crsname     <- "+proj=stere +lat_0=90 +lat_ts=60 +lon_0=263 +k=1 +x_0=3475000 +y_0=7475000 +ellps=WGS84 +datum=WGS84 +units=m"

shapename   <- "SHP/cb_2018_us_state_500k"

# system(paste0("unzip ", shapename, ".zip"))

shapename   <- "cb_2018_us_state_500k"

#########################################################################

vars        <- c("tasmin", "tasmax")
years       <- c("1979", "1981", "1986", "1991", "1996", "2001")
cols        <- hcl.colors(200, "RdBu", rev = TRUE)

region      <- c("Texas", "Oklahoma", "Kansas", "Arkansas", "Louisiana", "Mississippi") # South
# region      <- c("Arizona", "New Mexico", "Utah", "Colorado") # Southwest
# region      <- c("Kansas", "Missouri", "Oklahoma", "Arkansas") # Custom

#########################################################################

day.init    <- as.Date("2001-01-01 06:00:00")
daytime     <- day.init + 181:211

datnames    <- paste0("DATA/", vars, "/", vars, "_CRCM_ncep_", years[6], "010106.nc")

tasmin.CRCM.ncep <- list()
for(i in 1:length(daytime)) {

  tasmin.CRCM.ncep[[i]] <- mappingRaster(rastername = datnames[1], tpoint = i, shapename = shapename, crsname = crsname, states = region)
}

tasmin.ranges <- list()
for(i in 1:length(daytime)) {

  tasmin.ranges[[i]] <- range(rasterToPoints(tasmin.CRCM.ncep[[i]]$crop1)[, 3])
}

range(unlist(tasmin.ranges))

##

tasmax.CRCM.ncep <- list()
for(i in 1:length(daytime)) {

  tasmax.CRCM.ncep[[i]] <- mappingRaster(rastername = datnames[2], tpoint = i, shapename = shapename, crsname = crsname, states = region)
}

tasmax.ranges <- list()
for(i in 1:length(daytime)) {

  tasmax.ranges[[i]] <- range(rasterToPoints(tasmax.CRCM.ncep[[i]]$crop1)[, 3])
}

range(unlist(tasmax.ranges))

#########################################################################

library(dplyr)
library(ggplot2)
library(gridExtra)
library(rasterVis)

breaks    <- seq(247, 300, by = 1)
cols      <- hcl.colors(length(breaks)-1, "RdBu", rev = TRUE)

##

tasmin.plots <- list()
for(i in 1:length(daytime)) {

  tasmin.plots[[i]] <- levelplot(tasmin.CRCM.ncep[[i]]$crop1, margin = FALSE, xlab = NULL, ylab = NULL, scales = list(draw = FALSE), at = breaks, col.regions = cols, main = list(label = as.character(format(daytime[i], "%m/%d/%Y")), cex = 2)) + layer(sp.polygons(tasmin.CRCM.ncep[[i]]$shape, lwd = 1))
}

leg       <- tasmin.plots[[1]]$legend$right$args$key
legGrob   <- lattice::draw.colorkey(key = leg, vp = grid::viewport(height = 0.5))
for(i in 1:length(daytime)) {

  tasmin.plots[[i]]$legend <- list()
}

##

tasmax.plots <- list()
for(i in 1:length(daytime)) {

  tasmax.plots[[i]] <- levelplot(tasmax.CRCM.ncep[[i]]$crop1, margin = FALSE, xlab = NULL, ylab = NULL, scales = list(draw = FALSE), at = breaks, col.regions = cols, main = list(label = as.character(format(daytime[i], "%m/%d/%Y")), cex = 2), colorkey = list(labels = list(cex = 1.75), height = 1.6, width = 3)) + layer(sp.polygons(tasmax.CRCM.ncep[[i]]$shape, lwd = 1))
}

leg       <- tasmax.plots[[1]]$legend$right$args$key
legGrob   <- lattice::draw.colorkey(key = leg, vp = grid::viewport(height = 0.5))
for(i in 1:length(daytime)) {

  tasmax.plots[[i]]$legend <- list()
}

##

# joint.p <- grid.arrange(arrangeGrob(tasmin.plots[[1]], tasmax.plots[[1]], tasmin.plots[[2]], tasmax.plots[[2]], tasmin.plots[[3]], tasmax.plots[[3]], tasmin.plots[[4]], tasmax.plots[[4]],
#                                     tasmin.plots[[5]], tasmax.plots[[5]], tasmin.plots[[6]], tasmax.plots[[6]], tasmin.plots[[7]], tasmax.plots[[7]], tasmin.plots[[8]], tasmax.plots[[8]],
#                                     tasmin.plots[[9]], tasmax.plots[[9]], tasmin.plots[[10]], tasmax.plots[[10]], tasmin.plots[[11]], tasmax.plots[[11]], tasmin.plots[[12]], tasmax.plots[[12]],
#                                     tasmin.plots[[13]], tasmax.plots[[13]], tasmin.plots[[14]], tasmax.plots[[14]], tasmin.plots[[15]], tasmax.plots[[15]], tasmin.plots[[16]], tasmax.plots[[16]],
#                                     tasmin.plots[[17]], tasmax.plots[[17]], tasmin.plots[[18]], tasmax.plots[[18]], tasmin.plots[[19]], tasmax.plots[[19]], tasmin.plots[[20]], tasmax.plots[[20]],
#                                     tasmin.plots[[21]], tasmax.plots[[21]], tasmin.plots[[22]], tasmax.plots[[22]], tasmin.plots[[23]], tasmax.plots[[23]], tasmin.plots[[24]], tasmax.plots[[24]],
#                                     tasmin.plots[[25]], tasmax.plots[[25]], tasmin.plots[[26]], tasmax.plots[[26]], tasmin.plots[[27]], tasmax.plots[[27]], tasmin.plots[[28]], tasmax.plots[[28]],
#                                     tasmin.plots[[29]], tasmax.plots[[29]], tasmin.plots[[30]], tasmax.plots[[30]], tasmin.plots[[31]], tasmax.plots[[31]], ncol = 8), legGrob, ncol = 2, widths = c(10, 1))

joint.p <- grid.arrange(arrangeGrob(tasmin.plots[[1]], tasmin.plots[[15]], tasmin.plots[[31]],
                                    tasmax.plots[[1]], tasmax.plots[[15]], tasmax.plots[[31]], ncol = 3), legGrob, ncol = 2, widths = c(10, 1))


ggsave(file="weather_plot.pdf", joint.p, width = 15.2, height = 10)

#########################################################################

df.tasmin     <- data.frame(x1 = double(), x2 = double(), z = double(), t = double())
for(i in 1:length(daytime)) {

  tempdat       <- as.data.frame(rasterCleaning(tasmin.CRCM.ncep[[i]]$crop1)$dat)
  colnames(tempdat) <- c("x1", "x2", "z")

  tempdat$t     <- as.numeric( daytime[i] - min(daytime) ) / as.numeric( max(daytime) - min(daytime) )

  df.tasmin     <- rbind(df.tasmin, tempdat)
}

##

df.tasmax     <- data.frame(x1 = double(), x2 = double(), z = double(), t = double())
for(i in 1:length(daytime)) {

  tempdat       <- as.data.frame(rasterCleaning(tasmax.CRCM.ncep[[i]]$crop1)$dat)
  colnames(tempdat) <- c("x1", "x2", "z")

  tempdat$t     <- as.numeric( daytime[i] - min(daytime) ) / as.numeric( max(daytime) - min(daytime) )

  df.tasmax     <- rbind(df.tasmax, tempdat)
}

df.joint      <- rbind(data.frame(df.tasmin, d = 0), data.frame(df.tasmax, d = 1))

#########################################################################

save(df.joint, file = "DATA/joint_CRCM_NCEP_10112021.RData")
