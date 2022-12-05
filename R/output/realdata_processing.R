####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script is for data pre-processing. This script generates Figure 10 (weather_plot.pdf).
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
daytime     <- day.init + 151:242

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

breaks    <- seq(246, 312, by = 1)
cols      <- hcl.colors(length(breaks)-1, "RdBu", rev = TRUE)

##

tasmin.plots <- list()
for(i in 1:length(daytime)) {

  tasmin.plots[[i]] <- levelplot(tasmin.CRCM.ncep[[i]]$crop1, margin = FALSE, xlab = NULL, ylab = NULL, scales = list(draw = FALSE), at = breaks, col.regions = cols, main = list(label = as.character(format(daytime[i], "%m/%d/%Y")), cex = 2)) + latticeExtra::layer(sp.polygons(tasmin.CRCM.ncep[[i]]$shape, lwd = 1))
}

leg       <- tasmin.plots[[1]]$legend$right$args$key
legGrob   <- lattice::draw.colorkey(key = leg, vp = grid::viewport(height = 0.5))
for(i in 1:length(daytime)) {

  tasmin.plots[[i]]$legend <- list()
}

##

tasmax.plots <- list()
for(i in 1:length(daytime)) {

  tasmax.plots[[i]] <- levelplot(tasmax.CRCM.ncep[[i]]$crop1, margin = FALSE, xlab = NULL, ylab = NULL, scales = list(draw = FALSE), at = breaks, col.regions = cols, main = list(label = as.character(format(daytime[i], "%m/%d/%Y")), cex = 2), colorkey = list(labels = list(cex = 1.75), height = 1.6, width = 3)) + latticeExtra::layer(sp.polygons(tasmax.CRCM.ncep[[i]]$shape, lwd = 1))
}

leg       <- tasmax.plots[[1]]$legend$right$args$key
legGrob   <- lattice::draw.colorkey(key = leg, vp = grid::viewport(height = 0.5))
for(i in 1:length(daytime)) {

  tasmax.plots[[i]]$legend <- list()
}

##

joint.p <- grid.arrange(arrangeGrob(tasmin.plots[[1]], tasmin.plots[[31]], tasmin.plots[[62]],
                                    tasmax.plots[[1]], tasmax.plots[[31]], tasmax.plots[[62]], ncol = 3), legGrob, ncol = 2, widths = c(10, 1))


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
