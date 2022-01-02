####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes several R functions to perform real data analysis.
###
###   Contents:
###
####################################################################################

#' @title Project raster (rastername) to shape (shapename) using crs (crsname)
#'
#' @param rastername raster
#' @param tpoint timepoint
#' @param shapename shape
#' @param crsname crs
#' @param states states of interest
#'
#' @return list
#'
#' @export
#'
#' @examples
#' 1 + 1
mappingRaster <- function(rastername, tpoint, shapename, crsname = NULL, states = "all")
{
  ### Raster ###

  dat         <- raster(rastername, band = tpoint) # ; crs(dat) ; plot(dat) ; proj4string(dat)

  if( !is.null(crsname) ) {

    crs(dat)    <- crsname # ; crs(dat) ; plot(dat) ; dat ; image(dat)
  }

  ### Shape ###

  spdf        <- readOGR(dsn = getwd(), layer = shapename, verbose = TRUE) # ; names(spdf)

  if(states %in% c("all", "USA", "us", "full", "every")) {

    # CAUTION: hard-coded !!!
    spdf        <- subset(spdf, NAME %in% spdf$NAME[-c(14, 28, 38, 39, 43, 45, 46)])

  } else {

    spdf        <- subset(spdf, NAME %in% states)
  }

  spdf        <- spTransform(spdf, crs(dat)) # ; proj4string(my_spdf) ; plot(my_spdf)

  xp1         <- mask(crop(dat, extent(spdf)), spdf)
  xp2         <- crop(dat, extent(spdf))

  ### return ###

  return( list(crop1 = xp1, crop2 = xp2, data = dat, shape = spdf) )
}

#' @title Cleaning raster file
#'
#' @param rasterfile raster
#'
#' @return list
#'
#' @export
#'
#' @examples
#' 1 + 1
rasterCleaning <- function(rasterfile)
{
  rasterfile  <- rasterToPoints(rasterfile)

  min.x1      <- min(rasterfile[, 1])
  min.x2      <- min(rasterfile[, 2])
  max.x1      <- max(rasterfile[, 1])
  max.x2      <- max(rasterfile[, 2])

  maxlen      <- max(max.x1 - min.x1, max.x2 - min.x2)

  rasterfile[, 1] <- (rasterfile[, 1] - min.x1) / maxlen
  rasterfile[, 2] <- (rasterfile[, 2] - min.x2) / maxlen

  mean.y      <- mean(rasterfile[, 3])
  sd.y        <- sd(rasterfile[, 3])

  # rasterfile[, 3] <- (rasterfile[, 3] - mean.y) / sd.y

  return( list(dat = rasterfile, min.x1 = min.x1, min.x2 = min.x2, max.x1 = max.x1, max.x2 = max.x2, mean.y = mean.y, sd.y = sd.y) )
}

#' Spliting a dataset into training and test datasets
#'
#' @param df.joint dataset
#' @param method spacewise, timewise, or blockwise
#' @param num.locs number of locations
#' @param num.time number of time points
#' @param size.blck size of block
#'
#' @return list
#'
#' @export
#'
#' @examples
#' 1 + 1
splitData <- function(df.joint, method, num.locs = 12, num.time = 31, size.blck = 3^2)
{
  if(method %in% c("spacewise")) {

    locs.all        <- df.joint %>% distinct(x1, x2)
    n.space         <- locs.all %>% nrow() # = 852
    locs.all$idx    <- seq(n.space)
    df.joint        <- df.joint %>% left_join(locs.all, by = c("x1", "x2"))

    idx.test        <- sample(seq(n.space), as.integer(n.space / 3)) %>% sort()
    df.test         <- df.joint %>% filter(idx %in% idx.test) # 284 * 62 = 17608

    idx.train       <- seq(n.space)[-idx.test]
    df.train        <- df.joint %>% filter(idx %in% idx.train)

    df.joint        <- df.joint %>% as_tibble() %>% select(-idx)
    df.test         <- df.test %>% as_tibble() %>% select(-idx)
    df.train        <- df.train %>% as_tibble() %>% select(-idx)

  } else if(method %in% c("timewise")) {

    dts.all         <- df.joint %>% distinct(t)
    n.dts           <- dts.all %>% nrow()
    dts.all$idx     <- seq(n.dts)
    df.joint        <- df.joint %>% left_join(dts.all, by = c("t"))

    idx.test        <- seq(from = n.dts - as.integer(n.dts / 3) + 1, to = n.dts) # sample(seq(n.dts), as.integer(n.dts / 3)) %>% sort()
    df.test         <- df.joint %>% filter(idx %in% idx.test) # 852 * 2 * 10 = 17040

    idx.train       <- seq(n.dts)[-idx.test]
    df.train        <- df.joint %>% filter(idx %in% idx.train)

    df.joint        <- df.joint %>% select(-idx)
    df.test         <- df.test %>% select(-idx)
    df.train        <- df.train %>% select(-idx)

  } else if(method %in% c("blockwise")) {

    locs.all        <- df.joint %>% distinct(x1, x2)
    n.space         <- locs.all %>% nrow() # = 852
    distmat         <- dist(locs.all, method = "maximum") %>% as.matrix()
    # distmat[upper.tri(distmat)] <- NA

    distordmat      <- apply(distmat, 1, function(x) order(x, decreasing = FALSE))

    # par(mfrow = c(1, 2))
    #
    # plot(locs.all, cex = 0.5, xlab = "", ylab = "")
    #
    # pivot <- 400 ; size <- 3^2
    # plot(locs.all[-distordmat[1:size, pivot], ], xlim = range(locs.all[, 1]), ylim = range(locs.all[, 2]), xlab = "", ylab = "", cex = 0.5)
    # points(locs.all[distordmat[1, pivot], ], cex = 0.5, col = "red", pch = 16)
    # points(locs.all[distordmat[2:size, pivot], ], cex = 0.5, col = "green", pch = 16)
    #
    # pivot <- 200 ; size <- 3^2
    # plot(locs.all[-distordmat[1:size, pivot], ], xlim = range(locs.all[, 1]), ylim = range(locs.all[, 2]), xlab = "", ylab = "", cex = 0.5)
    # points(locs.all[distordmat[1, pivot], ], cex = 0.5, col = "red", pch = 16)
    # points(locs.all[distordmat[2:size, pivot], ], cex = 0.5, col = "green", pch = 16)
    #
    # pivot <- 1 ; size <- 3^2
    # plot(locs.all[-distordmat[1:size, pivot], ], xlim = range(locs.all[, 1]), ylim = range(locs.all[, 2]), xlab = "", ylab = "", cex = 0.5)
    # points(locs.all[distordmat[1, pivot], ], cex = 0.5, col = "red", pch = 16)
    # points(locs.all[distordmat[2:size, pivot], ], cex = 0.5, col = "green", pch = 16)
    #
    # pivot <- 454 ; size <- 3^2
    # plot(locs.all[-distordmat[1:size, pivot], ], xlim = range(locs.all[, 1]), ylim = range(locs.all[, 2]), xlab = "", ylab = "", cex = 0.5)
    # points(locs.all[distordmat[1, pivot], ], cex = 0.5, col = "red", pch = 16)
    # points(locs.all[distordmat[2:size, pivot], ], cex = 0.5, col = "green", pch = 16)
    #
    # pivot <- 852 ; size <- 3^2
    # plot(locs.all[-distordmat[1:size, pivot], ], xlim = range(locs.all[, 1]), ylim = range(locs.all[, 2]), xlab = "", ylab = "", cex = 0.5)
    # points(locs.all[distordmat[1, pivot], ], cex = 0.5, col = "red", pch = 16)
    # points(locs.all[distordmat[2:size, pivot], ], cex = 0.5, col = "green", pch = 16)
    #
    # par(mfrow = c(1, 1))

    dts.all         <- df.joint %>% distinct(t)
    n.dts           <- dts.all %>% nrow()

    locs.all$idx.s  <- seq(n.space)
    df.joint        <- df.joint %>% left_join(locs.all, by = c("x1", "x2"))

    dts.all$idx.t   <- seq(n.dts)
    df.joint        <- df.joint %>% left_join(dts.all, by = c("t"))

    ###

    idx.dts         <- dts.all$idx.t %>% sample(num.time) %>% sort()

    df.var1         <- df.joint %>% filter(d == 0)
    df.var1$split   <- "train"
    for(i in idx.dts) {

      idx.locs        <- locs.all$idx.s %>% sample(num.locs) %>% sort()
      idx.blck        <- c()
      for(j in idx.locs) idx.blck <- c(idx.blck, distordmat[seq(size.blck), j])

      # plot(locs.all[-idx.blck, 1:2], xlim = range(locs.all[, 1]), ylim = range(locs.all[, 2]), xlab = "", ylab = "", cex = 0.5)
      # for(j in idx.locs) points(locs.all[distordmat[seq(3^2), j], c("x1", "x2")], cex = 0.5, col = i + 1, pch = 16)

      if(i == 1) {

        df.var1$split[df.var1$idx.s %in% idx.blck & df.var1$idx.t %in% c(i, i + 1)] <- "test"

      } else if(i == n.dts) {

        df.var1$split[df.var1$idx.s %in% idx.blck & df.var1$idx.t %in% c(i - 1, i)] <- "test"

      } else {

        df.var1$split[df.var1$idx.s %in% idx.blck & df.var1$idx.t %in% c(i - 1, i, i + 1)] <- "test"

      }
    }

    df.var1.test    <- df.var1 %>% filter(split == "test")
    df.var1.train   <- df.var1 %>% filter(split == "train")

    # plot(df.var1.train[, c("x1", "x2")], xlim = range(locs.all[, 1]), ylim = range(locs.all[, 2]), xlab = "", ylab = "", cex = 0.5, col = "blue", pch = 16)
    # points(df.var1.test[, c("x1", "x2")], cex = 0.5, col = "red", pch = 16)

    ###

    idx.dts         <- dts.all$idx.t %>% sample(num.time) %>% sort()

    df.var2         <- df.joint %>% filter(d == 1)
    df.var2$split   <- "train"
    for(i in idx.dts) {

      idx.locs        <- locs.all$idx.s %>% sample(num.locs) %>% sort()
      idx.blck        <- c()
      for(j in idx.locs) idx.blck <- c(idx.blck, distordmat[seq(size.blck), j])

      # plot(locs.all[-idx.blck, 1:2], xlim = range(locs.all[, 1]), ylim = range(locs.all[, 2]), xlab = "", ylab = "", cex = 0.5)
      # for(j in idx.locs) points(locs.all[distordmat[seq(3^2), j], c("x1", "x2")], cex = 0.5, col = i + 1, pch = 16)

      if(i == 1) {

        df.var2$split[df.var2$idx.s %in% idx.blck & df.var2$idx.t %in% c(i, i + 1)] <- "test"

      } else if(i == n.dts) {

        df.var2$split[df.var2$idx.s %in% idx.blck & df.var2$idx.t %in% c(i - 1, i)] <- "test"

      } else {

        df.var2$split[df.var2$idx.s %in% idx.blck & df.var2$idx.t %in% c(i - 1, i, i + 1)] <- "test"

      }
    }

    df.var2.test    <- df.var2 %>% filter(split == "test")
    df.var2.train   <- df.var2 %>% filter(split == "train")

    # plot(df.var2.train[, c("x1", "x2")], xlim = range(locs.all[, 1]), ylim = range(locs.all[, 2]), xlab = "", ylab = "", cex = 0.5, col = "blue", pch = 16)
    # points(df.var2.test[, c("x1", "x2")], cex = 0.5, col = "red", pch = 16)

    ##

    df.joint        <- df.joint %>% select(-idx.s, -idx.t)
    df.test         <- rbind(df.var1.test, df.var2.test) %>% select(-idx.s, -idx.t, -split) # < 2 (n.process) * 7 (num.locs) * 5^2 (size.blck) * 7 (num.time) * 3 (length.time)
    df.train        <- rbind(df.var1.train, df.var2.train) %>% select(-idx.s, -idx.t, -split)

  }

  return(list(df.joint = df.joint, df.train = df.train, df.test = df.test))
}
