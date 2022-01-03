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

###################################################################################

.order_TIME <- function(locs)
{
  n         <- unlist(lapply(locs, nrow))
  p         <- length(locs)

  for(i in 1:p) rownames(locs[[i]]) <- rep(i, n[i])

  locsall   <- as.matrix(do.call(rbind, locs))
  ord.all   <- correlationVecchia::order_time(locsall)

  ord       <- list()
  for(i in 1:p) ord[[i]] <- rank(ord.all[rownames(locsall) == i])

  return( list(ord = ord, ord.all = ord.all) )
}

.order_SEMM <- function(locs)
{
  n         <- unlist(lapply(locs, nrow))
  n.cumsum  <- cumsum(n)
  p         <- length(locs)

  ord       <- list()
  for(i in 1:p) {
    ord[[i]]  <- GPvecchia::order_maxmin_exact( as.matrix(locs[[i]]) )

    if(i > 1) {
      ord.all   <- c(ord.all, ord[[i]] + n.cumsum[i-1])
    } else {
      ord.all   <- ord[[i]]
    }
  }

  return( list(ord = ord, ord.all = ord.all) )
}

.find_ordered_TNN <- function(locsord, m)
{
  cond.sets   <- matrix(NA, nrow = nrow(locsord), ncol = m + 1)
  for(i in 1:nrow(locsord)) {
    ind                 <- seq(from = 1, to = min(i, m + 1), by = 1)
    cond.sets[i, ind]   <- seq(from = i, by = -1, length.out = length(ind))
  }

  return(cond.sets)
}

.find_ordered_DENN <- function(locsord, m, n)
{
  p           <- length(n)
  n.cumsum    <- cumsum(n)

  cond.sets   <- matrix(NA, nrow = nrow(locsord), ncol = m + 1)
  cond.sets[seq(n[1]), seq(m+1)] <- GpGp::find_ordered_nn(locs = locsord[seq(n[1]), , drop = FALSE], m = m)
  for(i in 2:p) {
    for(j in seq(from = n.cumsum[[i-1]] + 1, to = n.cumsum[[i]], by = 1)) {
      history       <- locsord[seq(j), , drop = FALSE]
      distance      <- as.numeric(fields::rdist(x1 = locsord[j, , drop = FALSE], x2 = history))
      ord.distance  <- order(distance)

      siz                 <- min(m + 1, j)
      siz.each            <- rep(NA, times = i)
      siz.each[i]         <- min( ceiling((m + 1) / p), j - n.cumsum[[i-1]] )
      siz.each[seq(i-1)]  <- table(rep_len(seq(i-1), length.out = siz - siz.each[i]))

      cond.ind <- list()
      for(k in 1:i) {
        if(k == 1) {
          ind <- seq(n.cumsum[[1]])
        } else if(k == p){
          ind <- seq(from = n.cumsum[[k-1]] + 1, to = j, by = 1)
        } else {
          ind <- seq(from = n.cumsum[[k-1]] + 1, to = n.cumsum[[k]], by = 1)
        }

        cond.ind[[k]] <- order(distance[ind])[seq(siz.each[k])]
        if(k != 1) {
          cond.ind[[k]] <- cond.ind[[k]] + as.numeric(n.cumsum[k-1])
        }
      }

      cond.sets[j, ] <- unlist(cond.ind)[order(rank(distance[unlist(cond.ind)], ties.method = "last"))]
    }
  }

  return(cond.sets)
}

.find_ordered_JENN <- function(locsord, m)
{
  cond.sets   <- GpGp::find_ordered_nn(locs = locsord, m = m)

  return(cond.sets)
}

.find_ordered_SENN <- function(locsord, m, n)
{
  p         <- length(n)
  n.cumsum  <- cumsum(n)

  cond.each <- list()
  for(i in 1:p) {
    ind             <- seq(from = ifelse(i == 1, 1, n.cumsum[i-1] + 1), to = n.cumsum[i], by = 1)
    cond.each[[i]]  <- GpGp::find_ordered_nn(locs = locsord[ind, , drop = FALSE], m = m)

    if(i > 1) {
      cond.sets       <- rbind(cond.sets, cond.each[[i]] + n.cumsum[i-1])
    } else {
      cond.sets       <- cond.each[[i]]
    }
  }

  return(cond.sets)
}

.find_ordered_scaling <- function(locsord, m, scales = NULL)
{
  if(is.null(scales)) {

    cond.sets   <- GpGp::find_ordered_nn(locsord, m)

  } else {

    cond.sets   <- GpGp::find_ordered_nn(t(t(locsord)*scales), m)
  }

  return(cond.sets)
}

###################################################################################

#' @title Fitting parameters using baseline models for real-data analysis
#'
#' @param approx 1, 2, 3, 4, 5, 6, or 7
#' @param y data vector of length n
#' @param inputs nxd matrix of input coordinates
#' @param ms vector of conditioning-set sizes
#' @param trend options are 'zero' (no trend), 'intercept', 'linear' (incl intercept)
#' @param X nxp trend matrix (use if more complicated trend is desired)
#' @param nu smoothness parameter. 1.5,2.5,3.5,4.5 avoid bessel (faster). estimated if nu=NULL
#' @param nug nugget or noise variance. estimated if nug=NULL
#' @param scale scaling of inputs for ordering and conditioning. 'parms': by parameter estimates. 'ranges': to (0, 1). 'none': no scaling
#' @param var.ini initial value for GP variance parameter
#' @param ranges.ini initial values for d range parameters
#' @param select un-select input variables if estimated range parameter is above select (assuming standardized (0, 1) inputs)
#' @param print.level 0: no printing. 1: print outer loop. 2: print outer+inner loop
#' @param max.it maximum number of iterations for inner loop
#' @param tol.dec converged if dot product between the step and the gradient is less than \code{10^(-convtol)}
#' @param n.est subsample size for estimation
#' @param find.vcf logical
#' @param vcf.scorefun score function
#'
#' @return output of fit_scaled()
#'
#' @export
#'
#' @examples
#' 1 + 1
fit_scaled_bs_mulv <- function(approx,
                               y,inputs,ms=c(30),trend='pre',X,nu=3.5,nug=0,scale='parms',
                               var.ini,ranges.ini,select=Inf,print.level=0,max.it=32,tol.dec=4,
                               n.est=min(5e3,nrow(inputs)),find.vcf=FALSE,vcf.scorefun=lscore)
{
  ## inputs? (added)
  if(!is.list(inputs)) {

    message("The last column of the argument inputs is used to split inputs into groups.")

    inputs.lst  <- split(as.data.frame(inputs), inputs[, ncol(inputs)])
    inputs.lst  <- lapply(inputs.lst, function(x) as.matrix(x)[, -ncol(x), drop = FALSE])

  } else {

    inputs.lst  <- inputs
    inputs      <- as.matrix(do.call(rbind, inputs))
  }

  ## dimensions (modified)
  n=nrow(inputs)
  d=ncol(inputs)
  p=length(inputs.lst)

  ## specify trend covariates
  if(missing(X)) {
    if(trend=='zero'){
      X=as.matrix(sample(c(-1,1),n,replace=TRUE))
    } else if(trend=='intercept'){
      X=as.matrix(rep(1,n))
    } else if(trend=='linear'){
      X=cbind(rep(1,n),inputs)
    } else if(trend=='pre'){
      X=as.matrix(sample(c(-1,1),n,replace=TRUE))
      beta=mean(y)
      y=y-beta
    } else stop('invalid trend option specified')
  } else trend='X'

  ## default variance parameter
  if(missing(var.ini)) {
    cur.var=summary(stats::lm(y~X-1))$sigma^2
  } else cur.var=var.ini

  ## default range parameters
  input.ranges=apply(inputs,2,function(x) diff(range(x)))
  if(missing(ranges.ini)) cur.ranges=.2*input.ranges else cur.ranges=ranges.ini
  active=rep(TRUE,d)

  ## fixed nugget?
  if(is.null(nug)){
    fix.nug=FALSE; nug=.01*var(y)
  } else fix.nug=TRUE

  ## smoothness: fixed? bessel?
  if(is.null(nu)){
    covfun='matern_scaledim'
    cur.oth=c(3.5,nug)
    fix.nu=FALSE
  } else if(nu %in% (.5+(1:4))){
    covfun=paste0("matern",nu*10,"_scaledim")
    cur.oth=nug
    fix.nu=FALSE
  } else {
    covfun='matern_scaledim'
    cur.oth=c(nu,nug)
    fix.nu=TRUE
  }

  ## only use subsample for estimation? (modified)
  if(n.est<n){

    stop("n.est must be the same with n for now.")
  }

  ## decrease or remove m values larger than n
  ms=unique(ifelse(ms<length(y),ms,length(y)-1))


  ### for increasing m
  for(i.m in 1:length(ms)){

    m=ms[i.m]
    if(i.m<length(ms)){ tol=10^(-tol.dec-2) } else {tol=10^(-tol.dec)}

    ### increase maxit until convergence
    conv=FALSE
    maxit=2
    while(conv==FALSE & maxit<=max.it){

      if(print.level>0) {
        print(paste0('m=',m,', maxit=',maxit)); print(cur.ranges)}

      ## check for inactive input dims (large range params)
      active=(cur.ranges<input.ranges*select)
      if(sum(active,na.rm=TRUE)==0) stop('all inputs inactive. increase select?')
      cur.ranges[!active]=Inf

      ## specify how to scale input dimensions
      cur.ranges[!active]=Inf

      ## order and condition based on current params
      if(scale=='parms'){ scales=1/cur.ranges
      } else if(scale=='ranges'){ scales=1/input.ranges
      } else if(scale=='none'){ scales=1
      } else stop(paste0('invalid argument scale=',scale))

      ## order and condition based on current params (modified)
      if(approx == 1 | approx == "S-E-MM + D-E-NN") { # S-E-MM + D-E-NN

        ord.both    <- .order_SEMM(inputs.lst)

      } else if(approx == 2 | approx == "S-E-MM + J-E-NN") { # S-E-MM + J-E-NN

        ord.both    <- .order_SEMM(inputs.lst)

      } else if(approx == 3 | approx == "S-E-MM + S-E-NN") { # S-E-MM + S-E-NN

        ord.both    <- .order_SEMM(inputs.lst)

      } else if(approx == 4 | approx == "S-E-MM + C-NN") { # S-E-MM + C-NN

        ord.both    <- .order_SEMM(inputs.lst)

      } else if(approx == 5 | approx == "T-ord + T-NN") {

        ord.both    <- .order_TIME(inputs.lst)

      } else if(approx == 6 | approx == "T-ord + J-E-NN") {

        ord.both    <- .order_TIME(inputs.lst)

      } else if(approx == 7 | approx == "T-ord + C-NN") {

        ord.both    <- .order_TIME(inputs.lst)

      } else {

        stop("The argument approx must be 1, 2, 3, 4, 5, 6, or 7.")
      }

      inputs.ord  <- inputs[ord.both$ord.all , , drop = FALSE]

      for(i in seq(length(inputs.lst))) inputs.lst[[i]] <- inputs.lst[[i]][ord.both$ord[[i]], , drop = FALSE]

      y.ord       <- y[ord.both$ord.all]
      X.ord       <- X[ord.both$ord.all, , drop = FALSE]

      if(approx == 1 | approx == "S-E-MM + D-E-NN") { # S-E-MM + D-E-NN

        NNarray     <- .find_ordered_DENN(inputs.ord, m, unlist(lapply(inputs.lst, nrow)))

      } else if(approx == 2 | approx == "S-E-MM + J-E-NN") { # S-E-MM + J-E-NN

        NNarray     <- .find_ordered_JENN(inputs.ord, m)

      } else if(approx == 3 | approx == "S-E-MM + S-E-NN") { # S-E-MM + S-E-NN

        NNarray     <- .find_ordered_SENN(inputs.ord, m, unlist(lapply(inputs.lst, nrow)))

      } else if(approx == 4 | approx == "S-E-MM + C-NN") { # S-E-MM + C-NN

        NNarray     <- GpGp::find_ordered_nn(t(t(inputs.ord)*scales), m)

      } else if(approx == 5 | approx == "T-ord + T-NN") {

        NNarray     <- .find_ordered_TNN(inputs.ord, m)

      } else if(approx == 6 | approx == "T-ord + J-E-NN") {

        NNarray     <- .find_ordered_JENN(inputs.ord, m)

      } else if(approx == 7 | approx == "T-ord + C-NN") {

        NNarray     <- GpGp::find_ordered_nn(t(t(inputs.ord)*scales), m)

      } else {

        stop("The argument approx must be 1, 2, 3, 4, 5, 6, or 7.")
      }

      ## starting and fixed parameters
      cur.parms=c(cur.var,cur.ranges[active],cur.oth)
      fixed=NULL
      if(fix.nu) fixed=c(fixed,length(cur.parms)-1)
      if(fix.nug) fixed=c(fixed,length(cur.parms))

      ## fisher scoring
      fit=GpGp::fit_model(y.ord,inputs.ord[,active,drop=FALSE],X.ord,
                          NNarray=NNarray,m_seq=m,convtol=tol,
                          start_parms=cur.parms,max_iter=maxit,
                          covfun_name=covfun,silent=(print.level<2),
                          reorder=FALSE,fixed_parms=fixed)
      cur.var=fit$covparms[1]
      cur.ranges[active]=fit$covparms[1+(1:sum(active))]
      cur.oth=fit$covparms[-(1:(1+sum(active)))]
      conv=fit$conv
      maxit=maxit*2

    }
  }

  ### prepare fit object for subsequent prediction
  fit$covparms=c(cur.var,cur.ranges,cur.oth)
  fit$trend=trend
  if(n.est<n){
    fit$y=y.full
    fit$locs=inputs.full
    fit$X=X.full
  } else {
    fit$locs=inputs.ord
  }
  if(trend=='zero') {
    fit$X=as.matrix(rep(0,n))
  } else if(trend=='pre') {
    fit$betahat=beta
    fit$y=fit$y+beta
    fit$trend='intercept'
    fit$X=as.matrix(rep(1,n))
  }

  ### find variance correction factor, if requested
  if(find.vcf){
    fit$vcf=fit_vcf(fit,scale=scale,scorefun=vcf.scorefun)
  } else fit$vcf=1

  return(fit)
}

#' @title Fitting parameters using cvecchia for real-data analysis
#'
#' @param y data vector of length n
#' @param inputs nxd matrix of input coordinates
#' @param ms vector of conditioning-set sizes
#' @param trend options are 'zero' (no trend), 'intercept', 'linear' (incl intercept)
#' @param X nxp trend matrix (use if more complicated trend is desired)
#' @param nu smoothness parameter. 1.5,2.5,3.5,4.5 avoid bessel (faster). estimated if nu=NULL
#' @param nug nugget or noise variance. estimated if nug=NULL
#' @param scale scaling of inputs for ordering and conditioning. 'parms': by parameter estimates. 'ranges': to (0, 1). 'none': no scaling
#' @param var.ini initial value for GP variance parameter
#' @param ranges.ini initial values for d range parameters
#' @param select un-select input variables if estimated range parameter is above select (assuming standardized (0, 1) inputs)
#' @param print.level 0: no printing. 1: print outer loop. 2: print outer+inner loop
#' @param max.it maximum number of iterations for inner loop
#' @param tol.dec converged if dot product between the step and the gradient is less than \code{10^(-convtol)}
#' @param n.est subsample size for estimation
#' @param find.vcf logical
#' @param vcf.scorefun score function
#'
#' @return output of fit_scaled()
#'
#' @export
#'
#' @examples
#' 1 + 1
fit_scaled_cv_mulv <- function(y,inputs,ms=c(30),trend='pre',X,nu=3.5,nug=0,scale='parms',
                               var.ini,ranges.ini,select=Inf,print.level=0,max.it=32,tol.dec=4,
                               n.est=min(5e3,nrow(inputs)),find.vcf=FALSE,vcf.scorefun=lscore)
{
  ## inputs? (added)
  if(!is.list(inputs)) {

    message("The last column of the argument inputs is used to split inputs into groups.")

    inputs.lst  <- split(as.data.frame(inputs), inputs[, ncol(inputs)])
    inputs.lst  <- lapply(inputs.lst, function(x) as.matrix(x)[, -ncol(x), drop = FALSE])

  } else {

    inputs.lst  <- inputs
    inputs      <- as.matrix(do.call(rbind, inputs))
  }

  ## dimensions (modified)
  n=nrow(inputs)
  d=ncol(inputs)
  p=length(inputs.lst)

  ## specify trend covariates
  if(missing(X)) {
    if(trend=='zero'){
      X=as.matrix(sample(c(-1,1),n,replace=TRUE))
    } else if(trend=='intercept'){
      X=as.matrix(rep(1,n))
    } else if(trend=='linear'){
      X=cbind(rep(1,n),inputs)
    } else if(trend=='pre'){
      X=as.matrix(sample(c(-1,1),n,replace=TRUE))
      beta=mean(y)
      y=y-beta
    } else stop('invalid trend option specified')
  } else trend='X'

  ## default variance parameter
  if(missing(var.ini)) {
    cur.var=summary(stats::lm(y~X-1))$sigma^2
  } else cur.var=var.ini

  ## default range parameters
  input.ranges=apply(inputs,2,function(x) diff(range(x)))
  if(missing(ranges.ini)) cur.ranges=.2*input.ranges else cur.ranges=ranges.ini
  active=rep(TRUE,d)

  ## fixed nugget?
  if(is.null(nug)){
    fix.nug=FALSE; nug=.01*var(y)
  } else fix.nug=TRUE

  ## smoothness: fixed? bessel?
  if(is.null(nu)){
    covfun='matern_scaledim'
    cur.oth=c(3.5,nug)
    fix.nu=FALSE
  } else if(nu %in% (.5+(1:4))){
    covfun=paste0("matern",nu*10,"_scaledim")
    cur.oth=nug
    fix.nu=FALSE
  } else {
    covfun='matern_scaledim'
    cur.oth=c(nu,nug)
    fix.nu=TRUE
  }

  ## only use subsample for estimation? (modified)
  if(n.est<n){

    stop("n.est must be the same with n for now.")
  }

  ## decrease or remove m values larger than n
  ms=unique(ifelse(ms<length(y),ms,length(y)-1))


  ### for increasing m
  for(i.m in 1:length(ms)){

    m=ms[i.m]
    if(i.m<length(ms)){ tol=10^(-tol.dec-2) } else {tol=10^(-tol.dec)}

    ### increase maxit until convergence
    conv=FALSE
    maxit=2
    while(conv==FALSE & maxit<=max.it){

      if(print.level>0) {
        print(paste0('m=',m,', maxit=',maxit)); print(cur.ranges)}

      ## check for inactive input dims (large range params)
      active=(cur.ranges<input.ranges*select)
      if(sum(active,na.rm=TRUE)==0) stop('all inputs inactive. increase select?')
      cur.ranges[!active]=Inf

      ## specify how to scale input dimensions
      cur.ranges[!active]=Inf

      ## order and condition based on current params
      if(scale=='parms'){ scales=1/cur.ranges
      } else if(scale=='ranges'){ scales=1/input.ranges
      } else if(scale=='none'){ scales=1
      } else stop(paste0('invalid argument scale=',scale))

      ## order and condition based on current params
      ord=GPvecchia::order_maxmin_exact(t(t(inputs)*scales))
      inputs.ord=inputs[ord,,drop=FALSE]
      y.ord=y[ord]
      X.ord=X[ord,,drop=FALSE]
      NNarray=GpGp::find_ordered_nn(t(t(inputs.ord)*scales),m)

      ## starting and fixed parameters
      cur.parms=c(cur.var,cur.ranges[active],cur.oth)
      fixed=NULL
      if(fix.nu) fixed=c(fixed,length(cur.parms)-1)
      if(fix.nug) fixed=c(fixed,length(cur.parms))

      ## fisher scoring
      fit=GpGp::fit_model(y.ord,inputs.ord[,active,drop=FALSE],X.ord,
                          NNarray=NNarray,m_seq=m,convtol=tol,
                          start_parms=cur.parms,max_iter=maxit,
                          covfun_name=covfun,silent=(print.level<2),
                          reorder=FALSE,fixed_parms=fixed)
      cur.var=fit$covparms[1]
      cur.ranges[active]=fit$covparms[1+(1:sum(active))]
      cur.oth=fit$covparms[-(1:(1+sum(active)))]
      conv=fit$conv
      maxit=maxit*2

    }
  }


  ### prepare fit object for subsequent prediction
  fit$covparms=c(cur.var,cur.ranges,cur.oth)
  fit$trend=trend
  if(n.est<n){
    fit$y=y.full
    fit$locs=inputs.full
    fit$X=X.full
  } else {
    fit$locs=inputs.ord
  }
  if(trend=='zero') {
    fit$X=as.matrix(rep(0,n))
  } else if(trend=='pre') {
    fit$betahat=beta
    fit$y=fit$y+beta
    fit$trend='intercept'
    fit$X=as.matrix(rep(1,n))
  }

  ### find variance correction factor, if requested
  if(find.vcf){
    fit$vcf=fit_vcf(fit,scale=scale,scorefun=vcf.scorefun)
  } else fit$vcf=1

  return(fit)
}

###################################################################################

#' @title prediction using baseline models, using output from fit_scaled_bs_mulv()
#'
#' @param approx 2, 4, 5, 6, or 7
#' @param ns_obs ns_obs
#' @param ns_pred ns_pred
#' @param fit object returned from fit_scaled_bs_mulv()
#' @param locs_pred ns_pred x d matrix of test/prediction inputs/locations
#' @param m conditioning-set size (larger is more accurate but slower)
#' @param joint logical
#' @param nsims number of simulations
#' @param predvar logical
#' @param X_pred X_pred
#' @param scale scaling of inputs for ordering and conditioning. 'parms': by parameter estimates. 'ranges': to (0, 1). 'none': no scaling
#' @param tol numerical tolerance
#'
#' @return Vector of length n.p (\code{n.sims=0}, \code{predvar=FALSE}) or list with entries \code{means} and/or \code{vars} and/or \code{samples}
#'
#' @export
#'
#' @examples
#' 1 + 1
predictions_bs_mulv <- function(approx, ns_obs, ns_pred,
                                fit,locs_pred,m=100,joint=TRUE,nsims=0,
                                predvar=FALSE,X_pred,scale='parms',tol=1e-8)
{
  ## locs
  idx           <- rep(seq(length(ns_obs)) - 1, ns_obs)
  locs_obs      <- fit$locs
  locs_obs_lst  <- split(as.data.frame(locs_obs), idx)
  locs_obs_lst  <- lapply(locs_obs_lst, function(x) as.matrix(x))

  idx           <- rep(seq(length(ns_pred)) - 1, ns_pred)
  locs_pred_lst <- split(as.data.frame(locs_pred), idx)
  locs_pred_lst <- lapply(locs_pred_lst, function(x) as.matrix(x))

  n_obs         <- nrow(locs_obs)
  n_pred        <- nrow(locs_pred)

  ## basic
  y_obs = fit$y
  X_obs = fit$X
  beta = fit$betahat
  covparms = fit$covparms
  covfun_name = fit$covfun_name
  if(is.null(fit$vcf)) vcf=1 else vcf=fit$vcf

  # specify trend if missing
  if(missing(X_pred)){
    if(fit$trend=='zero'){
      X_pred=as.matrix(rep(0,n_pred))
    } else if(fit$trend=='intercept'){
      X_pred=as.matrix(rep(1,n_pred))
    } else if(fit$trend=='linear'){
      X_pred=cbind(rep(1,n_pred),locs_pred)
    } else stop('X_pred must be specified')
  }

  # specify how to scale input dimensions
  if(scale=='parms'){ scales=1/covparms[1+(1:ncol(locs_obs))]
  } else if(scale=='ranges'){ scales=1/apply(locs_obs,2,function(x) diff(range(x)))
  } else stop(paste0('invalid argument scale=',scale))

  ###
  if(joint){  # joint predictions

    # get orderings (modified)
    if(approx == 1 | approx == "S-E-MM + D-E-NN") { # S-E-MM + D-E-NN

      ord.both.1      <- .order_SEMM(locs_obs_lst)
      ord.both.2      <- .order_SEMM(locs_pred_lst)

    } else if(approx == 2 | approx == "S-E-MM + J-E-NN") { # S-E-MM + J-E-NN

      ord.both.1      <- .order_SEMM(locs_obs_lst)
      ord.both.2      <- .order_SEMM(locs_pred_lst)

    } else if(approx == 3 | approx == "S-E-MM + S-E-NN") { # S-E-MM + S-E-NN

      ord.both.1      <- .order_SEMM(locs_obs_lst)
      ord.both.2      <- .order_SEMM(locs_pred_lst)

    } else if(approx == 4 | approx == "S-E-MM + C-NN") { # S-E-MM + C-NN

      ord.both.1      <- .order_SEMM(locs_obs_lst)
      ord.both.2      <- .order_SEMM(locs_pred_lst)

    } else if(approx == 5 | approx == "T-ord + T-NN") {

      ord.both.1      <- .order_TIME(locs_obs_lst)
      ord.both.2      <- .order_TIME(locs_pred_lst)

    } else if(approx == 6 | approx == "T-ord + J-E-NN") {

      ord.both.1      <- .order_TIME(locs_obs_lst)
      ord.both.2      <- .order_TIME(locs_pred_lst)

    } else if(approx == 7 | approx == "T-ord + C-NN") {

      ord.both.1      <- .order_TIME(locs_obs_lst)
      ord.both.2      <- .order_TIME(locs_pred_lst)

    } else {

      stop("The argument approx must be 1, 2, 3, 4, 5, 6, or 7.")
    }

    ord1            <- ord.both.1[[2]]
    ord2            <- ord.both.2[[2]]

    # reorder stuff
    X_obs <- as.matrix(X_obs)
    X_pred <- as.matrix(X_pred)
    Xord_obs  <- X_obs[ord1,,drop=FALSE]
    Xord_pred <- X_pred[ord2,,drop=FALSE]
    y  <- y_obs[ord1] - Xord_obs %*% beta

    # put all coordinates together
    locs_all      <- rbind( locs_obs[ord1,,drop=FALSE], locs_pred[ord2,,drop=FALSE] )
    idx           <- rep(rep(seq(length(ns_obs)) - 1, length(ns_obs)), c(ns_obs, ns_pred))
    locs_all_lst  <- split(as.data.frame(locs_all), idx)
    locs_all_lst  <- lapply(locs_all_lst, function(x) as.matrix(x))

    # get nearest neighbor array (modified)
    sm = if (n_pred<1e5) 2 else 1.5

    if(approx == 1 | approx == "S-E-MM + D-E-NN") { # S-E-MM + D-E-NN

      NNarray     <- 1 ; stop("NOT YET!")

    } else if(approx == 2 | approx == "S-E-MM + J-E-NN") { # S-E-MM + J-E-NN

      NNarray     <- find_ordered_nn_pred(locs_all, m, fix.first = n_obs, searchmult = sm)

    } else if(approx == 3 | approx == "S-E-MM + S-E-NN") { # S-E-MM + S-E-NN

      NNarray     <- 1 ; stop("NOT YET!")

    } else if(approx == 4 | approx == "S-E-MM + C-NN") { # S-E-MM + C-NN

      NNarray     <- find_ordered_nn_pred(t(t(locs_all) * scales), m, fix.first = n_obs, searchmult = sm)

    } else if(approx == 5 | approx == "T-ord + T-NN") {

      NNarray     <- find_ordered_nn_pred(locs_all[, ncol(locs_all), drop = FALSE], m, fix.first = n_obs, searchmult = sm)

    } else if(approx == 6 | approx == "T-ord + J-E-NN") {

      NNarray     <- find_ordered_nn_pred(locs_all, m, fix.first = n_obs, searchmult = sm)

    } else if(approx == 7 | approx == "T-ord + C-NN") {

      NNarray     <- find_ordered_nn_pred(t(t(locs_all) * scales), m, fix.first = n_obs, searchmult = sm)

    } else {

      stop("The argument approx must be 1, 2, 3, 4, 5, 6, or 7.")
    }

    NNarray_pred=NNarray_all[-(1:n_obs),-1]

    means=numeric(length=n_pred)
    if(nsims>0) samples=array(dim=c(n_pred,nsims))

    # make predictions sequentially
    for(i in 1:n_pred){

      # NN conditioning sets
      NN=sort(NNarray_pred[i,])
      NN_obs=NN[NN<=n_obs]
      NN_pred=NN[NN>n_obs]-n_obs

      # (co-)variances (modified)
      K=get(covfun_name)(covparms, as.matrix(locs_all[c(NN,i+n_obs),]))
      cl=t(chol(K + tol * diag(nrow(K))))

      # prediction
      y.NN=y[NN_obs]
      means[i]=cl[m+1,1:m]%*%forwardsolve(cl[1:m,1:m],c(y.NN,means[NN_pred]))
      if(nsims>0){ # conditional simulation
        pred.var=cl[m+1,m+1]^2*vcf
        for(s in 1:nsims){
          pm=cl[m+1,1:m]%*%forwardsolve(cl[1:m,1:m],c(y.NN,samples[NN_pred,s]))
          samples[i,s]=stats::rnorm(1,pm,sqrt(pred.var))
        }
      }

    }

    # add (prior) mean and return to original ordering
    means[ord2] = means + c(Xord_pred %*% beta)
    if(nsims==0){
      preds=means
    } else {
      samples[ord2,] = samples + c(Xord_pred %*% beta)
      preds=list(means=means,samples=samples)
    }

  } else {  # separate predictions

    if(nsims>0) stop('cannot produce joint samples when joint=FALSE')

    y  = y_obs - X_obs %*% beta

    # find the NNs (modified)
    m=min(m,nrow(locs_obs))

    if(approx == 1 | approx == "S-E-MM + D-E-NN") { # S-E-MM + D-E-NN

      NNarray     <- 1 ; stop("NOT YET!")

    } else if(approx == 2 | approx == "S-E-MM + J-E-NN") { # S-E-MM + J-E-NN

      NNarray     <- FNN::get.knnx(locs_obs, locs_pred, m)$nn.index

    } else if(approx == 3 | approx == "S-E-MM + S-E-NN") { # S-E-MM + S-E-NN

      NNarray     <- 1 ; stop("NOT YET!")

    } else if(approx == 4 | approx == "S-E-MM + C-NN") { # S-E-MM + C-NN

      NNarray     <- FNN::get.knnx(t(t(locs_obs)*scales), t(t(locs_pred)*scales), m)$nn.index

    } else if(approx == 5 | approx == "T-ord + T-NN") {

      NNarray     <- FNN::get.knnx(locs_obs[, ncol(locs_obs), drop = FALSE], locs_pred[, ncol(locs_pred), drop = FALSE], m)$nn.index

    } else if(approx == 6 | approx == "T-ord + J-E-NN") {

      NNarray     <- FNN::get.knnx(locs_obs, locs_pred, m)$nn.index

    } else if(approx == 7 | approx == "T-ord + C-NN") {

      NNarray     <- FNN::get.knnx(t(t(locs_obs)*scales), t(t(locs_pred)*scales), m)$nn.index

    } else {

      stop("The argument approx must be 1, 2, 3, 4, 5, 6, or 7.")
    }

    means=vars=numeric(length=n_pred)
    for(i in 1:n_pred){

      # NN conditioning sets
      NN=NNarray[i,]

      # (co-)variances (modified)
      K=get(covfun_name)(covparms, as.matrix(rbind(locs_obs[NN,],locs_pred[i,])))
      cl=t(chol(K + tol * diag(nrow(K))))

      # prediction
      means[i]=cl[m+1,1:m]%*%forwardsolve(cl[1:m,1:m],y[NN])
      vars[i]=cl[m+1,m+1]^2*vcf

    }
    means=means+c(X_pred %*% beta)

    if(predvar==FALSE){
      preds=means
    } else {
      preds=list(means=means,vars=vars)
    }

  }

  return(preds)
}

#' @title prediction using cvecchia, using output from fit_scaled_cv_mulv()
#'
#' @param ns_obs ns_obs
#' @param ns_pred ns_pred
#' @param fit object returned from fit_scaled_bs_mulv()
#' @param locs_pred ns_pred x d matrix of test/prediction inputs/locations
#' @param m conditioning-set size (larger is more accurate but slower)
#' @param joint logical
#' @param nsims number of simulations
#' @param predvar logical
#' @param X_pred X_pred
#' @param scale scaling of inputs for ordering and conditioning. 'parms': by parameter estimates. 'ranges': to (0, 1). 'none': no scaling
#' @param tol numerical tolerance
#'
#' @return Vector of length n.p (\code{n.sims=0}, \code{predvar=FALSE}) or list with entries \code{means} and/or \code{vars} and/or \code{samples}
#'
#' @export
#'
#' @examples
#' 1 + 1
predictions_cv_mulv <- function(ns_obs, ns_pred,
                                fit,locs_pred,m=100,joint=TRUE,nsims=0,
                                predvar=FALSE,X_pred,scale='parms',tol=1e-8)
{
  ## locs
  idx           <- rep(seq(length(ns_obs)) - 1, ns_obs)
  locs_obs      <- fit$locs
  locs_obs_lst  <- split(as.data.frame(locs_obs), idx)
  locs_obs_lst  <- lapply(locs_obs_lst, function(x) as.matrix(x))

  idx           <- rep(seq(length(ns_pred)) - 1, ns_pred)
  locs_pred_lst <- split(as.data.frame(locs_pred), idx)
  locs_pred_lst <- lapply(locs_pred_lst, function(x) as.matrix(x))

  n_obs         <- nrow(locs_obs)
  n_pred        <- nrow(locs_pred)

  ## basic
  y_obs = fit$y
  X_obs = fit$X
  beta = fit$betahat
  covparms = fit$covparms
  covfun_name = fit$covfun_name
  if(is.null(fit$vcf)) vcf=1 else vcf=fit$vcf

  # specify trend if missing
  if(missing(X_pred)){
    if(fit$trend=='zero'){
      X_pred=as.matrix(rep(0,n_pred))
    } else if(fit$trend=='intercept'){
      X_pred=as.matrix(rep(1,n_pred))
    } else if(fit$trend=='linear'){
      X_pred=cbind(rep(1,n_pred),locs_pred)
    } else stop('X_pred must be specified')
  }

  # specify how to scale input dimensions
  if(scale=='parms'){ scales=1/covparms[1+(1:ncol(locs_obs))]
  } else if(scale=='ranges'){ scales=1/apply(locs_obs,2,function(x) diff(range(x)))
  } else stop(paste0('invalid argument scale=',scale))

  ###
  if(joint){  # joint predictions

    # get orderings (modified)
    temp=order_maxmin_pred(t(t(locs_obs)*scales),t(t(locs_pred)*scales))
    ord1=temp$ord
    ord2=temp$ord_pred

    # reorder stuff
    X_obs <- as.matrix(X_obs)
    X_pred <- as.matrix(X_pred)
    Xord_obs  <- X_obs[ord1,,drop=FALSE]
    Xord_pred <- X_pred[ord2,,drop=FALSE]
    y  <- y_obs[ord1] - Xord_obs %*% beta

    # put all coordinates together
    locs_all      <- rbind( locs_obs[ord1,,drop=FALSE], locs_pred[ord2,,drop=FALSE] )
    idx           <- rep(rep(seq(length(ns_obs)) - 1, length(ns_obs)), c(ns_obs, ns_pred))
    locs_all_lst  <- split(as.data.frame(locs_all), idx)
    locs_all_lst  <- lapply(locs_all_lst, function(x) as.matrix(x))

    # get nearest neighbor array (modified)
    sm = if (n_pred<1e5) 2 else 1.5
    NNarray_all <- find_ordered_nn_pred(t(t(locs_all)*scales),m,
                                        fix.first=n_obs,searchmult=sm)
    NNarray_pred=NNarray_all[-(1:n_obs),-1]

    means=numeric(length=n_pred)
    if(nsims>0) samples=array(dim=c(n_pred,nsims))

    # make predictions sequentially
    for(i in 1:n_pred){

      # NN conditioning sets
      NN=sort(NNarray_pred[i,])
      NN_obs=NN[NN<=n_obs]
      NN_pred=NN[NN>n_obs]-n_obs

      # (co-)variances (modified)
      K=get(covfun_name)(covparms, as.matrix(locs_all[c(NN,i+n_obs),]))
      cl=t(chol(K + tol * diag(nrow(K))))

      # prediction
      y.NN=y[NN_obs]
      means[i]=cl[m+1,1:m]%*%forwardsolve(cl[1:m,1:m],c(y.NN,means[NN_pred]))
      if(nsims>0){ # conditional simulation
        pred.var=cl[m+1,m+1]^2*vcf
        for(s in 1:nsims){
          pm=cl[m+1,1:m]%*%forwardsolve(cl[1:m,1:m],c(y.NN,samples[NN_pred,s]))
          samples[i,s]=stats::rnorm(1,pm,sqrt(pred.var))
        }
      }

    }

    # add (prior) mean and return to original ordering
    means[ord2] = means + c(Xord_pred %*% beta)
    if(nsims==0){
      preds=means
    } else {
      samples[ord2,] = samples + c(Xord_pred %*% beta)
      preds=list(means=means,samples=samples)
    }

  } else {  # separate predictions

    if(nsims>0) stop('cannot produce joint samples when joint=FALSE')

    y  = y_obs - X_obs %*% beta

    # find the NNs (modified)
    m=min(m,nrow(locs_obs))
    NNarray=FNN::get.knnx(t(t(locs_obs)*scales),
                          t(t(locs_pred)*scales),m)$nn.index


    means=vars=numeric(length=n_pred)
    for(i in 1:n_pred){

      # NN conditioning sets
      NN=NNarray[i,]

      # (co-)variances (modified)
      K=get(covfun_name)(covparms, as.matrix(rbind(locs_obs[NN,],locs_pred[i,])))
      cl=t(chol(K + tol * diag(nrow(K))))

      # prediction
      means[i]=cl[m+1,1:m]%*%forwardsolve(cl[1:m,1:m],y[NN])
      vars[i]=cl[m+1,m+1]^2*vcf

    }
    means=means+c(X_pred %*% beta)

    if(predvar==FALSE){
      preds=means
    } else {
      preds=list(means=means,vars=vars)
    }

  }

  return(preds)
}
