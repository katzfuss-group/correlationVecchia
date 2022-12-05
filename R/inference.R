####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes R functions to estimate parameters using the Fisher-scoring approach.
###
###   Contents:
###     fit_scaled / fit_scaled_spacetime_range
###
####################################################################################

#' @title Fitting parameters using scaled Vecchia, assuming matern covariance (just copied from scaledVecchia)
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
#'
#' @return Object containing fit information
#' @export
#'
#' @examples
#' \dontrun{
#' inputs=matrix(runif(40),ncol=2)
#' y=sin(rowSums(inputs*5))
#' fit=fit_scaled(y,inputs)
#' summary.GpGp_fit(fit)
#' }
fit_scaled <- function(y, inputs, ms = c(30), trend = 'zero', X, nu = 3.5, nug = 0, scale = 'parms', var.ini, ranges.ini, select = Inf, print.level = 0, max.it = 32, tol.dec = 4, n.est = min(5e3, nrow(inputs)))
{
  ## dimensions
  n = nrow(inputs)
  d = ncol(inputs)

  ## specify trend covariates
  if(missing(X)) {
    if(trend == 'zero'){
      X = as.matrix(sample(c(-1, 1), n, replace = TRUE))
    } else if(trend == 'intercept'){
      X = as.matrix(rep(1, n))
    } else if(trend == 'linear'){
      X = cbind(rep(1, n), inputs)
    } else stop('invalid trend option specified')
  } else trend = 'X'

  ## default variance parameter
  if(missing(var.ini)) {
    cur.var = summary(stats::lm(y~X-1))$sigma^2
  } else cur.var = var.ini

  ## default range parameters
  input.ranges = apply(inputs, 2, function(x) diff(range(x)))
  if(missing(ranges.ini)) cur.ranges = .2 * input.ranges else cur.ranges = ranges.ini
  active = rep(TRUE, d)

  ## fixed nugget?
  if(is.null(nug)){
    fix.nug = FALSE; nug = .01 * var(y)
  } else fix.nug = TRUE

  ## smoothness: fixed? bessel?
  if(is.null(nu)){
    covfun = 'matern_scaledim'
    cur.oth = c(3.5, nug)
    fix.nu = FALSE
  } else if(nu %in% (.5+(1:4))){
    covfun = paste0("matern", nu * 10, "_scaledim")
    cur.oth = nug
    fix.nu = FALSE
  } else {
    covfun = 'matern_scaledim'
    cur.oth = c(nu, nug)
    fix.nu = TRUE
  }

  ## only use subsample for estimation?
  if(n.est < n){
    ind.est = sample(1:n, n.est)
    y.full = y; inputs.full = inputs; X.full = X
    y = y[ind.est]; inputs = inputs[ind.est, , drop=FALSE]; X = X[ind.est, , drop=FALSE]
  }

  ## decrease or remove m values larger than n
  ms = unique(ifelse(ms < length(y), ms, length(y)-1))

  ### for increasing m
  for(i.m in 1:length(ms)){

    m = ms[i.m]
    if(i.m < length(ms)){ tol = 10^(-tol.dec-2) } else {tol = 10^(-tol.dec)}

    ### increase maxit until convergence
    conv = FALSE
    maxit = 2
    while(conv == FALSE & maxit <= max.it){

      if(print.level > 0) {
        print(paste0('m = ', m, ', maxit = ', maxit)); print(cur.ranges)}

      ## check for inactive input dims (large range params)
      active = (cur.ranges < input.ranges * select)
      if(sum(active, na.rm = TRUE) == 0) stop('all inputs inactive. increase select?')
      cur.ranges[!active] = Inf

      ## specify how to scale input dimensions
      cur.ranges[!active] = Inf

      ## order and condition based on current params
      if(scale == 'parms'){ scales = 1/cur.ranges
      } else if(scale == 'ranges'){ scales = 1/input.ranges
      } else if(scale == 'none'){ scales = 1
      } else stop(paste0('invalid argument scale = ', scale))

      ## order and condition based on current params
      ord = GPvecchia::order_maxmin_exact(t(t(inputs)*scales))
      inputs.ord = inputs[ord, , drop = FALSE]
      y.ord = y[ord]
      X.ord = X[ord, , drop = FALSE]
      NNarray = GpGp::find_ordered_nn(t(t(inputs.ord)*scales), m)

      ## starting and fixed parameters
      cur.parms = c(cur.var, cur.ranges[active], cur.oth)
      fixed = NULL
      if(fix.nu) fixed = c(fixed, length(cur.parms) - 1)
      if(fix.nug) fixed = c(fixed, length(cur.parms))

      ## fisher scoring
      fit = GpGp::fit_model(y.ord, inputs.ord[, active], X.ord, NNarray = NNarray, m_seq = m, convtol = tol, start_parms = cur.parms, max_iter = maxit, covfun_name = covfun, silent = (print.level < 2), reorder = FALSE, fixed_parms = fixed)
      cur.var = fit$covparms[1]
      cur.ranges[active] = fit$covparms[1+(1:sum(active))]
      cur.oth = fit$covparms[-(1:(1+sum(active)))]
      conv = fit$conv
      maxit = maxit * 2

    }
  }

  ### prepare fit object for subsequent prediction
  fit$covparms=c(cur.var, cur.ranges, cur.oth)
  fit$trend = trend
  if(n.est < n){
    fit$y = y.full
    fit$locs = inputs.full
    fit$X = X.full
  } else {
    fit$locs = inputs.ord
  }

  if(trend == 'zero') fit$X = as.matrix(rep(0, n))
  return(fit)
}

####################################################################################

#' @title Fitting range parameters using scaled Vecchia, assuming spacetime matern covariance (based on scaledVecchia)
#'
#' @param method b1, b2, b3, cc, or gp
#' @param y data vector of length n
#' @param inputs nxd matrix of input coordinates
#' @param ms vector of conditioning-set sizes
#' @param sig2 GP variance parameter
#' @param srange.ini initial value for spatial range parameter
#' @param trange.ini initial value for temporal range parameter
#' @param nu smoothness parameter. 1.5,2.5,3.5,4.5 avoid bessel (faster). estimated if nu=NULL
#' @param nug nugget or noise variance. estimated if nug=NULL
#' @param print.level 0: no printing. 1: print outer loop. 2: print outer+inner loop
#' @param max.it maximum number of iterations for inner loop
#' @param tol.dec converged if dot product between the step and the gradient is less than \code{10^(-convtol)}
#'
#' @return Object containing fit information, including for use in predictions_scaled()
#'
#' @export
#'
#' @examples
#' 1 + 1
fit_scaled_spacetime_range <- function(method, y, inputs, ms = NULL, sig2 = 1, srange.ini, trange.ini, nu = 0.5, nug = 0, print.level = 0, max.it = 100, tol.dec = 4)
{
  if( method %in% c("b1", "T-ord + T-NN") ) {

    if(is.null(ms)) stop("If the method is not the exact GP, the argument ms must be provided!")

    output <- .fit_scaled_spacetime_baseline_range(approx = 1, y = y, inputs = inputs, ms = ms, sig2 = sig2, srange.ini = srange.ini, trange.ini = trange.ini, nu = nu, nug = nug, print.level = print.level, max.it = max.it, tol.dec = tol.dec)

  } else if( method %in% c("b2", "T-ord + E-NN") ) {

    if(is.null(ms)) stop("If the method is not the exact GP, the argument ms must be provided!")

    output <- .fit_scaled_spacetime_baseline_range(approx = 2, y = y, inputs = inputs, ms = ms, sig2 = sig2, srange.ini = srange.ini, trange.ini = trange.ini, nu = nu, nug = nug, print.level = print.level, max.it = max.it, tol.dec = tol.dec)

  } else if( method %in% c("b3", "T-ord + C-NN") ) {

    if(is.null(ms)) stop("If the method is not the exact GP, the argument ms must be provided!")

    output <- .fit_scaled_spacetime_baseline_range(approx = 3, y = y, inputs = inputs, ms = ms, sig2 = sig2, srange.ini = srange.ini, trange.ini = trange.ini, nu = nu, nug = nug, print.level = print.level, max.it = max.it, tol.dec = tol.dec)

  } else if( method %in% c("cc", "C-MM + C-NN") ) {

    if(is.null(ms)) stop("If the method is not the exact GP, the argument ms must be provided!")

    output <- .fit_scaled_spacetime_cvecchia_range(y = y, inputs = inputs, ms = ms, sig2 = sig2, srange.ini = srange.ini, trange.ini = trange.ini, nu = nu, nug = nug, print.level = print.level, max.it = max.it, tol.dec = tol.dec)

  } else if( method %in% c("gp", "Exact GP") ) {

    output <- .fit_scaled_spacetime_exactgp_range(y = y, inputs = inputs, sig2 = sig2, srange.ini = srange.ini, trange.ini = trange.ini, nu = nu, nug = nug, print.level = print.level, max.it = max.it, tol.dec = tol.dec)

  } else {

    stop("Check the argument method!")
  }

  return(output)
}

.fit_scaled_spacetime_baseline_range <- function(approx, y, inputs, ms, sig2, srange.ini, trange.ini, nu, nug, print.level, max.it, tol.dec)
{
  ### dimensions
  n = nrow(inputs)
  d = ncol(inputs)

  # X = as.matrix(sample(c(-1, 1), n, replace = TRUE))
  X = as.matrix(rep(1, n))

  fix.nug     <- TRUE
  fix.nu      <- TRUE

  cur.var     <- sig2
  cur.ranges  <- c(srange.ini, trange.ini) ; cur.srange <- srange.ini ; cur.trange <- trange.ini
  cur.oth     <- c(nu, nug)

  ### decrease or remove m values larger than n
  ms = unique(ifelse(ms < length(y), ms, length(y)-1))

  ### for increasing m
  for(i.m in 1:length(ms)) {

    m = ms[i.m]
    if(i.m < length(ms)){ tol = 10^(-tol.dec-2) } else {tol = 10^(-tol.dec)}

    ### increase maxit until convergence
    conv = FALSE
    maxit = 2
    while(conv == FALSE & maxit <= max.it){

      if(print.level > 0) {
        print(paste0('m = ', m, ', maxit = ', maxit)); print(cur.ranges)
      }

      ## define scale parameters using range parameters
      scales = 1/c(rep(cur.srange, d-1), cur.trange)

      ## order and condition based on current params
      ord = order_time(locs = inputs, coordinate = NULL) # timely ordering
      inputs.ord = inputs[ord, , drop = FALSE]
      y.ord = y[ord]
      X.ord = X[ord, , drop = FALSE]

      if(approx == 1 | approx == "T-ord + T-NN") { # T-ord + T-NN
        NNarray = matrix(NA, nrow = n, ncol = m + 1)
        for(i in 1:n) {
          ind               <- seq(from = 1, to = min(i, m + 1), by = 1)
          NNarray[i, ind]   <- seq(from = i, by = -1, length.out = length(ind))
        }
      } else if(approx == 2 | approx == "T-ord + E-NN") { # T-ord + E-NN
        NNarray = GpGp::find_ordered_nn(inputs.ord, m)
      } else if(approx == 3 | approx == "T-ord + C-NN") { # T-ord + C-NN
        NNarray = GpGp::find_ordered_nn(t(t(inputs.ord)*scales), m)
      } else {
        stop("The argument baseline must be 1, 2, or 3.")
      }

      ## starting and fixed parameters
      cur.parms     <- c(sig2, cur.ranges, cur.oth)
      ind.parms     <- c(rep("var", length(sig2)), rep("range", length(c(cur.ranges))), rep("oth", length(cur.oth)))

      fixed         <- c(1, 4, 5)

      ## fisher scoring
      fit = GpGp::fit_model(y.ord, inputs.ord, X.ord, NNarray = NNarray, m_seq = m, convtol = tol, start_parms = cur.parms, max_iter = maxit, covfun_name = "matern_spacetime", silent = (print.level < 2), reorder = FALSE, fixed_parms = fixed)
      cur.var = fit$covparms[ind.parms == "var"]
      cur.ranges = fit$covparms[ind.parms == "range"] ; cur.srange <- cur.ranges[1] ; cur.trange <- cur.ranges[2]
      cur.oth = fit$covparms[ind.parms == "oth"]
      conv = fit$conv
      maxit = maxit * 2
    }
  }

  ### prepare fit object for subsequent prediction
  fit$covparms=c(cur.var, cur.ranges, cur.oth)
  fit$locs = inputs.ord
  return(fit)
}

.fit_scaled_spacetime_cvecchia_range <- function(y, inputs, ms, sig2, srange.ini, trange.ini, nu, nug, print.level, max.it, tol.dec)
{
  ### dimensions
  n = nrow(inputs)
  d = ncol(inputs)

  # X = as.matrix(sample(c(-1, 1), n, replace = TRUE))
  X = as.matrix(rep(1, n))

  fix.nug     <- TRUE
  fix.nu      <- TRUE

  cur.var     <- sig2
  cur.ranges  <- c(srange.ini, trange.ini) ; cur.srange <- srange.ini ; cur.trange <- trange.ini
  cur.oth     <- c(nu, nug)

  ### decrease or remove m values larger than n
  ms = unique(ifelse(ms < length(y), ms, length(y)-1))

  ### for increasing m
  for(i.m in 1:length(ms)) {

    m = ms[i.m]
    if(i.m < length(ms)){ tol = 10^(-tol.dec-2) } else {tol = 10^(-tol.dec)}

    ### increase maxit until convergence
    conv = FALSE
    maxit = 2
    while(conv == FALSE & maxit <= max.it){

      if(print.level > 0) {
        print(paste0('m = ', m, ', maxit = ', maxit)); print(cur.ranges)
      }

      ## define scale parameters using range parameters
      scales = 1/c(rep(cur.srange, d-1), cur.trange)

      ## order and condition based on current params
      ord = GPvecchia::order_maxmin_exact(t(t(inputs)*scales))
      inputs.ord = inputs[ord, , drop = FALSE]
      y.ord = y[ord]
      X.ord = X[ord, , drop = FALSE]
      NNarray = GpGp::find_ordered_nn(t(t(inputs.ord)*scales), m)

      ## starting and fixed parameters
      cur.parms     <- c(sig2, cur.ranges, cur.oth)
      ind.parms     <- c(rep("var", length(sig2)), rep("range", length(c(cur.ranges))), rep("oth", length(cur.oth)))

      fixed         <- c(1, 4, 5)

      ## fisher scoring
      fit = GpGp::fit_model(y = y.ord, locs = inputs.ord, X = X.ord, NNarray = NNarray, m_seq = m, convtol = tol, start_parms = cur.parms, max_iter = maxit, covfun_name = "matern_spacetime", silent = (print.level < 2), reorder = FALSE, fixed_parms = fixed)
      cur.var = fit$covparms[ind.parms == "var"]
      cur.ranges = fit$covparms[ind.parms == "range"] ; cur.srange <- cur.ranges[1] ; cur.trange <- cur.ranges[2]
      cur.oth = fit$covparms[ind.parms == "oth"]
      conv = fit$conv
      maxit = maxit * 2
    }
  }

  ### prepare fit object for subsequent prediction
  fit$covparms=c(cur.var, cur.ranges, cur.oth)
  fit$locs = inputs.ord
  return(fit)
}

.fit_scaled_spacetime_exactgp_range <- function(y, inputs, sig2, srange.ini, trange.ini, nu, nug, print.level, max.it, tol.dec)
{
  ### dimensions
  n = nrow(inputs)
  d = ncol(inputs)
  ms = n - 1

  # X = as.matrix(sample(c(-1, 1), n, replace = TRUE))
  X = as.matrix(rep(1, n))

  fix.nug     <- TRUE
  fix.nu      <- TRUE

  cur.var     <- sig2
  cur.ranges  <- c(srange.ini, trange.ini) ; cur.srange <- srange.ini ; cur.trange <- trange.ini
  cur.oth     <- c(nu, nug)

  ### decrease or remove m values larger than n
  ms = unique(ifelse(ms < length(y), ms, length(y)-1))

  ### for increasing m
  for(i.m in 1:length(ms)) {

    m = ms[i.m]
    if(i.m < length(ms)){ tol = 10^(-tol.dec-2) } else {tol = 10^(-tol.dec)}

    ### increase maxit until convergence
    conv = FALSE
    maxit = 2
    while(conv == FALSE & maxit <= max.it){

      if(print.level > 0) {
        print(paste0('m = ', m, ', maxit = ', maxit)); print(cur.ranges)
      }

      ## define scale parameters using range parameters
      scales = 1/c(rep(cur.srange, d-1), cur.trange)

      ## order and condition based on current params
      ord = order_time(locs = inputs, coordinate = NULL) # timely ordering
      inputs.ord = inputs[ord, , drop = FALSE]
      y.ord = y[ord]
      X.ord = X[ord, , drop = FALSE]

      NNarray = matrix(NA, nrow = n, ncol = m + 1)
      for(i in 1:n) {
        ind               <- seq(from = 1, to = min(i, m + 1), by = 1)
        NNarray[i, ind]   <- seq(from = i, by = -1, length.out = length(ind))
      }

      ## starting and fixed parameters
      cur.parms     <- c(sig2, cur.ranges, cur.oth)
      ind.parms     <- c(rep("var", length(sig2)), rep("range", length(c(cur.ranges))), rep("oth", length(cur.oth)))

      fixed         <- c(1, 4, 5)

      ## fisher scoring
      fit = GpGp::fit_model(y.ord, inputs.ord, X.ord, NNarray = NNarray, m_seq = m, convtol = tol, start_parms = cur.parms, max_iter = maxit, covfun_name = "matern_spacetime", silent = (print.level < 2), reorder = FALSE, fixed_parms = fixed)
      cur.var = fit$covparms[ind.parms == "var"]
      cur.ranges = fit$covparms[ind.parms == "range"] ; cur.srange <- cur.ranges[1] ; cur.trange <- cur.ranges[2]
      cur.oth = fit$covparms[ind.parms == "oth"]
      conv = fit$conv
      maxit = maxit * 2
    }
  }

  ### prepare fit object for subsequent prediction
  fit$covparms=c(cur.var, cur.ranges, cur.oth)
  fit$locs = inputs.ord
  return(fit)
}
