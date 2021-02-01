####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview:
###
###   Contents:
###
####################################################################################



#' Computing the approximate posterior of the range parameter for the mutivariate GP
#'
#' @param m m
#' @param approx b1, b2, b3, b4, or cc 
#' @param N N
#' @param limit extent of the search 
#' @param z observations the the GP
#' @param locs location matrix
#' @param covparms covariance parameters
#' @param nugget nugget (noise)
#' @param meanlog prior mean of log(range parameter)
#' @param sdlog prior standard deviation of log(range parameter)
#' @param logprior prior function
#'
#' @return list
#' 
#' @import RandomFields
#' @import mvtnorm
#' @import LaplacesDemon
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' n               <- 30^2
#' covparms        <- c(1, 0.1, 0.1)
#' process         <- generate_gp_space(nsim = 1, n = n, d = 2, p = 2, method.locs = 'random', covmodel = cov_bivariate_expo_latDim, covparms = covparms)
#' locs            <- rbind(process$sim$sim1$locs$locs1, process$sim$sim1$locs$locs2)
#' nugget          <- 0.00
#' 
#' z               <- process$sim$sim1$y + nugget * rnorm(2 * n)
#' par(mfrow = c(1, 2))
#' fields::quilt.plot(locs[, 1][1:n], locs[, 2][1:n], z[1:n], main = "process 1")
#' fields::quilt.plot(locs[, 1][1:n + n], locs[, 2][1:n + n], z[1:n + n], main = "process 2")
#' par(mfrow = c(1, 1))
#' 
#' # covmat          <- process$sim$sim1$covmat + nugget * diag(2 * n)
#' 
#' out1 <- posterior_multivariate(m = 10, N = 150, z = z, locs = process$sim$sim1$locs, covparms = covparms, nugget = nugget)
#' out2 <- posterior_multivariate(m = 30, N = 150, z = z, locs = process$sim$sim1$locs, covparms = covparms, nugget = nugget)
#' 
#' vis_posterior_multivariate(out1 = out1, out2 = out2)
#' }
posterior_multivariate <- function(m = 30, approx = "cc", N = 150, limit = c(0.088, 0.118), z, locs, covparms, nugget, meanlog = NULL, sdlog = NULL, logprior = NULL)
{
  ### prior information
  
  if(is.null(meanlog)) meanlog <- log(covparms[2])
  
  if(is.null(sdlog)) sdlog <- 0.6
  
  if(is.null(logprior)) logprior <- function(alpha) dlnorm(alpha, meanlog, sdlog, log = TRUE)
  
  ### posterior computation
  
  locsbind        <- rbind(locs$locs1, locs$locs2)
  alpha.grid      <- seq(limit[1], limit[2], length = N)
  logpost         <- array(dim = c(N, 2))
  for(i in 1:N) {
    
    print(i)
    alpha           <- alpha.grid[i]
    
    Sigma.true      <- cov_bivariate_expo_latDim(locs = locs, covparms = c(covparms[1], alpha, covparms[3])) + nugget * diag(2 * n)
    
    ## exact
    loglik          <- dmvnorm(z, mean = rep(0, 2 * n), sigma = Sigma.true, log = TRUE)
    logpost[i, 1]   <- loglik + logprior(alpha)
    
    ## Vecchia
    if(approx == "b1") {
      vecchia.approx  <- correlationVecchia::baseline_1_multivariate_specify(locs = locs, m = m)
    } else if(approx == "b2") {
      vecchia.approx  <- correlationVecchia::baseline_2_multivariate_specify(locs = locs, m = m)
    } else if(approx == "b3") {
      vecchia.approx  <- correlationVecchia::baseline_3_multivariate_specify(locs = locs, m = m)
    } else if(approx == "b4") {
      vecchia.approx  <- correlationVecchia::baseline_4_multivariate_specify(locs = locs, m = m, covmodel = cov_bivariate_expo_latDim, covparms = c(covparms[1], alpha, covparms[3]))
    } else {
      vecchia.approx  <- correlationVecchia::corrvecchia_specify_knownCovparms_2(locs = locsbind, m = m, covmodel = Sigma.true, covparms = c(covparms[1], alpha, covparms[3]))
    }
    
    Sigma.ord       <- Sigma.true[vecchia.approx$ord, vecchia.approx$ord]
    U               <- GPvecchia::createU(vecchia.approx, covparms = c(1), nuggets = nugget, covmodel = Sigma.ord)$U
    revord          <- order(vecchia.approx$ord)
    Sigma.vecchia   <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
    
    loglik          <- mvtnorm::dmvnorm(z, mean = rep(0, 2 * n), sigma = Sigma.vecchia, log = TRUE)
    logpost[i, 2]   <- loglik + logprior(alpha)
  }
  
  post            <- exp(logpost - mean(logpost))
  post.norm       <- post
  for(j in 1:2) post.norm[, j] <- post[, j] / sum(post[, j])
  
  ### return
  
  return(list(alpha = alpha.grid, post = post.norm))
}



#' Visualizing the approximate posteriors of the range parameter for the mutivariate GP
#'
#' @param out1 out1
#' @param out2 out2
#' @param color color
#' @param legend legend
#' @param ylimit ylimit
#' @param alpha alpha
#' @param size.point size.point
#' @param size.line size.line
#' @param size.legend size.legend
#' @param size.lab size.lab
#' @param size.text size.text
#' @param size.margin size.margin
#'
#' @return ggplot
#' 
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import RColorBrewer
#' @import gridExtra
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' n               <- 30^2
#' covparms        <- c(1, 0.1, 0.1)
#' process         <- generate_gp_space(nsim = 1, n = n, d = 2, p = 2, method.locs = 'random', covmodel = cov_bivariate_expo_latDim, covparms = covparms)
#' locs            <- rbind(process$sim$sim1$locs$locs1, process$sim$sim1$locs$locs2)
#' nugget          <- 0.00
#' 
#' z               <- process$sim$sim1$y + nugget * rnorm(2 * n)
#' par(mfrow = c(1, 2))
#' fields::quilt.plot(locs[, 1][1:n], locs[, 2][1:n], z[1:n], main = "process 1")
#' fields::quilt.plot(locs[, 1][1:n + n], locs[, 2][1:n + n], z[1:n + n], main = "process 2")
#' par(mfrow = c(1, 1))
#' 
#' # covmat          <- process$sim$sim1$covmat + nugget * diag(2 * n)
#' 
#' out1 <- posterior_multivariate(m = 10, N = 150, z = z, locs = process$sim$sim1$locs, covparms = covparms, nugget = nugget)
#' out2 <- posterior_multivariate(m = 30, N = 150, z = z, locs = process$sim$sim1$locs, covparms = covparms, nugget = nugget)
#' 
#' vis_posterior_multivariate(out1 = out1, out2 = out2)
#' }
vis_posterior_multivariate <- function(out1, out2, color = c("#377EB8", "#E41A1C", "gray1"), legend = c("m = 10", "m = 30", "True density"), ylimit = c(0, 0.025), alpha = 0.7, size.point = 4, size.line = 1, size.legend = 14, size.lab = 14, size.text = 12, size.margin = c(5.5, 20, 5.5, 5.5), tol = 1e-8)
{
  if(length(out1$alpha) != sum(out1$alpha == out2$alpha)) stop("Check alpha!")
  
  if(sum((out1$post[, 1] - out2$post[, 1])^2) > tol) stop("Check true density!")
  
  output <- data.frame(alpha = out1$alpha, trueval = out1$post[, 1], approx1 = out1$post[, 2], approx2 = out2$post[, 2])
  output <- output %>% tidyr::gather(key = "method", value = "value", trueval, approx1, approx2)
  
  pplot <- ggplot(output, aes(x = alpha, y = value, group = method, colour = method)) + geom_line(lwd = 0.8) +
    xlab("range") + ylab('value') + 
    scale_color_manual(values = color, labels = legend, guide = guide_legend(nrow = 1, byrow = TRUE)) + 
    #    ylim(ylimit) +
    coord_cartesian(ylim = ylimit) + 
    theme(axis.title.x = element_text(size = size.lab), 
          axis.text.x = element_text(size = size.text), 
          axis.title.y = element_text(size = size.lab), 
          axis.text.y = element_text(size = size.text), 
          legend.title = element_blank(), 
          legend.text = element_text(size = size.legend, margin = margin(r = 25, unit = 'pt')), 
          legend.position="top",
          legend.direction = 'horizontal', 
          #        legend.spacing.x = unit(15, 'pt'),
          plot.margin = unit(size.margin, "pt")) # t, r, b, l
  
  return(pplot)
}


# n               <- 30^2
# covparms        <- c(1, 0.1, 0.1)
# process         <- generate_gp_space(nsim = 1, n = n, d = 2, p = 2, method.locs = 'random', covmodel = cov_bivariate_expo_latDim, covparms = covparms)
# locs            <- rbind(process$sim$sim1$locs$locs1, process$sim$sim1$locs$locs2)
# nugget          <- 0.00
# 
# z               <- process$sim$sim1$y + nugget * rnorm(2 * n)
# par(mfrow = c(1, 2))
# fields::quilt.plot(locs[, 1][1:n], locs[, 2][1:n], z[1:n], main = "process 1")
# fields::quilt.plot(locs[, 1][1:n + n], locs[, 2][1:n + n], z[1:n + n], main = "process 2")
# par(mfrow = c(1, 1))
# 
# # covmat          <- process$sim$sim1$covmat + nugget * diag(2 * n)
# 
# out1 <- posterior_multivariate(m = 10, N = 50, z = z, locs = process$sim$sim1$locs, covparms = covparms, nugget = nugget)
# out2 <- posterior_multivariate(m = 30, N = 50, z = z, locs = process$sim$sim1$locs, covparms = covparms, nugget = nugget)
# 
# vis_posterior_multivariate(out1 = out1, out2 = out2)



#' Computing the approximate posterior of the parameters for cov_bivariate_expo_latDim
#'
#' @param z process
#' @param locs locs
#' @param covparms covariance parameters
#' @param nugget nugget (noise)
#' @param m the size of conditioning sets
#' @param approx b1, b2, b3, b4, or cc
#' @param target range or distance
#' @param N the number of grid points
#' @param xlim extent of the search
#' @param sdlog prior standard deviation of log(parameter)
#'
#' @return list
#' 
#' @import RandomFields
#' @import mvtnorm
#' @import LaplacesDemon
#'  
#' @export 
#'
#' @examples
#' \dontrun{
#' n               <- 20^2
#' covparms        <- c(1, 0.1, 0.1)
#' process         <- generate_gp_space(nsim = 1, n = n, d = 2, p = 2, method.locs = 'random', covmodel = cov_bivariate_expo_latDim, covparms = covparms)
#' 
#' locs            <- process$sim$sim1$locs
#' locsall         <- rbind(locs$locs1, locs$locs2)
#' 
#' nugget          <- 0.0
#' z               <- process$sim$sim1$y + nugget * rnorm(2 * n)
#' 
#' par(mfrow = c(1, 2))
#' fields::quilt.plot(locsall, 1][1:n], locsall[, 2][1:n], z[1:n], main = "component 1")
#' fields::quilt.plot(locsall[, 1][1:n + n], locsall[, 2][1:n + n], z[1:n + n], main = "component 2")
#' par(mfrow = c(1, 1))
#' 
#' out.range       <- posterior_for_cov_bivariate_expo_latDim(z = z, locs = locs, covparms = covparms, target = "range", xlim = c(0.065, 0.14))
#' out.distance    <- posterior_for_cov_bivariate_expo_latDim(z = z, locs = locs, covparms = covparms, target = "distance", xlim = c(0.065, 0.14))
#' 
#' par(mfrow = c(1, 2))
#' matplot(out.range$alpha, out.range$post.norm, type = 'l', xlim = range(out.range$alpha), lty = 1, lwd = 2, xlab = 'range', ylab = 'density', main = "posterior distribution for range parameter")
#' legend('topright', c('exact', 'approx'), col = 1:2, lty = 1, lwd = 2)
#' 
#' matplot(out.distance$alpha, out.distance$post.norm, type = 'l', xlim = range(out.distance$alpha), lty = 1, lwd = 2, xlab = 'distance in latent dimension', ylab = 'density', main = "posterior distribution for distance in latent dimension")
#' legend('topright', c('exact', 'approx'), col = 1:2, lty = 1, lwd = 2)
#' par(mfrow = c(1, 1))
#' }
posterior_for_cov_bivariate_expo_latDim <- function(z, locs, covparms, nugget = 0, m = 30, approx = "cc", target = "range", N = 100, xlim = c(0.05, 0.15), sdlog = 0.6)
{
  if(target == "range") {
    
    output <- .posterior_for_range_expo_latDim(z = z, locs = locs, covparms = covparms, nugget = nugget, m = m, approx = approx, N = N, xlim = xlim, sdlog = sdlog)
    
  } else if(target == "distance" | target == "latentDistance") {
    
    output <- .posterior_for_distance_expo_latDim(z = z, locs = locs, covparms = covparms, nugget = nugget, m = m, approx = approx, N = N, xlim = xlim, sdlog = sdlog)
    
  } else {
    
    stop("Check the target!")
    
  }
  
  return(output)
}

.posterior_for_range_expo_latDim <- function(z, locs, covparms, nugget, m, approx, N, xlim, sdlog)
{
  ### prior information
  
  n <- nrow(locs[[1]])
  
  meanlog   <- log(covparms[2])
  logprior  <- function(alpha) dlnorm(alpha, meanlog, sdlog, log = TRUE)
  
  ### posterior computation
  
  locsbind        <- rbind(locs$locs1, locs$locs2)
  alpha.grid      <- seq(xlim[1], xlim[2], length = N)
  logpost         <- array(dim = c(N, 2))
  for(i in 1:N) {
    
    print(i)
    alpha           <- alpha.grid[i]
    
    Sigma.true      <- cov_bivariate_expo_latDim(locs = locs, covparms = c(covparms[1], alpha, covparms[3])) + nugget * diag(2 * n)
    
    ## exact
    loglik          <- dmvnorm(z, mean = rep(0, 2 * n), sigma = Sigma.true, log = TRUE)
    logpost[i, 1]   <- loglik + logprior(alpha)
    
    ## Vecchia
    if(approx == "b1") {
      vecchia.approx  <- correlationVecchia::baseline_1_multivariate_specify(locs = locs, m = m)
    } else if(approx == "b2") {
      vecchia.approx  <- correlationVecchia::baseline_2_multivariate_specify(locs = locs, m = m)
    } else if(approx == "b3") {
      vecchia.approx  <- correlationVecchia::baseline_3_multivariate_specify(locs = locs, m = m)
    } else if(approx == "b4") {
      vecchia.approx  <- correlationVecchia::baseline_4_multivariate_specify(locs = locs, m = m, covmodel = cov_bivariate_expo_latDim, covparms = c(covparms[1], alpha, covparms[3]))
    } else {
      vecchia.approx  <- correlationVecchia::corrvecchia_specify_knownCovparms_2(locs = locsbind, m = m, covmodel = Sigma.true, covparms = c(covparms[1], alpha, covparms[3]))
    }
    
    Sigma.ord       <- Sigma.true[vecchia.approx$ord, vecchia.approx$ord]
    U               <- GPvecchia::createU(vecchia.approx, covparms = c(1), nuggets = nugget, covmodel = Sigma.ord)$U
    revord          <- order(vecchia.approx$ord)
    Sigma.vecchia   <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
    
    loglik          <- mvtnorm::dmvnorm(z, mean = rep(0, 2 * n), sigma = Sigma.vecchia, log = TRUE)
    logpost[i, 2]   <- loglik + logprior(alpha)
  }
  
  post            <- logpost - mean(logpost)
  A               <- mean(post[, 1])
  B               <- mean(post[, 2])
  
  post.norm       <- post
  post.norm[, 1]  <- exp(post[, 1] - A) / sum(exp(post[, 1] - A))
  post.norm[, 2]  <- exp(post[, 2] - B) / sum(exp(post[, 2] - B))
  
  post            <- exp(post)
  
  # post            <- exp(logpost - mean(logpost))
  # post.norm       <- post
  # for(j in 1:2) post.norm[, j] <- post[, j] / sum(post[, j])
  
  return(list(alpha = alpha.grid, logpost = logpost, post = post, post.norm = post.norm))
}

.posterior_for_distance_expo_latDim <- function(z, locs, covparms, nugget, m, approx, N, xlim, sdlog)
{
  ### prior information
  
  n <- nrow(locs[[1]])
  
  meanlog   <- log(covparms[3])
  logprior  <- function(alpha) dlnorm(alpha, meanlog, sdlog, log = TRUE)
  
  ### posterior computation
  
  locsbind        <- rbind(locs$locs1, locs$locs2)
  alpha.grid      <- seq(xlim[1], xlim[2], length = N)
  logpost         <- array(dim = c(N, 2))
  for(i in 1:N) {
    
    print(i)
    alpha           <- alpha.grid[i]
    
    Sigma.true      <- cov_bivariate_expo_latDim(locs = locs, covparms = c(covparms[1], covparms[2], alpha)) + nugget * diag(2 * n)
    
    ## exact
    loglik          <- dmvnorm(z, mean = rep(0, 2 * n), sigma = Sigma.true, log = TRUE)
    logpost[i, 1]   <- loglik + logprior(alpha)
    
    ## Vecchia
    if(approx == "b1") {
      vecchia.approx  <- correlationVecchia::baseline_1_multivariate_specify(locs = locs, m = m)
    } else if(approx == "b2") {
      vecchia.approx  <- correlationVecchia::baseline_2_multivariate_specify(locs = locs, m = m)
    } else if(approx == "b3") {
      vecchia.approx  <- correlationVecchia::baseline_3_multivariate_specify(locs = locs, m = m)
    } else if(approx == "b4") {
      vecchia.approx  <- correlationVecchia::baseline_4_multivariate_specify(locs = locs, m = m, covmodel = cov_bivariate_expo_latDim, covparms = c(covparms[1], covparms[2], alpha))
    } else {
      vecchia.approx  <- correlationVecchia::corrvecchia_specify_knownCovparms_2(locs = locsbind, m = m, covmodel = Sigma.true, covparms = c(covparms[1], covparms[2], alpha))
    }
    
    Sigma.ord       <- Sigma.true[vecchia.approx$ord, vecchia.approx$ord]
    U               <- GPvecchia::createU(vecchia.approx, covparms = c(1), nuggets = nugget, covmodel = Sigma.ord)$U
    revord          <- order(vecchia.approx$ord)
    Sigma.vecchia   <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
    
    loglik          <- mvtnorm::dmvnorm(z, mean = rep(0, 2 * n), sigma = Sigma.vecchia, log = TRUE)
    logpost[i, 2]   <- loglik + logprior(alpha)
  }
  
  post            <- logpost - mean(logpost)
  A               <- mean(post[, 1])
  B               <- mean(post[, 2])
  
  post.norm       <- post
  post.norm[, 1]  <- exp(post[, 1] - A) / sum(exp(post[, 1] - A))
  post.norm[, 2]  <- exp(post[, 2] - B) / sum(exp(post[, 2] - B))
  
  post            <- exp(post)
  
  # post            <- exp(logpost - mean(logpost))
  # post.norm       <- post
  # for(j in 1:2) post.norm[, j] <- post[, j] / sum(post[, j])
  
  return(list(alpha = alpha.grid, logpost = logpost, post = post, post.norm = post.norm))
}



#' Computing the approximate posterior of the (spatial or temporal) range parameter for cov_matern_spacetime
#'
#' @param z process
#' @param locs locs
#' @param covparms covariance parameters
#' @param nugget nugget (noise)
#' @param m the size of conditioning sets
#' @param approx b1, b2, b3, or cc
#' @param target spatialRange or temporalRange
#' @param N the number of grid points
#' @param xlim extent of the search
#' @param sdlog prior standard deviation of log(range parameter)
#'
#' @return list
#' 
#' @import RandomFields
#' @import mvtnorm
#' @import LaplacesDemon
#'  
#' @export 
#'
#' @examples
#' \dontrun{
#'   n               <- 20^2
#'   covparms        <- c(1, 0.1, 0.2, 0.5)
#'   process         <- generate_gp_spacetime(nsim = 2, n = n, d = 2, t.len = 1, method.locs = "satellite", covmodel = cov_matern_spacetime, covparms = covparms)
#'   locs            <- process$sim$sim1$locs
#'   
#'   z               <- process$sim$sim1$y
#'   fields::quilt.plot(locs[, 1], locs[, 2], z, main = "process")
#'   
#'   covmat          <- process$sim$sim1$covmat
#'   
#'   
#'   par(mfrow = c(1, 2))
#'   
#'   meanlog <- log(covparms[2])
#'   sdlog <- 0.6
#'   logprior <- function(alpha) dlnorm(alpha, meanlog, sdlog, log = TRUE)
#'   
#'   plot(seq(0.01, 1, 0.01), exp(logprior(seq(0.01, 1, 0.01))), type = "l", xlab = "range parameter", ylab = "density", main = "prior distribution")
#'   
#'   meanlog <- log(covparms[3])
#'   sdlog <- 0.6
#'   logprior <- function(alpha) dlnorm(alpha, meanlog, sdlog, log = TRUE)
#'   
#'   plot(seq(0.01, 1, 0.01), exp(logprior(seq(0.01, 1, 0.01))), type = "l", xlab = "range parameter", ylab = "density", main = "prior distribution")
#'   
#'   par(mfrow = c(1, 1)) ; rm(meanlog, sdlog, logprior)
#'   
#'   
#'   out.spatialRange <- posterior_for_cov_matern_spacetime(z = z, locs = locs, covparms = covparms, target = "spatialRange", xlim = c(0.065, 0.14))
#'   out.temporalRange <- posterior_for_cov_matern_spacetime(z = z, locs = locs, covparms = covparms, target = "temporalRange", xlim = c(0.05, 0.75))
#'   
#'   par(mfrow = c(1, 2))
#'   matplot(out.spatialRange$alpha, out.spatialRange$post.norm, type = 'l', xlim = range(out.spatialRange$alpha), lty = 1, lwd = 2, xlab = 'range', ylab = 'density', main = "posterior distribution for spatial range parameter")
#'   legend('topright', c('exact', 'approx'), col = 1:2, lty = 1, lwd = 2)
#'   
#'   matplot(out.temporalRange$alpha, out.temporalRange$post.norm, type = 'l', xlim = range(out.temporalRange$alpha), lty = 1, lwd = 2, xlab = 'range', ylab = 'density', main = "posterior distribution for temporal range parameter")
#'   legend('topright', c('exact', 'approx'), col = 1:2, lty = 1, lwd = 2)
#'   par(mfrow = c(1, 1))
#' }
posterior_for_cov_matern_spacetime <- function(z, locs, covparms, nugget = 0, m = 30, approx = "cc", target = "spatialRange", N = 100, xlim = c(0.05, 0.15), sdlog = 0.6)
{
  if(target == "spatialRange" | target == "spatial") {
    
    output <- .posterior_for_spatialRange(z = z, locs = locs, covparms = covparms, nugget = nugget, m = m, approx = approx, N = N, xlim = xlim, sdlog = sdlog)
    
  } else if(target == "temporalRange" | target == "temporal") {
    
    output <- .posterior_for_temporalRange(z = z, locs = locs, covparms = covparms, nugget = nugget, m = m, approx = approx, N = N, xlim = xlim, sdlog = sdlog)
    
  } else {
    
    stop("Check the target!")
    
  }
  
  return(output)
}

.posterior_for_spatialRange <- function(z, locs, covparms, nugget, m, approx, N, xlim, sdlog)
{
  n               <- nrow(locs)
  
  meanlog         <- log(covparms[2])
  logprior        <- function(alpha) dlnorm(alpha, meanlog, sdlog, log = TRUE)
  
  alpha.grid      <- seq(xlim[1], xlim[2], length = N)
  logpost         <- array(dim = c(N, 2))
  for(i in 1:N) {
    
    print(i)
    alpha           <- alpha.grid[i]
    
    Sigma.true      <- cov_matern_spacetime(locs = locs, covparms = c(covparms[1], alpha, covparms[3], covparms[4])) + nugget * diag(n)
    
    ## exact
    loglik          <- dmvnorm(z, mean = rep(0, n), sigma = Sigma.true, log = TRUE)
    logpost[i, 1]   <- loglik + logprior(alpha)
    
    ## Vecchia
    if(approx == "b1") {
      vecchia.approx  <- correlationVecchia::baseline_1_spacetime_specify(locs = locs, m = m)
    } else if(approx == "b2") {
      vecchia.approx  <- correlationVecchia::baseline_2_spacetime_specify(locs = locs, m = m)
    } else if(approx == "b3") {
      vecchia.approx  <- correlationVecchia::baseline_3_spacetime_specify(locs = locs, m = m, covmodel = cov_matern_spacetime, covparms = c(covparms[1], alpha, covparms[3], covparms[4]))
    } else {
      vecchia.approx  <- correlationVecchia::corrvecchia_specify_knownCovparms_2(locs = locs, m = m, covmodel = Sigma.true, covparms = c(covparms[1], alpha, covparms[3], covparms[4]))
    }
    
    Sigma.ord       <- Sigma.true[vecchia.approx$ord, vecchia.approx$ord]
    U               <- GPvecchia::createU(vecchia.approx, covparms = c(1), nuggets = nugget, covmodel = Sigma.ord)$U
    revord          <- order(vecchia.approx$ord)
    Sigma.vecchia   <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
    
    loglik          <- mvtnorm::dmvnorm(z, mean = rep(0, n), sigma = Sigma.vecchia, log = TRUE)
    logpost[i, 2]   <- loglik + logprior(alpha)
  }
  
  post            <- logpost - mean(logpost)
  A               <- mean(post[, 1])
  B               <- mean(post[, 2])
  
  post.norm       <- post
  post.norm[, 1]  <- exp(post[, 1] - A) / sum(exp(post[, 1] - A))
  post.norm[, 2]  <- exp(post[, 2] - B) / sum(exp(post[, 2] - B))
  
  post            <- exp(post)
  
  # post            <- exp(logpost - mean(logpost))
  # post.norm       <- post
  # for(j in 1:2) post.norm[, j] <- post[, j] / sum(post[, j])
  
  return(list(alpha = alpha.grid, logpost = logpost, post = post, post.norm = post.norm))
}

.posterior_for_temporalRange <- function(z, locs, covparms, nugget, m, approx, N, xlim, sdlog)
{
  n               <- nrow(locs)
  
  meanlog         <- log(covparms[3])
  logprior        <- function(alpha) dlnorm(alpha, meanlog, sdlog, log = TRUE)
  
  alpha.grid      <- seq(xlim[1], xlim[2], length = N)
  logpost         <- array(dim = c(N, 2))
  for(i in 1:N) {
    
    print(i)
    alpha           <- alpha.grid[i]
    
    Sigma.true      <- cov_matern_spacetime(locs = locs, covparms = c(covparms[1], covparms[2], alpha, covparms[4])) + nugget * diag(n)
    
    ## exact
    loglik          <- mvtnorm::dmvnorm(z, mean = rep(0, n), sigma = Sigma.true, log = TRUE)
    logpost[i, 1]   <- loglik + logprior(alpha)
    
    ## Vecchia
    if(approx == "b1") {
      vecchia.approx  <- correlationVecchia::baseline_1_spacetime_specify(locs = locs, m = m)
    } else if(approx == "b2") {
      vecchia.approx  <- correlationVecchia::baseline_2_spacetime_specify(locs = locs, m = m)
    } else if(approx == "b3") {
      vecchia.approx  <- correlationVecchia::baseline_3_spacetime_specify(locs = locs, m = m, covmodel = cov_matern_spacetime, covparms = c(covparms[1], covparms[2], alpha, covparms[4]))
    } else {
      vecchia.approx  <- correlationVecchia::corrvecchia_specify_knownCovparms_2(locs = locs, m = m, covmodel = Sigma.true, covparms = c(covparms[1], covparms[2], alpha, covparms[4]))
    }
    
    Sigma.ord       <- Sigma.true[vecchia.approx$ord, vecchia.approx$ord]
    U               <- GPvecchia::createU(vecchia.approx, covparms = c(1), nuggets = nugget, covmodel = Sigma.ord)$U
    revord          <- order(vecchia.approx$ord)
    Sigma.vecchia   <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
    
    loglik          <- mvtnorm::dmvnorm(z, mean = rep(0, n), sigma = Sigma.vecchia, log = TRUE)
    logpost[i, 2]   <- loglik + logprior(alpha)
  }
  
  post            <- logpost - mean(logpost)
  A               <- mean(post[, 1])
  B               <- mean(post[, 2])
  
  post.norm       <- post
  post.norm[, 1]  <- exp(post[, 1] - A) / sum(exp(post[, 1] - A))
  post.norm[, 2]  <- exp(post[, 2] - B) / sum(exp(post[, 2] - B))
  
  post            <- exp(post)
  
  # post            <- exp(logpost - mean(logpost))
  # post.norm       <- post
  # for(j in 1:2) post.norm[, j] <- post[, j] / sum(post[, j])
  
  return(list(alpha = alpha.grid, logpost = logpost, post = post, post.norm = post.norm))
}
