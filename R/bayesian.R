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
    
    loglik          <- dmvnorm(z, mean = rep(0, 2 * n), sigma = Sigma.vecchia, log = TRUE)
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
