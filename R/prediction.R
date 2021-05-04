####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: 
###
###   Contents:
###       FTN - 
###
####################################################################################



#' Wrapper function for Matern covariance with scaled dimensions (GpGp package)
#'
#' @param locs A matrix of locations
#' @param covparms A numerical vector with covariance parameters = (variance, range_1, ..., range_d, smoothness)
#'
#' @return A matrix with \code{n} rows and \code{n} column
#' 
#' @export
#'
#' @examples
#' 1 + 1
cov_matern_scaledim <- function(locs, covparms)
{
  return(GpGp::matern_scaledim(c(covparms, 0), locs)) # 0 <- nugget
}



#' Simple function for calculating logarithmic prediction score
#'
#' @param mu.pred predicted mean vector
#' @param var.pred predicted variance vector
#' @param z observations of the process of interest
#'
#' @return A numeric value
#' 
#' @export
#'
#' @examples
#' 1 + 1
logscore <- function(mu.pred, var.pred, z)
{
  S <- dnorm(x = z, mean = mu.pred, sd = sqrt(var.pred), log = TRUE)
  
  return(S)
}



#' Prediction function of baseline approximations (for cov_matern_spacetime)
#'
#' @param z input
#' @param locs input 
#' @param locs.pred input
#' @param m input
#' @param coordinate input 
#' @param covparms input
#' @param nuggets input
#' @param method input
#' @param predcond.method input 
#' @param var.exact input
#' @param return.values input
#'
#' @return output
#' 
#' @export
#'
#' @examples
#' 1 + 1
prediction_baseline_for_spacetime <- function(z, locs, locs.pred, m, coordinate = NULL, covparms, nuggets, method = "b1", predcond.method = "general", var.exact, return.values = "all")
{
  ### main
  if(method == "b1") {
    
    ord.obsv          <- order_time(locs = locs, coordinate = coordinate)
    ord.pred          <- order_time(locs = locs.pred, coordinate = coordinate)
    ord               <- c(ord.obsv, ord.pred + nrow(locs))
    
    locsord           <- rbind(locs, locs.pred)[ord, , drop = FALSE]
    obs               <- c(rep(TRUE, nrow(locs)), rep(FALSE, nrow(locs.pred)))
    
    ord.z             <- ord[seq(nrow(locs))]
    ord.pred          <- "obspred"
    cond.yz           <- "y"
    
    cond.sets   <- matrix(NA, nrow = nrow(locsord), ncol = m + 1)
    for(i in 1:nrow(locsord)) {
      ind                 <- seq(from = 1, to = min(i, m + 1), by = 1)
      cond.sets[i, ind]   <- seq(from = i, by = -1, length.out = length(ind)) 
    }
    
    Cond        <- matrix(NA, nrow(cond.sets), ncol(cond.sets)); Cond[!is.na(cond.sets)] <- TRUE
    U.prep      <- GPvecchia:::U_sparsity(locsord, cond.sets, obs, Cond)
    
    vecchia.approx    <- list(locsord = locsord, obs = obs, ord = ord, ord.z = ord.z, ord.pred = predcond.method, U.prep = U.prep, cond.yz = cond.yz, ordering = "maxmin", conditioning = "NN", ic0 = FALSE)
    
    covmat            <- correlationVecchia:::.correlation(locs = locsord, covmodel = cov_matern_spacetime, covparms = covparms, abs.corr = FALSE)
    vecchia.predict   <- vecchia_prediction(z = z, vecchia.approx = vecchia.approx, covparms = covparms[1], nuggets = nuggets, var.exact = var.exact, covmodel = covmat, return.values = return.values)
    
  } else if(method == "b2") {
    
    ord.obsv          <- order_time(locs = locs, coordinate = coordinate)
    ord.pred          <- order_time(locs = locs.pred, coordinate = coordinate)
    ord               <- c(ord.obsv, ord.pred + nrow(locs))
    
    locsord           <- rbind(locs, locs.pred)[ord, , drop = FALSE]
    obs               <- c(rep(TRUE, nrow(locs)), rep(FALSE, nrow(locs.pred)))
    
    ord.z             <- ord[seq(nrow(locs))]
    ord.pred          <- "obspred"
    cond.yz           <- "y"
    
    d.spacetime <- fields::rdist(x1 = locsord[, 1:2, drop = FALSE]) + 1 * fields::rdist(x1 = locsord[, 3, drop = FALSE])
    
    cond.sets   <- conditioning_nn(m = m, d = d.spacetime)
    
    Cond        <- matrix(NA, nrow(cond.sets), ncol(cond.sets)); Cond[!is.na(cond.sets)] <- TRUE
    U.prep      <- GPvecchia:::U_sparsity(locsord, cond.sets, obs, Cond)
    
    vecchia.approx    <- list(locsord = locsord, obs = obs, ord = ord, ord.z = ord.z, ord.pred = predcond.method, U.prep = U.prep, cond.yz = cond.yz, ordering = "maxmin", conditioning = "NN", ic0 = FALSE)
    
    covmat            <- correlationVecchia:::.correlation(locs = locsord, covmodel = cov_matern_spacetime, covparms = covparms, abs.corr = FALSE)
    vecchia.predict   <- vecchia_prediction(z = z, vecchia.approx = vecchia.approx, covparms = covparms[1], nuggets = nuggets, var.exact = var.exact, covmodel = covmat, return.values = return.values)
    
  } else if(method == "b3") {
    
    ord.obsv          <- order_time(locs = locs, coordinate = coordinate)
    ord.pred          <- order_time(locs = locs.pred, coordinate = coordinate)
    ord               <- c(ord.obsv, ord.pred + nrow(locs))
    
    locsord           <- rbind(locs, locs.pred)[ord, , drop = FALSE]
    obs               <- c(rep(TRUE, nrow(locs)), rep(FALSE, nrow(locs.pred)))
    
    ord.z             <- ord[seq(nrow(locs))]
    ord.pred          <- "obspred"
    cond.yz           <- "y"
    
    covmat            <- correlationVecchia:::.correlation(locs = locsord, covmodel = cov_matern_spacetime, covparms = covparms, abs.corr = FALSE)
    cond.sets         <- conditioning_nn(m = m, d = covparms[1] - covmat)
    Cond              <- matrix(NA, nrow(cond.sets), ncol(cond.sets)); Cond[!is.na(cond.sets)] <- TRUE
    U.prep            <- GPvecchia:::U_sparsity(locsord, cond.sets, obs, Cond)
    
    vecchia.approx    <- list(locsord = locsord, obs = obs, ord = ord, ord.z = ord.z, ord.pred = predcond.method, U.prep = U.prep, cond.yz = cond.yz, ordering = "maxmin", conditioning = "NN", ic0 = FALSE)
    
    vecchia.predict   <- vecchia_prediction(z = z, vecchia.approx = vecchia.approx, covparms = covparms[1], nuggets = nuggets, var.exact = var.exact, covmodel = covmat, return.values = return.values)
    
  } else {
    
    stop("Check the argument method!")
    
  }
  
  ### return
  return(list(approx = vecchia.approx, predict = vecchia.predict))
}



#' Prediction function of Vecchia approximations (for cov_matern_spacetime)
#'
#' @param z input
#' @param locs input
#' @param locs.pred input
#' @param m input
#' @param initial.pt input 
#' @param covparms input
#' @param nuggets input
#' @param method input
#' @param predcond.method input 
#' @param var.exact input
#' @param return.values input
#'
#' @return output
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' n               <- 30^2
#' n.pred          <- 10^2
#' m               <- 10
#' 
#' locs            <- matrix(runif(3 * n), n, 3)
#' locs.pred       <- matrix(runif(3 * n.pred), n.pred, 3)
#' 
#' nugget          <- 0.0
#' covparms        <- c(1, 0.1, 0.1, 0.5)
#' covmat          <- cov_matern_spacetime(locs = rbind(locs, locs.pred), covparms = covparms) # + nugget * diag(nrow(locs) + nrow(locs.pred))
#' 
#' z.all            <- as.numeric(t(chol(x = covmat)) %*% as.matrix(rnorm(nrow(locs) + nrow(locs.pred))))
#' z               <- z.all[seq(nrow(locs))]
#' z.pred          <- z.all[seq(nrow(locs.pred)) + nrow(locs)]
#' 
#' out.euclidean   <- prediction_corrvecchia_for_spacetime(z = z, locs = locs, locs.pred = locs.pred, m = m, initial.pt = NULL, covparms = covparms, nuggets = nugget, method = "euclidean", predcond.method = "general", var.exact = TRUE, return.values = "all")
#' out.corrvecchia <- prediction_corrvecchia_for_spacetime(z = z, locs = locs, locs.pred = locs.pred, m = m, initial.pt = NULL, covparms = covparms, nuggets = nugget, method = "correlation", predcond.method = "general", var.exact = TRUE, return.values = "all")
#' 
#' # MSPE
#' 
#' mean((out.corrvecchia$predict$mu.pred - z.pred)^2)
#' mean((out.euclidean$predict$mu.pred - z.pred)^2)
#' 
#' # logscore
#' 
#' mean(logscore(out.corrvecchia$predict$mu.pred, out.corrvecchia$predict$var.pred, z.pred))
#' mean(logscore(out.euclidean$predict$mu.pred, out.euclidean$predict$var.pred, z.pred))
#' }
prediction_corrvecchia_for_spacetime <- function(z, locs, locs.pred, m, initial.pt = NULL, covparms, nuggets, method = "correlation", predcond.method = "general", var.exact, return.values = "all")
{
  ### initial input
  if(is.null(initial.pt)) {
    cen               <- t(as.matrix(colMeans(locs)))
    initial.pt        <- which.min(rowSums((locs - matrix(as.numeric(cen), nrow = nrow(locs), ncol = ncol(locs), byrow = T))^2))
  }
  
  ### main
  if(method == "euclidean") {
    
    ord               <- predSortSparse_Rcpp(locs, locs.pred, 2.01, initial.pt-1, "euclidean", "cov_matern_spacetime_cpp", covparms)$P + 1
    locsord           <- rbind(locs, locs.pred)[ord, , drop = FALSE]
    obs               <- c(rep(TRUE, nrow(locs)), rep(FALSE, nrow(locs.pred)))
    
    ord.z             <- ord[seq(nrow(locs))]
    ord.pred          <- "obspred"
    cond.yz           <- "y"
    
    cond.sets         <- GpGp::find_ordered_nn(locs = locsord, m = m)
    Cond              <- matrix(NA, nrow(cond.sets), ncol(cond.sets)); Cond[!is.na(cond.sets)] <- TRUE
    U.prep            <- GPvecchia:::U_sparsity(locsord, cond.sets, obs, Cond)
    
    vecchia.approx    <- list(locsord = locsord, obs = obs, ord = ord, ord.z = ord.z, ord.pred = predcond.method, U.prep = U.prep, cond.yz = cond.yz, ordering = "maxmin", conditioning = "NN", ic0 = FALSE)
    
    covmat            <- correlationVecchia:::.correlation(locs = locsord, covmodel = cov_matern_spacetime, covparms = covparms, abs.corr = FALSE)
    vecchia.predict   <- vecchia_prediction(z = z, vecchia.approx = vecchia.approx, covparms = covparms[1], nuggets = nuggets, var.exact = var.exact, covmodel = covmat, return.values = return.values)
    
  } else if(method == "correlation") {
    
    ord               <- predSortSparse_Rcpp(locs, locs.pred, 2.01, initial.pt-1, "correlation", "cov_matern_spacetime_cpp", covparms)$P + 1
    locsord           <- rbind(locs, locs.pred)[ord, , drop = FALSE]
    obs               <- c(rep(TRUE, nrow(locs)), rep(FALSE, nrow(locs.pred)))
    
    ord.z             <- ord[seq(nrow(locs))]
    ord.pred          <- "obspred"
    cond.yz           <- "y"
    
    covmat            <- correlationVecchia:::.correlation(locs = locsord, covmodel = cov_matern_spacetime, covparms = covparms, abs.corr = FALSE)
    cond.sets         <- conditioning_nn(m = m, d = covparms[1] - covmat)
    Cond              <- matrix(NA, nrow(cond.sets), ncol(cond.sets)); Cond[!is.na(cond.sets)] <- TRUE
    U.prep            <- GPvecchia:::U_sparsity(locsord, cond.sets, obs, Cond)
    
    vecchia.approx    <- list(locsord = locsord, obs = obs, ord = ord, ord.z = ord.z, ord.pred = predcond.method, U.prep = U.prep, cond.yz = cond.yz, ordering = "maxmin", conditioning = "NN", ic0 = FALSE)
    
    vecchia.predict   <- vecchia_prediction(z = z, vecchia.approx = vecchia.approx, covparms = covparms[1], nuggets = nuggets, var.exact = var.exact, covmodel = covmat, return.values = return.values)
    
  } else {
    
    stop("Check the argument method!")
    
  }
  
  ### return
  return(list(approx = vecchia.approx, predict = vecchia.predict))
}



#' Prediction function of Vecchia approximations (for cov_matern_scaledim)
#'
#' @param z input
#' @param locs input
#' @param locs.pred input
#' @param m input
#' @param initial.pt input
#' @param covparms input
#' @param nuggets input
#' @param method input
#' @param predcond.method input
#' @param var.exact input
#' @param return.values input
#'
#' @return output
#' 
#' @export
#'
#' @examples
#' 1 + 1
prediction_corrvecchia_scaledim <- function(z, locs, locs.pred, m, initial.pt = NULL, covparms, nuggets, method = "correlation", predcond.method = "general", var.exact, return.values = "all")
{
  ### initial input
  if(is.null(initial.pt)) {
    cen               <- t(as.matrix(colMeans(locs)))
    initial.pt        <- which.min(rowSums((locs - matrix(as.numeric(cen), nrow = nrow(locs), ncol = ncol(locs), byrow = T))^2))
  }
  
  ### main
  if(method == "euclidean") {
    
    sig2 <- covparms[1] ; nu <- covparms[length(covparms)-1] ; ranges <- covparms[seq(from = 2, to = length(covparms)-2)]
    
    ord               <- predSortSparse_Rcpp(locs, locs.pred, 2.01, initial.pt-1, "euclidean", "cov_matern_scaledim_cpp", covparms)$P + 1
    locsord           <- rbind(locs, locs.pred)[ord, , drop = FALSE]
    obs               <- c(rep(TRUE, nrow(locs)), rep(FALSE, nrow(locs.pred)))
    
    ord.z             <- ord[seq(nrow(locs))]
    ord.pred          <- "obspred"
    cond.yz           <- "y"
    
    cond.sets         <- GpGp::find_ordered_nn(locs = locsord, m = m)
    Cond              <- matrix(NA, nrow(cond.sets), ncol(cond.sets)); Cond[!is.na(cond.sets)] <- TRUE
    U.prep            <- GPvecchia:::U_sparsity(locsord, cond.sets, obs, Cond)
    
    vecchia.approx    <- list(locsord = locsord, obs = obs, ord = ord, ord.z = ord.z, ord.pred = predcond.method, U.prep = U.prep, cond.yz = cond.yz, ordering = "maxmin", conditioning = "NN", ic0 = FALSE)
    
    covmat            <- cov_matern_scaledim(locsord, covparms)
    vecchia.predict   <- vecchia_prediction(z = z, vecchia.approx = vecchia.approx, covparms = covparms[1], nuggets = nuggets, var.exact = var.exact, covmodel = covmat, return.values = return.values)
    
  } else if(method == "correlation") {
    
    sig2 <- covparms[1] ; nu <- covparms[length(covparms)-1] ; ranges <- covparms[seq(from = 2, to = length(covparms)-2)]
    
    ord               <- predSortSparse_Rcpp(locs, locs.pred, 2.01, initial.pt-1, "correlation", "cov_matern_scaledim_cpp", covparms)$P + 1
    locsord           <- rbind(locs, locs.pred)[ord, , drop = FALSE]
    obs               <- c(rep(TRUE, nrow(locs)), rep(FALSE, nrow(locs.pred)))
    
    ord.z             <- ord[seq(nrow(locs))]
    ord.pred          <- "obspred"
    cond.yz           <- "y"
    
    covmat            <- cov_matern_scaledim(locsord, covparms)
    cond.sets         <- conditioning_nn(m = m, d = covparms[1] - covmat)
    Cond              <- matrix(NA, nrow(cond.sets), ncol(cond.sets)); Cond[!is.na(cond.sets)] <- TRUE
    U.prep            <- GPvecchia:::U_sparsity(locsord, cond.sets, obs, Cond)
    
    vecchia.approx    <- list(locsord = locsord, obs = obs, ord = ord, ord.z = ord.z, ord.pred = predcond.method, U.prep = U.prep, cond.yz = cond.yz, ordering = "maxmin", conditioning = "NN", ic0 = FALSE)
    
    vecchia.predict   <- vecchia_prediction(z = z, vecchia.approx = vecchia.approx, covparms = covparms[1], nuggets = nuggets, var.exact = var.exact, covmodel = covmat, return.values = return.values)
    
  } else {
    
    stop("Check the argument method!")
    
  }
  
  ### return
  return(list(approx = vecchia.approx, predict = vecchia.predict))
}



#' Visualization function for prediction
#'
#' @param vdat1 input
#' @param vdat2 input
#' @param vdat3 input
#' @param ylab_vdat1 input
#' @param ylab_vdat2 input
#' @param ylab_vdat3 input
#' @param legend input
#' @param color input
#' @param shape input
#' @param alpha input
#' @param size.point input
#' @param size.line input
#' @param size.legend input
#' @param size.lab input
#' @param size.text input
#' @param size.margin input
#'
#' @return output
#' 
#' @export
#'
#' @examples
#' 1 + 1
vis_prediction <- function(vdat1, vdat2, vdat3, ylab_vdat1 = "log10( logscore )", ylab_vdat2 = "log10( logscore )", ylab_vdat3 = "log10( logscore )", legend, color, shape, alpha = 0.7, size.point = 4, size.line = 1, size.legend = 14, size.lab = 14, size.text = 12, size.margin = c(5.5, 20, 5.5, 5.5))
{
  # library(ggplot2) ; library(dplyr) ; library(RColorBrewer) ; library(gridExtra)
  
  # n.space   <- 10
  # legend    <- paste0( legend, paste0(rep(" ", n.space), collapse = "") )
  
  ### Gathering vdat1 ###
  
  vars1     <- vdat1 %>% select(-index, -starts_with("approx")) %>% colnames()
  if(length(vars1) > 1) stop("Please check the argument vdat1.")
  
  xlabel1   <- vdat1 %>% pull(vars1) %>% unique() %>% sort()
  vis1      <- vdat1 %>% tidyr::gather(key = "approx", value = "logscore", -index, -vars1)
  
  ### plot vis1 ###
  
  plot1     <- vis1 %>% ggplot(aes(x = get(vars1), y = logscore, colour = approx, shape = approx)) + 
    geom_point(size = size.point) + geom_line(size = size.line, alpha = alpha) + 
    scale_x_continuous(name = vars1, limits = range(xlabel1), breaks = xlabel1) +
    scale_colour_manual(name = NULL, values = color, labels = legend, guide = guide_legend(nrow = 1, byrow = TRUE)) + 
    scale_shape_manual(name = NULL, values = shape, labels = legend, guide = guide_legend(nrow = 1, byrow = TRUE)) +
    xlab(vars1) + ylab(ylab_vdat1) + 
    theme(axis.title.x = element_text(size = size.lab), 
          axis.text.x = element_text(size = size.text), 
          axis.title.y = element_text(size = size.lab), 
          axis.text.y = element_text(size = size.text), 
          legend.title = element_blank(), 
          legend.text = element_text(size = size.legend, margin = margin(r = 25, unit = 'pt')), 
          legend.direction = 'horizontal', legend.spacing.x = unit(15, 'pt'),
          plot.margin = unit(size.margin, "pt")) # t, r, b, l
  
  ### Gathering vdat2 ###
  
  vars2     <- vdat2 %>% select(-index, -starts_with("approx")) %>% colnames()
  if(length(vars2) > 1) stop("Please check the argument vdat12.")
  
  xlabel2   <- vdat2 %>% pull(vars2) %>% unique() %>% sort()
  vis2      <- vdat2 %>% tidyr::gather(key = "approx", value = "logscore", -index, -vars2)
  
  ### plot vis1 ###
  
  plot2     <- vis2 %>% ggplot(aes(x = get(vars2), y = logscore, colour = approx, shape = approx)) + 
    geom_point(size = size.point) + geom_line(size = size.line, alpha = alpha) + 
    scale_x_continuous(name = vars2, limits = range(xlabel2), breaks = xlabel2) +
    scale_colour_manual(name = NULL, values = color, labels = legend, guide = guide_legend(nrow = 1, byrow = TRUE)) + 
    scale_shape_manual(name = NULL, values = shape, labels = legend, guide = guide_legend(nrow = 1, byrow = TRUE)) +
    xlab(vars2) + ylab(ylab_vdat2) + 
    theme(axis.title.x = element_text(size = size.lab), 
          axis.text.x = element_text(size = size.text), 
          axis.title.y = element_text(size = size.lab), 
          axis.text.y = element_text(size = size.text), 
          legend.title = element_blank(), 
          legend.text = element_text(size = size.legend, margin = margin(r = 25, unit = 'pt')), 
          legend.direction = 'horizontal', legend.spacing.x = unit(15, 'pt'),
          plot.margin = unit(size.margin, "pt")) # t, r, b, l
  
  ### Gathering vdat1 ###
  
  vars3     <- vdat3 %>% select(-index, -starts_with("approx")) %>% colnames()
  if(length(vars3) > 1) stop("Please check the argument vdat3.")
  
  xlabel3   <- vdat3 %>% pull(vars3) %>% unique() %>% sort()
  vis3      <- vdat3 %>% tidyr::gather(key = "approx", value = "logscore", -index, -vars3)
  
  ### plot vis1 ###
  
  plot3     <- vis3 %>% ggplot(aes(x = get(vars3), y = logscore, colour = approx, shape = approx)) + 
    geom_point(size = size.point) + geom_line(size = size.line, alpha = alpha) + 
    scale_x_continuous(name = vars3, limits = range(xlabel3), breaks = xlabel3) +
    scale_colour_manual(name = NULL, values = color, labels = legend, guide = guide_legend(nrow = 1, byrow = TRUE)) + 
    scale_shape_manual(name = NULL, values = shape, labels = legend, guide = guide_legend(nrow = 1, byrow = TRUE)) +
    xlab(vars3) + ylab(ylab_vdat3) + 
    theme(axis.title.x = element_text(size = size.lab), 
          axis.text.x = element_text(size = size.text), 
          axis.title.y = element_text(size = size.lab), 
          axis.text.y = element_text(size = size.text), 
          legend.title = element_blank(), 
          legend.text = element_text(size = size.legend, margin = margin(r = 25, unit = 'pt')), 
          legend.direction = 'horizontal', legend.spacing.x = unit(15, 'pt'),
          plot.margin = unit(size.margin, "pt")) # t, r, b, l
  
  ### Creating the legend ###
  
  tmp       <- ggplot_gtable(ggplot_build(plot1))
  leg       <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  mylegend  <- tmp$grobs[[leg]]
  
  ### Merge the two plots ###
  
  result    <- grid.arrange(mylegend, arrangeGrob(plot1 + theme(legend.position="none"), plot2 + theme(legend.position="none"), plot3 + theme(legend.position="none"), nrow = 1), nrow = 2, heights = c(1, 10))
  
  ### return ###
  
  return(result)
}



### CAUTION: The functions below are basically identical to the corresponding functions of GPvecchia. I need to work on this. This is actually redundant.

vecchia_prediction <- function(z, vecchia.approx, covparms, nuggets, var.exact, covmodel = 'matern', return.values = 'all') 
{
  # remove NAs in data and U
  GPvecchia:::removeNAs()
  
  # create the U matrix
  U.obj         <- GPvecchia::createU(vecchia.approx,covparms,nuggets,covmodel)
  U.obj$ic0     <- FALSE
  
  # compute cholesky V for posterior inference
  V.ord         <- U2V(U.obj)
  
  if(length(U.obj$zero.nugg) > 0) warning('Rows/cols of V have been removed for data with zero noise')
  
  # compute the posterior mean
  vecchia.mean  <- vecchia_mean(z,U.obj,V.ord)
  
  # return what is requested
  return.list   <- list(mu.pred = vecchia.mean$mu.pred, mu.obs = vecchia.mean$mu.obs, var.pred = NULL, var.obs = NULL, V.ord = NULL, U.obj = NULL)
  
  if(return.values == 'meanmat' | return.values == 'all') {
    return.list$V.ord <- V.ord ; return.list$U.obj <- U.obj
  }
  
  if(return.values == 'meanvar' | return.values == 'all') {
    
    # compute posterior variances
    if(missing(var.exact)) var.exact <- (sum(!vecchia.approx$obs)<4*1e4)
    
    vars.vecchia <- vecchia_var(U.obj, V.ord, exact = var.exact)
    
    return.list$var.pred  <- vars.vecchia$vars.pred
    return.list$var.obs   <- vars.vecchia$vars.obs
    
  }
  
  return(return.list)
}

U2V <- function(U.obj)
{
  ### when changing this function make sure it returns a dtCMatrix!
  ### Otherwise solve in the parent function will be very slow
  
  U.y           <- U.obj$U[U.obj$latent, ]
  
  if(U.obj$cond.yz == 'zy') {
    
    V.ord         <- GPvecchia:::revMat(U.y[ ,U.obj$latent, drop = FALSE])
    
  } else if(U.obj$ord.pred != 'obspred') {
    
    W             <- Matrix::tcrossprod(U.y)
    W.rev         <- GPvecchia:::revMat(W)
    
    if(U.obj$ic0) {
      V.ord         <- Matrix::t(ichol(W.rev))
    } else {
      V.ord         <- Matrix::t(Matrix::chol(W.rev))
    }
    
    V.ord         <- as(V.ord, 'dtCMatrix')
    
    
  } else {  # for obspred ordering
    
    last.obs      <- max(which(!U.obj$latent))
    latents.before <- sum(U.obj$latent[1:last.obs])
    latents.after <- sum(U.obj$latent[-(1:last.obs)])
    
    # pred columns are unchanged
    V.pr          <- revMat(U.y[ ,(last.obs+1):ncol(U.y), drop = FALSE])
    
    # have to compute cholesky for obs block
    U.oo          <- U.y[1:latents.before, 1:last.obs]
    A             <- Matrix::tcrossprod(U.oo)
    A.rev         <- GPvecchia:::revMat(A)
    
    if(U.obj$ic0) { 
      V.oor         <- Matrix::t(ichol(A.rev))
    } else {
      V.oor         <- Matrix::t(Matrix::chol(A.rev))
    } 
    
    # combine the blocks into one matrix
    zeromat.sparse <- Matrix::sparseMatrix(c(), c(), dims = c(latents.after, latents.before))
    V.or          <- rbind(zeromat.sparse, V.oor)
    
    V.ord         <- methods::as(cbind(V.pr, V.or), 'dtCMatrix')
    
  }
  
  return(V.ord)
}

ichol <- function(M, S = NULL)
{
  if(!is(M, "sparseMatrix")) warning("Passing a dense matrix")
  
  if(!is(M, "CsparseMatrix") || !Matrix::isTriangular(M)){
    M               <- as(Matrix::triu(M), "CsparseMatrix")
  }
  
  if(!is.null(S)){
    
    if(!is(S, "sparseMatrix")) S <- as(Matrix::triu(S), "CsparseMatrix")
    
    p <- S@p ; i <- S@i
    
  } else {
    p <- M@p ; i <- M@i
  }
  
  vals            <- ic0(p, i, M@x)
  Msp             <- Matrix::sparseMatrix(i = i, p = p, x = vals, index1 = FALSE)
  Msp@x           <- vals
  
  return(Msp)
}

vecchia_mean <- function(z, U.obj, V.ord)
{
  U               <- U.obj$U
  
  # compute entire posterior mean vector
  z.ord           <- z[U.obj$ord.z]
  z1              <- Matrix::crossprod(U[!U.obj$latent, ], z.ord)
  z2              <- as.numeric(U[U.obj$latent, ] %*% z1)
  
  temp            <- Matrix::solve(V.ord, rev(z2))
  mu.rev          <- -Matrix::solve(Matrix::t(V.ord), temp)
  mu.ord          <- rev(mu.rev)
  
  # for zero nugget, observations are posterior means
  if(length(U.obj$zero.nugg) > 0){
    obs.zero        <- z.ord[U.obj$zero.nugg$inds.z]
    mu.ord          <- c(mu.ord, obs.zero)
  }
  
  # extract obs and pred parts; return to original ordering
  orig.order      <- order(U.obj$ord)
  mu              <- mu.ord[orig.order]
  obs.orig        <- U.obj$obs[orig.order]
  mu.obs          <- mu[obs.orig]
  mu.pred         <- mu[!obs.orig]
  
  return(list(mu.obs = mu.obs, mu.pred = mu.pred))
}

vecchia_lincomb <- function(H, U.obj, V.ord, cov.mat = FALSE) 
{
  ord             <- U.obj$ord
  
  if(length(U.obj$zero.nugg)>0) ord <- order(order(ord[1:(length(ord)-length(U.obj$zero.nugg$inds.U))]))
  
  H.tt            <- Matrix::t(H[, rev(ord), drop = FALSE])
  temp            <- Matrix::solve(V.ord, H.tt)
  
  if(cov.mat) {
    
    lincomb.cov     <- as.matrix(Matrix::t(temp) %*% temp)
    return(lincomb.cov)
    
  } else {
    
    lincomb.vars    <- as.numeric(Matrix::t(temp * temp) %*% rep(1, ncol(H)))
    return(lincomb.vars)
    
  }
}

SelInv <- function(cholmat)
{
  n.all           <- nrow(cholmat)
  sparseinv::Takahashi_Davis(Q = Matrix::sparseMatrix(c(), c(), dims = c(n.all, 1)), cholQp = cholmat, P = Matrix::sparseMatrix(i = 1:n.all, j = 1:n.all, x = 1))
}

vecchia_var <- function(U.obj, V.ord, exact = FALSE)
{
  # compute selected inverse and extract variances
  inv.sparse      <- SelInv(V.ord)
  vars.ord        <- rev(Matrix::diag(inv.sparse))
  
  # for zero nugget, add zero variances
  if(length(U.obj$zero.nugg) > 0) vars.ord <- c(vars.ord, rep(0, length(U.obj$zero.nugg$inds.z)))
  
  
  # return to original ordering
  orig.order      <- order(U.obj$ord)
  vars            <- vars.ord[orig.order]
  
  # extract obs and pred variances
  obs.orig        <- U.obj$obs[orig.order]
  vars.obs        <- vars[obs.orig]
  vars.pred       <- vars[!obs.orig]
  
  # if exact variances are desired:
  if(exact & U.obj$ord.pred == 'obspred'){
    
    n               <- sum(U.obj$obs)
    n.p             <- sum(!U.obj$obs)
    H               <- Matrix::sparseMatrix(i = 1:(n + n.p), j = 1:(n + n.p), x = 1)
    
    if(U.obj$cond.yz != 'zy') {
      
      if(sum(!obs.orig)>0) {
        vars.pred       <- vecchia_lincomb(H[(n + 1):(n + n.p), ], U.obj, V.ord)
      } else {
        vars.pred       <- c()
      }
      
    } else {
      
      vars.exact      <- vecchia_lincomb(H, U.obj, V.ord)
      vars.obs        <- vars.exact[1:n]
      vars.pred       <- (if(n.p>0) vars.exact[n + (1:n.p)] else c())
      
    }
    
  }
  
  return(list(vars.obs = vars.obs, vars.pred = vars.pred))
}

V2covmat <- function(preds)
{
  # compute joint covariance matrix
  Sigma.ord       <- Matrix::solve(as.matrix(GPvecchia:::revMat(preds$V.ord%*%Matrix::t(preds$V.ord))))
  
  # for zero nugget, add zero rows/columns
  if(length(preds$U.obj$zero.nugg)>0) {
    k               <- nrow(Sigma.ord)
    l               <- length(preds$U.obj$zero.nugg$inds.z)
    Sigma.ord       <- rbind(cbind(Sigma.ord, matrix(0, nrow = k, ncol = l)), matrix(0, nrow = l, ncol = k + l))
  }
  
  # return to original ordering
  orig.order      <- order(preds$U.obj$ord)
  Sigma           <- Sigma.ord[orig.order, orig.order]
  
  # extract parts corresponding to obs and pred locs
  obs.orig        <- preds$U.obj$obs[orig.order]
  Sigma.obs       <- Sigma[obs.orig, obs.orig]
  Sigma.pred      <- Sigma[!obs.orig, !obs.orig]
  
  return(list(Sigma.obs = Sigma.obs, Sigma.pred = Sigma.pred))
}
