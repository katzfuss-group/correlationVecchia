####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: 
###
####################################################################################

rm(list = ls())
set.seed(12122020)

### load packages
library(correlationVecchia)
library(GPvecchia)
library(RandomFields)
library(mvtnorm)
library(LaplacesDemon)

### setting #############################################################################################

n               <- 30^2
covparms        <- c(1, 0.1, 0.1)
process         <- generate_gp_space(nsim = 1, n = n, d = 2, p = 2, method.locs = 'random', covmodel = cov_bivariate_expo_latDim, covparms = covparms)
locs            <- rbind(process$sim$sim1$locs$locs1, process$sim$sim1$locs$locs2)

nugget          <- 0.00

N               <- 150
limit.range     <- c(0.088, 0.118)

### simulated data #############################################################################################

z               <- process$sim$sim1$y + nugget * rnorm(2 * n)
par(mfrow = c(1, 2))
fields::quilt.plot(locs[, 1][1:n], locs[, 2][1:n], z[1:n], main = "process 1")
fields::quilt.plot(locs[, 1][1:n + n], locs[, 2][1:n + n], z[1:n + n], main = "process 2")
par(mfrow = c(1, 1))

covmat          <- process$sim$sim1$covmat + nugget * diag(2 * n)

### prior distribution #############################################################################################

meanlog         <- log(covparms[2])
sdlog           <- 0.6
logprior        <- function(alpha) dlnorm(alpha, meanlog, sdlog, log = TRUE)

plot(seq(0.01, 1, 0.01), exp(logprior(seq(0.01, 1, 0.01))), type = "l", xlab = "x", ylab = "density", main = "prior distribution")

### posterior over grid #############################################################################################

postDist <- function(m = 30, approx = "cc", N = 200, limit = c(0.088, 0.118), locs, covparms, nugget)
{
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
  
  return(list(alpha = alpha.grid, post = post.norm))
}

### temp #############################################################################################

temp <- postDist(m = 30, approx = "cc", N = N, limit = limit.range, locs = process$sim$sim1$locs, covparms = covparms, nugget = nugget)

matplot(temp$alpha, temp$post, type = 'l', xlim = range(temp$alpha), lty = 1, lwd = 2, xlab = 'range', ylab = 'posterior density')
legend('topright', c('exact', 'approx'), col = 1:2, lty = 1, lwd = 2)

### output #############################################################################################

Sys.time()

# out.m10.b2 <- postDist(m = 10, approx = "b2", N = N, limit = limit.range, locs = process$sim$sim1$locs, covparms = covparms, nugget = nugget)
# out.m30.b2 <- postDist(m = 30, approx = "b2", N = N, limit = limit.range, locs = process$sim$sim1$locs, covparms = covparms, nugget = nugget)
# 
# save(n, covparms, process, nugget, N, limit.range, meanlog, sdlog, logprior, postDist, out.m10.b2, out.m30.b2, file = "R/SCSG2021/bayesian_b2.RData")
# Sys.time()
# 
# out.m10.cc <- postDist(m = 10, approx = "cc", N = N, limit = limit.range, locs = process$sim$sim1$locs, covparms = covparms, nugget = nugget)
# out.m30.cc <- postDist(m = 30, approx = "cc", N = N, limit = limit.range, locs = process$sim$sim1$locs, covparms = covparms, nugget = nugget)
# 
# save(n, covparms, process, nugget, N, limit.range, meanlog, sdlog, logprior, postDist, out.m10.cc, out.m30.cc, file = "R/SCSG2021/bayesian_cc.RData")
# Sys.time()
# 
# out.m10.b1 <- postDist(m = 10, approx = "b1", N = N, limit = limit.range, locs = process$sim$sim1$locs, covparms = covparms, nugget = nugget)
# out.m30.b1 <- postDist(m = 30, approx = "b1", N = N, limit = limit.range, locs = process$sim$sim1$locs, covparms = covparms, nugget = nugget)
# 
# save(n, covparms, process, nugget, N, limit.range, meanlog, sdlog, logprior, postDist, out.m10.b1, out.m30.b1, file = "R/SCSG2021/bayesian_b1.RData")
# Sys.time()
# 
# save(n, covparms, process, nugget, N, limit.range, meanlog, sdlog, logprior, postDist, out.m10.b2, out.m30.b2, out.m10.cc, out.m30.cc, out.m10.b1, out.m30.b1, file = "R/SCSG2021/bayesian.RData")

### visualization #############################################################################################

rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)

load("C:/Users/kmjst/Dropbox/4_PhD_Personal/10_2020W/2_Research/70_12152020/SCSG_final/bayesian/bayesian_b1.RData")
load("C:/Users/kmjst/Dropbox/4_PhD_Personal/10_2020W/2_Research/70_12152020/SCSG_final/bayesian/bayesian_b2.RData")
load("C:/Users/kmjst/Dropbox/4_PhD_Personal/10_2020W/2_Research/70_12152020/SCSG_final/bayesian/bayesian_b3.RData")
load("C:/Users/kmjst/Dropbox/4_PhD_Personal/10_2020W/2_Research/70_12152020/SCSG_final/bayesian/bayesian_b4.RData")
load("C:/Users/kmjst/Dropbox/4_PhD_Personal/10_2020W/2_Research/70_12152020/SCSG_final/bayesian/bayesian_cc.RData")

### alpha
# sum(out.m10.b1$alpha == out.m30.b1$alpha)
# sum(out.m10.b2$alpha == out.m30.b2$alpha)
# sum(out.m10.b3$alpha == out.m30.b3$alpha)
# 
# sum(out.m10.b1$alpha == out.m10.b2$alpha)
# sum(out.m10.b1$alpha == out.m10.b3$alpha)

### true density
# sum((out.m10.b1$post[, 1] - out.m30.b1$post[, 1])^2) 
# sum((out.m10.b2$post[, 1] - out.m30.b2$post[, 1])^2) 
# sum((out.m10.b3$post[, 1] - out.m30.b3$post[, 1])^2) 
# 
# sum((out.m10.b1$post[, 1] - out.m10.b2$post[, 1])^2) 
# sum((out.m10.b1$post[, 1] - out.m10.b3$post[, 1])^2) 

### plot
# matplot(out.m10.b1$alpha, out.m10.b1$post, type = "l")
# matplot(out.m30.b1$alpha, out.m30.b1$post, type = "l")
# 
# matplot(out.m10.b2$alpha, out.m10.b2$post, type = "l")
# matplot(out.m30.b2$alpha, out.m30.b2$post, type = "l")
# 
# matplot(out.m10.b3$alpha, out.m10.b3$post, type = "l")
# matplot(out.m30.b3$alpha, out.m30.b3$post, type = "l")

outFunc <- function(out1, out2, tol = 1e-8)
{
  if(length(out1$alpha) != sum(out1$alpha == out2$alpha)) stop("Check alpha!")
  
  if(sum((out1$post[, 1] - out2$post[, 1])^2) > tol) stop("Check true density!")
  
  output <- data.frame(alpha = out1$alpha, trueval = out1$post[, 1], approx1 = out1$post[, 2], approx2 = out2$post[, 2])
  output <- output %>% tidyr::gather(key = "method", value = "value", trueval, approx1, approx2)
  return(output)
}

out.b1 <- outFunc(out.m10.b1, out.m30.b1)
out.b2 <- outFunc(out.m10.b2, out.m30.b2)
out.b3 <- outFunc(out.m10.b3, out.m30.b3)
out.b4 <- outFunc(out.m10.b4, out.m30.b4)
out.cc <- outFunc(out.m10.cc, out.m30.cc)

# par(mfrow = c(1, 3))
# matplot(out.b1$alpha, select(out.b1, -alpha), type = "l")
# matplot(out.b2$alpha, select(out.b2, -alpha), type = "l")
# matplot(out.b3$alpha, select(out.b3, -alpha), type = "l")
# par(mfrow = c(1, 1))

vis_posterior <- function(output, color = c("#377EB8", "#E41A1C", "gray1"), legend = c("m = 10", "m = 30", "True density"), ylimit = c(0, 0.025), alpha = 0.7, size.point = 4, size.line = 1, size.legend = 14, size.lab = 14, size.text = 12, size.margin = c(5.5, 20, 5.5, 5.5))
{
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

ggsave("R/SCSG2021/posterior_b1.pdf", vis_posterior(out.b1), width = 5.0, height = 5.7)
ggsave("R/SCSG2021/posterior_b2.pdf", vis_posterior(out.b2), width = 5.0, height = 5.7)
ggsave("R/SCSG2021/posterior_b3.pdf", vis_posterior(out.b3), width = 5.0, height = 5.7)
ggsave("R/SCSG2021/posterior_b4.pdf", vis_posterior(out.b4), width = 5.0, height = 5.7)
ggsave("R/SCSG2021/posterior_cc.pdf", vis_posterior(out.cc), width = 5.0, height = 5.7)

tempdat.true <- out.b1 %>% filter(method == "trueval") %>% select(-method) %>% cbind(data.frame(method = rep("trueval", N)))
tempdat.b1 <- out.b1 %>% filter(method == "approx2") %>% select(-method) %>% cbind(data.frame(method = rep("baseline1", N)))
tempdat.b2 <- out.b2 %>% filter(method == "approx2") %>% select(-method) %>% cbind(data.frame(method = rep("baseline2", N)))
tempdat.b3 <- out.b3 %>% filter(method == "approx2") %>% select(-method) %>% cbind(data.frame(method = rep("baseline3", N)))
tempdat.b4 <- out.b4 %>% filter(method == "approx2") %>% select(-method) %>% cbind(data.frame(method = rep("baseline4", N)))
tempdat.cc <- out.cc %>% filter(method == "approx2") %>% select(-method) %>% cbind(data.frame(method = rep("cvecchia", N)))

visout <- rbind(tempdat.true, tempdat.b1, tempdat.b2, tempdat.b3, tempdat.b4, tempdat.cc)

# color = c("#984EA3", "#4DAF4A", "#377EB8", "#FF7F00", "#E41A1C")
# c("S-E-MM + HH-E-NN", "S-E-MM + J-E-NN", "S-E-MM + S-E-NN", "S-E-MM + C-NN", "C-MM + C-NN")
vis_posterior(visout, color = c("gray1", "#984EA3", "#4DAF4A", "#377EB8", "#FF7F00", "#E41A1C"), legend = c("True density", "S-E-MM + HH-E-NN", "S-E-MM + J-E-NN", "S-E-MM + S-E-NN", "S-E-MM + C-NN", "C-MM + C-NN"))

ggsave("R/SCSG2021/posterior_b1_cpt.pdf", vis_posterior(out.b1 %>% filter(method != "approx1"), color = c("#984EA3", "gray1"), legend = c("S-E-MM + HH-E-NN", "true"), ylim = c(0, 0.03)), width = 5, height = 5.7)
ggsave("R/SCSG2021/posterior_b2_cpt.pdf", vis_posterior(out.b2 %>% filter(method != "approx1"), color = c("#4DAF4A", "gray1"), legend = c("S-E-MM + J-E-NN", "true"), ylim = c(0, 0.03)), width = 5, height = 5.7)
ggsave("R/SCSG2021/posterior_b3_cpt.pdf", vis_posterior(out.b3 %>% filter(method != "approx1"), color = c("#377EB8", "gray1"), legend = c("S-E-MM + S-E-NN", "true"), ylim = c(0, 0.03)), width = 5, height = 5.7)
ggsave("R/SCSG2021/posterior_b4_cpt.pdf", vis_posterior(out.b4 %>% filter(method != "approx1"), color = c("#FF7F00", "gray1"), legend = c("S-E-MM + C-NN", "true"), ylim = c(0, 0.03)), width = 5, height = 5.7)
ggsave("R/SCSG2021/posterior_cc_cpt.pdf", vis_posterior(out.cc %>% filter(method != "approx1"), color = c("#E41A1C", "gray1"), legend = c("C-MM + C-NN", "true"), ylim = c(0, 0.03)), width = 5, height = 5.7)

