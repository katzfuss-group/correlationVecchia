####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: 
###
####################################################################################

rm(list = ls())
set.seed(12122020)

library(correlationVecchia)
library(dplyr)
library(foreach)
library(ggplot2)
library(gridExtra)

### setting #############################################################################################

nsim    <- 10
cand.m  <- c(3, 5, 10, 15, 20, 25, 30, 35, 40, 45)

### visualization #############################################################################################

vis_SCSG <- function(vdat1, legend, color, shape, alpha = 0.7, size.point = 4, size.line = 1, size.legend = 14, size.lab = 14, size.text = 12, size.margin = c(5.5, 20, 5.5, 5.5))
{
  vars1     <- vdat1 %>% select(-index, -starts_with("approx")) %>% colnames()
  if(length(vars1) > 1) stop("Please check the argument vdat1.")
  
  xlabel1   <- vdat1 %>% pull(vars1) %>% unique() %>% sort()
  vis1      <- vdat1 %>% tidyr::gather(key = "approx", value = "kldiv", -index, -vars1)
  
  plot1 <- vis1 %>% ggplot(aes(x = get(vars1), y = log10(kldiv), col = approx, shape = approx)) + 
    geom_point(size = size.point) + 
    geom_line(size = size.line, alpha = alpha) + 
    scale_x_continuous(name = vars1, limits = range(xlabel1), breaks = xlabel1) +
    scale_color_manual(values = color, labels = legend, guide = guide_legend(nrow = 2, byrow = TRUE)) + 
    scale_shape_manual(values = shape, labels = legend, guide = guide_legend(nrow = 2, byrow = TRUE)) +
    xlab(vars1) + ylab('log10(KL)') + 
    theme(axis.title.x = element_text(size = size.lab), 
          axis.text.x = element_text(size = size.text), 
          axis.title.y = element_text(size = size.lab), 
          axis.text.y = element_text(size = size.text), 
          legend.title = element_blank(), 
          legend.position="top",
          legend.text = element_text(size = size.legend, margin = margin(r = 25, unit = 'pt')), 
          legend.direction = 'horizontal', 
          legend.spacing.x = unit(15, 'pt'),
          plot.margin = unit(size.margin, "pt"))
  
  return(plot1)
}

### anisotropic case #############################################################################################

cand.a  <- c(10, 20)

# legend  <- c("E-Maxmin + E-NN", "E-Maxmin + C-NN", "C-Maxmin + E-NN", "C-Maxmin + C-NN", "X-Coord + E-NN", "Y-Coord + E-NN")
# color   <- c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "#FF7F00", "#FFFF33")
# shape   <- c(18, 15, 17, 16, 8, 13)

out <- parallel_simulate_anisotropic_knownCovparms(cand.m = cand.m, cand.a = cand.a, nsim = nsim, n = 30^2, d = 2, corr.dist = "sqrt(1-abs(rho))", covmodel = cov_expo_aniso, covparms = c(1, 0.1), method.locs = "random", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)

save(out, file = "R/SCSG2021/aniso.RData")

vdat1   <- out$vars %>% left_join(out$kldiv, by = "index") %>% filter(a == 10) %>% select(-a)
vis     <- vis_SCSG(vdat1 = vdat1, legend = c("E-MM + E-NN", "E-MM + C-NN", "C-MM + E-NN", "C-MM + C-NN", "X-ord + E-NN", "Y-ord + E-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "#FF7F00", "#FFFF33"), shape = c(18, 15, 17, 16, 8, 13))

ggplot2::ggsave("R/SCSG2021/aniso.pdf", vis, width = 9.0, height = 5.7)

Sys.time()

### multivariate case #############################################################################################

cand.d  <- c(0.2, 0.4)

# legend  <- c("Baseline 1", "Baseline 2", "Baseline 3", "Baseline 4", "C-Maxmin + C-NN") 
# color   <- c("#984EA3", "#4DAF4A", "#377EB8", "#FF7F00", "#E41A1C")
# shape   <- c(18, 15, 17, 8, 16)

out1 <- parallel_simulate_multivariate_knownCovparms(cand.m = cand.m, cand.d = cand.d, nsim = nsim, n = 20^2, d = 2, p = 2, corr.dist = "sqrt(1-abs(rho))", covmodel = cov_bivariate_expo_latDim, covparms = c(1, 0.1), method.locs = "random", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)

save(out1, file = "R/SCSG2021/multi.RData")

vdat1     <- out1$vars %>% left_join(out1$kldiv, by = "index") %>% filter(d == 0.4) %>% select(-d)
vis_ran   <- vis_SCSG(vdat1 = vdat1, legend = c("S-E-MM + HH-E-NN", "S-E-MM + J-E-NN", "S-E-MM + S-E-NN", "S-E-MM + C-NN", "C-MM + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#FF7F00", "#E41A1C"), shape = c(18, 15, 17, 8, 16))
# vis_ran   <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("Baseline 1", "Baseline 2", "Baseline 3", "Baseline 4", "C-Maxmin + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#FF7F00", "#E41A1C"), shape = c(18, 15, 17, 8, 16))

ggplot2::ggsave("R/SCSG2021/multi.pdf", vis_ran, width = 9.0, height = 5.7)

Sys.time()

### spacetime case #############################################################################################

covparms  <- c(1, 0.1, 0.2, 0.5, 0) # var, range_space, range_time, nu, nugget
cov_spacetime <- function(locs, covparms) GpGp::matern_spacetime(covparms = covparms, locs = locs)

# legend  <- c("Baseline 1", "Baseline 2", "Baseline 3", "C-Maxmin + C-NN") 
# color   <- c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C")
# shape   <- c(18, 15, 17, 16)

gen22 <- generate_gp_spacetime(nsim = 2, n = 25, d = 2, t.len = 36, method.locs = "space.random.time.grid", covmodel = cov_spacetime, covparms = covparms, method.modify = "eigen-I", pivot = FALSE, tol = 1e-3)

gen24 <- generate_gp_spacetime(nsim = 2, n = 900, d = 2, t.len = 1, method.locs = "satellite", covmodel = cov_spacetime, covparms = covparms)

fields::quilt.plot(gen22$sim$sim1$locs[, 1], gen22$sim$sim1$locs[, 2], gen22$sim$sim1$y)
plot(gen22$sim$sim1$locs[, 3], gen22$sim$sim1$y)

fields::quilt.plot(gen24$sim$sim2$locs[, 1], gen24$sim$sim2$locs[, 2], gen24$sim$sim2$y)
plot(gen24$sim$sim1$locs[, 3], gen24$sim$sim1$y, type = 'o', col = as.factor(gen24$sim$sim1$locs[, 2]))

out22 <- parallel_simulate_spacetime_knownCovparms(cand.m = cand.m, nsim = nsim, n = 25, d = 2, t.len = 30, corr.dist = "sqrt(1-abs(rho))", covmodel = cov_spacetime, covparms = covparms, method.locs = "space.random.time.grid", method.modify = "eigen-I", pivot = FALSE, tol = 1e-4, ncores = NULL)

out24 <- parallel_simulate_spacetime_knownCovparms(cand.m = cand.m, nsim = nsim, n = 30^2, d = 2, t.len = 1, corr.dist = "sqrt(1-abs(rho))", covmodel = cov_spacetime, covparms = covparms, method.locs = "satellite", method.modify = "eigen-I", pivot = FALSE, tol = 1e-4, ncores = NULL)

save(out22, out24, file = "R/SCSG2021/spti.RData")

vdat1     <- out22$vars %>% left_join(out22$kldiv, by = "index")
vdat2     <- out24$vars %>% left_join(out24$kldiv, by = "index")

vis_spti1 <- vis_SCSG(vdat1 = vdat1, legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), shape = c(18, 15, 17, 16))
vis_spti2 <- vis_SCSG(vdat1 = vdat2, legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), shape = c(18, 15, 17, 16))

ggplot2::ggsave("R/SCSG2021/spti1.pdf", vis_spti1, width = 9.0, height = 5.7)
ggplot2::ggsave("R/SCSG2021/spti2.pdf", vis_spti2, width = 9.0, height = 5.7)

Sys.time()

