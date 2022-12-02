####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview:
###
####################################################################################

rm(list = ls())

#########################################################################

set.seed(10112021)

library(dplyr) ; library(ggplot2) ; library(gridExtra) ; library(grid)

#########################################################################

load("DATA/joint_CRCM_NCEP_10152021_prediction.RData")

visout.left   <- mspe.m.varying$joint
visout.right  <- table.runtime[, c(2, 1, 3, 5, 4, 6), drop = FALSE] # CAUTION: I changed the order of approximations!!!

rm(ns_obs, ns_pred, z, z.pred, locs, locs.pred, X.train, X.pred)
rm(output.joint.m.fixed, output.joint.m.varying, output.marginal.m.fixed, output.marginal.m.varying)
rm(mspe.m.fixed, mspe.m.varying, mls.m.fixed, mls.m.varying, table.runtime)

#########################################################################

xticks1     <- ms

vdat1       <- data.frame(index = seq(length(ms)), m = ms, approx_1 = visout.left$b3, approx_2 = visout.left$b2, approx_3 = visout.left$b5, approx_4 = visout.left$b8, approx_5 = visout.left$b7, approx_6 = visout.left$cc) # CAUTION: I changed the order of approximations!!!
vdat1       <- vdat1 %>% select(-index) %>% tidyr::gather(key = "approx", value = "mspe", -m)

#########################################################################

xticks2     <- c(0, 5000, 10000, 15000, 22000)
colnames(visout.right) <- paste0("approx_", 1:6)

vdat2       <- visout.right %>% tidyr::gather(key = "approx", value = "runtime")
vdat2       <- vdat2 %>% mutate(m = rep(ms, 6)) %>% left_join(vdat1, by = c("approx", "m")) %>% select(-m)

####################################################################################

# b2 = S-E-MM + J-E-NN
# b3 = S-E-MM + S-E-NN
# b5 = T-ord + T-NN
# b7 = T-ord + J-C-NN
# b8 = T-ord + S-C-NN
# cc = C-MM + C-NN

legs        <- c(approx_1 = "S-E-MM + S-E-NN", approx_2 = "S-E-MM + J-E-NN", approx_3 = "T-ord + T-NN", approx_4 = "T-ord + S-C-NN", approx_5 = "T-ord + J-C-NN", approx_6 = "C-MM + C-NN")
cols        <- c(approx_1 = "#4DAF4A", approx_2 = "#FF7F00", approx_3 = "#984EA3", approx_4 = "#FFFF33", approx_5 = "#377EB8", approx_6 = "#E41A1C")
shps        <- c(approx_1 = 15, approx_2 = 8, approx_3 = 18, approx_4 = 13, approx_5 = 17, approx_6 = 16)
# ftn         <- log10
xlabs       <- c("increasing m", "runtime (s)")
ylabs       <- c("RMSPE", "")
xlims       <- list(xl1 = NULL, xl2 = NULL)
ylims       <- list(yl1 = NULL, yl2 = NULL)
alpha       <- 0.7
size.point  <- 5
size.tick   <- 1
size.line   <- 1
size.legend <- 18
size.lab    <- 18
size.text   <- 16
size.margin <- c(5.5, 20, 5.5, 5.5)

####################################################################################

ftn1 <- function(x) sqrt(x)
ftn2 <- function(x) log(x)

# plot vis1

plot1     <- vdat1 %>% ggplot(aes(x = m, y = ftn1(mspe), col = approx, shape = approx)) +
  geom_point(size = size.point) +
  geom_line(size = size.line, alpha = alpha) +
  scale_color_manual(values = cols, labels = legs, guide = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_shape_manual(values = shps, labels = legs, guide = guide_legend(nrow = 2, byrow = TRUE)) +
  xlab(xlabs[1]) + ylab(ylabs[1]) + coord_cartesian(xlim = xlims[[1]], ylim = ylims[[1]]) +
  scale_y_continuous(trans='log10', breaks = c(1, 2, 3, 4, 5, 6, 7)) + 
  scale_x_continuous(breaks = xticks1) +
  theme(axis.title.x = element_text(size = size.lab),
        axis.text.x = element_text(size = size.text),
        axis.title.y = element_text(size = size.lab),
        axis.text.y = element_text(size = size.text),
        legend.title = element_blank(),
        legend.text = element_text(size = size.legend, margin = margin(r = 25, unit = 'pt')),
        legend.direction = 'horizontal',
        legend.spacing.x = unit(15, 'pt'),
        plot.margin = unit(c(5.5, 20, 5.5, 5.5), "pt")) # t, r, b, l



plot2     <- vdat2 %>% ggplot(aes(x = runtime, y = ftn1(mspe), col = approx, shape = approx)) +
  geom_point(size = size.point) +
  geom_line(size = size.line, alpha = alpha) +
  # scale_x_continuous(name = vars1, limits = range(xlabel1), breaks = xlabel1) +
  scale_color_manual(values = cols, labels = legs, guide = guide_legend(nrow = 1, byrow = TRUE)) +
  scale_shape_manual(values = shps, labels = legs, guide = guide_legend(nrow = 1, byrow = TRUE)) +
  xlab(xlabs[2]) + ylab(ylabs[2]) + coord_cartesian(xlim = xlims[[2]], ylim = ylims[[2]]) +
  scale_y_continuous(trans='log10', breaks = c(1, 2, 3, 4, 5, 6, 7)) +
  scale_x_continuous(breaks = xticks2) +
  theme(axis.title.x = element_text(size = size.lab),
        axis.text.x = element_text(size = size.text),
        axis.title.y = element_text(size = size.lab),
        axis.text.y = element_text(size = size.text),
        legend.title = element_blank(),
        legend.text = element_text(size = size.legend, margin = margin(r = 25, unit = 'pt')),
        legend.direction = 'horizontal',
        legend.spacing.x = unit(15, 'pt'),
        plot.margin = unit(c(5.5, 20, 5.5, 5.5), "pt")) # t, r, b, l

tmp       <- ggplot_gtable(ggplot_build(plot1))
leg       <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
mylegend  <- tmp$grobs[[leg]]

result    <- grid.arrange(mylegend, arrangeGrob(plot1 + theme(legend.position="none"), plot2 + theme(legend.position="none"), nrow=1), nrow=2, heights=c(1.5, 9.5))

####################################################################################

ggplot2::ggsave("realdata_performance.pdf", result, width = 15.2, height = 5.7)
