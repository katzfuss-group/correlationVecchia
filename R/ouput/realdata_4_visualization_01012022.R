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

xticks <- ms

vdat1 <- data.frame(index = seq(length(ms)), m = ms, approx_1 = mspe.b5, approx_2 = mspe.b2, approx_3 = mspe.cc)
vdat2 <- data.frame(index = seq(length(ms)), m = ms, approx_1 = mls.b5, approx_2 = mls.b2, approx_3 = mls.cc)

vdat1      <- vdat1 %>% select(-index) %>% tidyr::gather(key = "approx", value = "kldiv", -m)
vdat2      <- vdat2 %>% select(-index) %>% tidyr::gather(key = "approx", value = "kldiv", -m)

####################################################################################

legs        <- c("T-ord + T-NN", "S-E-MM + J-E-NN", "C-MM + C-NN")
cols        <- c("#984EA3", "#4DAF4A", "#E41A1C")
shps        <- c(18, 15, 16)
# ftn         <- log10
xlabs       <- c("increasing m", "increasing m")
ylabs       <- c("RMSPE", "logarithmic score")
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

# plot vis1
ftn1 <- function(x) sqrt(x)

plot1     <- vdat1 %>% ggplot(aes(x = m, y = ftn1(kldiv), col = approx, shape = approx)) +
  geom_point(size = size.point) +
  geom_line(size = size.line, alpha = alpha) +
  # scale_x_continuous(name = vars1, limits = range(xlabel1), breaks = xlabel1) +
  scale_color_manual(values = cols, labels = legs, guide = guide_legend(nrow = 1, byrow = TRUE)) +
  scale_shape_manual(values = shps, labels = legs, guide = guide_legend(nrow = 1, byrow = TRUE)) +
  xlab(xlabs[1]) + ylab(ylabs[1]) + coord_cartesian(xlim = xlims[[1]], ylim = ylims[[1]]) +
  scale_y_continuous(trans='log10', breaks = c(1, 2, 3, 4, 5, 6)) + # ylim(0.75, 7.25) +
  scale_x_continuous(breaks = xticks) +
  theme(axis.title.x = element_text(size = size.lab),
        axis.text.x = element_text(size = size.text),
        axis.title.y = element_text(size = size.lab),
        axis.text.y = element_text(size = size.text),
        legend.title = element_blank(),
        legend.text = element_text(size = size.legend, margin = margin(r = 25, unit = 'pt')),
        legend.direction = 'horizontal',
        legend.spacing.x = unit(15, 'pt'),
        plot.margin = unit(c(5.5, 20, 5.5, 5.5), "pt")) # t, r, b, l

ftn2 <- function(x) -x

plot2     <- vdat2 %>% ggplot(aes(x = m, y = ftn2(kldiv), col = approx, shape = approx)) +
  geom_point(size = size.point) +
  geom_line(size = size.line, alpha = alpha) +
  # scale_x_continuous(name = vars1, limits = range(xlabel1), breaks = xlabel1) +
  scale_color_manual(values = cols, labels = legs, guide = guide_legend(nrow = 1, byrow = TRUE)) +
  scale_shape_manual(values = shps, labels = legs, guide = guide_legend(nrow = 1, byrow = TRUE)) +
  xlab(xlabs[2]) + ylab(ylabs[2]) + coord_cartesian(xlim = xlims[[2]], ylim = ylims[[2]]) +
  scale_y_continuous(trans='log10', breaks = c(1.5, 2.0, 2.5, 3.0, 3.5, 4.0)) +
  scale_x_continuous(breaks = xticks) +
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

result    <- grid.arrange(mylegend, arrangeGrob(plot1 + theme(legend.position="none"), plot2 + theme(legend.position="none"), nrow=1), nrow=2, heights=c(1, 10))

####################################################################################

ggplot2::ggsave("realdata_performance.pdf", result, width = 15.2, height = 5.7)
