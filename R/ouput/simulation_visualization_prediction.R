####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview:
###
####################################################################################

rm(list = ls())

library(dplyr) ; library(ggplot2) ; library(gridExtra) ; library(grid)

####################################################################################

load("simout_prediction_allrandom_05242021.RData")

rm(nsim, n, n.pred, m, t, nugget, covparms)

load("simout_prediction_monitoring_05242021.RData")

rm(nsim, n, n.pred, m, t, nugget, covparms)

load("simout_prediction_satellite_05242021.RData")

rm(nsim, n, n.pred, m, t, nugget, covparms)

####################################################################################

vdat1 <- data.frame(index = seq(length(output.sptm.pred.random$setting$m)), m = output.sptm.pred.random$setting$m, approx_1 = output.sptm.pred.random$output$logscore[[1]], approx_2 = output.sptm.pred.random$output$logscore[[2]], approx_3 = output.sptm.pred.random$output$logscore[[3]], approx_4 = output.sptm.pred.random$output$logscore[[4]])
vdat2 <- data.frame(index = seq(length(output.sptm.pred.monitoring$setting$m)), m = output.sptm.pred.monitoring$setting$m, approx_1 = output.sptm.pred.monitoring$output$logscore[[1]], approx_2 = output.sptm.pred.monitoring$output$logscore[[2]], approx_3 = output.sptm.pred.monitoring$output$logscore[[3]], approx_4 = output.sptm.pred.monitoring$output$logscore[[4]])
vdat3 <- data.frame(index = seq(length(output.sptm.pred.satellite$setting$m)), m = output.sptm.pred.satellite$setting$m, approx_1 = output.sptm.pred.satellite$output$logscore[[1]], approx_2 = output.sptm.pred.satellite$output$logscore[[2]], approx_3 = output.sptm.pred.satellite$output$logscore[[3]], approx_4 = output.sptm.pred.satellite$output$logscore[[4]])

vdat1   <- vdat1 %>% filter(m != 3)
vdat2   <- vdat2 %>% filter(m != 3)
vdat3   <- vdat3 %>% filter(m != 3)

vdat1      <- vdat1 %>% select(-index) %>% tidyr::gather(key = "approx", value = "kldiv", -m)
vdat2      <- vdat2 %>% select(-index) %>% tidyr::gather(key = "approx", value = "kldiv", -m)
vdat3      <- vdat3 %>% select(-index) %>% tidyr::gather(key = "approx", value = "kldiv", -m)

vdat1$kldiv <- -vdat1$kldiv
vdat2$kldiv <- -vdat2$kldiv
vdat3$kldiv <- -vdat3$kldiv

rm(output.sptm.pred.random, output.sptm.pred.monitoring, output.sptm.pred.satellite)

vdat1$scen <- "random"
vdat2$scen <- "station"
vdat3$scen <- "satellite"

####################################################################################

legs        <- c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN")
cols        <- c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C")
shps        <- c(18, 15, 17, 16)
ftn         <- log10
xlabs       <- c("increasing m", "increasing m", "increasing m")
ylabs       <- c("log10(logscore)", "log10(logscore)", "log10(logscore)")
xlims       <- list(xl1 = NULL, xl2 = NULL, xl3 = NULL)
ylims       <- list(yl1 = NULL, yl2 = NULL, yl3 = NULL)
alpha       <- 0.7
size.point  <- 5
size.tick   <- 1
size.line   <- 1
size.legend <- 18
size.lab    <- 18
size.text   <- 16
size.margin <- c(5.5, 20, 5.5, 5.5)

####################################################################################

p1 <- vdat1 %>% ggplot(aes(x = m, y = ftn(kldiv), col = approx, shape = approx)) +
  geom_point(size = size.point) + geom_line(size = size.line, alpha = alpha) +
  scale_color_manual(values = cols, labels = legs, guide = guide_legend(nrow = 1, byrow = TRUE)) +
  scale_shape_manual(values = shps, labels = legs, guide = guide_legend(nrow = 1, byrow = TRUE)) +
  xlab(NULL) + ylab(NULL) + coord_cartesian(xlim = xlims[[1]], ylim = ylims[[1]]) +
  theme(axis.title.x = element_text(size = size.lab),
        axis.text.x = element_text(size = size.text),
        axis.title.y = element_text(size = size.lab),
        axis.text.y = element_text(size = size.text),
        strip.text = element_text(size = size.lab),
        legend.title = element_blank(),
        legend.text = element_text(size = size.legend, margin = margin(r = 25, unit = 'pt')),
        legend.direction = 'horizontal',
        legend.spacing.x = unit(15, 'pt'),
        plot.margin = unit(size.margin, "pt")) +
  facet_wrap(.~scen, scales = "free")

p2 <- vdat2 %>% ggplot(aes(x = m, y = ftn(kldiv), col = approx, shape = approx)) +
  geom_point(size = size.point) + geom_line(size = size.line, alpha = alpha) +
  scale_color_manual(values = cols, labels = legs, guide = guide_legend(nrow = 1, byrow = TRUE)) +
  scale_shape_manual(values = shps, labels = legs, guide = guide_legend(nrow = 1, byrow = TRUE)) +
  xlab(NULL) + ylab(NULL) + coord_cartesian(xlim = xlims[[2]], ylim = ylims[[2]]) +
  theme(axis.title.x = element_text(size = size.lab),
        axis.text.x = element_text(size = size.text),
        axis.title.y = element_text(size = size.lab),
        axis.text.y = element_text(size = size.text),
        strip.text = element_text(size = size.lab),
        legend.title = element_blank(),
        legend.text = element_text(size = size.legend, margin = margin(r = 25, unit = 'pt')),
        legend.direction = 'horizontal',
        legend.spacing.x = unit(15, 'pt'),
        plot.margin = unit(size.margin, "pt")) +
  facet_wrap(.~scen, scales = "free")

p3 <- vdat3 %>% ggplot(aes(x = m, y = ftn(kldiv), col = approx, shape = approx)) +
  geom_point(size = size.point) + geom_line(size = size.line, alpha = alpha) +
  scale_color_manual(values = cols, labels = legs, guide = guide_legend(nrow = 1, byrow = TRUE)) +
  scale_shape_manual(values = shps, labels = legs, guide = guide_legend(nrow = 1, byrow = TRUE)) +
  xlab(NULL) + ylab(NULL) + coord_cartesian(xlim = xlims[[3]], ylim = ylims[[3]]) +
  theme(axis.title.x = element_text(size = size.lab),
        axis.text.x = element_text(size = size.text),
        axis.title.y = element_text(size = size.lab),
        axis.text.y = element_text(size = size.text),
        strip.text = element_text(size = size.lab),
        legend.title = element_blank(),
        legend.text = element_text(size = size.legend, margin = margin(r = 25, unit = 'pt')),
        legend.direction = 'horizontal',
        legend.spacing.x = unit(15, 'pt'),
        plot.margin = unit(size.margin, "pt")) +
  facet_wrap(.~scen, scales = "free")

####################################################################################

tmp       <- ggplot_gtable(ggplot_build(p1))
leg       <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
mylegend  <- tmp$grobs[[leg]]

### Merge the two plots ###

result    <- grid.arrange(mylegend, arrangeGrob(p1 + theme(legend.position="none"),
                                                p2 + theme(legend.position="none"),
                                                p3 + theme(legend.position="none"),
                                                nrow=1), nrow=2, heights=c(1, 10),
                          left = textGrob(ylabs[1], rot = 90, gp = gpar(fontsize = size.lab)),
                          bottom = textGrob(xlabs[1], gp = gpar(fontsize = size.lab)))

####################################################################################

ggplot2::ggsave("simulation_spacetime_prediction.pdf", result, width = 15.2, height = 5.7)
