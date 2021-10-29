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

load("C:/Users/kmjst/Dropbox/4_PhD_Personal/12_2021Su/1_Research/06292021/output/knownCovparms/simout_wave_05242021.RData")

out     <- output.wave.dc.d1
vdat1   <- out$vars %>% left_join(out$kldiv, by = "index") %>% filter(period == 0.25) %>% select(-period)
vdat2   <- out$vars %>% left_join(out$kldiv, by = "index") %>% filter(m == 30) %>% select(-m)

out     <- output.wave.dc.d2
vdat3   <- out$vars %>% left_join(out$kldiv, by = "index") %>% filter(period == 0.25) %>% select(-period)
vdat4   <- out$vars %>% left_join(out$kldiv, by = "index") %>% filter(m == 30) %>% select(-m)

vdat1      <- vdat1 %>% select(-index) %>% tidyr::gather(key = "approx", value = "kldiv", -m)
vdat2      <- vdat2 %>% select(-index) %>% tidyr::gather(key = "approx", value = "kldiv", -period)
vdat3      <- vdat3 %>% select(-index) %>% tidyr::gather(key = "approx", value = "kldiv", -m)
vdat4      <- vdat4 %>% select(-index) %>% tidyr::gather(key = "approx", value = "kldiv", -period)

####################################################################################

vdat1$scen <- "dimension of domain = 1"
vdat2$scen <- "dimension of domain = 1"
vdat3$scen <- "dimension of domain = 2"
vdat4$scen <- "dimension of domain = 2"

####################################################################################

legs        <- c("E-MM + E-NN", "E-MM + C-NN", "C-MM + E-NN", "C-MM + C-NN   ")
cols        <- c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C")
shps        <- c(18, 15, 17, 16)
ftn         <- log10
xlabs       <- c("increasing m", "increasing period", "increasing m", "increasing period")
ylabs       <- c("log10(KL)", "log10(KL)", "log10(KL)", "log10(KL)")
xlims       <- list(xl1 = NULL, xl2 = NULL, xl3 = NULL, xl4 = NULL)
ylims       <- list(yl1 = ftn(range(c(vdat1$kldiv, vdat2$kldiv, vdat3$kldiv, vdat4$kldiv))), yl2 = ftn(range(c(vdat1$kldiv, vdat2$kldiv, vdat3$kldiv, vdat4$kldiv))), yl3 = ftn(range(c(vdat1$kldiv, vdat2$kldiv, vdat3$kldiv, vdat4$kldiv))), yl4 = ftn(range(c(vdat1$kldiv, vdat2$kldiv, vdat3$kldiv, vdat4$kldiv))))
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
  xlab(xlabs[1]) + ylab(NULL) + coord_cartesian(xlim = xlims[[1]], ylim = ylims[[1]]) +
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

p2 <- vdat2 %>% ggplot(aes(x = period, y = ftn(kldiv), col = approx, shape = approx)) +
  geom_point(size = size.point) + geom_line(size = size.line, alpha = alpha) +
  scale_color_manual(values = cols, labels = legs, guide = guide_legend(nrow = 1, byrow = TRUE)) +
  scale_shape_manual(values = shps, labels = legs, guide = guide_legend(nrow = 1, byrow = TRUE)) +
  xlab(xlabs[2]) + ylab(NULL) + coord_cartesian(xlim = xlims[[2]], ylim = ylims[[2]]) +
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
  xlab(xlabs[3]) + ylab(NULL) + coord_cartesian(xlim = xlims[[3]], ylim = ylims[[3]]) +
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

p4 <- vdat4 %>% ggplot(aes(x = period, y = ftn(kldiv), col = approx, shape = approx)) +
  geom_point(size = size.point) + geom_line(size = size.line, alpha = alpha) +
  scale_color_manual(values = cols, labels = legs, guide = guide_legend(nrow = 1, byrow = TRUE)) +
  scale_shape_manual(values = shps, labels = legs, guide = guide_legend(nrow = 1, byrow = TRUE)) +
  xlab(xlabs[4]) + ylab(NULL) + coord_cartesian(xlim = xlims[[4]], ylim = ylims[[4]]) +
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
                                                p4 + theme(legend.position="none"),
                                                nrow=2), nrow=2, heights=c(1, 10),
                          left = textGrob(ylabs[1], rot = 90, gp = gpar(fontsize = size.lab)))

####################################################################################

ggplot2::ggsave("simulation_supp_periodic_cosine.pdf", result, width = 15.2, height = 10)

####################################################################################

rm(list = ls())

library(dplyr) ; library(ggplot2) ; library(gridExtra) ; library(grid)

####################################################################################

load("C:/Users/kmjst/Dropbox/4_PhD_Personal/12_2021Su/1_Research/06292021/output/knownCovparms/simout_wave_05242021.RData")

out     <- output.wave.ds.d1
vdat1   <- out$vars %>% left_join(out$kldiv, by = "index") %>% filter(period == 0.25) %>% select(-period)
vdat2   <- out$vars %>% left_join(out$kldiv, by = "index") %>% filter(m == 30) %>% select(-m)

out     <- output.wave.ds.d2
vdat3   <- out$vars %>% left_join(out$kldiv, by = "index") %>% filter(period == 0.25) %>% select(-period)
vdat4   <- out$vars %>% left_join(out$kldiv, by = "index") %>% filter(m == 30) %>% select(-m)

vdat1      <- vdat1 %>% select(-index) %>% tidyr::gather(key = "approx", value = "kldiv", -m)
vdat2      <- vdat2 %>% select(-index) %>% tidyr::gather(key = "approx", value = "kldiv", -period)
vdat3      <- vdat3 %>% select(-index) %>% tidyr::gather(key = "approx", value = "kldiv", -m)
vdat4      <- vdat4 %>% select(-index) %>% tidyr::gather(key = "approx", value = "kldiv", -period)

####################################################################################

vdat1$scen <- "dimension of domain = 1"
vdat2$scen <- "dimension of domain = 1"
vdat3$scen <- "dimension of domain = 2"
vdat4$scen <- "dimension of domain = 2"

####################################################################################

legs        <- c("E-MM + E-NN", "E-MM + C-NN", "C-MM + E-NN", "C-MM + C-NN   ")
cols        <- c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C")
shps        <- c(18, 15, 17, 16)
ftn         <- log10
xlabs       <- c("increasing m", "increasing period", "increasing m", "increasing period")
ylabs       <- c("log10(KL)", "log10(KL)", "log10(KL)", "log10(KL)")
xlims       <- list(xl1 = NULL, xl2 = NULL, xl3 = NULL, xl4 = NULL)
ylims       <- list(yl1 = ftn(range(c(vdat1$kldiv, vdat2$kldiv, vdat3$kldiv, vdat4$kldiv))), yl2 = ftn(range(c(vdat1$kldiv, vdat2$kldiv, vdat3$kldiv, vdat4$kldiv))), yl3 = ftn(range(c(vdat1$kldiv, vdat2$kldiv, vdat3$kldiv, vdat4$kldiv))), yl4 = ftn(range(c(vdat1$kldiv, vdat2$kldiv, vdat3$kldiv, vdat4$kldiv))))
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
  xlab(xlabs[1]) + ylab(NULL) + coord_cartesian(xlim = xlims[[1]], ylim = ylims[[1]]) +
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

p2 <- vdat2 %>% ggplot(aes(x = period, y = ftn(kldiv), col = approx, shape = approx)) +
  geom_point(size = size.point) + geom_line(size = size.line, alpha = alpha) +
  scale_color_manual(values = cols, labels = legs, guide = guide_legend(nrow = 1, byrow = TRUE)) +
  scale_shape_manual(values = shps, labels = legs, guide = guide_legend(nrow = 1, byrow = TRUE)) +
  xlab(xlabs[2]) + ylab(NULL) + coord_cartesian(xlim = xlims[[2]], ylim = ylims[[2]]) +
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
  xlab(xlabs[3]) + ylab(NULL) + coord_cartesian(xlim = xlims[[3]], ylim = ylims[[3]]) +
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

p4 <- vdat4 %>% ggplot(aes(x = period, y = ftn(kldiv), col = approx, shape = approx)) +
  geom_point(size = size.point) + geom_line(size = size.line, alpha = alpha) +
  scale_color_manual(values = cols, labels = legs, guide = guide_legend(nrow = 1, byrow = TRUE)) +
  scale_shape_manual(values = shps, labels = legs, guide = guide_legend(nrow = 1, byrow = TRUE)) +
  xlab(xlabs[4]) + ylab(NULL) + coord_cartesian(xlim = xlims[[4]], ylim = ylims[[4]]) +
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
                                                p4 + theme(legend.position="none"),
                                                nrow=2), nrow=2, heights=c(1, 10),
                          left = textGrob(ylabs[1], rot = 90, gp = gpar(fontsize = size.lab)))

####################################################################################

ggplot2::ggsave("simulation_supp_periodic_sine.pdf", result, width = 15.2, height = 10)

####################################################################################

rm(list = ls())

library(dplyr) ; library(ggplot2) ; library(gridExtra) ; library(grid)

####################################################################################

load("C:/Users/kmjst/Dropbox/4_PhD_Personal/12_2021Su/1_Research/06292021/output/knownCovparms/simout_wave_05242021.RData")

out     <- output.wave.bj.d1
vdat1   <- out$vars %>% left_join(out$kldiv, by = "index") %>% filter(period == 0.25) %>% select(-period)
vdat2   <- out$vars %>% left_join(out$kldiv, by = "index") %>% filter(m == 30) %>% select(-m)

out     <- output.wave.bj.d2
vdat3   <- out$vars %>% left_join(out$kldiv, by = "index") %>% filter(period == 0.25) %>% select(-period)
vdat4   <- out$vars %>% left_join(out$kldiv, by = "index") %>% filter(m == 30) %>% select(-m)

vdat1      <- vdat1 %>% select(-index) %>% tidyr::gather(key = "approx", value = "kldiv", -m)
vdat2      <- vdat2 %>% select(-index) %>% tidyr::gather(key = "approx", value = "kldiv", -period)
vdat3      <- vdat3 %>% select(-index) %>% tidyr::gather(key = "approx", value = "kldiv", -m)
vdat4      <- vdat4 %>% select(-index) %>% tidyr::gather(key = "approx", value = "kldiv", -period)

####################################################################################

vdat1$scen <- "dimension of domain = 1"
vdat2$scen <- "dimension of domain = 1"
vdat3$scen <- "dimension of domain = 2"
vdat4$scen <- "dimension of domain = 2"

####################################################################################

legs        <- c("E-MM + E-NN", "E-MM + C-NN", "C-MM + E-NN", "C-MM + C-NN   ")
cols        <- c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C")
shps        <- c(18, 15, 17, 16)
ftn         <- log10
xlabs       <- c("increasing m", "increasing period", "increasing m", "increasing period")
ylabs       <- c("log10(KL)", "log10(KL)", "log10(KL)", "log10(KL)")
xlims       <- list(xl1 = NULL, xl2 = NULL, xl3 = NULL, xl4 = NULL)
ylims       <- list(yl1 = ftn(range(c(vdat1$kldiv, vdat2$kldiv, vdat3$kldiv, vdat4$kldiv))), yl2 = ftn(range(c(vdat1$kldiv, vdat2$kldiv, vdat3$kldiv, vdat4$kldiv))), yl3 = ftn(range(c(vdat1$kldiv, vdat2$kldiv, vdat3$kldiv, vdat4$kldiv))), yl4 = ftn(range(c(vdat1$kldiv, vdat2$kldiv, vdat3$kldiv, vdat4$kldiv))))
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
  xlab(xlabs[1]) + ylab(NULL) + coord_cartesian(xlim = xlims[[1]], ylim = ylims[[1]]) +
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

p2 <- vdat2 %>% ggplot(aes(x = period, y = ftn(kldiv), col = approx, shape = approx)) +
  geom_point(size = size.point) + geom_line(size = size.line, alpha = alpha) +
  scale_color_manual(values = cols, labels = legs, guide = guide_legend(nrow = 1, byrow = TRUE)) +
  scale_shape_manual(values = shps, labels = legs, guide = guide_legend(nrow = 1, byrow = TRUE)) +
  xlab(xlabs[2]) + ylab(NULL) + coord_cartesian(xlim = xlims[[2]], ylim = ylims[[2]]) +
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
  xlab(xlabs[3]) + ylab(NULL) + coord_cartesian(xlim = xlims[[3]], ylim = ylims[[3]]) +
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

p4 <- vdat4 %>% ggplot(aes(x = period, y = ftn(kldiv), col = approx, shape = approx)) +
  geom_point(size = size.point) + geom_line(size = size.line, alpha = alpha) +
  scale_color_manual(values = cols, labels = legs, guide = guide_legend(nrow = 1, byrow = TRUE)) +
  scale_shape_manual(values = shps, labels = legs, guide = guide_legend(nrow = 1, byrow = TRUE)) +
  xlab(xlabs[4]) + ylab(NULL) + coord_cartesian(xlim = xlims[[4]], ylim = ylims[[4]]) +
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
                                                p4 + theme(legend.position="none"),
                                                nrow=2), nrow=2, heights=c(1, 10),
                          left = textGrob(ylabs[1], rot = 90, gp = gpar(fontsize = size.lab)))

####################################################################################

ggplot2::ggsave("simulation_supp_periodic_bessel.pdf", result, width = 15.2, height = 10)
