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

load("C:/Users/kmjst/Dropbox/4_PhD_Personal/12_2021Su/1_Research/06292021/output/knownCovparms/simout_biv_05242021.RData")

out1      <- output.biv.random
vdat1     <- out1$vars %>% left_join(out1$kldiv, by = "index") %>% filter(d == 0.4) %>% select(-d)
vdat2     <- out1$vars %>% left_join(out1$kldiv, by = "index") %>% filter(m == 20) %>% select(-m)

rm(out1, output.biv.random, output.biv.overlap, cand.m, cand.d, nsim)

load("C:/Users/kmjst/Dropbox/4_PhD_Personal/12_2021Su/1_Research/06292021/output/knownCovparms/simout_triv_05242021.RData")

out1      <- output.triv.random
vdat3     <- out1$vars %>% left_join(out1$kldiv, by = "index") %>% filter(d == 0.4) %>% select(-d)
vdat4     <- out1$vars %>% left_join(out1$kldiv, by = "index") %>% filter(m == 20) %>% select(-m)

rm(out1, output.triv.random, output.triv.overlap, cand.m, cand.d, nsim)

####################################################################################

vdat1      <- vdat1 %>% select(-index) %>% tidyr::gather(key = "approx", value = "kldiv", -m)
vdat2      <- vdat2 %>% select(-index) %>% tidyr::gather(key = "approx", value = "kldiv", -d)
vdat3      <- vdat3 %>% select(-index) %>% tidyr::gather(key = "approx", value = "kldiv", -m)
vdat4      <- vdat4 %>% select(-index) %>% tidyr::gather(key = "approx", value = "kldiv", -d)

vdat1$scen <- "bivariate"
vdat2$scen <- "bivariate"
vdat3$scen <- "trivariate"
vdat4$scen <- "trivariate"

####################################################################################

legs        <- c("S-E-MM + D-E-NN", "S-E-MM + J-E-NN", "S-E-MM + S-E-NN", "S-E-MM + C-NN", "C-MM + C-NN   ")
cols        <- c("#984EA3", "#4DAF4A", "#377EB8", "#FF7F00", "#E41A1C")
shps        <- c(18, 15, 17, 8, 16)
ftn         <- log10
xlabs       <- c("increasing m", "increasing distance in latent dimension", "increasing m", "increasing distance in latent dimension")
ylabs       <- c("log10(KL)", "log10(KL)", "log10(KL)", "log10(KL)")
xlims       <- list(xl1 = NULL, xl2 = NULL, xl3 = NULL, xl4 = NULL)
ylims       <- list(yl1 = ftn(range(c(vdat1$kldiv, vdat2$kldiv))), yl2 = ftn(range(c(vdat1$kldiv, vdat2$kldiv))), yl3 = ftn(range(c(vdat3$kldiv, vdat4$kldiv))), yl4 = ftn(range(c(vdat3$kldiv, vdat4$kldiv))))
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

p2 <- vdat2 %>% ggplot(aes(x = d, y = ftn(kldiv), col = approx, shape = approx)) +
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

p4 <- vdat4 %>% ggplot(aes(x = d, y = ftn(kldiv), col = approx, shape = approx)) +
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

ggplot2::ggsave("simulation_multivariate.pdf", result, width = 15.2, height = 10)
