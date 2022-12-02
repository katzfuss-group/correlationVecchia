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

load("simout_aniso_05242021.RData")

out     <- output.aniso
vdat1   <- out$vars %>% left_join(out$kldiv, by = "index") %>% filter(a == 10) %>% select(-a)
vdat2   <- out$vars %>% left_join(out$kldiv, by = "index") %>% filter(m == 30) %>% select(-m)

xticks1 <- cand.m
xticks2 <- cand.a

rm(out, output.aniso, cand.a, cand.m, nsim)

load("simout_nonst_05242021.RData")

out1    <- output.nonst1
out2    <- output.nonst4
vdat3   <- out1$vars %>% left_join(out1$kldiv, by = "index")
vdat4   <- out2$vars %>% left_join(out2$kldiv, by = "index")

xticks3 <- cand.m
xticks4 <- cand.m

rm(out1, out2, output.nonst1, output.nonst2, output.nonst3, output.nonst4, nsim)

####################################################################################

vdat1      <- vdat1 %>% select(-index) %>% tidyr::gather(key = "approx", value = "kldiv", -m)
vdat2      <- vdat2 %>% select(-index) %>% tidyr::gather(key = "approx", value = "kldiv", -a)
vdat3      <- vdat3 %>% select(-index) %>% tidyr::gather(key = "approx", value = "kldiv", -m)
vdat4      <- vdat4 %>% select(-index) %>% tidyr::gather(key = "approx", value = "kldiv", -m)

vdat1$scen <- "anisotropic"
vdat2$scen <- "anisotropic"
vdat3$scen <- "varying smoothness"
vdat4$scen <- "varying rotation"

####################################################################################

legs        <- c("E-MM + E-NN", "E-MM + C-NN", "C-MM + E-NN", "C-MM + C-NN", "X-ord + E-NN", "Y-ord + E-NN   ")
cols        <- c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "#FF7F00", "#FFFF33")
shps        <- c(18, 15, 17, 16, 8, 13)
ftn         <- log10
xlabs       <- c("increasing m", "increasing a", "increasing m", "increasing m")
ylabs       <- c("log10(KL)", "log10(KL)", "log10(KL)", "log10(KL)")
xlims       <- list(xl1 = NULL, xl2 = NULL, xl3 = NULL, xl4 = NULL)
ylims       <- list(yl1 = ftn(range(c(vdat1$kldiv, vdat2$kldiv))), yl2 = ftn(range(c(vdat1$kldiv, vdat2$kldiv))), yl3 = NULL, yl4 = NULL)
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
  scale_x_continuous(breaks = xticks1) +
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

p2 <- vdat2 %>% ggplot(aes(x = a, y = ftn(kldiv), col = approx, shape = approx)) +
  geom_point(size = size.point) + geom_line(size = size.line, alpha = alpha) +
  scale_color_manual(values = cols, labels = legs, guide = guide_legend(nrow = 1, byrow = TRUE)) +
  scale_shape_manual(values = shps, labels = legs, guide = guide_legend(nrow = 1, byrow = TRUE)) +
  xlab(xlabs[2]) + ylab(NULL) + coord_cartesian(xlim = xlims[[2]], ylim = ylims[[2]]) +
  scale_x_continuous(breaks = xticks2) +
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
  scale_x_continuous(breaks = xticks3) +
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

p4 <- vdat4 %>% ggplot(aes(x = m, y = ftn(kldiv), col = approx, shape = approx)) +
  geom_point(size = size.point) + geom_line(size = size.line, alpha = alpha) +
  scale_color_manual(values = cols, labels = legs, guide = guide_legend(nrow = 1, byrow = TRUE)) +
  scale_shape_manual(values = shps, labels = legs, guide = guide_legend(nrow = 1, byrow = TRUE)) +
  xlab(xlabs[4]) + ylab(NULL) + coord_cartesian(xlim = xlims[[4]], ylim = ylims[[4]]) +
  scale_x_continuous(breaks = xticks4) +
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

ggplot2::ggsave("simulation_anisononst.pdf", result, width = 15.2, height = 10)
