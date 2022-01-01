####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview:
###
####################################################################################

library(correlationVecchia) ; library(microbenchmark)

library(dplyr) ; library(ggplot2) ; library(gridExtra) ; library(grid)

####################################################################################

rm(list = ls())

###

load(file = "comparison_m_and_rho_1.RData")

###

vdat1       <- data.frame(x1 = locs[, 1], x2 = locs[, 2])
vdat2       <- data.frame(index = seq(length(ratio.m)), mtype = ratio.m, rhotype = ratio.rho)
vdat3       <- data.frame(index = seq(length(sizes.rho)), mtype = rep(m, times = length(sizes.rho)), rhotype = sizes.rho)

vdat3[, "mtype"][seq(m)] <- seq(m) - 1

# vdat3[, "mtype"]    <- vdat3[, "mtype"] - 1
vdat3[, "rhotype"]  <- vdat3[, "rhotype"] - 1

vdat2       <- vdat2 %>% tidyr::gather(key = "approx", value = "kldiv", -index)
vdat3       <- vdat3 %>% tidyr::gather(key = "approx", value = "kldiv", -index)

###

# legs        <- c(paste0(expression(m),"-type EVecchia"), paste0(expression(rho),"-type EVecchia"))
legs        <- c("EVecchia with size-based NN", "EVecchia with radius-based NN")
cols        <- c("#E41A1C", "#377EB8")
ftn         <- function(x) log(x)
alpha       <- 0.7
size.point  <- 1.5
size.tick   <- 1
size.line   <- 1
size.legend <- 18
size.lab    <- 18
size.text   <- 16
size.margin <- c(5.5, 20, 5.5, 5.5)

###

p1 <- vdat1 %>% ggplot(aes(x = x1, y = x2)) +
  geom_point(size = 0.4) + xlab("x1") + ylab("x2") +
  theme(axis.title.x = element_text(size = size.lab),
        axis.text.x = element_text(size = size.text),
        axis.title.y = element_text(size = size.lab),
        axis.text.y = element_text(size = size.text),
        strip.text = element_text(size = size.lab),
        legend.title = element_blank(),
        legend.text = element_text(size = size.legend, margin = margin(r = 25, unit = 'pt')),
        legend.direction = 'horizontal',
        legend.spacing.x = unit(15, 'pt'),
        plot.margin = unit(size.margin, "pt"))

shps        <- c(15, 16)

p2 <- vdat2 %>% ggplot(aes(x = index, y = ftn(kldiv), col = approx, shape = approx)) +
  geom_point(size = size.point) + geom_line(size = size.line, alpha = alpha) +
  scale_x_continuous(breaks = c(1, 100, 200, 300, 400)) +
  # scale_y_continuous(trans='log10') +
  scale_color_manual(values = cols, labels = legs, guide = guide_legend(nrow = 1, byrow = TRUE)) +
  scale_shape_manual(values = shps, labels = legs, guide = guide_legend(nrow = 1, byrow = TRUE)) +
  xlab("index") + ylab("log( variance ratio )") + coord_cartesian(xlim = c(1, nrow(locs))) +
  theme(axis.title.x = element_text(size = size.lab),
        axis.text.x = element_text(size = size.text),
        axis.title.y = element_text(size = size.lab),
        axis.text.y = element_text(size = size.text),
        strip.text = element_text(size = size.lab),
        legend.title = element_blank(),
        legend.text = element_text(size = size.legend, margin = margin(r = 25, unit = 'pt')),
        legend.direction = 'horizontal',
        legend.spacing.x = unit(15, 'pt'),
        plot.margin = unit(size.margin, "pt"))

p3 <- vdat3 %>% ggplot(aes(x = index, y = kldiv, col = approx, shape = approx)) +
  geom_point(size = size.point) + geom_line(size = size.line, alpha = alpha) +
  scale_x_continuous(breaks = c(1, 100, 200, 300, 400)) +
  scale_color_manual(values = cols, labels = legs, guide = guide_legend(nrow = 1, byrow = TRUE)) +
  scale_shape_manual(values = shps, labels = legs, guide = guide_legend(nrow = 1, byrow = TRUE)) +
  xlab("index") + ylab("size of conditioning set") + coord_cartesian(xlim = c(1, nrow(locs))) +
  theme(axis.title.x = element_text(size = size.lab),
        axis.text.x = element_text(size = size.text),
        axis.title.y = element_text(size = size.lab),
        axis.text.y = element_text(size = size.text),
        strip.text = element_text(size = size.lab),
        legend.title = element_blank(),
        legend.text = element_text(size = size.legend, margin = margin(r = 25, unit = 'pt')),
        legend.direction = 'horizontal',
        legend.spacing.x = unit(15, 'pt'),
        plot.margin = unit(size.margin, "pt"))

###

tmp       <- ggplot_gtable(ggplot_build(p3))
leg       <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
mylegend  <- tmp$grobs[[leg]]

### Merge the two plots ###

result    <- grid.arrange(mylegend, arrangeGrob(p1 + theme(legend.position="none"),
                                                p3 + theme(legend.position="none"),
                                                p2 + theme(legend.position="none"),
                                                nrow=1), nrow=2, heights=c(1, 10))
###

ggplot2::ggsave("simulation_supp_m_and_rho_1.pdf", result, width = 15.2, height = 5.7)

####################################################################################

load(file = "comparison_m_and_rho_2.RData")

###

vdat1       <- data.frame(x1 = locs[, 1], x2 = locs[, 2])
vdat2       <- data.frame(index = seq(length(ratio.m)), mtype = ratio.m, rhotype = ratio.rho)
vdat3       <- data.frame(index = seq(length(sizes.rho)), mtype = rep(m, times = length(sizes.rho)), rhotype = sizes.rho)

vdat3[, "mtype"][seq(m)] <- seq(m) - 1

# vdat3[, "mtype"]    <- vdat3[, "mtype"] - 1
vdat3[, "rhotype"]  <- vdat3[, "rhotype"] - 1

vdat2       <- vdat2 %>% tidyr::gather(key = "approx", value = "kldiv", -index)
vdat3       <- vdat3 %>% tidyr::gather(key = "approx", value = "kldiv", -index)

###

# legs        <- c(paste0(expression(m),"-type EVecchia"), paste0(expression(rho),"-type EVecchia"))
legs        <- c("EVecchia with size-based NN", "EVecchia with radius-based NN")
cols        <- c("#E41A1C", "#377EB8")
ftn         <- function(x) log(x)
alpha       <- 0.7
size.point  <- 1.5
size.tick   <- 1
size.line   <- 1
size.legend <- 18
size.lab    <- 18
size.text   <- 16
size.margin <- c(5.5, 20, 5.5, 5.5)

###

p1 <- vdat1 %>% ggplot(aes(x = x1, y = x2)) +
  geom_point(size = 0.4) + xlab("x1") + ylab("x2") +
  theme(axis.title.x = element_text(size = size.lab),
        axis.text.x = element_text(size = size.text),
        axis.title.y = element_text(size = size.lab),
        axis.text.y = element_text(size = size.text),
        strip.text = element_text(size = size.lab),
        legend.title = element_blank(),
        legend.text = element_text(size = size.legend, margin = margin(r = 25, unit = 'pt')),
        legend.direction = 'horizontal',
        legend.spacing.x = unit(15, 'pt'),
        plot.margin = unit(size.margin, "pt"))

shps        <- c(15, 16)

p2 <- vdat2 %>% ggplot(aes(x = index, y = ftn(kldiv), col = approx, shape = approx)) +
  geom_point(size = size.point) + geom_line(size = size.line, alpha = alpha) +
  scale_x_continuous(breaks = c(1, 100, 200, 300, 400)) +
  # scale_y_continuous(trans='log10') +
  scale_color_manual(values = cols, labels = legs, guide = guide_legend(nrow = 1, byrow = TRUE)) +
  scale_shape_manual(values = shps, labels = legs, guide = guide_legend(nrow = 1, byrow = TRUE)) +
  xlab("index") + ylab("log( variance ratio )") + coord_cartesian(xlim = c(1, nrow(locs))) +
  theme(axis.title.x = element_text(size = size.lab),
        axis.text.x = element_text(size = size.text),
        axis.title.y = element_text(size = size.lab),
        axis.text.y = element_text(size = size.text),
        strip.text = element_text(size = size.lab),
        legend.title = element_blank(),
        legend.text = element_text(size = size.legend, margin = margin(r = 25, unit = 'pt')),
        legend.direction = 'horizontal',
        legend.spacing.x = unit(15, 'pt'),
        plot.margin = unit(size.margin, "pt"))

p3 <- vdat3 %>% ggplot(aes(x = index, y = kldiv, col = approx, shape = approx)) +
  geom_point(size = size.point) + geom_line(size = size.line, alpha = alpha) +
  scale_x_continuous(breaks = c(1, 100, 200, 300, 400)) +
  scale_color_manual(values = cols, labels = legs, guide = guide_legend(nrow = 1, byrow = TRUE)) +
  scale_shape_manual(values = shps, labels = legs, guide = guide_legend(nrow = 1, byrow = TRUE)) +
  xlab("index") + ylab("size of conditioning set") + coord_cartesian(xlim = c(1, nrow(locs))) +
  theme(axis.title.x = element_text(size = size.lab),
        axis.text.x = element_text(size = size.text),
        axis.title.y = element_text(size = size.lab),
        axis.text.y = element_text(size = size.text),
        strip.text = element_text(size = size.lab),
        legend.title = element_blank(),
        legend.text = element_text(size = size.legend, margin = margin(r = 25, unit = 'pt')),
        legend.direction = 'horizontal',
        legend.spacing.x = unit(15, 'pt'),
        plot.margin = unit(size.margin, "pt"))

###

tmp       <- ggplot_gtable(ggplot_build(p3))
leg       <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
mylegend  <- tmp$grobs[[leg]]

### Merge the two plots ###

result    <- grid.arrange(mylegend, arrangeGrob(p1 + theme(legend.position="none"),
                                                p3 + theme(legend.position="none"),
                                                p2 + theme(legend.position="none"),
                                                nrow=1), nrow=2, heights=c(1, 10))
###

ggplot2::ggsave("simulation_supp_m_and_rho_2.pdf", result, width = 15.2, height = 5.7)
