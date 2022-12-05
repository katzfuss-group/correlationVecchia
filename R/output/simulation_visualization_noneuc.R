####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script generates Figure 6 (simulation_noneuc.pdf).
###
####################################################################################

rm(list = ls())

library(dplyr) ; library(ggplot2) ; library(gridExtra) ; library(grid)

####################################################################################

load("simout_noneuc_10012021.RData")

xticks      <- ms

vdat1       <- matrix(NA, nrow = length(ms), ncol = 4 + length(kls3))
vdat1[, 1]  <- seq(length(ms))
vdat1[, 2]  <- ms
vdat1[, 3]  <- kls1
vdat1[, 4]  <- kls2

for(i in 1:length(kls3)) vdat1[, i + 4] <- kls3[[i]]

vdat1       <- as.data.frame(vdat1)
colnames(vdat1) <- c("index", "m", paste0("approx_", seq(length(kls3)+2)))

vdat1[1:5, 1:10]

vdat2       <- vdat1[, c("index", "m", "approx_3", "approx_2", "approx_1")]

vdat2

vdat3       <- vdat1 %>% select(index, m, approx_1)
vdat4       <- vdat1 %>% select(index, m, approx_2)
vdat5       <- vdat1 %>% select(-approx_1, -approx_2)

vdat1      <- vdat1 %>% select(-index) %>% tidyr::gather(key = "approx", value = "kldiv", -m)
vdat2      <- vdat2 %>% select(-index) %>% tidyr::gather(key = "approx", value = "kldiv", -m)
vdat3      <- vdat3 %>% select(-index) %>% tidyr::gather(key = "approx", value = "kldiv", -m)
vdat4      <- vdat4 %>% select(-index) %>% tidyr::gather(key = "approx", value = "kldiv", -m)
vdat5      <- vdat5 %>% select(-index) %>% tidyr::gather(key = "approx", value = "kldiv", -m)

####################################################################################

legs        <- c("R-ord + R-N", "L-ord + C-NN", "C-MM + C-NN")
cols        <- c("#377EB8", "#4DAF4A", "#E41A1C")
shps        <- c(15, 17, 16)
ftn         <- log10
alpha       <- 0.7
size.point  <- 5
size.tick   <- 1
size.line   <- 1
size.legend <- 18
size.lab    <- 18
size.text   <- 16
size.margin <- c(5.5, 20, 5.5, 5.5)

####################################################################################

p2 <- vdat2 %>% ggplot(aes(x = m, y = ftn(kldiv), col = approx, shape = approx)) +
  geom_point(size = size.point) + geom_line(size = size.line, alpha = alpha) +
  scale_color_manual(values = cols, labels = legs, guide = guide_legend(nrow = 1, byrow = TRUE)) +
  scale_shape_manual(values = shps, labels = legs, guide = guide_legend(nrow = 1, byrow = TRUE)) +
  xlab("increasing m") + ylab("log10(KL)") + # coord_cartesian(xlim = xlims[[1]], ylim = ylims[[1]]) +
  scale_x_continuous(breaks = xticks) +
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

tmp       <- ggplot_gtable(ggplot_build(p2))
leg       <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
mylegend  <- tmp$grobs[[leg]]

####################################################################################

p1 <- ggplot() +
  geom_point(data = vdat3, aes(x = m, y = ftn(kldiv)), colour = cols[3], shape = shps[3], size = size.point) +
  geom_line(data = vdat3, aes(x = m, y = ftn(kldiv)), colour = cols[3], size = size.line, alpha = alpha) +
  geom_point(data = vdat4, aes(x = m, y = ftn(kldiv)), colour = cols[2], shape = shps[2], size = size.point) +
  geom_line(data = vdat4, aes(x = m, y = ftn(kldiv)), colour = cols[2], size = size.line, alpha = alpha) +
  geom_point(data = vdat5, aes(x = m, y = ftn(kldiv), group = approx), colour = cols[1], shape = shps[1], size = size.point) +
  geom_line(data = vdat5, aes(x = m, y = ftn(kldiv), group = approx), colour = cols[1], size = size.line, alpha = alpha) +
  xlab("increasing m") + ylab("log10(KL)") + # coord_cartesian(xlim = xlims[[1]], ylim = ylims[[1]]) +
  scale_x_continuous(breaks = xticks) +
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

####################################################################################

result    <- grid.arrange(mylegend, p1 + theme(legend.position="none"), nrow=2, heights=c(1, 10))

ggplot2::ggsave("simulation_noneuc.pdf", result, width = 7.6, height = 5)

####################################################################################
