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

load("C:/Users/kmjst/Dropbox/4_PhD_Personal/12_2021Su/1_Research/06292021/output/fisher/simout_fisher_random_05242021.RData")

xticks1 <- cand.m

rm(cand.m, covparms, nsim)

load("C:/Users/kmjst/Dropbox/4_PhD_Personal/12_2021Su/1_Research/06292021/output/fisher/simout_fisher_monitoring_05242021.RData")

xticks2 <- cand.m

rm(cand.m, covparms, nsim)

load("C:/Users/kmjst/Dropbox/4_PhD_Personal/12_2021Su/1_Research/06292021/output/fisher/simout_fisher_satellite_05242021.RData")

xticks3 <- cand.m

rm(cand.m, covparms, nsim)

####################################################################################

vis_arrange_fisher <- function(vdat1, vdat2, vdat3, legend, color, shape, ftn = log10, xticks, xlim = list(xl1 = NULL, xl2 = NULL, xl3 = NULL), ylim = list(yl1 = NULL, yl2 = NULL, yl3 = NULL), xlab = c("xl1", "xl2", "xl3"), ylab = c("yl1", "yl2", "yl3"), alpha = 0.7, size.point = 4, size.line = 1, size.legend = 18, size.lab = 18, size.text = 16, size.margin = c(5.5, 20, 5.5, 5.5))
{
  size.margin.plot1 = size.margin # + c(0, 100, 0, 100)

  ### Visualize the first data table ###

  # Gathering vdat1
  vars1     <- vdat1 %>% select(-index, -starts_with("approx")) %>% colnames()
  if(length(vars1) > 1) stop("Please check the argument vdat1.")

  xlabel1   <- vdat1 %>% pull(all_of(vars1)) %>% unique() %>% sort()
  vis1      <- vdat1 %>% tidyr::gather(key = "approx", value = "kldiv", -index, -all_of(vars1))

  # plot vis1
  plot1     <- vis1 %>% ggplot(aes(x = get(vars1), y = ftn(kldiv), col = approx, shape = approx)) +
    geom_point(size = size.point) +
    geom_line(size = size.line, alpha = alpha) +
    # scale_x_continuous(name = vars1, limits = range(xlabel1), breaks = xlabel1) +
    scale_color_manual(values = color, labels = legend, guide = guide_legend(nrow = 1, byrow = TRUE)) +
    scale_shape_manual(values = shape, labels = legend, guide = guide_legend(nrow = 1, byrow = TRUE)) +
    xlab(NULL) + ylab(ylab[1]) + coord_cartesian(xlim = xlim[[1]], ylim = ylim[[1]]) +
    scale_x_continuous(breaks = xticks) +
    theme(axis.title.x = element_text(size = size.lab),
          axis.text.x = element_text(size = size.text),
          axis.title.y = element_text(size = size.lab),
          axis.text.y = element_text(size = size.text),
          legend.title = element_blank(),
          legend.text = element_text(size = size.legend, margin = margin(r = 25, unit = 'pt')),
          legend.direction = 'horizontal',
          legend.spacing.x = unit(15, 'pt'),
          plot.margin = unit(size.margin.plot1, "pt")) # t, r, b, l

  ### Visualize the second data table ###

  # Gathering vdat2
  vars2     <- vdat2 %>% select(-index, -starts_with("approx")) %>% colnames()
  if(length(vars2) > 1) stop("Please check the argument vdat2.")

  xlabel2   <- vdat2 %>% pull(all_of(vars2)) %>% unique() %>% sort()
  vis2      <- vdat2 %>% tidyr::gather(key = "approx", value = "kldiv", -index, -all_of(vars2))

  # plot vis2
  plot2     <- vis2 %>% ggplot(aes(x = get(vars2), y = ftn(kldiv), col = approx, shape = approx)) +
    geom_point(size = size.point) +
    geom_line(size = size.line, alpha = alpha) +
    # scale_x_continuous(name = vars2, limits = range(xlabel2), breaks = xlabel2) +
    scale_color_manual(values = color, labels = legend) +
    scale_shape_manual(values = shape, labels = legend) +
    xlab(NULL) + ylab(ylab[2]) + coord_cartesian(xlim = xlim[[2]], ylim = ylim[[2]]) +
    scale_x_continuous(breaks = xticks) +
    theme(axis.title.x = element_text(size = size.lab),
          axis.text.x = element_text(size = size.text),
          axis.title.y = element_text(size = size.lab-4),
          axis.text.y = element_text(size = size.text),
          legend.title = element_blank(),
          legend.text = element_text(size = size.legend),
          legend.direction = 'horizontal',
          plot.margin = unit(size.margin, "pt")) # t, r, b, l

  ### Visualize the third data table ###

  # Gathering vdat3
  vars3     <- vdat3 %>% select(-index, -starts_with("approx")) %>% colnames()
  if(length(vars3) > 1) stop("Please check the argument vdat3.")

  xlabel3   <- vdat3 %>% pull(all_of(vars3)) %>% unique() %>% sort()
  vis3      <- vdat3 %>% tidyr::gather(key = "approx", value = "kldiv", -index, -all_of(vars3))

  # plot vis3
  plot3     <- vis3 %>% ggplot(aes(x = get(vars3), y = ftn(kldiv), col = approx, shape = approx)) +
    geom_point(size = size.point) +
    geom_line(size = size.line, alpha = alpha) +
    # scale_x_continuous(name = vars3, limits = range(xlabel3), breaks = xlabel3) +
    scale_color_manual(values = color, labels = legend) +
    scale_shape_manual(values = shape, labels = legend) +
    xlab(NULL) + ylab(ylab[3]) + coord_cartesian(xlim = xlim[[3]], ylim = ylim[[3]]) +
    scale_x_continuous(breaks = xticks) +
    theme(axis.title.x = element_text(size = size.lab),
          axis.text.x = element_text(size = size.text),
          axis.title.y = element_text(size = size.lab-4),
          axis.text.y = element_text(size = size.text),
          legend.title = element_blank(),
          legend.text = element_text(size = size.legend),
          legend.direction = 'horizontal',
          plot.margin = unit(size.margin, "pt")) # t, r, b, l

  ### Creating the legend ###

  tmp       <- ggplot_gtable(ggplot_build(plot1))
  leg       <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  mylegend  <- tmp$grobs[[leg]]

  ### Merge the three plots ###

  result    <- grid.arrange(mylegend, arrangeGrob(plot1 + theme(legend.position="none"),
                                                  plot2 + theme(legend.position="none"),
                                                  plot3 + theme(legend.position="none"),
                                                  layout_matrix = rbind(c(1, 2), c(1, 3))),
                            nrow=2, heights=c(1, 10),
                            bottom = textGrob(xlab[1], gp = gpar(fontsize = size.lab)))

  return( result )
}
####################################################################################

vdat1   <- output.monitoring$kls.average ; vdat1[, !(names(vdat1) %in% c("index", "m"))] <- log10(vdat1[, !(names(vdat1) %in% c("index", "m"))])
vdat2   <- data.frame(index = seq(length(output.monitoring$kls.average$m)), m = output.monitoring$kls.average$m, approx_1 = sqrt(output.monitoring$msd.srange[[1]]), approx_2 = sqrt(output.monitoring$msd.srange[[2]]), approx_3 = sqrt(output.monitoring$msd.srange[[3]]), approx_4 = sqrt(output.monitoring$msd.srange[[4]]))
vdat3   <- data.frame(index = seq(length(output.monitoring$kls.average$m)), m = output.monitoring$kls.average$m, approx_1 = sqrt(output.monitoring$msd.trange[[1]]), approx_2 = sqrt(output.monitoring$msd.trange[[2]]), approx_3 = sqrt(output.monitoring$msd.trange[[3]]), approx_4 = sqrt(output.monitoring$msd.trange[[4]]))

vdat1   <- vdat1 %>% filter(m != 3)
vdat2   <- vdat2 %>% filter(m != 3)
vdat3   <- vdat3 %>% filter(m != 3)

ylim1   <- vdat1 %>% select(approx_1, approx_2, approx_3, approx_4, approx_5) %>% range()
ylim2   <- vdat2 %>% select(approx_1, approx_2, approx_3, approx_4) %>% range()
ylim3   <- vdat3 %>% select(approx_2, approx_3, approx_4) %>% range()

vis.mon <- vis_arrange_fisher(vdat1 = vdat1, vdat2 = vdat2, vdat3 = vdat3, ftn = function(x) x, xticks = xticks1[-1], xlim = list(xl1 = NULL, xl2 = NULL, xl3 = NULL), ylim = list(yl1 = ylim1, yl2 = ylim2, yl3 = ylim3), xlab = c("increasing m", "increasing m", "increasing m"), ylab = c("log10(KL)", "RMSD for spatial range", "RMSD for temporal range"), legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "Exact GP"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "gray30"), shape = c(18, 15, 17, 16, NA))

ggplot2::ggsave("simulation_spacetime_fisher_mon.pdf", vis.mon, width = 12.0, height = 5.7)

vdat1   <- output.satellite$kls.average ; vdat1[, !(names(vdat1) %in% c("index", "m"))] <- log10(vdat1[, !(names(vdat1) %in% c("index", "m"))])
vdat2   <- data.frame(index = seq(length(output.satellite$kls.average$m)), m = output.satellite$kls.average$m, approx_1 = sqrt(output.satellite$msd.srange[[1]]), approx_2 = sqrt(output.satellite$msd.srange[[2]]), approx_3 = sqrt(output.satellite$msd.srange[[3]]), approx_4 = sqrt(output.satellite$msd.srange[[4]]))
vdat3   <- data.frame(index = seq(length(output.satellite$kls.average$m)), m = output.satellite$kls.average$m, approx_1 = sqrt(output.satellite$msd.trange[[1]]), approx_2 = sqrt(output.satellite$msd.trange[[2]]), approx_3 = sqrt(output.satellite$msd.trange[[3]]), approx_4 = sqrt(output.satellite$msd.trange[[4]]))

# vdat2   <- data.frame(index = seq(length(output.satellite$kls.average$m)), m = output.satellite$kls.average$m, approx_1 = sqrt(output.satellite$mse.srange[[1]]), approx_2 = sqrt(output.satellite$mse.srange[[2]]), approx_3 = sqrt(output.satellite$mse.srange[[3]]), approx_4 = sqrt(output.satellite$mse.srange[[4]]))
# vdat3   <- data.frame(index = seq(length(output.satellite$kls.average$m)), m = output.satellite$kls.average$m, approx_1 = sqrt(output.satellite$mse.trange[[1]]), approx_2 = sqrt(output.satellite$mse.trange[[2]]), approx_3 = sqrt(output.satellite$mse.trange[[3]]), approx_4 = sqrt(output.satellite$mse.trange[[4]]))

vdat1   <- vdat1 %>% filter(m != 3)
vdat2   <- vdat2 %>% filter(m != 3)
vdat3   <- vdat3 %>% filter(m != 3)

ylim1   <- vdat1 %>% select(approx_1, approx_2, approx_3, approx_4, approx_5) %>% range()
ylim2   <- vdat2 %>% select(approx_1, approx_2, approx_3, approx_4) %>% range()
ylim3   <- vdat3 %>% select(approx_2, approx_3, approx_4) %>% range()

vis.sat <- vis_arrange_fisher(vdat1 = vdat1, vdat2 = vdat2, vdat3 = vdat3, ftn = function(x) x, xticks = xticks2[-1], xlim = list(xl1 = NULL, xl2 = NULL, xl3 = NULL), ylim = list(yl1 = ylim1, yl2 = ylim2, yl3 = ylim3), xlab = c("increasing m", "increasing m", "increasing m"), ylab = c("log10(KL)", "RMSD for spatial range", "RMSD for temporal range"), legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "Exact GP"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "gray30"), shape = c(18, 15, 17, 16, NA))

ggplot2::ggsave("simulation_spacetime_fisher_sat.pdf", vis.sat, width = 12.0, height = 5.7)

vdat1   <- output.random$kls.average ; vdat1[, !(names(vdat1) %in% c("index", "m"))] <- log10(vdat1[, !(names(vdat1) %in% c("index", "m"))])
vdat2   <- data.frame(index = seq(length(output.random$kls.average$m)), m = output.random$kls.average$m, approx_1 = sqrt(output.random$msd.srange[[1]]), approx_2 = sqrt(output.random$msd.srange[[2]]), approx_3 = sqrt(output.random$msd.srange[[3]]), approx_4 = sqrt(output.random$msd.srange[[4]]))
vdat3   <- data.frame(index = seq(length(output.random$kls.average$m)), m = output.random$kls.average$m, approx_1 = sqrt(output.random$msd.trange[[1]]), approx_2 = sqrt(output.random$msd.trange[[2]]), approx_3 = sqrt(output.random$msd.trange[[3]]), approx_4 = sqrt(output.random$msd.trange[[4]]))

# vdat2   <- data.frame(index = seq(length(output.random$kls.average$m)), m = output.random$kls.average$m, approx_1 = sqrt(output.random$mse.srange[[1]]), approx_2 = sqrt(output.random$mse.srange[[2]]), approx_3 = sqrt(output.random$mse.srange[[3]]), approx_4 = sqrt(output.random$mse.srange[[4]]))
# vdat3   <- data.frame(index = seq(length(output.random$kls.average$m)), m = output.random$kls.average$m, approx_1 = sqrt(output.random$mse.trange[[1]]), approx_2 = sqrt(output.random$mse.trange[[2]]), approx_3 = sqrt(output.random$mse.trange[[3]]), approx_4 = sqrt(output.random$mse.trange[[4]]))

vdat1   <- vdat1 %>% filter(m != 3)
vdat2   <- vdat2 %>% filter(m != 3)
vdat3   <- vdat3 %>% filter(m != 3)

ylim1   <- vdat1 %>% select(approx_1, approx_2, approx_3, approx_4, approx_5) %>% range()
ylim2   <- vdat2 %>% select(approx_1, approx_2, approx_3, approx_4) %>% range()
ylim3   <- vdat3 %>% select(approx_2, approx_3, approx_4) %>% range()

ylim3[2]<- 0.2

vis.ran <- vis_arrange_fisher(vdat1 = vdat1, vdat2 = vdat2, vdat3 = vdat3, ftn = function(x) x, xticks = xticks3[-1], xlim = list(xl1 = NULL, xl2 = NULL, xl3 = NULL), ylim = list(yl1 = ylim1, yl2 = ylim2, yl3 = ylim3), xlab = c("increasing m", "increasing m", "increasing m"), ylab = c("log10(KL)", "RMSD for spatial range", "RMSD for temporal range"), legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "Exact GP"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "gray30"), shape = c(18, 15, 17, 16, NA))

ggplot2::ggsave("simulation_spacetime_fisher_ran.pdf", vis.ran, width = 12.0, height = 5.7)

