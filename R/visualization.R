####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes several R functions to visualize simulation results.
###
###   Contents:
###
####################################################################################

# library(dplyr) ; library(tidyr) ; library(ggplot2) ; library(gridExtra) ; library(RColorBrewer)
# 
# load("C:/Users/kmjst/Desktop/temp/simout_aniso_05242021.RData")
# 
# out     <- output.aniso
# vdat1   <- out$vars %>% left_join(out$kldiv, by = "index") %>% filter(a == 10) %>% select(-a)
# vdat2   <- out$vars %>% left_join(out$kldiv, by = "index") %>% filter(m == 20) %>% select(-m)
# vdat3   <- out$vars %>% left_join(out$kldiv, by = "index") %>% filter(m == 40) %>% select(-m)
# 
# vis     <- vis_arrange_ena(vdat1 = vdat1, ftn = log10, xlab = c("m"), ylab = "log10(KL)", xlim = NULL, ylim = NULL, legend = c("E-MM + E-NN", "E-MM + C-NN", "C-MM + E-NN", "C-MM + C-NN", "X-ord + E-NN", "Y-ord + E-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "#FF7F00", "#FFFF33"), shape = c(18, 15, 17, 16, 8, 13))
# 
# ggplot2::ggsave("example_ena.pdf", vis, width = 15.2, height = 5.7)
# 
# vis     <- vis_arrange_dio(vdat1 = vdat1, vdat2 = vdat2, ftn = log10, xlab = c("m", "a"), ylab = c("log10(KL)", "log10(KL)"), xlim = list(xl1 = NULL, xl2 = NULL), ylim = list(yl1 = NULL, yl2 = NULL), legend = c("E-MM + E-NN", "E-MM + C-NN", "C-MM + E-NN", "C-MM + C-NN", "X-ord + E-NN", "Y-ord + E-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "#FF7F00", "#FFFF33"), shape = c(18, 15, 17, 16, 8, 13))
# 
# ggplot2::ggsave("example_dio.pdf", vis, width = 15.2, height = 5.7)
# 
# vis     <- vis_arrange_tria(vdat1 = vdat1, vdat2 = vdat2, vdat3 = vdat3, ftn = log10, xlab = c("m", "a", "a"), ylab = c("log10(KL)", "log10(KL)", "log10(KL)"), xlim = list(xl1 = NULL, xl2 = NULL, xl3 = NULL), ylim = list(yl1 = NULL, yl2 = NULL, yl3 = NULL), legend = c("E-MM + E-NN", "E-MM + C-NN", "C-MM + E-NN", "C-MM + C-NN", "X-ord + E-NN", "Y-ord + E-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "#FF7F00", "#FFFF33"), shape = c(18, 15, 17, 16, 8, 13))
# 
# ggplot2::ggsave("example_tria.pdf", vis, width = 15.2, height = 5.7)

#' @title Visualization
#'
#' @param vdat1 dataset consisting of the columns: index, variable (for example, m), approx_1, approx_2, ...
#' @param legend legend
#' @param color color
#' @param shape shape 
#' @param ftn scale function
#' @param xlim xlim
#' @param ylim ylim
#' @param xlab xlab
#' @param ylab ylab
#' @param alpha alpha
#' @param size.point size.point
#' @param size.line size.line
#' @param size.legend size.legend
#' @param size.lab size.lab
#' @param size.text size.text
#' @param size.margin size.margin 
#'
#' @return output of grid.arrange(). Use plot(vis_arrange(...)) to plot it
#' 
#' @import ggplot2
#' @import dplyr
#' @import RColorBrewer
#' @importFrom gridExtra arrangeGrob
#' @importFrom gridExtra grid.arrange
#' 
#' @export
#'
#' @examples
#' 1 + 1
vis_arrange_ena <- function(vdat1, legend, color, shape, ftn = log10, xlim = NULL, ylim = NULL, xlab = "m", ylab = "log10(KL)", alpha = 0.7, size.point = 4, size.line = 1, size.legend = 14, size.lab = 14, size.text = 12, size.margin = c(5.5, 20, 5.5, 5.5))
{
  ### Visualize the data table ###
  
  # Gathering vdat1
  vars1     <- vdat1 %>% select(-index, -starts_with("approx")) %>% colnames()
  if(length(vars1) > 1) stop("Please check the argument vdat1.")
  
  xlabel1   <- vdat1 %>% pull(all_of(vars1)) %>% unique() %>% sort()
  vis1      <- vdat1 %>% tidyr::gather(key = "approx", value = "kldiv", -index, -all_of(vars1))
  
  # plot vis1
  plot1     <- vis1 %>% ggplot(aes(x = get(vars1), y = ftn(kldiv), col = approx, shape = approx)) + 
    geom_point(size = size.point) + 
    geom_line(size = size.line, alpha = alpha) + 
    scale_x_continuous(name = vars1, limits = range(xlabel1), breaks = xlabel1) +
    scale_color_manual(values = color, labels = legend, guide = guide_legend(nrow = 1, byrow = TRUE)) + 
    scale_shape_manual(values = shape, labels = legend, guide = guide_legend(nrow = 1, byrow = TRUE)) +
    xlab(xlab) + ylab(ylab) + coord_cartesian(xlim = xlim, ylim = ylim) +
    theme(axis.title.x = element_text(size = size.lab), 
          axis.text.x = element_text(size = size.text), 
          axis.title.y = element_text(size = size.lab), 
          axis.text.y = element_text(size = size.text), 
          legend.title = element_blank(), 
          legend.text = element_text(size = size.legend, margin = margin(r = 25, unit = 'pt')), 
          legend.direction = 'horizontal', 
          legend.spacing.x = unit(15, 'pt'),
          plot.margin = unit(size.margin, "pt")) # t, r, b, l
  
  ### Creating the legend ###
  
  tmp       <- ggplot_gtable(ggplot_build(plot1))
  leg       <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  mylegend  <- tmp$grobs[[leg]]
  
  ### return ###
  
  result    <- grid.arrange(mylegend, plot1 + theme(legend.position="none"), nrow=2, heights=c(1, 10))
  
  return( result )
}

#' @title Visualization
#'
#' @param vdat1 First dataset consisting of the columns: index, variable (for example, m), approx_1, approx_2, ...
#' @param vdat2 Second dataset consisting of the columns: index, variable (for example, a), approx_1, approx_2, ...
#' @param legend legend
#' @param color color
#' @param shape shape 
#' @param ftn scale function
#' @param xlim xlim
#' @param ylim ylim
#' @param xlab xlab
#' @param ylab ylab
#' @param alpha alpha
#' @param size.point size.point
#' @param size.line size.line
#' @param size.legend size.legend
#' @param size.lab size.lab
#' @param size.text size.text
#' @param size.margin size.margin 
#'
#' @return output of grid.arrange(). Use plot(vis_arrange(...)) to plot it
#' 
#' @import ggplot2
#' @import dplyr
#' @import RColorBrewer
#' @importFrom gridExtra arrangeGrob
#' @importFrom gridExtra grid.arrange
#' 
#' @export
#'
#' @examples
#' 1 + 1
vis_arrange_dio <- function(vdat1, vdat2, legend, color, shape, ftn = log10, xlim = c("m", "m"), ylim = c("log10(KL)", "log10(KL)"), xlab = list(xl1 = NULL, xl2 = NULL), ylab = list(yl1 = NULL, yl2 = NULL), alpha = 0.7, size.point = 2, size.line = 1, size.legend = 14, size.lab = 14, size.text = 12, size.margin = c(5.5, 20, 5.5, 5.5))
{
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
    scale_x_continuous(name = vars1, limits = range(xlabel1), breaks = xlabel1) +
    scale_color_manual(values = color, labels = legend, guide = guide_legend(nrow = 1, byrow = TRUE)) + 
    scale_shape_manual(values = shape, labels = legend, guide = guide_legend(nrow = 1, byrow = TRUE)) +
    xlab(xlab[1]) + ylab(ylab[1]) + coord_cartesian(xlim = xlim[[1]], ylim = ylim[[1]]) +
    theme(axis.title.x = element_text(size = size.lab), 
          axis.text.x = element_text(size = size.text), 
          axis.title.y = element_text(size = size.lab), 
          axis.text.y = element_text(size = size.text), 
          legend.title = element_blank(), 
          legend.text = element_text(size = size.legend, margin = margin(r = 25, unit = 'pt')), 
          legend.direction = 'horizontal', 
          legend.spacing.x = unit(15, 'pt'),
          plot.margin = unit(size.margin, "pt")) # t, r, b, l
  
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
    scale_x_continuous(name = vars2, limits = range(xlabel2), breaks = xlabel2) +
    scale_color_manual(values = color, labels = legend) + 
    scale_shape_manual(values = shape, labels = legend) +
    xlab(xlab[2]) + ylab(ylab[2]) + coord_cartesian(xlim = xlim[[2]], ylim = ylim[[2]]) +
    theme(axis.title.x = element_text(size = size.lab), 
          axis.text.x = element_text(size = size.text), 
          axis.title.y = element_text(size = size.lab), 
          axis.text.y = element_text(size = size.text), 
          legend.title = element_blank(), 
          legend.text = element_text(size = size.legend), 
          legend.direction = 'horizontal', 
          plot.margin = unit(size.margin, "pt")) # t, r, b, l
  
  ### Creating the legend ###
  
  tmp       <- ggplot_gtable(ggplot_build(plot1))
  leg       <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  mylegend  <- tmp$grobs[[leg]]
  
  ### Merge the two plots ###
  
  result    <- grid.arrange(mylegend, arrangeGrob(plot1 + theme(legend.position="none"), plot2 + theme(legend.position="none"), nrow=1), nrow=2, heights=c(1, 10))
  
  return( result )
}

#' @title Visualization
#'
#' @param vdat1 First dataset consisting of the columns: index, variable (for example, m), approx_1, approx_2, ...
#' @param vdat2 Second dataset consisting of the columns: index, variable (for example, m), approx_1, approx_2, ...
#' @param vdat3 Third dataset consisting of the columns: index, variable (for example, m), approx_1, approx_2, ...
#' @param legend legend
#' @param color color
#' @param shape shape 
#' @param ftn scale function
#' @param xlim xlim
#' @param ylim ylim
#' @param xlab xlab
#' @param ylab ylab
#' @param alpha alpha
#' @param size.point size.point
#' @param size.line size.line
#' @param size.legend size.legend
#' @param size.lab size.lab
#' @param size.text size.text
#' @param size.margin size.margin 
#'
#' @return output of grid.arrange(). Use plot(vis_arrange(...)) to plot it
#' 
#' @import ggplot2
#' @import dplyr
#' @import RColorBrewer
#' @importFrom gridExtra arrangeGrob
#' @importFrom gridExtra grid.arrange
#' 
#' @export
#'
#' @examples
#' 1 + 1
vis_arrange_tria <- function(vdat1, vdat2, vdat3, legend, color, shape, ftn = log10, xlim = list(xl1 = NULL, xl2 = NULL, xl3 = NULL), ylim = list(yl1 = NULL, yl2 = NULL, yl3 = NULL), xlab = c("xl1", "xl2", "xl3"), ylab = c("yl1", "yl2", "yl3"), alpha = 0.7, size.point = 4, size.line = 1, size.legend = 14, size.lab = 14, size.text = 12, size.margin = c(5.5, 20, 5.5, 5.5))
{
  size.margin.plot1 = size.margin + c(0, 100, 0, 100)
  
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
    scale_x_continuous(name = vars1, limits = range(xlabel1), breaks = xlabel1) +
    scale_color_manual(values = color, labels = legend, guide = guide_legend(nrow = 1, byrow = TRUE)) + 
    scale_shape_manual(values = shape, labels = legend, guide = guide_legend(nrow = 1, byrow = TRUE)) +
    xlab(xlab[1]) + ylab(ylab[1]) + coord_cartesian(xlim = xlim[[1]], ylim = ylim[[1]]) +
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
    scale_x_continuous(name = vars2, limits = range(xlabel2), breaks = xlabel2) +
    scale_color_manual(values = color, labels = legend) + 
    scale_shape_manual(values = shape, labels = legend) +
    xlab(xlab[2]) + ylab(ylab[2]) + coord_cartesian(xlim = xlim[[2]], ylim = ylim[[2]]) +
    theme(axis.title.x = element_text(size = size.lab), 
          axis.text.x = element_text(size = size.text), 
          axis.title.y = element_text(size = size.lab), 
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
    scale_x_continuous(name = vars3, limits = range(xlabel3), breaks = xlabel3) +
    scale_color_manual(values = color, labels = legend) + 
    scale_shape_manual(values = shape, labels = legend) +
    xlab(xlab[3]) + ylab(ylab[3]) + coord_cartesian(xlim = xlim[[3]], ylim = ylim[[3]]) +
    theme(axis.title.x = element_text(size = size.lab), 
          axis.text.x = element_text(size = size.text), 
          axis.title.y = element_text(size = size.lab), 
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
  
  result    <- grid.arrange(mylegend, arrangeGrob(plot1 + theme(legend.position="none"), plot2 + theme(legend.position="none"), plot3 + theme(legend.position="none"), layout_matrix = rbind(c(1, 1), c(2, 3))), nrow=2, heights=c(1, 10))
  
  return( result )
}

#' @title Visualization
#'
#' @param vdat1 First dataset consisting of the columns: index, variable (for example, m), approx_1, approx_2, ...
#' @param vdat2 Second dataset consisting of the columns: index, variable (for example, m), approx_1, approx_2, ...
#' @param vdat3 Third dataset consisting of the columns: index, variable (for example, m), approx_1, approx_2, ...
#' @param legend legend
#' @param color color
#' @param shape shape 
#' @param ftn scale function
#' @param xlim xlim
#' @param ylim ylim
#' @param xlab xlab
#' @param ylab ylab
#' @param alpha alpha
#' @param size.point size.point
#' @param size.line size.line
#' @param size.legend size.legend
#' @param size.lab size.lab
#' @param size.text size.text
#' @param size.margin size.margin 
#'
#' @return output of grid.arrange(). Use plot(vis_arrange(...)) to plot it
#' 
#' @import ggplot2
#' @import dplyr
#' @import RColorBrewer
#' @importFrom gridExtra arrangeGrob
#' @importFrom gridExtra grid.arrange
#'  
#' @export
#'
#' @examples
#' 1 + 1
vis_arrange_tria_simple <- function(vdat1, vdat2, vdat3, legend, color, shape, ftn = log10, xlim = list(xl1 = NULL, xl2 = NULL, xl3 = NULL), ylim = list(yl1 = NULL, yl2 = NULL, yl3 = NULL), xlab = c("xl1", "xl2", "xl3"), ylab = c("yl1", "yl2", "yl3"), alpha = 0.7, size.point = 4, size.line = 1, size.legend = 14, size.lab = 14, size.text = 12, size.margin = c(5.5, 20, 5.5, 5.5))
{
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
    scale_x_continuous(name = vars1, limits = range(xlabel1), breaks = xlabel1) +
    scale_color_manual(values = color, labels = legend, guide = guide_legend(nrow = 1, byrow = TRUE)) + 
    scale_shape_manual(values = shape, labels = legend, guide = guide_legend(nrow = 1, byrow = TRUE)) +
    xlab(xlab[1]) + ylab(ylab[1]) + coord_cartesian(xlim = xlim[[1]], ylim = ylim[[1]]) +
    theme(axis.title.x = element_text(size = size.lab), 
          axis.text.x = element_text(size = size.text), 
          axis.title.y = element_text(size = size.lab), 
          axis.text.y = element_text(size = size.text), 
          legend.title = element_blank(), 
          legend.text = element_text(size = size.legend, margin = margin(r = 25, unit = 'pt')), 
          legend.direction = 'horizontal', 
          legend.spacing.x = unit(15, 'pt'),
          plot.margin = unit(size.margin, "pt")) # t, r, b, l
  
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
    scale_x_continuous(name = vars2, limits = range(xlabel2), breaks = xlabel2) +
    scale_color_manual(values = color, labels = legend) + 
    scale_shape_manual(values = shape, labels = legend) +
    xlab(xlab[2]) + ylab(ylab[2]) + coord_cartesian(xlim = xlim[[2]], ylim = ylim[[2]]) +
    theme(axis.title.x = element_text(size = size.lab), 
          axis.text.x = element_text(size = size.text), 
          axis.title.y = element_text(size = size.lab), 
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
    scale_x_continuous(name = vars3, limits = range(xlabel3), breaks = xlabel3) +
    scale_color_manual(values = color, labels = legend) + 
    scale_shape_manual(values = shape, labels = legend) +
    xlab(xlab[3]) + ylab(ylab[3]) + coord_cartesian(xlim = xlim[[3]], ylim = ylim[[3]]) +
    theme(axis.title.x = element_text(size = size.lab), 
          axis.text.x = element_text(size = size.text), 
          axis.title.y = element_text(size = size.lab), 
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
  
  result    <- grid.arrange(mylegend, arrangeGrob(plot1 + theme(legend.position="none"), plot2 + theme(legend.position="none"), plot3 + theme(legend.position="none"), nrow = 1), nrow = 2, heights = c(1, 10))
  
  return( result )
}
