####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes several R functions to visualize simulation results.
###
###   Contents:
###     vis_arrange_ena / vis_arrange_dio / vis_arrange_tria
###     vis_arrange_tria_simple / vis_arrange_posterior
###     vis_arrange_tetra / vis_arrange_posterior_ic0
###
###   Note: For generating figures in the manuscript, we used modified versions of these functions. Please take a look into the output folder.
###
####################################################################################

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

#' @title Visualization
#'
#' @param output output from posterior_... functions
#' @param legend legend
#' @param color color
#' @param xlim xlim
#' @param ylim ylim
#' @param xlab xlab
#' @param ylab ylab
#' @param alpha alpha
#' @param size.line size.line
#' @param size.legend size.legend
#' @param size.lab size.lab
#' @param size.text size.text
#' @param size.margin size.margin
#'
#' @return plot
#'
#' @export
vis_arrange_posterior <- function(output, legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "True"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "black"), xlim = NULL, ylim = NULL, xlab = "x", ylab = "y", alpha = 0.7, size.line = 1, size.legend = 14, size.lab = 14, size.text = 12, size.margin = c(5.5, 20, 5.5, 5.5))
{
  candids   <- output$setting$candid
  ms        <- unique(candids$m)
  approxs   <- unique(candids$approx)

  idx.list <- list()
  for(i in 1:length(ms)) idx.list[[i]] <- which(candids$m == ms[i])

  visdat.ms <- list()
  for(i in 1:length(ms)) {

    idx             <- idx.list[[i]]
    visdat.ms[[i]]  <- data.frame(gp = "tr", x = output$simout[[ idx[1] ]]$alpha, y = output$simout[[ idx[1] ]]$post.norm[, 1])

    for(j in 1:length(approxs)) {

      tempdat         <- data.frame(gp = approxs[j], x = output$simout[[ idx[j] ]]$alpha, y = output$simout[[ idx[j] ]]$post.norm[, 2])
      visdat.ms[[i]]  <- rbind(visdat.ms[[i]], tempdat)
    }

    visdat.ms[[i]]  <- visdat.ms[[i]] %>% filter(gp != "ex")
  }

  plots <- list()
  for(i in 1:length(visdat.ms)) {

    plots[[i]] <- visdat.ms[[i]] %>%
      ggplot(aes(x = x, y = y, col = gp)) +
      geom_line(size = size.line, alpha = alpha) +
      scale_color_manual(values = color, labels = legend) +
      xlab(xlab) + ylab(ylab) +
      coord_cartesian(xlim = xlim, ylim = ylim) +
      theme(axis.title.x = element_text(size = size.lab),
            axis.text.x = element_text(size = size.text),
            axis.title.y = element_text(size = size.lab),
            axis.text.y = element_text(size = size.text),
            legend.title = element_blank(),
            legend.text = element_text(size = size.legend, margin = margin(r = 25, unit = 'pt')),
            legend.direction = 'horizontal',
            legend.spacing.x = unit(15, 'pt'),
            plot.margin = unit(size.margin, "pt"))
  }

  tmp       <- ggplot_gtable(ggplot_build(plots[[1]]))
  leg       <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  mylegend  <- tmp$grobs[[leg]]

  ### CAUTION: hard-coded
  result    <- grid.arrange(mylegend, arrangeGrob(plots[[1]] + theme(legend.position="none"),
                                                  plots[[2]] + theme(legend.position="none"),
                                                  plots[[3]] + theme(legend.position="none"),
                                                  plots[[4]] + theme(legend.position="none"),
                                                  plots[[5]] + theme(legend.position="none"), nrow = 1),
                            nrow = 2, heights = c(1, 10))

  return( result )
}

#' @title Visualization
#'
#' @param vdat1 vdat1
#' @param vdat2 vdat2
#' @param vdat3 vdat3
#' @param vdat4 vdat4
#' @param legend legend
#' @param color color
#' @param shape shape
#' @param ftn ftn
#' @param xlabs xlabs
#' @param ylabs ylabs
#' @param xlims xlims
#' @param ylims ylims
#' @param alpha alpha
#' @param size.point size.point
#' @param size.line size.line
#' @param size.legend size.legend
#' @param size.lab size.lab
#' @param size.text size.text
#' @param size.margin size.margin
#'
#' @return plot
#'
#' @export
vis_arrange_tetra <- function(vdat1, vdat2, vdat3, vdat4, legend, color, shape, ftn = log10, xlabs = c("m", "m", "m", "m"), ylabs = c("log10(KL)", "log10(KL)", "log10(KL)", "log10(KL)"), xlims = list(xl1 = NULL, xl2 = NULL, xl3 = NULL, xl4 = NULL), ylims = list(yl1 = NULL, yl2 = NULL, yl3 = NULL, yl4 = NULL), alpha = 0.7, size.point = 2, size.line = 1, size.legend = 14, size.lab = 14, size.text = 12, size.margin = c(5.5, 20, 5.5, 5.5))
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
    # scale_x_continuous(name = vars1, limits = range(xlabel1), breaks = xlabel1) +
    scale_color_manual(values = color, labels = legend, guide = guide_legend(nrow = 1, byrow = TRUE)) +
    scale_shape_manual(values = shape, labels = legend, guide = guide_legend(nrow = 1, byrow = TRUE)) +
    xlab(xlabs[1]) + ylab(ylabs[1]) + coord_cartesian(xlim = xlims[[1]], ylim = ylims[[1]]) +
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
    # scale_x_continuous(name = vars2, limits = range(xlabel2), breaks = xlabel2) +
    scale_color_manual(values = color, labels = legend) +
    scale_shape_manual(values = shape, labels = legend) +
    xlab(xlabs[2]) + ylab(ylabs[2]) + coord_cartesian(xlim = xlims[[2]], ylim = ylims[[2]]) +
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
    # scale_x_continuous(name = vars1, limits = range(xlabel3), breaks = xlabel3) +
    scale_color_manual(values = color, labels = legend, guide = guide_legend(nrow = 1, byrow = TRUE)) +
    scale_shape_manual(values = shape, labels = legend, guide = guide_legend(nrow = 1, byrow = TRUE)) +
    xlab(xlabs[3]) + ylab(ylabs[3]) + coord_cartesian(xlim = xlims[[3]], ylim = ylims[[3]]) +
    theme(axis.title.x = element_text(size = size.lab),
          axis.text.x = element_text(size = size.text),
          axis.title.y = element_text(size = size.lab),
          axis.text.y = element_text(size = size.text),
          legend.title = element_blank(),
          legend.text = element_text(size = size.legend, margin = margin(r = 25, unit = 'pt')),
          legend.direction = 'horizontal',
          legend.spacing.x = unit(15, 'pt'),
          plot.margin = unit(size.margin, "pt")) # t, r, b, l

  ### Visualize the fourth data table ###

  # Gathering vdat4
  vars4     <- vdat4 %>% select(-index, -starts_with("approx")) %>% colnames()
  if(length(vars4) > 1) stop("Please check the argument vdat4.")

  xlabel4   <- vdat4 %>% pull(all_of(vars4)) %>% unique() %>% sort()
  vis4      <- vdat4 %>% tidyr::gather(key = "approx", value = "kldiv", -index, -all_of(vars4))

  # plot vis4
  plot4     <- vis4 %>% ggplot(aes(x = get(vars4), y = ftn(kldiv), col = approx, shape = approx)) +
    geom_point(size = size.point) +
    geom_line(size = size.line, alpha = alpha) +
    # scale_x_continuous(name = vars1, limits = range(xlabel4), breaks = xlabel4) +
    scale_color_manual(values = color, labels = legend, guide = guide_legend(nrow = 1, byrow = TRUE)) +
    scale_shape_manual(values = shape, labels = legend, guide = guide_legend(nrow = 1, byrow = TRUE)) +
    xlab(xlabs[4]) + ylab(ylabs[4]) + coord_cartesian(xlim = xlims[[4]], ylim = ylims[[4]]) +
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

  ### Merge the two plots ###

  result    <- grid.arrange(mylegend, arrangeGrob(plot1 + theme(legend.position="none"),
                                                  plot2 + theme(legend.position="none"),
                                                  plot3 + theme(legend.position="none"),
                                                  plot4 + theme(legend.position="none"),
                                                  nrow=2), nrow=2, heights=c(1, 10))

  return( result )
}

#' @title Visualization
#'
#' @param output output from posterior_... functions
#' @param legend legend
#' @param color color
#' @param xlim xlim
#' @param ylim ylim
#' @param xlab xlab
#' @param ylab ylab
#' @param alpha alpha
#' @param size.line size.line
#' @param size.legend size.legend
#' @param size.lab size.lab
#' @param size.text size.text
#' @param size.margin size.margin
#'
#' @return plot
#'
#' @export
vis_arrange_posterior_ic0 <- function(output, legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "True"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "black"), xlim = NULL, ylim = NULL, xlab = "x", ylab = "y", alpha = 0.7, size.line = 1, size.legend = 14, size.lab = 14, size.text = 12, size.margin = c(5.5, 20, 5.5, 5.5))
{
  candids   <- output$setting$candid
  ms        <- unique(candids$m)
  approxs   <- unique(candids$approx)

  idx.list <- list()
  for(i in 1:length(ms)) idx.list[[i]] <- which(candids$m == ms[i])

  visdat.ms <- list()
  for(i in 1:length(ms)) {

    idx             <- idx.list[[i]]
    visdat.ms[[i]]  <- data.frame(gp = "tr", x = output$simout.ic0[[ idx[1] ]]$alpha, y = output$simout.ic0[[ idx[1] ]]$post.norm[, 1])

    for(j in 1:length(approxs)) {

      tempdat         <- data.frame(gp = approxs[j], x = output$simout.ic0[[ idx[j] ]]$alpha, y = output$simout.ic0[[ idx[j] ]]$post.norm[, 2])
      visdat.ms[[i]]  <- rbind(visdat.ms[[i]], tempdat)
    }

    visdat.ms[[i]]  <- visdat.ms[[i]] %>% filter(gp != "ex")
  }

  plots <- list()
  for(i in 1:length(visdat.ms)) {

    plots[[i]] <- visdat.ms[[i]] %>%
      ggplot(aes(x = x, y = y, col = gp)) +
      geom_line(size = size.line, alpha = alpha) +
      scale_color_manual(values = color, labels = legend) +
      xlab(xlab) + ylab(ylab) +
      coord_cartesian(xlim = xlim, ylim = ylim) +
      theme(axis.title.x = element_text(size = size.lab),
            axis.text.x = element_text(size = size.text),
            axis.title.y = element_text(size = size.lab),
            axis.text.y = element_text(size = size.text),
            legend.title = element_blank(),
            legend.text = element_text(size = size.legend, margin = margin(r = 25, unit = 'pt')),
            legend.direction = 'horizontal',
            legend.spacing.x = unit(15, 'pt'),
            plot.margin = unit(size.margin, "pt"))
  }

  tmp       <- ggplot_gtable(ggplot_build(plots[[1]]))
  leg       <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  mylegend  <- tmp$grobs[[leg]]

  ### CAUTION: hard-coded
  result    <- grid.arrange(mylegend, arrangeGrob(plots[[1]] + theme(legend.position="none"),
                                                  plots[[2]] + theme(legend.position="none"),
                                                  plots[[3]] + theme(legend.position="none"),
                                                  plots[[4]] + theme(legend.position="none"),
                                                  plots[[5]] + theme(legend.position="none"), nrow = 1),
                            nrow = 2, heights = c(1, 10))

  return( result )
}
