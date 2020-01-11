####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes several R functions to visualize simulation results.
###
####################################################################################



vis_arrange <- function(vdat1, vdat2 = NULL, legend = c("E-Maxmin + E-NN", "E-Maxmin + C-NN", "C-Maxmin + E-NN", "C-Maxmin + C-NN", "X-Coord + E-NN", "Y-Coord + E-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "#FF7F00", "#FFFF33"), shape = c(18, 15, 17, 16, 8, 13), alpha = 0.7, size.point = 4, size.line = 1, size.legend = 14, size.lab = 14, size.text = 12, size.margin = c(5.5, 20, 5.5, 5.5))
{
  n.space   <- 10
  legend    <- paste0( legend, paste0(rep(" ", n.space), collapse = "") )
  
  if(is.null(vdat2)) {
    
    stop("not yet!")
    
  } else {
    
    ### Visualize the first data table ###
    
    # Gathering vdat1
    vars1     <- vdat1 %>% select(-index, -starts_with("approx")) %>% colnames()
    if(length(vars1) > 1) stop("Please check the argument vdat1.")
    
    xlabel1   <- vdat1 %>% pull(vars1) %>% unique() %>% sort()
    vis1      <- vdat1 %>% tidyr::gather(key = "approx", value = "kldiv", -index, -vars1)
    
    # plot vis1
    plot1     <- vis1 %>% ggplot(aes(x = get(vars1), y = log10(kldiv), col = approx, shape = approx)) + 
      geom_point(size = size.point) + 
      geom_line(size = size.line, alpha = alpha) + 
      scale_x_continuous(name = vars1, limits = range(xlabel1), breaks = xlabel1) +
      scale_color_manual(values = color, labels = legend) + 
      scale_shape_manual(values = shape, labels = legend) +
      xlab(vars1) + ylab('log10(KL)') + 
      theme(axis.title.x = element_text(size = size.lab), 
            axis.text.x = element_text(size = size.text), 
            axis.title.y = element_text(size = size.lab), 
            axis.text.y = element_text(size = size.text), 
            legend.title = element_blank(), 
            legend.text = element_text(size = size.legend), 
            legend.direction = 'horizontal', 
            plot.margin = unit(size.margin, "pt")) # t, r, b, l
    
    ### Visualize the second data table ###
    
    # Gathering vdat2
    vars2     <- vdat2 %>% select(-index, -starts_with("approx")) %>% colnames()
    if(length(vars2) > 1) stop("Please check the argument vdat2.")
    
    xlabel2   <- vdat2 %>% pull(vars2) %>% unique() %>% sort()
    vis2      <- vdat2 %>% tidyr::gather(key = "approx", value = "kldiv", -index, -vars2)
    
    # plot vis2
    plot2     <- vis2 %>% ggplot(aes(x = get(vars2), y = log10(kldiv), col = approx, shape = approx)) + 
      geom_point(size = size.point) + 
      geom_line(size = size.line, alpha = alpha) + 
      scale_x_continuous(name = vars2, limits = range(xlabel2), breaks = xlabel2) +
      scale_color_manual(values = color, labels = legend) + 
      scale_shape_manual(values = shape, labels = legend) +
      xlab(vars2) + ylab('log10(KL)') + 
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
    
  }
  
  return(result)
}