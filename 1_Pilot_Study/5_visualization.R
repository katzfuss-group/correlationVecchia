##########################################################################
#####
##### Author: Myeongjong Kang (kmj.stat@gmail.com)
#####
##### Description: 
#####
##########################################################################

# library(ggplot2) ; library(gridExtra) ; library(RColorBrewer)

vis_arrange <- function(vdat1, vdat2, combined.legend, color.pal = brewer.pal(4, "Set1"), alpha.value = 0.7, size.legend = 28, size.lab = 28, size.text = 18){
  
  xlabel1 <- sort(unique(vdat1$m))
  plot1   <- ggplot(vdat1, aes(x=m, y = log10(KL), col = method)) + 
             geom_point(aes(shape = method), size = 3) + 
             geom_line(size = 1, alpha = alpha.value) +
             ylab('log10(KL)') + 
             scale_x_discrete(name = 'm', limits=xlabel1, labels=as.character(xlabel1)) +
             scale_color_manual(values = color.pal, labels = combined.legend) +
             scale_shape_manual(values = c(15, 16, 17, 18), labels = combined.legend) +
             theme(axis.title.x = element_text(size = size.lab), 
                   axis.text.x = element_text(size = size.text),
                   axis.title.y = element_text(size = size.lab), 
                   axis.text.y = element_text(size = size.text),
                   legend.title = element_blank(),
                   legend.text=element_text(size = size.legend),
                   legend.direction = 'horizontal',
                   plot.margin = unit(c(5.5, 20, 5.5, 5.5), "pt")) # t, r, b, l
  
  xlabel2 <- sort(unique(vdat2$scale))
  plot2   <- ggplot(vdat2, aes(x=scale, y = log10(KL), col = method)) + 
             geom_point(aes(shape = method), size = 2) + 
             geom_line(size = 1, alpha = alpha.value) +
             ylab('log10(KL)') + 
             scale_x_discrete(name = 'a', limits=xlabel2, labels=as.character(xlabel2)) +
             scale_color_manual(values = color.pal, labels = combined.legend) +
             scale_shape_manual(values = c(15, 16, 17, 18), labels = combined.legend) +
             theme(axis.title.x = element_text(size = size.lab), 
                   axis.text.x = element_text(size = size.text),
                   axis.title.y = element_text(size = size.lab), 
                   axis.text.y = element_text(size = size.text),
                   plot.margin = unit(c(5.5, 5.5, 5.5, 20), "pt"))
  
  tmp <- ggplot_gtable(ggplot_build(plot1))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  mylegend <- tmp$grobs[[leg]]
  
  grid.arrange(mylegend, 
               arrangeGrob(plot1 + theme(legend.position="none"), 
                           plot2 + theme(legend.position="none"),
                           nrow=1), 
               nrow=2,heights=c(1, 10))
}


# size <- 28
# p1 <- ggplot(vis.dat1, aes(x=m, y = log10(KL), col = method)) + 
#   geom_point(aes(shape = method), size = 5) + geom_line(size = 2, alpha = 0.7) +
#   ylab('log10(KL)') + scale_x_discrete(name = 'm', limits=cand.m, labels=as.character(cand.m)) +
#   scale_color_manual(values = brewer.pal(4, "Set1"), labels = kls.legend) +
#   scale_shape_manual(values = c(15, 16, 17, 18), labels = kls.legend) +
#   theme(axis.title.x = element_text(size = size), 
#         axis.text.x = element_text(size = size-4),
#         axis.title.y = element_text(size = size), 
#         axis.text.y = element_text(size = size-4),
#         legend.position = 'top',
#         legend.title = element_blank(),
#         legend.direction = 'horizontal',
#         legend.text=element_text(size = size)) +
#   guides(color=guide_legend(nrow=2,byrow=FALSE))
# print(p)
# 
# size <- 28
# p2 <- ggplot(vis.dat2, aes(x=scale, y = log10(KL), col = method)) + 
#   geom_point(aes(shape = method), size = 5) + geom_line(size = 2, alpha = 0.7) +
#   ylab('log10(KL)') + scale_x_discrete(name = 'a', limits=cand.scale, labels=as.character(cand.scale)) +
#   scale_color_manual(values = brewer.pal(4, "Set1"), labels = kls.legend) +
#   scale_shape_manual(values = c(15, 16, 17, 18), labels = kls.legend) +
#   theme(axis.title.x = element_text(size = size), 
#         axis.text.x = element_text(size = size-4),
#         axis.title.y = element_text(size = size), 
#         axis.text.y = element_text(size = size-4),
#         legend.position = 'top',
#         legend.title = element_blank(),
#         legend.direction = 'horizontal',
#         legend.text=element_text(size = size)) +
#   guides(color=guide_legend(nrow=2,byrow=FALSE))
# print(p)