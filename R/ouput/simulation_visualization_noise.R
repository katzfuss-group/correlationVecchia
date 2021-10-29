####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview:
###
####################################################################################

rm(list = ls())

library(dplyr) ; library(ggplot2) ; library(gridExtra) ; library(grid) ; library(ggh4x)

####################################################################################

load("C:/Users/kmjst/Dropbox/4_PhD_Personal/12_2021Su/1_Research/06292021/output/noise/simout_posterior_allrandom_05242021.RData")

rm(covparms, d, ms, n, nugget, t)

load("C:/Users/kmjst/Dropbox/4_PhD_Personal/12_2021Su/1_Research/06292021/output/noise/simout_posterior_monitoring_05242021.RData")

rm(covparms, d, ms, n, nugget, t)

load("C:/Users/kmjst/Dropbox/4_PhD_Personal/12_2021Su/1_Research/06292021/output/noise/simout_posterior_satellite_05242021.RData")

rm(covparms, d, ms, n, nugget, t)

####################################################################################

plots_post <- function(output.s, output.t, label, legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "True"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "black"), alpha = 0.7, size.line = 1, size.legend = 18, size.lab = 18, size.text = 16, size.margin = c(5.5, 20, 5.5, 5.5))
{
  candids   <- output.s$setting$candid
  ms        <- unique(candids$m)
  approxs   <- unique(candids$approx)

  idx.list <- list()
  for(i in 1:length(ms)) idx.list[[i]] <- which(candids$m == ms[i])

  visdat.s.str <- list()
  visdat.s.ic0 <- list()
  for(i in 1:length(ms)) {

    idx             <- idx.list[[i]]
    visdat.s.str[[i]]  <- data.frame(gp = "tr", x = output.s$simout[[ idx[1] ]]$alpha, y = output.s$simout[[ idx[1] ]]$post.norm[, 1])
    visdat.s.ic0[[i]]  <- data.frame(gp = "tr", x = output.s$simout.ic0[[ idx[1] ]]$alpha, y = output.s$simout.ic0[[ idx[1] ]]$post.norm[, 1])

    for(j in 1:length(approxs)) {

      tempdat         <- data.frame(gp = approxs[j], x = output.s$simout[[ idx[j] ]]$alpha, y = output.s$simout[[ idx[j] ]]$post.norm[, 2])
      visdat.s.str[[i]]  <- rbind(visdat.s.str[[i]], tempdat)

      tempdat         <- data.frame(gp = approxs[j], x = output.s$simout.ic0[[ idx[j] ]]$alpha, y = output.s$simout.ic0[[ idx[j] ]]$post.norm[, 2])
      visdat.s.ic0[[i]]  <- rbind(visdat.s.ic0[[i]], tempdat)
    }

    visdat.s.str[[i]]  <- visdat.s.str[[i]] %>% filter(gp != "ex")
    visdat.s.ic0[[i]]  <- visdat.s.ic0[[i]] %>% filter(gp != "ex")
  }

  visdat.t.str <- list()
  visdat.t.ic0 <- list()
  for(i in 1:length(ms)) {

    idx             <- idx.list[[i]]
    visdat.t.str[[i]]  <- data.frame(gp = "tr", x = output.t$simout[[ idx[1] ]]$alpha, y = output.t$simout[[ idx[1] ]]$post.norm[, 1])
    visdat.t.ic0[[i]]  <- data.frame(gp = "tr", x = output.t$simout.ic0[[ idx[1] ]]$alpha, y = output.t$simout.ic0[[ idx[1] ]]$post.norm[, 1])

    for(j in 1:length(approxs)) {

      tempdat         <- data.frame(gp = approxs[j], x = output.t$simout[[ idx[j] ]]$alpha, y = output.t$simout[[ idx[j] ]]$post.norm[, 2])
      visdat.t.str[[i]]  <- rbind(visdat.t.str[[i]], tempdat)

      tempdat         <- data.frame(gp = approxs[j], x = output.t$simout.ic0[[ idx[j] ]]$alpha, y = output.t$simout.ic0[[ idx[j] ]]$post.norm[, 2])
      visdat.t.ic0[[i]]  <- rbind(visdat.t.ic0[[i]], tempdat)
    }

    visdat.t.str[[i]]  <- visdat.t.str[[i]] %>% filter(gp != "ex")
    visdat.t.ic0[[i]]  <- visdat.t.ic0[[i]] %>% filter(gp != "ex")
  }

  for(i in 1:length(ms)) {

    visdat.s.str[[i]]$m <- ms[i]
    visdat.s.ic0[[i]]$m <- ms[i]
    visdat.t.str[[i]]$m <- ms[i]
    visdat.t.ic0[[i]]$m <- ms[i]

    visdat.s.str[[i]]$target <- "spatial range"
    visdat.s.ic0[[i]]$target <- "spatial range"

    visdat.t.str[[i]]$target <- "temporal range"
    visdat.t.ic0[[i]]$target <- "temporal range"

    visdat.s.str[[i]]$method <- "naive"
    visdat.t.str[[i]]$method <- "naive"

    visdat.s.ic0[[i]]$method <- "ic0-based"
    visdat.t.ic0[[i]]$method <- "ic0-based"
  }

  visdat <- bind_rows(visdat.s.str, visdat.t.str, visdat.s.ic0, visdat.t.ic0)
  visdat$scen <- label

  visdat$method <- factor(visdat$method, levels = c("naive", "ic0-based"))

  result <- visdat %>% ggplot(aes(x = x, y = y, col = gp)) +
    geom_line(size = size.line, alpha = alpha) +
    scale_color_manual(values = color, labels = legend) +
    xlab(NULL) + ylab(NULL) +
    coord_cartesian(xlim = NULL, ylim = NULL) +
    theme(axis.title.x = element_text(size = size.lab),
          axis.text.x = element_text(size = size.text - 8, angle = 90),
          axis.title.y = element_text(size = size.lab),
          axis.text.y = element_text(size = size.text - 8, angle = 90),
          strip.text = element_text(size = size.lab),
          legend.title = element_blank(),
          legend.text = element_text(size = size.legend, margin = margin(r = 25, unit = 'pt')),
          legend.direction = 'horizontal',
          legend.spacing.x = unit(15, 'pt'),
          plot.margin = unit(size.margin, "pt")) +
    facet_nested(method ~ scen + target + m, scales = "free_x")

  return( result )
}

####################################################################################

ps_random <- plots_post(output.s = output.sptm.posterior.random.srange, output.t = output.sptm.posterior.random.trange, label = "random")
ps_monitoring <- plots_post(output.s = output.sptm.posterior.monitoring.srange, output.t = output.sptm.posterior.monitoring.trange, label = "station")
ps_satellite <- plots_post(output.s = output.sptm.posterior.satellite.srange, output.t = output.sptm.posterior.satellite.trange, label = "satellite")

tmp       <- ggplot_gtable(ggplot_build(ps_random))
leg       <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
mylegend  <- tmp$grobs[[leg]]

rm(tmp, leg,
   output.sptm.posterior.random.srange, output.sptm.posterior.random.trange,
   output.sptm.posterior.monitoring.srange, output.sptm.posterior.monitoring.trange,
   output.sptm.posterior.satellite.srange, output.sptm.posterior.satellite.trange)

####################################################################################

result    <- grid.arrange(mylegend, arrangeGrob(ps_random + theme(legend.position="none"),
                                                ps_monitoring + theme(legend.position="none"),
                                                ps_satellite + theme(legend.position="none"),
                                                nrow=3), nrow=2, heights=c(1, 20))

####################################################################################

ggplot2::ggsave("simulation_spacetime_noise.pdf", result, width = 12, height = 14)

