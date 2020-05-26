library(igraph)
library(dplyr)
library(ggplot2)
library(ggnetwork)

set.seed(20200524)

n     <- 30
m     <- as.integer(sqrt(n))

# locs      <- matrix(runif(2 * n), nrow = n, ncol = 2)
# locs      <- as.matrix(expand.grid(seq(from = 0.2, to = 0.8, length.out = m), seq(from = 0.1, to = 0.9, length.out = m)))
# locs      <- locs + matrix(0.15 * runif(2 * m^2, min = 0, max = 1), nrow = m^2, ncol = 2)
# n         <- nrow(locs)

locs      <- as.matrix(expand.grid(seq(from = 0.2, to = 0.8, length.out = m), seq(from = 0.1, to = 0.9, length.out = m)))
noise     <- matrix(0.15 * runif(2 * m^2, min = 0, max = 1), nrow = m^2, ncol = 2)
noise     <- noise - mean(noise)
locs      <- locs + noise
n         <- nrow(locs)

#####

# edmat     <- fields::rdist(x1 = locs, x2 = NULL)
# eord      <- order_maxmin_euclidean(locs = locs, initial.pt = NULL)
# locseord  <- locs[eord, , drop = FALSE]
# econd     <- GpGp::find_ordered_nn(locs = locseord, m = m)
# 
# elist <- list()
# econd <- econd[, -1, drop = FALSE]
# for(i in 1:n) {
#   elist[[i]] <- sort(which(econd == i, arr.ind = TRUE)[, 1])
# }
# 
# econd2 <- matrix(NA, nrow = 1, ncol = 2)
# for(i in 1:n) {
#   if(length(elist[[i]]) != 0) {
#     econd2 <- rbind(econd2, cbind(i, elist[[i]]))
#   }
# } ; econd2 <- econd2[-1, ]
# 
# adjmat.firsttry <- get.adjacency(graph.edgelist(econd2, directed = TRUE)) %>% as.matrix()

#####

edmat     <- fields::rdist(x1 = locs, x2 = NULL)
eord      <- correlationVecchia::order_maxmin_euclidean(locs = locs, initial.pt = NULL)
locseord  <- locs[eord, , drop = FALSE]
econd     <- GpGp::find_ordered_nn(locs = locseord, m = m)

econd.table   <- cbind(as.vector(econd[, -1, drop = FALSE]), rep(econd[, 1], times = m)) %>% na.omit()
adjmat        <- get.adjacency(graph.edgelist(econd.table, directed = TRUE)) %>% as.matrix()

pivind        <- as.integer(n * 0.5)
cond.piv      <- which(adjmat[, pivind] == 1)

adjmat[, pivind] <- adjmat[, pivind] * 10 

net           <- graph.adjacency(adjmat, mode = "directed", weighted = TRUE, diag = FALSE)

V(net)$gp <- rep("lightsteelblue2", times = n)
V(net)$gp[cond.piv] <- "gold"
V(net)$gp[pivind] <- "tomato"

plot.igraph(net, vertex.label = seq(n), layout = locseord, 
            edge.width = E(net)$weight / 3, edge.color =  c("gray80", "black")[(E(net)$weight==10)+1], edge.arrow.size = 0.,
            vertex.size = 10, vertex.color = V(net)$gp)
title(main="E-MM + E-NN")
axis(1, labels=FALSE, tick=TRUE)
axis(2, labels=FALSE, tick=TRUE)

#####

ggnet <- ggnetwork(net, layout = locseord)

ggnet$vertex.names <- NA
ggnet$vertex.names[which(is.na(ggnet$weight))] <- 1:n

ggnet$edge.names <- "ordinary"
ggnet$edge.names[which(ggnet$weight == 10)] <- "important"

p.eucl <- ggplot(data = ggnet) + 
  geom_edges(aes(x = x, y = y, xend = xend, yend = yend, color = edge.names), arrow = arrow(length = unit(8, "pt"), type = "closed")) + 
  geom_nodes(aes(x = x, y = y, color = gp), size = 8) +
  geom_nodetext(aes(x = x, y = y, label = vertex.names)) +
  scale_color_manual(values = c("gold", "black", "lightsteelblue2", "grey80", "tomato")) + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), legend.position = "none", axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12)) +
  xlab("") + ylab("")

p.eucl

#####

rho       <- correlationVecchia::cov_expo_aniso(locs, covparms = c(1, 1, 10))
cord      <- correlationVecchia::order_maxmin_correlation_inverseDist(locs = locs, d.inv = rho, initial.pt = eord[1])
locscord  <- locs[cord, , drop = FALSE]
ccond     <- correlationVecchia::conditioning_nn(m = m, d = 1 - rho[cord, cord])

ccond.table   <- cbind(as.vector(ccond[, -1, drop = FALSE]), rep(ccond[, 1], times = m)) %>% na.omit()
adjmat        <- get.adjacency(graph.edgelist(ccond.table, directed = TRUE)) %>% as.matrix()

pivind        <- which(fields::rdist(x1 = locseord[pivind, , drop = FALSE], x2 = locscord) == 0)
cond.piv      <- which(adjmat[, pivind] == 1)

adjmat[, pivind] <- adjmat[, pivind] * 10 

net           <- graph.adjacency(adjmat, mode = "directed", weighted = TRUE, diag = FALSE)

V(net)$gp <- rep("lightsteelblue2", times = n)
V(net)$gp[cond.piv] <- "gold"
V(net)$gp[pivind] <- "tomato"

plot.igraph(net, vertex.label = seq(n), layout = locscord, 
            edge.width = E(net)$weight / 3, edge.color =  c("gray80", "black")[(E(net)$weight==10)+1], edge.arrow.size = 0.5,
            vertex.size = 10, vertex.color = V(net)$gp)
title(main="c-MM + c-NN")
axis(1, labels=FALSE, tick=TRUE)
axis(2, labels=FALSE, tick=TRUE)

#####

# edgelist <- get.edgelist(net)
# edges <- data.frame(locscord[edgelist[,1],], locscord[edgelist[,2],])
# colnames(edges) <- c("x","y","xend","yend")
# edges$type <- as.character(E(net)$weight)
# 
# nodes <- as.data.frame(locscord)
# colnames(nodes) <- c("x", "y")
# nodes$type <- V(net)$gp
# nodes$ord <- 1:n
# 
# ggplot() + geom_segment(aes(x = x, y = y, xend = xend, yend = yend), data = edges, size = 0.5, colour="grey", arrow = arrow(type = "closed", length = unit(0.2, "cm"))) + geom_point(aes(x = x, y = y), data = nodes)
# 
# ggplot() + geom_edges(aes(x = x, y = y, xend = xend, yend = yend, linetype = type), data = edges, color = "grey50", arrow = arrow(length = unit(8, "pt"), type = "closed")) + geom_nodes(aes(x = x, y = y, color = type), data = nodes, size = 6) + geom_nodetext(aes(x = x, y = y, label = ord), data = nodes)

#####

ggnet <- ggnetwork(net, layout = locscord)

ggnet$vertex.names <- NA
ggnet$vertex.names[which(is.na(ggnet$weight))] <- 1:n

ggnet$edge.names <- "ordinary"
ggnet$edge.names[which(ggnet$weight == 10)] <- "important"

p.corr <- ggplot(data = ggnet) + 
  geom_edges(aes(x = x, y = y, xend = xend, yend = yend, color = edge.names), arrow = arrow(length = unit(8, "pt"), type = "closed")) + 
  geom_nodes(aes(x = x, y = y, color = gp), size = 8) +
  geom_nodetext(aes(x = x, y = y, label = vertex.names)) +
  scale_color_manual(values = c("gold", "black", "lightsteelblue2", "grey80", "tomato")) + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), legend.position = "none", axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12)) +
  xlab("") + ylab("")

p.corr

#####

ggplot2::ggsave("Euclidean_MN.pdf", p.eucl, width = 7.6, height = 7.6)
ggplot2::ggsave("correlation_MN.pdf", p.corr, width = 7.6, height = 7.6)
