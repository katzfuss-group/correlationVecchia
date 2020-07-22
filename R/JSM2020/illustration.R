####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: 
###
####################################################################################

library(dplyr)
library(ggplot2)
library(correlationVecchia)
library(gridExtra)

load(file = "R/JSM2020/aniso.RData")

out$kldiv <- out$kldiv[, c(1, 2, 5, 7)]

vdat1   <- out$vars %>% left_join(out$kldiv, by = "index") %>% filter(a == 10) %>% select(-a)
vdat2   <- out$vars %>% left_join(out$kldiv, by = "index") %>% filter(m == 30) %>% select(-m)

colnames(vdat1)[2] <- "m"
colnames(vdat2)[2] <- "a"

vis     <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("Euclidean-based Vecchia approximation", "Correlation-based Vecchia approximation", "Coordinate-based Vecchia approximation = Vecchia's original approach"), color = c("#377EB8", "#E41A1C", "#984EA3"), shape = c(15, 17, 16))


load(file = "R/JSM2020/aniso.RData")

out$kldiv <- out$kldiv[, c(1, 2, 5, 6)]

vdat1   <- out$vars %>% left_join(out$kldiv, by = "index") %>% filter(a == 10) %>% select(-a)
vdat2   <- out$vars %>% left_join(out$kldiv, by = "index") %>% filter(m == 30) %>% select(-m)

colnames(vdat1)[2] <- "size of conditioning sets (m)"
colnames(vdat2)[2] <- "degree of anisotropy (a)"

vis     <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("Euclidean-based Vecchia", "Correlation-based Vecchia", "Coordinate-based Vecchia = Vecchia's original approach"), color = c("#377EB8", "#E41A1C", "#984EA3"), shape = c(15, 17, 16))

