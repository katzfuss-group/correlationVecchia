<!-- README.md is generated from README.Rmd. Please edit that file -->

# correlationVecchia

This repository (or R package) includes multiple source code files for
correlation-based Vecchia with MJ

## Installation

    # install.packages("devtools")
    devtools::install_github("katzfuss-group/correlationVecchia")

## Example

    locs    <- matrix(runif(400), 200, 2)
    m       <- 10
    cormat  <- cov_expo_iso(locs = locs, covparms = c(1, 0.1))

    ord     <- GPvecchia::order_maxmin_exact(locs)
    all.equal(ord, order_maxmin_euclidean(locs = locs, initial.pt = ord[1]))
    all.equal(ord, order_maxmin_correlation(locs = locs, d.inv = cormat, initial.pt = ord[1]))

    locsord <- locs[ord, , drop = FALSE]
    corord  <- cormat[ord, ord]

    all.equal(GpGp::find_ordered_nn_brute(locs = locsord, m = m),
              conditioning_m_Rcpp(m = m, d = 1 - corord) + 1)
    all.equal(GpGp::find_ordered_nn(locs = locsord, m = m),
              conditioning_m_Rcpp(m = m, d = 1 - corord) + 1)

    ### Example

    n             <- 20^2
    m             <- 20
    locs          <- matrix(runif(n * 2, 0, 1), n, 2)
    covparms      <- c(1, 0.1, 10)

    # true cov matrix
    covmat <- cov_expo_aniso(locs, covparms)

    # Visualize the process
    y <- as.numeric(t(chol(covmat)) %*% rnorm(n))
    fields::quilt.plot(locs[,1], locs[,2], y)

    out01 <- cvecchia_m_specify(locs = locs, m = m, initial.pt = NULL,
                                covmodel = covmat, covparms = covparms,
                                abs.corr = FALSE, coordinate = c(1),
                                ordering = "coord", ordering.method = "euc",
                                conditioning = "NN", conditioning.method = "euc")

    out02 <- cvecchia_m_specify(locs = locs, m = m, initial.pt = NULL,
                                covmodel = covmat, covparms = covparms,
                                abs.corr = FALSE, coordinate = c(2),
                                ordering = "coord", ordering.method = "euc",
                                conditioning = "NN", conditioning.method = "euc")

    out03 <- cvecchia_m_specify(locs = locs, m = m, initial.pt = NULL,
                                covmodel = covmat, covparms = covparms,
                                abs.corr = FALSE, coordinate = NULL,
                                ordering = "MM", ordering.method = "euc",
                                conditioning = "NN", conditioning.method = "euc")

    out04 <- cvecchia_m_specify(locs = locs, m = m, initial.pt = NULL,
                                covmodel = covmat, covparms = covparms,
                                abs.corr = FALSE, coordinate = NULL,
                                ordering = "MM", ordering.method = "euc",
                                conditioning = "NN", conditioning.method = "cor")

    out05 <- cvecchia_m_specify(locs = locs, m = m, initial.pt = NULL,
                                covmodel = covmat, covparms = covparms,
                                abs.corr = FALSE, coordinate = NULL,
                                ordering = "MM", ordering.method = "cor",
                                conditioning = "NN", conditioning.method = "euc")

    out06 <- cvecchia_m_specify(locs = locs, m = m, initial.pt = NULL,
                                covmodel = covmat, covparms = covparms,
                                abs.corr = FALSE, coordinate = NULL,
                                ordering = "MM", ordering.method = "cor",
                                conditioning = "NN", conditioning.method = "cor")

    kls.coord.x    <- performance(vecchia.approx = out01, locs = locs,
                                  covmodel = cov_expo_aniso, covparms = covparms)
    kls.coord.y    <- performance(vecchia.approx = out02, locs = locs,
                                  covmodel = cov_expo_aniso, covparms = covparms)
    kls.euc.euc    <- performance(vecchia.approx = out03, locs = locs,
                                  covmodel = cov_expo_aniso, covparms = covparms)
    kls.euc.cor    <- performance(vecchia.approx = out04, locs = locs,
                                  covmodel = cov_expo_aniso, covparms = covparms)
    kls.cor.euc    <- performance(vecchia.approx = out05, locs = locs,
                                  covmodel = cov_expo_aniso, covparms = covparms)
    kls.cor.cor    <- performance(vecchia.approx = out06, locs = locs,
                                  covmodel = cov_expo_aniso, covparms = covparms)

    barplot(log10(c(kls.coord.x, kls.coord.y, kls.euc.euc,
                    kls.euc.cor, kls.cor.euc, kls.cor.cor)),
            names.arg = c("X-ord + E-NN", "Y-ord + E-NN",
                          "E-MM + E-NN", "E-MM + C-NN",
                          "C-MM + E-NN", "C-MM + C-NN"),
            main = "Vecchia Approximations", ylab = "log10-scale KL divergence")

## Reproducing results

These R files reproduce the figures and results from Kang and Katzfuss
(2021): - Figure 1: simulation_visualization_ordering_and_conditioning.R
- Figure 2: simulation_knownCovparms_05242021.R -\>
simulation_visualization_anisononst.R - Figure 3:
simulation_knownCovparms_05242021.R -\>
simulation_visualization_multivariate.R - Figure 4:
simulation_visualization_spacetime_scenarios.R - Figure 5:
simulation_knownCovparms_05242021.R -\>
simulation_visualization_spacetime.R - Figure 6:
simulation_noneuc_09272021.R -\> simulation_visualization_noneuc.R -
Figure 7: simulation_fisherScoring_05242021.R -\>
simulation_visualization_fisher.R - Figure 8:
simulation_prediction_05242021.R -\>
simulation_visualization_prediction_2.R - Figure 9:
simulation_bayesianPosterior_nugget_05242021.R -\>
simulation_visualization_noise.R - Figures 10 and 11 and Table 1:
realdata_1\_processing_01012022.R -\> realdata_2\_estimation_01012022.R
-\> realdata_3\_performance_01012022.R -\>
realdata_4\_visualization_01012022.R

## Reference

Kang, M., & Katzfuss, M. (2021). Correlation-based sparse inverse
Cholesky factorization for fast Gaussian-process inference.
[*arXiv:2112.14591v1*](https://arxiv.org/abs/2112.14591).
