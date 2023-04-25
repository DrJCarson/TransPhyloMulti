#' Infer transmission tree given a phylogenetic tree
#' @param ptree Phylogenetic tree
#' @param w.shape Shape parameter of the Gamma distribution representing the generation time
#' @param w.scale Scale parameter of the Gamma distribution representing the generation time
#' @param ws.shape Shape parameter of the Gamma distribution representing the sampling time
#' @param ws.scale Scale parameter of the Gamma distribution representing the sampling time
#' @param w.mean Mean of the Gamma distribution representing the generation time
#' @param w.std Std of the Gamma distribution representing the generation time
#' @param ws.mean Mean of the Gamma distribution representing the sampling time
#' @param ws.std Std of the Gamma distribution representing the sampling time
#' @param mcmcIterations Number of MCMC iterations to run the algorithm for
#' @param thinning MCMC thinning interval between two sampled iterations
#' @param start_const Starting value of within-host population size
#' @param start_rate Starting value for the within-host population growth rate
#' @param startOff.r Starting value of parameter off.r
#' @param startOff.p Starting value of parameter off.p
#' @param startPi Starting value of sampling proportion pi
#' @param updateconst Whether of not to update the parameter lm_const
#' @param updaterate Whether of not to update the parameter lm_rate
#' @param updateOff.r Whether or not to update the parameter off.r
#' @param updateOff.p Whether or not to update the parameter off.p
#' @param updatePi Whether or not to update the parameter pi
#' @param qconst Proposal kernel range for parameter lm_const
#' @param qrate Proposal kernel range for parameter lm_rate
#' @param qOff.r Proposal kernel range for parameter Off.r
#' @param qOff.p Proposal kernel range for parameter Off.p
#' @param qPi Proposal kernel range for parameter pi
#' @param startCTree Optional combined tree to start from
#' @param updateTTree Whether or not to update the transmission tree
#' @param dateS Date when observations start
#' @param dateT Date when process stops (this can be Inf for fully simulated outbreaks)
#' @param delta Grid precision (smaller is better but slower)
#' @param verbose Whether or not to use verbose mode (default is false)
#' @export
inferTTreeM <- function(ptree, w.shape = 2, w.scale = 1, ws.shape = NA, ws.scale = NA,
                         w.mean = NA, w.std = NA, ws.mean = NA, ws.std = NA, mcmcIterations = 1000,
                         thinning = 1, tree_updates = 1, start_const = 2, start_rate = 2, startOff.r = 1, startOff.p = 0.5,
                         startPi = 0.5, updateconst = TRUE, updaterate = TRUE, updateOff.r = TRUE, updateOff.p = FALSE,
                         updatePi = TRUE, qconst = NA, qrate = NA, qcr = NA, qOff.r = NA, qOff.p = NA, qPi = NA,
                         bw = 0.8, rd = 1, startCTree = NA, updateTTree = TRUE, dateS = -Inf, dateT = Inf,
                        delta = NA, verbose = F) {

  ptree$ptree[, 1] <- ptree$ptree[, 1] + runif(nrow(ptree$ptree)) * 1e-10 #Ensure that all leaves have unique times

  if (dateT < max(ptree$ptree[, 1])) {

    stop('The parameter dateT cannot be smaller than the date of the last sample')

  }

  if (dateS > min(ptree$ptree[which(ptree$ptree[, 2] == 0), 1])) {

    stop('The parameter dateS cannot be greater than the date of the first sample')

  }

  for (i in (ceiling(nrow(ptree$ptree) / 2) + 1):nrow(ptree$ptree)) {

    for (j in 2:3) {

      if (ptree$ptree[ptree$ptree[i, j], 1] - ptree$ptree[i, 1] < 0) {

        stop("The phylogenetic tree contains negative branch lengths!")

      }

    }

  }

  if (!is.na(w.mean) && !is.na(w.std)) {

    w.shape <- w.mean ^ 2 /  w.std ^ 2
    w.scale <- w.std ^ 2 / w.mean

  }

  if (!is.na(ws.mean)&&!is.na(ws.std)) {

    ws.shape <- ws.mean ^ 2 / ws.std ^ 2
    ws.scale <- ws.std ^ 2 / ws.mean

  }

  if (is.na(ws.shape)) {

    ws.shape <- w.shape

  }

  if (is.na(ws.scale)) {

    ws.scale <- w.scale

  }

  if (is.na(qconst)) {

    qconst <- 0.5

  }

  if (is.na(qrate)) {

    qrate <- 0.5

  }

  if (is.na(qcr)) {

    qcr <- 0.5

  }

  if (is.na(qOff.r)) {

    qOff.r <- 0.5

  }

  if (is.na(qOff.p)) {

    qOff.p <- 0.5

  }

  if (is.na(qPi)) {

    qPi <- 0.5

  }

  if (is.na(delta)) {

    delta <- 0.001 * (max(ptree$ptree[, 1]) - min(ptree$ptree[, 1]))

  }

  lm_const <- start_const
  lm_rate <- start_rate
  off.r <- startOff.r
  off.p <- startOff.p
  pi <- startPi

  if (sum(is.na(startCTree))) {

    ctree <- init_ctree(ptree)

  } else {

    ctree <- startCTree

  }

  ttree <- extractTTreeM(ctree)

  record <- vector('list', mcmcIterations / thinning)

  grid <- seq(dateT, min(ttree$ttree[, 1]) - 0.5 * rd - delta, by = - delta)

  fn_list <- num_approx_disc(grid, delta, off.r, off.p, pi, w.shape, w.scale,
                             ws.shape, ws.scale, dateS, dateT)

  pTTree <- log_lik_ttree(ttree, grid, fn_list, off.r, off.p, pi, w.shape, w.scale,
                          ws.shape, ws.scale, dateS, dateT, delta, NA)

  pPTree <- log_lik_ptree_given_ctree(ctree, lm_const, lm_rate, NA)

  if (verbose == F) {

    pb <- utils::txtProgressBar(min = 0, max = mcmcIterations, style = 3)

  }

  for (i in 1:mcmcIterations) {

    if (i %% thinning == 0) {

      if (verbose == F) {

        utils::setTxtProgressBar(pb, i)

      }

      if (verbose==T) {

        message(sprintf('it = %d, lm_const = %f, lm_rate = %f, off.r = %f, off.p = %f, pi = %f, Prior = %e, Likelihood = %e, nind = %d', i, lm_const, lm_rate, off.r, off.p, pi, pTTree, pPTree, nrow(ttree$ttree)))

      }

      record[[i / thinning]]$ctree <- ctree
      record[[i / thinning]]$pTTree <- pTTree
      record[[i / thinning]]$pPTree <- pPTree
      record[[i / thinning]]$lm_const <- lm_const
      record[[i / thinning]]$lm_rate <- lm_rate
      record[[i / thinning]]$off.r <- off.r
      record[[i / thinning]]$off.p <- off.p
      record[[i / thinning]]$pi <- pi
      record[[i / thinning]]$w.shape <- w.shape
      record[[i / thinning]]$w.scale <- w.scale
      record[[i / thinning]]$ws.shape <- ws.shape
      record[[i / thinning]]$ws.scale <- ws.scale

      record[[i / thinning]]$source <- ctree$ctree[ctree$ctree[which(ctree$ctree[, 4] == 0), 2], 4]
      if (record[[i / thinning]]$source <= length(ctree$nam)) {

        record[[i / thinning]]$source <- ctree$nam[record[[i / thinning]]$source]

      } else {

        record[[i / thinning]]$source <- 'Unsampled'

      }

    }

    if (updateTTree) {

      if (verbose) {

        message("Proposing ttree update")

      }

      for (j in 1:tree_updates) {

        u <- runif(1)
        if (u < 1 / 3) {

          proptree <- add_transmission_9(ctree = ctree)
          prop_type <- 1

        } else if (u < 2 / 3) {

          proptree <- remove_transmission_9(ctree = ctree)
          prop_type <- 2

        } else {

          proptree <- remove_add_9(ctree = ctree)
          prop_type <- 3

        }


  #      if (proptree$is_valid == 1 & proptree$is_possible == 1) {

        prop_acc <- 0
        if (proptree$is_possible == 1) {

          ctree2 <- proptree$ctree

          ttree2 <- extractTTreeM(ctree2)

          pTTree_part <- log_lik_ttree(ttree, grid, fn_list, off.r, off.p, pi, w.shape, w.scale, ws.shape,
                                       ws.scale, dateS, dateT, delta, proptree$curr_hosts)

          pTTree_part2 <- log_lik_ttree(ttree2, grid, fn_list, off.r, off.p, pi, w.shape, w.scale, ws.shape,
                                        ws.scale, dateS, dateT, delta, proptree$prop_hosts)

          pPTree_part <- log_lik_ptree_given_ctree(ctree, lm_const, lm_rate, proptree$curr_hosts)

          pPTree_part2 <- log_lik_ptree_given_ctree(ctree2, lm_const, lm_rate, proptree$prop_hosts)

          #pTTree2 <- log_lik_ttree(ttree2, fn_list, off.r, off.p, pi, w.shape, w.scale, ws.shape,
          #                        ws.scale, dateS, dateT, delta)

          #pPTree2 <- log_lik_ptree_given_ctree(ctree2, lm_const, lm_rate)


          if (log(runif(1)) < (pTTree_part2 + pPTree_part2 + proptree$rev_density - pTTree_part - pPTree_part - proptree$prop_density)) {

            ctree <- ctree2
            ttree <- ttree2
  #          pTTree <- pTTree2
 #           pPTree <- pPTree2

            prop_acc <- 1

            if ((grid[length(grid)] - delta) > (min(ttree$ttree[, 1]) - 0.5 * rd - delta)) {

              grid <- seq(dateT, min(ttree$ttree[, 1]) - 0.5 * rd - delta, by = - delta)

              fn_list <- num_approx_disc(grid, delta, off.r, off.p, pi, w.shape, w.scale,
                                         ws.shape, ws.scale, dateS, dateT)

            }

          }

        }

      }

      if (i %% thinning == 0) {

        record[[i / thinning]]$prop_type <- prop_type
        record[[i / thinning]]$prop_acc <- prop_acc

      }

      pTTree <- log_lik_ttree(ttree, grid, fn_list, off.r, off.p, pi, w.shape, w.scale, ws.shape,
                              ws.scale, dateS, dateT, delta, NA)

      pPTree <- log_lik_ptree_given_ctree(ctree, lm_const, lm_rate, NA)


    }

    if (updateconst) {

      lm_const2 <- abs(lm_const + rnorm(1, mean = 0, sd = qconst))

      if (verbose) {

        message(sprintf("Proposing lm_const update %f->%f", lm_const , lm_const2))

      }

      pPTree2 <- log_lik_ptree_given_ctree(ctree, lm_const2, lm_rate, NA)

      if (log(runif(1)) < (pPTree2 - pPTree - lm_const2 + lm_const)) {

        lm_const <- lm_const2
        pPTree <- pPTree2

      }

    }

    if (updaterate) {

      lm_rate2 <- abs(lm_rate + rnorm(1, mean = 0, sd = qrate))

      if (verbose) {

        message(sprintf("Proposing lm_rate update %f->%f", lm_rate , lm_rate2))

      }

      pPTree2 <- log_lik_ptree_given_ctree(ctree, lm_const, lm_rate2, NA)

      if (log(runif(1)) < (pPTree2 - pPTree - lm_rate2 + lm_rate)) {

        lm_rate <- lm_rate2
        pPTree <- pPTree2

      }

    }

    if (updateconst & updaterate) {

      step <- rnorm(1, mean = 0, sd = qcr)

      lm_const2 <- lm_const + step
      lm_rate2 <- lm_rate - step

      if (lm_const2 > 0 & lm_rate2 > 0) {

        pPTree2 <- log_lik_ptree_given_ctree(ctree, lm_const2, lm_rate2, NA)

        if (log(runif(1)) < (pPTree2 - pPTree - lm_const2 + lm_const - lm_rate2 + lm_rate)) {

          lm_const <- lm_const2
          lm_rate <- lm_rate2
          pPTree <- pPTree2

        }

      }

    }

    if (updateOff.r) {

      #Metropolis update for off.r, assuming Exp(1) prior
      off.r2 <- abs(off.r + rnorm(1, mean = 0, sd = qOff.r))

      if (verbose) {

        message(sprintf("Proposing off.r update %f->%f", off.r, off.r2))

      }

      fn_list2 <- num_approx_disc(grid, delta, off.r2, off.p, pi, w.shape, w.scale,
                                 ws.shape, ws.scale, dateS, dateT)

      pTTree2 <- log_lik_ttree(ttree, grid, fn_list2, off.r2, off.p, pi, w.shape, w.scale, ws.shape,
                               ws.scale, dateS, dateT, delta, NA)

      if (log(runif(1)) < (pTTree2 - pTTree - off.r2 + off.r)) {

        off.r <- off.r2
        pTTree <- pTTree2
        fn_list <- fn_list2

      }

    }

    if (updateOff.p) {

      #Metropolis update for off.p, assuming Unif(0,1) prior
      off.p2 <- abs(off.p + (runif(1) - 0.5) * qOff.p)

      if (off.p2 > 1) {

        off.p2 <- 2 - off.p2

      }

      if (verbose) {

        message(sprintf("Proposing off.p update %f->%f", off.p, off.p2))

      }

      fn_list2 <- num_approx_disc(grid, delta, off.r, off.p2, pi, w.shape, w.scale,
                                 ws.shape, ws.scale, dateS, dateT)

      pTTree2 <- log_lik_ttree(ttree, grid, fn_list2, off.r, off.p2, pi, w.shape, w.scale, ws.shape,
                               ws.scale, dateS, dateT, delta, NA)

      if (log(runif(1)) < (pTTree2 - pTTree)) {

        off.p <- off.p2
        pTTree <- pTTree2
        fn_list <- fn_list2

      }

    }

    if (updatePi) {

      #Metropolis update for pi, assuming Unif(0,1) prior
      pi2 <- abs(pi + (runif(1) - 0.5) * qPi)

      if (pi2 > 1) {

        pi2 <- 2 - pi2

      }

      if (verbose) {

        message(sprintf("Proposing pi update %f->%f", pi, pi2))

      }

      fn_list2 <- num_approx_disc(grid, delta, off.r, off.p, pi2, w.shape, w.scale,
                                 ws.shape, ws.scale, dateS, dateT)

      pTTree2 <- log_lik_ttree(ttree, grid, fn_list2, off.r, off.p, pi2, w.shape, w.scale, ws.shape,
                               ws.scale, dateS, dateT, delta, NA)

      if (log(runif(1)) < (pTTree2 - pTTree)) {

        pi <- pi2
        pTTree <- pTTree2

        fn_list <- fn_list2

      }

    }

  }#End of main MCMC loop

  class(record)<-'resTransPhyloM'
  return(record)

}
