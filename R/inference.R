infereTTreeM <- function(ptree, w.shape = 2, w.scale = 1, ws.shape = NA, ws.scale = NA,
                         w.mean = NA, w.std = NA, ws.mean = NA, ws.std = NA, mcmcIterations = 1000,
                         thinning = 1, start_const = 2, start_rate = 2, startOff.r = 1, startOff.p = 0.5,
                         startPi = 0.5, updateconst = TRUE, updaterate = TRUE, updateOff.r = TRUE, updateOff.p = FALSE,
                         updatePi = TRUE, qconst = NA, qrate = NA, qOff.r = NA, qOff.p = NA, qPi = NA,
                         startCTree = NA, updateTTree = TRUE, dateS = -Inf, dateT = Inf, delta = 1 / 365, verbose = F) {

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

  if (is.na(delta_t)) {

    delta_t <- 0.001 * (max(ptree$ptree[, 1]) - min(ptree$ptree[, 1]))

  }

  if (is.na(qconst)) {

    qNeg <- 0.1

  }

  if (is.na(qrate)) {

    qNeg <- 0.1

  }

  if (is.na(qOff.r)) {

    qOff.r <- min(0.5, 5 / length(ptree$nam))

  }

  if (is.na(qOff.p)) {

    qOff.p <- min(0.1, 1 / length(ptree$nam))

  }

  if (is.na(qPi)) {

    qPi <- min(0.1, 1 / length(ptree$nam))

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

  ttree <- extractTTree(ctree)

  record <- vector('list', mcmcIterations / thinning)

  pTTree <- log_lik_ttree(ttree, off.r, off.p, pi, w.shape, w.scale, ws.shape,
                          ws.scale, dateS, dateT, delta)
  pPTree <- log_lik_ptree_given_ctree(ctree, lm_const, lm_rate)

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

      u <- runif(1)
      if (u < 1 / 3) {

        proptree <- add_transmission(ctree)

      } else if (u < 2 / 3) {

        proptree <- remove_transmission(ctree)

      } else {

        proptree <- remove_add_local(ctree)

      }

      ctree2 <- proptree$ctree

      ttree2 <- extractTTreeM(ctree2)

      pTTree2 <- log_lik_ttree(ttree2, off.r, off.p, pi, w.shape, w.scale, ws.shape,
                              ws.scale, dateS, dateT, delta)
      pPTree2 <- log_lik_ptree_given_ctree(ctree2, lm_const, lm_rate)

      if (log(runif(1)) < (pTTree2 + pPTree2 + proptree$prop_density - pTTree - pPTree - proptree$rev_density)) {

        ctree <- ctree2
        ttree <- ttree2
        pTTree <- pTTree2
        pPTree <- pPTree2

      }

    }

    if (updateconst) {

      lm_const2 <- abs(lm_const + rnorm(1, mean = 0, sd = qconst))

      if (verbose) {

        message(sprintf("Proposing lm_const update %f->%f", lm_const , lm_const2))

      }

      pPTree2 <- log_lik_ptree_given_ctree(ctree2, lm_const2, lm_rate)

      if (log(runif(1)) < (pPtree2 - pPtree - lm_const2 + lm_const)) {

        lm_const <- lm_const2
        pPtree <- pPTree2

      }

    }

    if (updaterate) {

      lm_rate2 <- abs(lm_rate + rnorm(1, mean = 0, sd = qrate))

      if (verbose) {

        message(sprintf("Proposing lm_rate update %f->%f", lm_rate , lm_rate2))

      }

      pPTree2 <- log_lik_ptree_given_ctree(ctree2, lm_const, lm_rate2)

      if (log(runif(1)) < (pPtree2 - pPtree - lm_rate2 + lm_rate)) {

        lm_rate <- lm_rate2
        pPtree <- pPTree2

      }

    }

    if (updateOff.r) {

      #Metropolis update for off.r, assuming Exp(1) prior
      off.r2 <- abs(off.r + rnorm(1, mean = 0, sd = qOff.r))

      if (verbose) {

        message(sprintf("Proposing off.r update %f->%f", off.r, off.r2))

      }

      pTTree2 <- log_lik_ttree(ttree2, off.r2, off.p, pi, w.shape, w.scale, ws.shape,
                               ws.scale, dateS, dateT, delta)

      if (log(runif(1)) < (pTTree2 - pTTree - off.r2 + off.r)) {

        off.r <- off.r2
        pTTree <- pTTree2

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

      pTTree2 <- log_lik_ttree(ttree2, off.r, off.p2, pi, w.shape, w.scale, ws.shape,
                               ws.scale, dateS, dateT, delta)

      if (log(runif(1)) < (pTTree2 - pTTree)) {

        off.p <- off.p2
        pTTree <- pTTree2

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

      pTTree2 <- log_lik_ttree(ttree2, off.r, off.p, pi2, w.shape, w.scale, ws.shape,
                               ws.scale, dateS, dateT, delta)

      if (log(runif(1)) < (pTTree2 - pTTree)) {

        pi <- pi2
        pTTree <- pTTree2

      }

    }

  }#End of main MCMC loop

  class(record)<-'resTransPhylo'
  return(record)

}
