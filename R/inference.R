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
#' @param mcmc.length Number of MCMC iterations to run the algorithm for
#' @param mcmc.thinning MCMC thinning interval between two sampled iterations
#' @param mcmc.tree.updates Number of transmission tree updates per parameter update
#' @param mcmc.cov Initial proposal covariance
#' @param init.r Starting value of offspring distribution parameter r
#' @param init.p Starting value of offspring distribution parameter p
#' @param init.pi Starting value of sampling proportion pi
#' @param init.kappa Starting value for the initial within-host population size
#' @param init.lambda Starting value for the within-host population growth rate
#' @param init.ctree Starting value for the combined tree
#' @param update.r Whether to update the parameter r
#' @param update.p Whether to update the parameter p
#' @param update.pi Whether to update the parameter pi
#' @param update.kappa Whether to update the parameter kappa
#' @param update.lambda Whether to update the parameter lambda
#' @param update.ctree Whether to update the transmission tree
#' @param r.shape Shape parameter for the Gamma prior of parameter r
#' @param r.scale Scale parameter for the Gamma prior of parameter r
#' @param p.shape1 Shape1 parameter for the Beta prior of parameter p
#' @param p.shape2 Shape2 parameter for the Beta prior of parameter p
#' @param pi.shape1 Shape1 parameter for the Beta prior of parameter pi
#' @param pi.shape2 Shape2 parameter for the Beta prior of parameter pi
#' @param kappa.shape Shape parameter for the Gamma prior of parameter kappa
#' @param kappa.scale Scale parameter for the Gamma prior of parameter kappa
#' @param lambda.shape Shape parameter for the Gamma prior of parameter lambda
#' @param lambda.scale Scale parameter for the Gamma prior of parameter lambda
#' @param obs.start Start date for observations
#' @param obs.end End date for observations
#' @param grid.delta Grid resolution for approximating exclusion probabilities
#' @param verbose Whether or not to use verbose mode (default is false)
#' @export
inferTTreeM <- function(ptree,
                        w.shape = 2,
                        w.scale = 1,
                        ws.shape = NA,
                        ws.scale = NA,
                        w.mean = NA,
                        w.std = NA,
                        ws.mean = NA,
                        ws.std = NA,
                        mcmc.length = 12000,
                        mcmc.thinning = 1,
                        mcmc.tree.updates = NA,
                        mcmc.cov = NA,
                        init.r = 2,
                        init.p = 0.5,
                        init.pi = 0.5,
                        init.kappa = 0.5,
                        init.lambda = 0.5,
                        init.ctree = NA,
                        update.r = T,
                        update.p = F,
                        update.pi = T,
                        update.kappa = T,
                        update.lambda = T,
                        update.ctree = T,
                        r.shape = 1,
                        r.scale = 1,
                        p.shape1 = 1,
                        p.shape2 = 1,
                        pi.shape1 = 1,
                        pi.shape2 = 1,
                        kappa.shape = 1,
                        kappa.scale = 1,
                        lambda.shape = 1,
                        lambda.scale = 1,
                        obs.start = -Inf,
                        obs.end = NA,
                        grid.delta = NA,
                        verbose = F) {

  # Ensure that all leaves have unique times
  ptree$ptree[, 1] <- ptree$ptree[, 1] + runif(nrow(ptree$ptree)) * 1e-10

  # Ensure branch lengths of ptree are positive
  for (i in (ceiling(nrow(ptree$ptree) / 2) + 1):nrow(ptree$ptree)) {

    for (j in 2:3) {

      if (ptree$ptree[ptree$ptree[i, j], 1] - ptree$ptree[i, 1] < 0) {

        stop("The phylogenetic tree contains negative branch lengths!")

      }

    }

  }

  # Determine vector of primary observation times
  prim_obs_times <- calc_prim_obs(ptree)

  # Ensure observation start time is consistent with observation times
  if (obs.start > min(ptree$ptree[which(ptree$ptree[, 2] == 0), 1])) {

    stop('The parameter obs.start cannot be later than any observation dates')

  }

  # If no observation end date is provided, set to approximately 00:00 today
  if (is.na(obs.end)) {

    obs.end <- 1970 + as.numeric(Sys.Date()) / 365.25

  }

  # Ensure observation end time is consistent with primary observation times
  if (obs.end < max(prim_obs_times)) {

    stop('The parameter obs.end cannot be earlier than any primary observation dates')

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

  if (is.na(mcmc.tree.updates)) {

    mcmc.tree.updates <- length(prim_obs_times)

  }

  if (is.na(mcmc.cov)) {

    mcmc.cov <- diag(c(0.5 * as.numeric(update.r),
                       0.25 * as.numeric(update.p),
                       0.25 * as.numeric(update.pi),
                       0.25 * as.numeric(update.kappa),
                       0.25 * as.numeric(update.lambda)))

  }

  if (is.na(init.ctree)) {

    ctree <- init_ctree(ptree)

  } else {

    ctree <- init.ctree

  }

  if (is.na(grid.delta)) {

    grid.delta <- 0.001 * (max(ptree$ptree[, 1]) - min(ptree$ptree[, 1]))

  }

  off.r <- init.r
  off.p <- init.p
  pi <- init.pi
  kappa <- init.kappa
  lambda <- init.lambda

  parms <- c(off.r, off.p, pi, kappa, lambda)

  ttree <- extractTTreeM(ctree)

  record <- vector('list', mcmc.length / mcmc.thinning)

  grid <- seq(obs.end, min(ttree$ttree[, 1]) - 0.5 * 1 - grid.delta, by = - grid.delta)

  fn_list <- num_approx_disc(grid, grid.delta, off.r, off.p, pi, w.shape, w.scale,
                             ws.shape, ws.scale, obs.start, obs.end)

  pTTree <- log_lik_ttree(ttree, grid, fn_list, off.r, off.p, pi, w.shape, w.scale,
                          ws.shape, ws.scale, obs.start, obs.end, grid.delta, NA)

  pPTree <- log_lik_ptree_given_ctree(ctree, kappa, lambda, NA)

  if (verbose == F) {

    pb <- utils::txtProgressBar(min = 0, max = mcmc.length, style = 3)

  }

  for (i in 1:mcmc.length) {

    if (update.ctree) {

      if (verbose) {

        message("Proposing ttree update")

      }

      tree_prop_count <- c("add" = 0,
                           "remove" = 0,
                           "move" = 0)

      tree_acc_count <- c("add" = 0,
                          "remove" = 0,
                          "move" = 0)

      tree_acc_rates <- c("add" = 0,
                          "remove" = 0,
                          "move" = 0)

      for (j in 1:mcmc.tree.updates) {

        u <- runif(1)
        if (u < 1 / 3) {

          proptree <- add_transmission(ctree = ctree)
          prop_type <- 1
          tree_prop_count$add <- tree_prop_count$add + 1

        } else if (u < 2 / 3) {

          proptree <- remove_transmission(ctree = ctree)
          prop_type <- 2
          tree_prop_count$remove <- tree_prop_count$remove + 1

        } else {

          proptree <- remove_add(ctree = ctree)
          prop_type <- 3
          tree_prop_count$move <- tree_prop_count$move + 1

        }

        if (proptree$is_possible == 1) {

          ctree2 <- proptree$ctree

          ttree2 <- extractTTreeM(ctree2)

          pTTree_part <- log_lik_ttree(ttree, grid, fn_list, off.r, off.p, pi, w.shape, w.scale, ws.shape,
                                       ws.scale, obs.start, obs.end, grid.delta, proptree$curr_hosts)

          pTTree_part2 <- log_lik_ttree(ttree2, grid, fn_list, off.r, off.p, pi, w.shape, w.scale, ws.shape,
                                        ws.scale, obs.start, obs.end, grid.delta, proptree$prop_hosts)

          pPTree_part <- log_lik_ptree_given_ctree(ctree, kappa, lambda, proptree$curr_hosts)

          pPTree_part2 <- log_lik_ptree_given_ctree(ctree2, kappa, lambda, proptree$prop_hosts)

          if (log(runif(1)) < (pTTree_part2 + pPTree_part2 + proptree$rev_density -
                               pTTree_part - pPTree_part - proptree$prop_density)) {

            ctree <- ctree2
            ttree <- ttree2

            pTTree <- pTTree + pTTree_part2 - pTTree_part
            pPTree <- pPTree + pPTree_part2 - pPTree_part

            if (prop_type == 1) {

              tree_acc_count$add <- tree_acc_count$add + 1

            } else if (prop_type == 2) {

              tree_acc_count$remove <- tree_acc_count$remove + 1

            } else {

              tree_acc_count$move <- tree_acc_count$move + 1

            }

            if (grid[length(grid)] > (min(ttree$ttree[, 1]) - 0.5 * 1)) {

              grid <- seq(obs.end, min(ttree$ttree[, 1]) - 0.5 * 1 - grid.delta, by = - grid.delta)

              fn_list <- num_approx_disc(grid, grid.delta, off.r, off.p, pi, w.shape, w.scale,
                                         ws.shape, ws.scale, obs.start, obs.end)

            }

          }

        }

      }

      if (i %% thinning == 0) {

        tree_acc_rate$add <- tree_add_count$add / tree_prop_count$add
        tree_acc_rate$remove <- tree_add_count$remove / tree_prop_count$remove
        tree_acc_rate$move <- tree_add_count$move / tree_prop_count$move

        record[[i / thinning]]$tree.prop.counts <- tree_prop_count
        record[[i / thinning]]$tree.acc.rates <- tree_acc_rates

      }

#      pTTree <- log_lik_ttree(ttree, grid, fn_list, off.r, off.p, pi, w.shape, w.scale, ws.shape, ws.scale, obs.start, obs.end, grid.delta, NA)

#      pPTree <- log_lik_ptree_given_ctree(ctree, kappa, lambda, NA)


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

    if (i %% thinning == 0) {

      if (verbose == F) {

        utils::setTxtProgressBar(pb, i)

      }

      if (verbose==T) {

        message(sprintf('it = %d, lm_const = %f, lm_rate = %f, off.r = %f, off.p = %f, pi = %f, Prior = %e, Likelihood = %e, nind = %d', i, lm_const, lm_rate, off.r, off.p, pi, pTTree, pPTree, nrow(ttree$ttree)))

      }

      rec_ctree <- trim_root(ctree)
      rec_ttree <- extractTTreeM(rec_ctree)

      rec_pTTree <- log_lik_ttree(rec_ttree, grid, fn_list, off.r, off.p, pi, w.shape, w.scale, ws.shape,
                                  ws.scale, dateS, dateT, delta, NA)
      rec_pPTree <- log_lik_ptree_given_ctree(rec_ctree, lm_const, lm_rate, NA)


      record[[i / thinning]]$ctree <- rec_ctree
      record[[i / thinning]]$pTTree <- rec_pTTree
      record[[i / thinning]]$pPTree <- rec_pPTree
      record[[i / thinning]]$lm_const <- lm_const
      record[[i / thinning]]$lm_rate <- lm_rate
      record[[i / thinning]]$off.r <- off.r
      record[[i / thinning]]$off.p <- off.p
      record[[i / thinning]]$pi <- pi
      record[[i / thinning]]$w.shape <- w.shape
      record[[i / thinning]]$w.scale <- w.scale
      record[[i / thinning]]$ws.shape <- ws.shape
      record[[i / thinning]]$ws.scale <- ws.scale

      record[[i / thinning]]$source <- rec_ctree$ctree[rec_ctree$ctree[which(rec_ctree$ctree[, 4] == 0), 2], 4]
      if (record[[i / thinning]]$source <= length(rec_ctree$nam)) {

        record[[i / thinning]]$source <- rec_ctree$nam[record[[i / thinning]]$source]

      } else {

        record[[i / thinning]]$source <- 'Unsampled'

      }

    }

  }#End of main MCMC loop

  class(record)<-'resTransPhyloM'
  return(record)

}
