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
#' @param mcmc.cov.tr Initial proposal covariance for transmission parameters
#' @param mcmc.cov.coa Initial proposal covariance for coalescent parameters
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
                        mcmc.cov.tr = NA,
                        mcmc.cov.coa = NA,
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

  if (is.na(sum(mcmc.cov.tr))) {

    mcmc.cov.tr <- diag(c(0.5 ^ 2 * as.numeric(update.r),
                          0.25 ^ 2 * as.numeric(update.p),
                          0.25 ^ 2 * as.numeric(update.pi)))

  }

  if (is.na(sum(mcmc.cov.coa))) {

    mcmc.cov.coa <- diag(c(0.1 ^ 2 * as.numeric(update.kappa),
                           0.1 ^ 2 * as.numeric(update.lambda)))

  }

  if (sum(is.na(init.ctree))) {

    ctree <- init_ctree(ptree)

  } else {

    ctree <- init.ctree

  }

  if (is.na(grid.delta)) {

    grid.delta <- 0.001 * (max(ptree$ptree[, 1]) - min(ptree$ptree[, 1]))

  }

  const.pi <- 3.14159265359

  parms.init.tr <- c(r = init.r,
                     p = init.p,
                     pi = init.pi)
  parms.init.coa <- c(kappa = init.kappa,
                      lambda = init.lambda)

  parms.curr.tr <- parms.init.tr
  parms.curr.coa <- parms.init.coa

  ss.a <- 0.234

  update.tr <- c(update.r, update.p, update.pi)
  update.coa <- c(update.kappa, update.lambda)

  ss.d.tr <- sum(as.numeric(update.tr))
  ss.d.coa <- sum(as.numeric(update.coa))

  ss.v0.tr <- ss.d.tr
  ss.v0.coa <- ss.d.coa

  #ss.f <- floor(0.5 * sqrt(2 * (1:mcmc.length)))
  ss.f <- floor(0.5 * (1:mcmc.length))

  ss.min <- 0.1

  ss.c.tr <- 2.38 ^ 2 / ss.d.tr
  ss.c.coa <- 2.38 ^ 2 / ss.d.coa

  ss.lamstart <- 1

  ss.lam.tr <- ss.lamstart
  ss.lam.coa <- ss.lamstart

  ss.nstart <- 5 / (ss.a * (1 - ss.a))

  ss.A <- -qnorm(ss.a / 2)

  ss.del.tr <- (1 - (1 / ss.d.tr)) * ((sqrt(2 * const.pi) * exp(ss.A ^ 2 / 2)) / (2 * ss.A)) + (1 / (ss.d.tr * ss.a * (1 - ss.a)))
  ss.del.coa <- (1 - (1 / ss.d.coa)) * ((sqrt(2 * const.pi) * exp(ss.A ^ 2 / 2)) / (2 * ss.A)) + (1 / (ss.d.coa * ss.a * (1 - ss.a)))

  ttree <- extractTTreeM(ctree)

  record <- vector('list', mcmc.length / mcmc.thinning)

  trace.r <- numeric(mcmc.length)
  trace.p <- numeric(mcmc.length)
  trace.pi <- numeric(mcmc.length)
  trace.kappa <- numeric(mcmc.length)
  trace.lambda <- numeric(mcmc.length)

  grid <- seq(obs.end, min(ttree$ttree[, 1]) - 0.5 * 1 - grid.delta, by = - grid.delta)

  fn_list <- num_approx_disc(grid, init.r, init.p, init.pi, w.shape, w.scale,
                             ws.shape, ws.scale, obs.start, obs.end)

  pTTree <- log_lik_ttree(ttree, grid, fn_list, init.r, init.p, init.pi, w.shape, w.scale,
                          ws.shape, ws.scale, obs.start, obs.end, grid.delta, NA)

  pPTree <- log_lik_ptree_given_ctree(ctree, init.kappa, init.lambda, NA)

  if (verbose == F) {

    pb <- utils::txtProgressBar(min = 0, max = mcmc.length, style = 3)

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

  for (i in 1:mcmc.length) {

    if (update.ctree) {

      if (verbose) {

        message("Proposing ttree update")

      }

      tree_prop_count["add"] <- 0
      tree_prop_count["remove"] <- 0
      tree_prop_count["move"] <- 0

      tree_acc_count["add"] <- 0
      tree_acc_count["remove"] <- 0
      tree_acc_count["move"] <- 0

      for (j in 1:mcmc.tree.updates) {

        u <- runif(1)
        if (u < 1 / 3) {

          proptree <- add_transmission(ctree = ctree)
          prop_type <- 1
          tree_prop_count["add"] <- tree_prop_count["add"] + 1

        } else if (u < 2 / 3) {

          proptree <- remove_transmission(ctree = ctree)
          prop_type <- 2
          tree_prop_count["remove"] <- tree_prop_count["remove"] + 1

        } else {

          proptree <- remove_add(ctree = ctree)
          prop_type <- 3
          tree_prop_count["move"] <- tree_prop_count["move"] + 1

        }

        if (proptree$is_possible == 1) {

          ctree2 <- proptree$ctree

          ttree2 <- extractTTreeM(ctree2)

          pTTree_part <- log_lik_ttree(ttree, grid, fn_list, parms.curr.tr[["r"]], parms.curr.tr[["p"]], parms.curr.tr[["pi"]], w.shape, w.scale, ws.shape,
                                       ws.scale, obs.start, obs.end, grid.delta, proptree$curr_hosts)

          pTTree_part2 <- log_lik_ttree(ttree2, grid, fn_list, parms.curr.tr[["r"]], parms.curr.tr[["p"]], parms.curr.tr[["pi"]], w.shape, w.scale, ws.shape,
                                        ws.scale, obs.start, obs.end, grid.delta, proptree$prop_hosts)

          pPTree_part <- log_lik_ptree_given_ctree(ctree, parms.curr.coa[["kappa"]], parms.curr.coa[["lambda"]], proptree$curr_hosts)

          pPTree_part2 <- log_lik_ptree_given_ctree(ctree2, parms.curr.coa[["kappa"]], parms.curr.coa[["lambda"]], proptree$prop_hosts)

          if (log(runif(1)) < (pTTree_part2 + pPTree_part2 + proptree$rev_density -
                               pTTree_part - pPTree_part - proptree$prop_density)) {

            ctree <- ctree2
            ttree <- ttree2

            pTTree <- pTTree + pTTree_part2 - pTTree_part
            pPTree <- pPTree + pPTree_part2 - pPTree_part

            if (prop_type == 1) {

              tree_acc_count["add"] <- tree_acc_count["add"] + 1

            } else if (prop_type == 2) {

              tree_acc_count["remove"] <- tree_acc_count["remove"] + 1

            } else {

              tree_acc_count["move"] <- tree_acc_count["move"] + 1

            }

            if (grid[length(grid)] > (min(ttree$ttree[, 1]) - 0.5 * 1)) {

              grid <- seq(obs.end, min(ttree$ttree[, 1]) - 0.5 * 1 - grid.delta, by = - grid.delta)

              fn_list <- num_approx_disc(grid, parms.curr.tr[["r"]], parms.curr.tr[["p"]], parms.curr.tr[["pi"]], w.shape, w.scale,
                                         ws.shape, ws.scale, obs.start, obs.end)

            }

          }

        }

      }

      if (i %% mcmc.thinning == 0) {

        tree_acc_rates["add"] <- tree_acc_count["add"] / tree_prop_count["add"]
        tree_acc_rates["remove"] <- tree_acc_count["remove"] / tree_prop_count["remove"]
        tree_acc_rates["move"] <- tree_acc_count["move"] / tree_prop_count["move"]

        record[[i / mcmc.thinning]]$tree.prop.counts <- tree_prop_count
        record[[i / mcmc.thinning]]$tree.acc.rates <- tree_acc_rates

      }

    }

    if (ss.d.tr > 0) {

      parms.prop.tr <- MASS::mvrnorm(1, mu = parms.curr.tr, Sigma = ss.lam.tr * ss.c.tr * mcmc.cov.tr)

      ss.u.tr <- log(runif(1))

      if (parms.prop.tr["r"] > 0 & parms.prop.tr["p"] > 0 & parms.prop.tr["p"] <= 1 &
          parms.prop.tr["pi"] > 0 & parms.prop.tr["pi"] <= 1) {

        fn_list2 <- num_approx_disc(grid, parms.prop.tr[["r"]], parms.prop.tr[["p"]], parms.prop.tr[["pi"]], w.shape, w.scale,
                                    ws.shape, ws.scale, obs.start, obs.end)

        pTTree2 <- log_lik_ttree(ttree, grid, fn_list2, parms.prop.tr[["r"]], parms.prop.tr[["p"]],
                                 parms.prop.tr[["pi"]], w.shape, w.scale, ws.shape,
                                 ws.scale, obs.start, obs.end, grid.delta, NA)

        ss.alpha.tr <- (pTTree2 - pTTree) +
          (dgamma(parms.prop.tr[["r"]], shape = r.shape, scale = r.scale, log = T) -
             dgamma(parms.curr.tr[["r"]], shape = r.shape, scale = r.scale, log = T)) +
          (dbeta(parms.prop.tr[["p"]], shape1 = p.shape1, shape2 = p.shape2, log = T) -
             dbeta(parms.curr.tr[["p"]], shape1 = p.shape1, shape2 = p.shape2, log = T)) +
          (dbeta(parms.prop.tr[["pi"]], shape1 = pi.shape1, shape2 = pi.shape2, log = T) -
             dbeta(parms.curr.tr[["pi"]], shape1 = pi.shape1, shape2 = pi.shape2, log = T))

      } else {

        ss.alpha.tr <- -Inf

      }

      if (ss.u.tr < ss.alpha.tr) {

        parms.curr.tr[which(update.tr)] <- parms.prop.tr[which(update.tr)]

        fn_list <- fn_list2
        pTTree <- pTTree2

      }

      trace.r[i] <- parms.curr.tr[["r"]]
      trace.p[i] <- parms.curr.tr[["p"]]
      trace.pi[i] <- parms.curr.tr[["pi"]]

      if (i == 1) {

        mcmc.mu.tr <- 0.5 * (parms.init.tr + parms.curr.tr)

        mcmc.cov.tr <- (1 / (ss.v0.tr + ss.d.tr + 3)) * (parms.init.tr %*% t(parms.init.tr) +
                                                  parms.curr.tr %*% t(parms.curr.tr) +
                                                  (ss.v0.tr + ss.d.tr + 1) * mcmc.cov.tr -
                                                  2 * mcmc.mu.tr %*% t(mcmc.mu.tr))

      } else if (ss.f[i] == ss.f[i - 1]) {

        mcmc.mu.new.tr <- ((i - ss.f[i]) / (i - ss.f[i] + 1)) * mcmc.mu.tr +
          (1 / (i - ss.f[i] + 1)) * parms.curr.tr

        mcmc.cov.tr <- (1 / (i - ss.f[i] + ss.v0.tr + ss.d.tr + 2)) *
          ((i - ss.f[i] + ss.v0.tr + ss.d.tr + 1) * mcmc.cov.tr +
             parms.curr.tr %*% t(parms.curr.tr) +
             (i - ss.f[i]) * mcmc.mu.tr %*% t(mcmc.mu.tr) -
             (i - ss.f[i] + 1) * mcmc.mu.new.tr %*% t(mcmc.mu.new.tr))

        mcmc.mu.tr <- mcmc.mu.new.tr

      } else {

        rem.el <- ss.f[i] - 1

        if (rem.el == 0) {

          parms.rem.tr <- parms.init.tr

        } else {

          parms.rem.tr <- c(r = trace.r[rem.el],
                            p = trace.p[rem.el],
                            pi = trace.pi[rem.el])

        }

        mcmc.mu.new.tr <- mcmc.mu.tr + (1 / (i - ss.f[i] + 1)) * (parms.curr.tr - parms.rem.tr)

        mcmc.cov.tr <- mcmc.cov.tr + (1 / (i - ss.f[i] + ss.v0.tr + ss.d.tr + 2)) *
          (parms.curr.tr %*% t(parms.curr.tr) - parms.rem.tr %*% t(parms.rem.tr) +
             (i - ss.f[i] + 1) * (mcmc.mu.tr %*% t(mcmc.mu.tr) - mcmc.mu.new.tr %*% t(mcmc.mu.new.tr)))

        mcmc.mu.tr <- mcmc.mu.new.tr

      }

      ss.lam.tr <- max(c(ss.min, ss.lam.tr * exp((ss.del.tr / (ss.nstart + i)) * (min(c(1, exp(ss.alpha.tr))) - ss.a))))

    }

    if (ss.d.coa > 0) {

      parms.prop.coa <- MASS::mvrnorm(1, mu = parms.curr.coa, Sigma = ss.lam.coa * ss.c.coa * mcmc.cov.coa)

      ss.u.coa <- log(runif(1))

      if (parms.prop.coa["kappa"] >= 0 & parms.prop.coa["lambda"] >= 0) {

        pPTree2 <- log_lik_ptree_given_ctree(ctree, parms.prop.coa[["kappa"]],
                                             parms.prop.coa[["lambda"]], NA)

        ss.alpha.coa <- (pPTree2 - pPTree) +
          (dgamma(parms.prop.coa[["kappa"]], shape = kappa.shape, scale = kappa.scale, log = T) -
             dgamma(parms.curr.coa[["kappa"]], shape = kappa.shape, scale = kappa.scale, log = T)) +
          (dgamma(parms.prop.coa[["lambda"]], shape = lambda.shape, scale = lambda.scale, log = T) -
             dgamma(parms.curr.coa[["lambda"]], shape = lambda.shape, scale = lambda.scale, log = T))

      } else {

        ss.alpha.coa <- -Inf

      }

      if (ss.u.coa < ss.alpha.coa) {

        parms.curr.coa[which(update.coa)] <- parms.prop.coa[which(update.coa)]

        pPTree <- pPTree2

      }

      trace.kappa[i] <- parms.curr.coa[["kappa"]]
      trace.lambda[i] <- parms.curr.coa[["lambda"]]

      if (i == 1) {

        mcmc.mu.coa <- 0.5 * (parms.init.coa + parms.curr.coa)

        mcmc.cov.coa <- (1 / (ss.v0.coa + ss.d.coa + 3)) * (parms.init.coa %*% t(parms.init.coa) +
                                                           parms.curr.coa %*% t(parms.curr.coa) +
                                                           (ss.v0.coa + ss.d.coa + 1) * mcmc.cov.coa -
                                                           2 * mcmc.mu.coa %*% t(mcmc.mu.coa))

      } else if (ss.f[i] == ss.f[i - 1]) {

        mcmc.mu.new.coa <- ((i - ss.f[i]) / (i - ss.f[i] + 1)) * mcmc.mu.coa +
          (1 / (i - ss.f[i] + 1)) * parms.curr.coa

        mcmc.cov.coa <- (1 / (i - ss.f[i] + ss.v0.coa + ss.d.coa + 2)) *
          ((i - ss.f[i] + ss.v0.coa + ss.d.coa + 1) * mcmc.cov.coa +
             parms.curr.coa %*% t(parms.curr.coa) +
             (i - ss.f[i]) * mcmc.mu.coa %*% t(mcmc.mu.coa) -
             (i - ss.f[i] + 1) * mcmc.mu.new.coa %*% t(mcmc.mu.new.coa))

        mcmc.mu.coa <- mcmc.mu.new.coa

      } else {

        rem.el <- ss.f[i] - 1

        if (rem.el == 0) {

          parms.rem.coa <- parms.init.coa

        } else {

          parms.rem.coa <- c(kappa = trace.kappa[rem.el],
                             lambda = trace.lambda[rem.el])

        }

        mcmc.mu.new.coa <- mcmc.mu.coa + (1 / (i - ss.f[i] + 1)) * (parms.curr.coa - parms.rem.coa)

        mcmc.cov.coa <- mcmc.cov.coa + (1 / (i - ss.f[i] + ss.v0.coa + ss.d.coa + 2)) *
          (parms.curr.coa %*% t(parms.curr.coa) - parms.rem.coa %*% t(parms.rem.coa) +
             (i - ss.f[i] + 1) * (mcmc.mu.coa %*% t(mcmc.mu.coa) - mcmc.mu.new.coa %*% t(mcmc.mu.new.coa)))

        mcmc.mu.coa <- mcmc.mu.new.coa

      }

      ss.lam.coa <- max(c(ss.min, ss.lam.coa * exp((ss.del.coa / (ss.nstart + i)) * (min(c(1, exp(ss.alpha.coa))) - ss.a))))

    }

    if (i %% mcmc.thinning == 0) {

      if (verbose == F) {

        utils::setTxtProgressBar(pb, i)

      }

      if (verbose==T) {

        message(sprintf('it = %d, r = %f, off.p = %f, pi = %f, kappa = %f, lambda = %f, prior = %e, likelihood = %e, nind = %d', i, parms.curr.tr["r"], parms.curr.tr["p"], parms.curr.tr["pi"], parms.curr.coa["kappa"], parms.curr.coa["lambda"], pTTree, pPTree, nrow(ttree$ttree)))

      }

      rec_ctree <- trim_root(ctree)

      record[[i / mcmc.thinning]]$ctree <- rec_ctree
      record[[i / mcmc.thinning]]$pTTree <- pTTree
      record[[i / mcmc.thinning]]$pPTree <- pPTree
      record[[i / mcmc.thinning]]$kappa <- parms.curr.coa[["kappa"]]
      record[[i / mcmc.thinning]]$lambda <- parms.curr.coa[["lambda"]]
      record[[i / mcmc.thinning]]$off.r <- parms.curr.tr[["r"]]
      record[[i / mcmc.thinning]]$off.p <- parms.curr.tr[["p"]]
      record[[i / mcmc.thinning]]$pi <- parms.curr.tr[["pi"]]
      record[[i / mcmc.thinning]]$w.shape <- w.shape
      record[[i / mcmc.thinning]]$w.scale <- w.scale
      record[[i / mcmc.thinning]]$ws.shape <- ws.shape
      record[[i / mcmc.thinning]]$ws.scale <- ws.scale
      record[[i / mcmc.thinning]]$mcmc.cov.tr <- mcmc.cov.tr
      record[[i / mcmc.thinning]]$ss.lam.tr <- ss.lam.tr
      record[[i / mcmc.thinning]]$mcmc.cov.coa <- mcmc.cov.coa
      record[[i / mcmc.thinning]]$ss.lam.coa <- ss.lam.coa

      record[[i / mcmc.thinning]]$source <- rec_ctree$ctree[rec_ctree$ctree[which(rec_ctree$ctree[, 4] == 0), 2], 4]
      if (record[[i / mcmc.thinning]]$source <= length(rec_ctree$nam)) {

        record[[i / mcmc.thinning]]$source <- rec_ctree$nam[record[[i / mcmc.thinning]]$source]

      } else {

        record[[i / mcmc.thinning]]$source <- 'Unsampled'

      }

    }

  }#End of main MCMC loop

  class(record)<-'resTransPhyloM'
  return(record)

}
