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

  if (sum(is.na(init.ctree))) {

    ctree <- init_ctree(ptree)

  } else {

    ctree <- init.ctree

  }

  if (is.na(grid.delta)) {

    grid.delta <- 0.001 * (max(ptree$ptree[, 1]) - min(ptree$ptree[, 1]))

  }

  const.pi <- 3.14159265359

  parms.init <- c(r = init.r,
                  p = init.p,
                  pi = init.pi,
                  kappa = init.kappa,
                  lambda = init.lambda)

  parms.curr <- parms.init

  ss.a <- 0.234
  update.ind <- c(update.r, update.p, update.pi, update.kappa, update.lambda)
  ss.d <- sum(as.numeric(update.ind))
  ss.v0 <- ss.d
  ss.f <- floor(0.5 * sqrt(2 * (1:mcmc.length)))
  ss.min <- 0.1

  ss.c <- 2.38 ^ 2 / ss.d
  ss.lamstart <- 1
  ss.lam <- ss.lamstart
  ss.nstart <- 5 / (ss.a * (1 - ss.a))
  ss.A <- -qnorm(ss.a / 2)
  ss.del <- (1 - (1 / ss.d)) * ((sqrt(2 * const.pi) * exp(ss.A ^ 2 / 2)) / (2 * ss.A)) + (1 / (ss.d * ss.a * (1 - ss.a)))

  ttree <- extractTTreeM(ctree)

  record <- vector('list', mcmc.length / mcmc.thinning)

  trace.r <- numeric(mcmc.length)
  trace.p <- numeric(mcmc.length)
  trace.pi <- numeric(mcmc.length)
  trace.kappa <- numeric(mcmc.length)
  trace.lambda <- numeric(mcmc.length)

  grid <- seq(obs.end, min(ttree$ttree[, 1]) - 0.5 * 1 - grid.delta, by = - grid.delta)

  fn_list <- num_approx_disc(grid, grid.delta, init.r, init.p, init.pi, w.shape, w.scale,
                             ws.shape, ws.scale, obs.start, obs.end)

  pTTree <- log_lik_ttree(ttree, grid, fn_list, init.r, init.p, init.pi, w.shape, w.scale,
                          ws.shape, ws.scale, obs.start, obs.end, grid.delta, NA)

  pPTree <- log_lik_ptree_given_ctree(ctree, init.kappa, init.lambda, NA)

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

          pTTree_part <- log_lik_ttree(ttree, grid, fn_list, parms.curr["r"], parms.curr["p"], parms.curr["pi"], w.shape, w.scale, ws.shape,
                                       ws.scale, obs.start, obs.end, grid.delta, proptree$curr_hosts)

          pTTree_part2 <- log_lik_ttree(ttree2, grid, fn_list, parms.curr["r"], parms.curr["p"], parms.curr["pi"], w.shape, w.scale, ws.shape,
                                        ws.scale, obs.start, obs.end, grid.delta, proptree$prop_hosts)

          pPTree_part <- log_lik_ptree_given_ctree(ctree, parms.curr["kappa"], parms.curr["lambda"], proptree$curr_hosts)

          pPTree_part2 <- log_lik_ptree_given_ctree(ctree2, parms.curr["kappa"], parms.curr["lambda"], proptree$prop_hosts)

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

              fn_list <- num_approx_disc(grid, grid.delta, parms.curr["r"], parms.curr["p"], parms.curr["pi"], w.shape, w.scale,
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

    if (ss.d > 0) {

      parms.prop <- MASS::mvrnorm(1, mu = parms.curr, Sigma = ss.lam * ss.c * mcmc.cov)

      ss.u <- log(runif(1))

      if (parms.prop["r"] > 0 & parms.prop["p"] > 0 & parms.prop["p"] <= 1 &
          parms.prop["pi"] > 0 & parms.prop["pi"] <= 1 & parms.prop["kappa"] > 0 &
          parms.prop["lambda"] > 0) {

        fn_list2 <- num_approx_disc(grid, grid.delta, parms.prop["r"], parms.prop["p"], parms.prop["pi"], w.shape, w.scale,
                                   ws.shape, ws.scale, obs.start, obs.end)

        pTTree2 <- log_lik_ttree(ttree, grid, fn_list2, parms.prop["r"], parms.prop["p"],
                                 parms.prop["pi"], w.shape, w.scale, ws.shape,
                                 ws.scale, obs.start, obs.end, grid.delta, NA)

        pPTree2 <- log_lik_ptree_given_ctree(ctree, parms.prop["kappa"],
                                             parms.prop["lambda"], NA)

        ss.alpha <- (pTTree2 - pTTree) + (pPTree2 - pPTree) +
          (dgamma(parms.prop["r"], shape = r.shape, scale = r.scale, log = T) -
             dgamma(parms.curr["r"], shape = r.shape, scale = r.scale, log = T)) +
          (dbeta(parms.prop["p"], shape1 = p.shape1, shape2 = p.shape2, log = T) -
             dbeta(parms.curr["p"], shape1 = p.shape1, shape2 = p.shape2, log = T)) +
          (dbeta(parms.prop["pi"], shape1 = pi.shape1, shape2 = pi.shape2, log = T) -
             dbeta(parms.curr["pi"], shape1 = pi.shape1, shape2 = pi.shape2, log = T)) +
          (dgamma(parms.prop["kappa"], shape = kappa.shape, scale = kappa.scale, log = T) -
             dgamma(parms.curr["kappa"], shape = kappa.shape, scale = kappa.scale, log = T)) +
          (dgamma(parms.prop["lambda"], shape = lambda.shape, scale = lambda.scale, log = T) -
             dgamma(parms.curr["lambda"], shape = lambda.shape, scale = lambda.scale, log = T))

      } else {

        ss.alpha <- -Inf

      }

      if (ss.u < ss.alpha) {

        parms.curr[which(update.ind)] <- parms.prop[which(update.ind)]

        fn_list <- fn_list2
        pTTree <- pTTree2
        pPTree <- pPTree2

      }

      trace.r[i] <- parms.curr["r"]
      trace.p[i] <- parms.curr["p"]
      trace.pi[i] <- parms.curr["pi"]
      trace.kappa[i] <- parms.curr["kappa"]
      trace.lambda[i] <- parms.curr["lambda"]

      if (i == 1) {

        mcmc.mu <- 0.5 * (parms.init + parms.curr)

        mcmc.cov <- (1 / (ss.v0 + ss.d + 3)) * (parms.init %*% t(parms.init) +
                                                  parms.curr %*% t(parms.curr) +
                                                  (ss.v0 + ss.d + 1) * mcmc.cov -
                                                  2 * mcmc.mu %*% t(mcmc.mu))

      } else if (ss.f[i] == ss.f[i - 1]) {

        mcmc.mu.new <- ((i - ss.f[i]) / (i - ss.f[i] + 1)) * mcmc.mu +
          (1 / (i - ss.f[i] + 1)) * parms.curr

        mcmc.cov <- (1 / (i - ss.f[i] + ss.v0 + ss.d + 2)) *
          ((i - ss.f[i] + ss.v0 + ss.d + 1) * mcmc.cov +
             parms.curr %*% t(parms.curr) +
             (i - ss.f[i]) * mcmc.mu %*% t(mcmc.mu) -
             (i - ss.f[i] + 1) * mcmc.mu.new %*% t(mcmc.mu.new))

        mcmc.mu <- mcmc.mu.new

      } else {

        rem.el <- ss.f[i] - 1

        if (rem.el == 0) {

          parms.rem <- parms.init

        } else {

          parms.rem <- c(r = trace.r[rem.el],
                         p = trace.p[rem.el],
                         pi = trace.pi[rem.el],
                         kappa = trace.kappa[rem.el],
                         lambda = trace.lambda[rem.el])

        }

        mcmc.mu.new <- mcmc.mu + (1 / (i - ss.f[i] + 1)) * (parms.curr - parms.rem)

        mcmc.cov <- mcmc.cov + (1 / (i - ss.f[i] + ss.v0 + ss.d + 2)) *
          (parms.curr %*% t(parms.curr) - parms.rem %*% t(parms.rem) +
             (i - ss.f[i] + 1) * (mcmc.mu %*% t(mcmc.mu) - mcmc.mu.new %*% t(mcmc.mu.new)))

        mcmc.mu <- mcmc.mu.new

      }

      ss.lam <- max(c(ss.min, ss.lam * exp((ss.del / (ss.nstart + i)) * (min(c(1, exp(ss.alpha))) - ss.a))))

#      if (abs(log(ss.lam) - log(ss.lamstart)) > log(3)) {

#        ss.lam <- ss.lamstart
#        ss.nstart <- (5 / (ss.a * (1 - ss.a))) - i

#      }

    }

    if (i %% mcmc.thinning == 0) {

      if (verbose == F) {

        utils::setTxtProgressBar(pb, i)

      }

      if (verbose==T) {

        message(sprintf('it = %d, r = %f, off.p = %f, pi = %f, kappa = %f, lambda = %f, prior = %e, likelihood = %e, nind = %d', i, parms.curr["r"], parms.curr["p"], parms.curr["pi"], parms.curr["kappa"], parms.curr["lambda"], pTTree, pPTree, nrow(ttree$ttree)))

      }

      rec_ctree <- trim_root(ctree)

      record[[i / mcmc.thinning]]$ctree <- rec_ctree
      record[[i / mcmc.thinning]]$pTTree <- pTTree
      record[[i / mcmc.thinning]]$pPTree <- pPTree
      record[[i / mcmc.thinning]]$lm_const <- parms.curr["kappa"]
      record[[i / mcmc.thinning]]$lm_rate <- parms.curr["lambda"]
      record[[i / mcmc.thinning]]$off.r <- parms.curr["r"]
      record[[i / mcmc.thinning]]$off.p <- parms.curr["p"]
      record[[i / mcmc.thinning]]$pi <- parms.curr["pi"]
      record[[i / mcmc.thinning]]$w.shape <- w.shape
      record[[i / mcmc.thinning]]$w.scale <- w.scale
      record[[i / mcmc.thinning]]$ws.shape <- ws.shape
      record[[i / mcmc.thinning]]$ws.scale <- ws.scale
      record[[i / mcmc.thinning]]$mcmc.cov <- mcmc.cov
      record[[i / mcmc.thinning]]$ss.lam <- ss.lam

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
