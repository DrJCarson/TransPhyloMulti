#' Calculate transmission tree likelihood
#'
#' @param ttree Transmission tree
#' @param grid Discrete grid over which to evaluate functions
#' @param fn_list Precalculated discrete approximations of exclusion probabilities
#' @param off.r Shape parameter for the number of offspring
#' @param off.p Probability parameter for the number of offspring
#' @param pi Probability of host being sampled
#' @param w.shape Shape parameter of generation time distribution
#' @param w.scale Scale parameter of generation time distribution
#' @param ws.shape Shape parameter of primary sampling time distribution
#' @param ws.scale Scale parameter of primary sampling time distribution
#' @param obs.start Start time of outbreak sampling
#' @param obs.end Stop time of outbreak sampling
#' @param grid.delta Discrete time step
#' @param hosts Hosts over which likelihood is calculated
log_lik_ttree <- function(ttree, grid, fn_list, off.r, off.p, pi, w.shape, w.scale, ws.shape,
                          ws.scale, obs.start, obs.end, grid.delta = 1 / 365, hosts = NA) {

  obs <- ttree$obs
  ttree <- ttree$ttree

  log_lik <- 0

  omega <- fn_list$omega
  omega_bar <- fn_list$omega_bar
  pit <- fn_list$pit
  gamma_prob <- fn_list$gamma_prob

  if (is.na(sum(hosts))) {

    hosts <- 1:nrow(ttree)

  }

  sum_lim <- qnbinom(1 - 1e-10, size = off.r, prob = off.p)

  for (i in hosts) {

    tidx <- 1 + floor((obs.end - ttree[i, 1]) / grid.delta)

    omega_int <- omega[tidx] + ((ttree[i, 1] - grid[tidx]) / (grid[tidx + 1] - grid[tidx])) * (omega[tidx + 1] - omega[tidx])
    omega_bar_int <- omega_bar[tidx] + ((ttree[i, 1] - grid[tidx]) / (grid[tidx + 1] - grid[tidx])) * (omega_bar[tidx + 1] - omega_bar[tidx])
    pit_int <- pit[tidx] + ((ttree[i, 1] - grid[tidx]) / (grid[tidx + 1] - grid[tidx])) * (pit[tidx + 1] - pit[tidx])

    if (ttree[i, 2] > 0) {

      log_lik <- log_lik + log(pi)

      obs_time <- min(obs[which(obs[, 2] == i), 1])

      log_lik <- log_lik + dgamma(obs_time - ttree[i, 1], shape = ws.shape, scale = ws.scale, log = T)

    } else {

      log_lik <- log_lik + log(1 - pit_int)

    }

    if (ttree[i, 3] == 0) {

      log_lik <- log_lik - log(1 - omega_int)

    }

    inc_off_idx <- which(ttree[, 3] == i)
    inc_off <- length(inc_off_idx)

    alpha_sum <- sum(dnbinom(inc_off:sum_lim, size = off.r, prob = off.p) *
                       choose(inc_off:sum_lim, inc_off) *
                       omega_bar_int ^ (0:(sum_lim - inc_off)))

    log_lik <- log_lik + log(alpha_sum)

    log_lik <- log_lik + lfactorial(inc_off)

    if (inc_off > 0) {

      for (j in 1:inc_off) {

        inf_host <- inc_off_idx[j]

        log_lik <- log_lik + dgamma(ttree[inf_host, 1] - ttree[i, 1], shape = w.shape, scale = w.scale, log = T)

      }

    }

  }

  return(log_lik)

}


#' Likelihood evaluation for the linear growth model
#'
#' @param infected_time Time at which individual was infected
#' @param final_time Lower bound of time period
#' @param start_time Upper bound of time period
#' @param kappa Initial population vale
#' @param lambda Growth rate
#' @param branch_combs Number of possible coalescence possibilities
#' @param coalescence Whether or not a coalescence occurs at final_time
log_likelihood_coalescence_linear <- function(infected_time, final_time,
                                              start_time, kappa, lambda,
                                              branch_combs, coalescence) {

  if (coalescence == 1) {

    if (lambda == 0) {

      log_likelihood_increment <- (-log(kappa) * kappa + branch_combs * final_time - branch_combs * start_time) / kappa

    } else {

      log_likelihood_increment <- (- log(lambda * (final_time - infected_time) +
                                          kappa) - (branch_combs / lambda) * (log(
                                            lambda * (start_time - infected_time) +
                                              kappa) - log(lambda * (final_time -
                                                infected_time) + kappa)))

    }

  } else {

    if (branch_combs > 0) {

      if (lambda == 0) {

        log_likelihood_increment <- (branch_combs * final_time - branch_combs * start_time) / kappa

      } else {

        log_likelihood_increment <- (- (branch_combs / lambda) * (log(lambda *
                                (start_time - infected_time) + kappa) - log(
                                lambda * (final_time - infected_time) + kappa)))

      }

    } else {

      log_likelihood_increment <- 0

    }

  }

  return(log_likelihood_increment)

}

#' Likelihood evaluation for a phylogenetic tree conditional on a transmission tree
#'
#' @param ctree Combined tree
#' @param kappa Initial pathogen population
#' @param lambda Pathogen growth rate
#' @param hosts Hosts over which likelihood is calculated
log_lik_ptree_given_ctree <- function(ctree, kappa, lambda, hosts = NA) {

  log_lik <- 0

  ctree <- ctree$ctree

  if (is.na(sum(hosts))) {

    hosts <- 1:max(ctree[, 4])

  }

  for (host in hosts) {

    host_rows <- which(ctree[, 4] == host)
    host_rows <- host_rows[order(ctree[host_rows, 1], decreasing = T)]

    inf_time <- ctree[which(ctree[, 2] == max(host_rows))[1], 1]

    lineages <- 0

    for (i in 1:length(host_rows)) {

      if (ctree[host_rows[i], 3] == 0) { # leaf

        lineages <- lineages + 1

      } else {

        lineages <- lineages - 1

      }

      t1 <- ctree[host_rows[i], 1]

      if (i < length(host_rows)) {

        t2 <- ctree[host_rows[i + 1], 1]

        if (ctree[host_rows[i + 1], 3] == 0) {

          is_coa <- 0

        } else {

          is_coa <- 1

        }

      } else{

        t2 <- inf_time

        is_coa <- 0

      }

      log_lik <- log_lik + log_likelihood_coalescence_linear(inf_time, t2, t1, kappa, lambda, choose(lineages, 2), is_coa)

    }

  }

  return(log_lik)

}
