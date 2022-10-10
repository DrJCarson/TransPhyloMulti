#' Calculate transmission tree likelihood
#'
#' @param ttree Transmission tree.
#' @param off.r Shape parameter for the number of offspring.
#' @param off.p Probability parameter for the number of offspring.
#' @param pi Probability of host being sampled.
#' @param w.shape Shape parameter of generation time distribution.
#' @param w.scale Scale parameter of generation time distribution.
#' @param ws.shape Shape parameter of primary sampling time distribution.
#' @param ws.scale Scale parameter of primary sampling time distribution.
#' @param dateStartOutbreak Outbreak start date.
#' @param dateS Start time of outbreak sampling.
#' @param dateT Stop time of outbreak sampling.
#' @param delta Discrete time step.
#'
#' @export
log_lik_ttree <- function(ttree, off.r, off.p, pi, w.shape, w.scale, ws.shape,
                          ws.scale, dateS, dateT, delta = 1 / 365) {

  obs <- ttree$obs
  ttree <- ttree$ttree

  log_lik <- 0

  grid <- seq(dateT, min(ttree[, 1]) - delta, by = - delta)

  fn_list <- num_approx_disc(grid, delta, off.r, off.p, pi, w.shape, w.scale,
                             ws.shape, ws.scale, dateS, dateT)

  omega <- fn_list$omega
  omega_bar <- fn_list$omega_bar
  phi <- fn_list$phi
  pit <- fn_list$pit
  gamma_prob <- fn_list$gamma_prob

  nhosts <- nrow(ttree)
  sum_lim <- qnbinom(0.999999, size = off.r, prob = off.p)

  for (i in 1:nhosts) {

    tidx <- 1 + round((dateT - ttree[i, 1]) / delta)

    if (ttree[i, 2] > 0) {

      log_lik <- log_lik + log(pi) - log(1 - omega[tidx])

      obs_time <- min(obs[which(obs[, 2] == i), 1])

      log_lik <- log_lik + dgamma(obs_time - ttree[i, 1], shape = ws.shape, scale = ws.scale, log = T)

    } else {

      log_lik <- log_lik + log(1 - pit[tidx]) - log(1 - omega[tidx])

    }

    inc_off_idx <- which(ttree[, 3] == i)
    inc_off <- length(inc_off_idx)

    alpha_sum <- sum(dnbinom(inc_off:sum_lim, size = off.r, prob = off.p) *
                       choose(inc_off:sum_lim, inc_off) *
                       omega_bar[tidx] ^ (0:(sum_lim - inc_off)))

    log_lik <- log_lik + log(alpha_sum)

    if (inc_off > 0) {

      for (j in 1:inc_off) {

        inf_host <- inc_off_idx[j]

        tidx2 <- 1 + round((dateT - ttree[inf_host, 1]) / delta)

        log_lik <- log_lik + log(1 - omega[tidx2])

        log_lik <- log_lik + dgamma(ttree[inf_host, 1] - ttree[i, 1], shape = w.shape, scale = w.scale, log = T)

      }

    }

  }

  return(log_lik)

}


#' Likelihood evaluation for the linear growth model
#'
#' @param infected_time Time at which individual was infected.
#' @param final_time Lower bound of time period.
#' @param start_time Upper bound of time period.
#' @param lm_const Initial population vale.
#' @param lm_rate Growth rate.
#' @param branch_combs Number of possible coalescence possibilities.
#' @param coalescence Whether or not a coalescence occurs at final_time.
#'
#' @export
log_likelihood_coalescence_linear <- function(infected_time, final_time,
                                              start_time, lm_const, lm_rate,
                                              branch_combs, coalescence) {

  if (coalescence == 1) {

    log_likelihood_increment <- - log(lm_rate * (final_time -
                                                   infected_time) + lm_const) -
      (branch_combs / lm_rate) * (log(lm_rate *
                                       (start_time - infected_time) +
                                       lm_const) - log(lm_rate * (final_time -
                                                                    infected_time) + lm_const))

  } else {

    if (branch_combs > 0) {

      log_likelihood_increment <- - (branch_combs / lm_rate) * (log(lm_rate *
                                                                     (start_time - infected_time) +
                                                                     lm_const) - log(lm_rate * (final_time -
                                                                                                  infected_time) + lm_const))

    } else {

      log_likelihood_increment <- 0

    }

  }

  return(log_likelihood_increment)

}

#' Likelihood evaluation for a phylogenetic tree conditional on a transmission tree
#'
#' @param ctree Combined tree.
#' @param lm_const Initial pathogen population.
#' @param lm_rate Pathogen growth rate
#'
#' @export
log_lik_ptree_given_ctree <- function(ctree, lm_const, lm_rate) {

  log_lik <- 0

  ctree <- sim$ctree

  for (host in 1:max(ctree[, 4])) {

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

      log_lik <- log_lik + log_likelihood_coalescence_linear(inf_time, t1, t2, lm_const, lm_rate, choose(lineages, 2), is_coa)

    }

  }

  return(log_lik)

}
