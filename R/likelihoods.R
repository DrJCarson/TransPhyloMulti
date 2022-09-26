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
                       omega_bar[t_idx] ^ (0:(sum_lim - inc_off)))

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
