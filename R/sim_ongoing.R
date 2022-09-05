#' Numeric approximations (discretised) of branching process
#'
#' @param grid Discrete grid over which to evaluate functions.
#' @param delta Discrete time step.
#' @param off.r Shape parameter for the number of offspring.
#' @param off.p Probability parameter for the number of offspring.
#' @param pi Probability of host being sampled.
#' @param w.shape Shape parameter of generation time distribution.
#' @param w.scale Scale parameter of generation time distribution.
#' @param ws.shape Shape parameter of primary sampling time distribution.
#' @param ws.scale Scale parameter of primary sampling time distribution.
#' @param dateS Start time of outbreak sampling.
#' @param dateT Stop time of outbreak sampling.
num_approx_disc <- function(grid,
                            delta,
                            off.r,
                            off.p,
                            pi,
                            w.shape,
                            w.scale,
                            ws.shape,
                            ws.scale,
                            dateS,
                            dateT) {

  grid_size <- length(grid)

  omega <- numeric(grid_size)
  omega_bar <- numeric(grid_size)
  phi <- numeric(grid_size)
  pit <- numeric(grid_size)

  omega[1] <- 1
  omega_bar[1] <- 1
  phi[1] <- 1

  pit <- pi * (pgamma(dateT - grid, shape = ws.shape, scale = ws.scale) - pgamma(dateS - grid, shape = ws.shape, scale = ws.scale))

  gamma_prob <- pgamma((0:(grid_size - 1)) * delta,  shape = w.shape, scale = w.scale) - pgamma((0:(grid_size - 1)) * delta - delta,  shape = w.shape, scale = w.scale)

  ft <- 1 - cumsum(gamma_prob)

  for (g in 2:grid_size) {

    omega_bar[g] <- ft[g] + sum(gamma_prob[g:2] * omega[1:(g - 1)])

    phi[g] <- ((1 - off.p) / (1 - off.p * omega_bar[g])) ^ off.r

    omega[g] <- (1 - pit[g]) * phi[g]

  }

  return(list(omega = omega,
              omega_bar = omega_bar,
              phi = phi,
              pit = pit,
              gamma_prob = gamma_prob))

}


#' Simulate an ongoing outbreak
#'
#' @param off.r Shape parameter for the number of offspring.
#' @param off.p Probability parameter for the number of offspring.
#' @param lm_const Initial pathogen population.
#' @param lm_rate Pathogen population linear growth rate.
#' @param pi Probability of host being sampled.
#' @param add.prob Probability of repeated samples.
#' @param add.sep Time between repeated samples.
#' @param w.shape Shape parameter of generation time distribution.
#' @param w.scale Scale parameter of generation time distribution.
#' @param ws.shape Shape parameter of primary sampling time distribution.
#' @param ws.scale Scale parameter of primary sampling time distribution.
#' @param w.mean Mean parameter of generation time distribution.
#' @param w.std Standard deviation of generation time distribution.
#' @param ws.mean Mean parameter of primary sampling time distribution.
#' @param ws.std Standard deviation of primary sampling time distribution.
#' @param dateStartOutbreak Outbreak start date.
#' @param dateS Start time of outbreak sampling.
#' @param dateT Stop time of outbreak sampling.
#' @param delta Discrete time step.
#'
#' @export
sim_ongoing <- function(off.r = 1,
                        off.p = 0.5,
                        lm_const = 0.2,
                        lm_rate = 0.2,
                        pi = 0.5,
                        add.prob = 0.1,
                        add.sep = 14 / 365,
                        w.shape = 2,
                        w.scale = 1,
                        ws.shape = NULL,
                        ws.scale = NULL,
                        w.mean = NULL,
                        w.std = NULL,
                        ws.mean = NULL,
                        ws.std = NULL,
                        dateStartOutbreak = 2000,
                        dateS = NULL,
                        dateT = 2010,
                        delta = 1 / 365) {

  if (!is.null(w.mean) && !is.null(w.std)) {

    w.shape <- w.mean ^ 2 / w.std ^ 2
    w.scale <- w.std ^ 2 / w.mean

  }

  if (!is.null(ws.mean) && !is.null(ws.std)) {

    ws.shape <- ws.mean ^ 2 / ws.std ^ 2
    ws.scale <- ws.std ^ 2 / ws.mean

  }

  if (is.null(ws.shape)) {

    ws.shape <- w.shape

  }

  if (is.null(ws.scale)) {

    ws.scale <- w.scale

  }

  if (is.null(dateS)) {

    dateS <- dateStartOutbreak

  }

  grid <- seq(dateT, dateStartOutbreak - delta, by = - delta)

  fn_list <- num_approx_disc(grid, delta, off.r, off.p, pi, w.shape, w.scale,
                             ws.shape, ws.scale, dateS, dateT)

  omega <- fn_list$omega
  omega_bar <- fn_list$omega_bar
  phi <- fn_list$phi
  pit <- fn_list$pit
  gamma_prob <- fn_list$gamma_prob

  repeat {

    ttree <- matrix(0,1,3)

    i <- 1
    n <- 1
    t_inf <- dateStartOutbreak

    ttree[i, 1] <- t_inf

    repeat {

      tidx <- 1 + round((dateT - t_inf) / delta)

      if (runif(1) < (pit[tidx] / (1 - omega[tidx]))) {

        ttree[i, 2] <- 1

      } else {

        ttree[i, 2] <- NA

      }

      repeat {

        offspring <- rnbinom(1, size = off.r, prob = off.p)

        if (offspring > 0) {

          offspring_inc <- rbinom(1, size = offspring, prob = 1 - omega_bar[tidx])

        } else {

          offspring_inc <- 0

        }


        if (!is.na(ttree[i, 2]) | offspring_inc > 0) {

          break

        }

      }

      if (offspring_inc > 0) {

        ttree <- rbind(ttree, matrix(0, offspring_inc, 3))

        gent <- (gamma_prob[tidx:2] * (1 - omega[1:(tidx - 1)])) / (1 - omega_bar[tidx])

        gent_sam <- sample(1:(tidx - 1), size = offspring_inc, prob = gent, replace = T)

        ttree[n + 1:offspring_inc, 1] <- grid[gent_sam]
        ttree[n + 1:offspring_inc, 3] <- i

        n <- n + offspring_inc

      }

      i <- i + 1

      if (i > nrow(ttree)) {

        break

      } else {

        t_inf <- ttree[i, 1]

      }

    }

    ord <- c(which(!is.na(ttree[,2])), which(is.na(ttree[,2])))
    invord = 1:length(ord)
    invord[ord] = 1:length(ord)
    ttree<-ttree[ord, , drop=FALSE]
    ttree[ttree[, 3] > 0, 3] = invord[ttree[ttree[, 3] > 0, 3]]

    obs <- matrix(0, 0, 2)

    for (h in which(ttree[, 2] == 1)) {

      t_inf <- ttree[h, 1]

      repeat {

        t_sam <- t_inf + rgamma(1, shape = ws.shape, scale = ws.scale)

        if (t_sam > dateS & t_sam < dateT) {

          break

        }

      }

      obs <- rbind(obs, c(t_sam, h))

      repeat {

        if ((t_sam + add.sep) < dateT & runif(1) < add.prob) {

          t_sam <- t_sam + add.sep

          obs <- rbind(obs, c(t_sam, h))

        } else {

          break

        }

      }

    }

    obs_df <- data.frame(time = obs[, 1], individual = obs[, 2])
    tr_df <- data.frame(time = ttree[, 1], infector = ttree[, 3], infected = 1:length(ttree[, 1]))
    ige_out <- sample_phylo_conditioned(observations = obs_df,
                                        transmissions = tr_df,
                                        growth_model = "linear",
                                        lm_const = lm_const,
                                        lm_rate = lm_rate,
                                        show_plot = F)

    if (length(which(is.na(ige_out$nodes$ancestor))) == 1) {

      break

    }

  }

  sams <- length(obs[, 1])

  ctree <- matrix(0, sams, 4)
  ctree[, 1] <- obs[, 1]
  ctree[, 4] <- obs[, 2]

  if (sams > 1) {

    ctree <- rbind(ctree, matrix(0, sams - 1, 4))

    for (r in (sams + 1):(2 * sams - 1)) {

      ctree[r, ] <- c(ige_out$nodes$time[r],
                      which(ige_out$nodes$ancestor == r),
                      ige_out$nodes$individual[r])

    }

    ord <- c(1:sams, sams + order(ctree[(sams + 1):(2 * sams - 1), 1], decreasing = T))
    invord = 1:length(ord)
    invord[ord] = 1:length(ord)

    ctree[, ] <- ctree[ord, ]
    ctree[which(ctree[, 2] > 0), 2] <- invord[ctree[which(ctree[, 2] > 0), 2]]
    ctree[which(ctree[, 3] > 0), 3] <- invord[ctree[which(ctree[, 3] > 0), 3]]

  }

  i <- 1
  j <- length(ctree[, 1]) + 1
  repeat {

    r <- which(ctree[, 2] == i | ctree[, 3] == i)

    if (length(r) == 0) {

      ti <- ctree[i, 4]

      if (ti > 0) {

        ctree <- rbind(ctree, rep(0, 4))

        ctree[j, ] <- c(ttree[ti, 1], i, 0, ttree[ti, 3])
#        ctree[r, which(ctree[r, 2:3] == i) + 1] <- j

        j <- j + 1

      }

    } else if (ctree[i, 4] != ctree[r, 4]) {

      ti <- ctree[i, 4]

      ctree <- rbind(ctree, rep(0, 4))

      ctree[j, ] <- c(ttree[ti, 1], i, 0, ttree[ti, 3])
      ctree[r, which(ctree[r, 2:3] == i) + 1] <- j

      j <- j + 1

    }

    if (i < length(ctree[, 1])) {

      i <- i + 1

    } else {

      break

    }

  }


  ord <- c(1:sams, sams + order(ctree[(sams + 1):length(ctree[, 1]), 1], decreasing = T))
  invord = 1:length(ord)
  invord[ord] = 1:length(ord)

  ctree[, ] <- ctree[ord, ]
  ctree[which(ctree[, 2] > 0), 2] <- invord[ctree[which(ctree[, 2] > 0), 2]]
  ctree[which(ctree[, 3] > 0), 3] <- invord[ctree[which(ctree[, 3] > 0), 3]]

  nam_host <- ctree[1:sams, 4]
  nam_num <- integer(sams)

  host_count <- integer(length(unique(nam_host)))

  for (i in 1:sams) {

    host_count[nam_host[i]] <- host_count[nam_host[i]] + 1

    nam_num[i] <- host_count[nam_host[i]]

  }

  out <- list(ctree = ctree, nam = paste(nam_host, ".", nam_num, sep = ""))
  class(out) <- 'ctree'

  return(out)

}


