#' Simulate an ongoing outbreak with specified observation end date
#'
#' @param off.r Shape parameter for the number of offspring
#' @param off.p Probability parameter for the number of offspring
#' @param kappa Initial pathogen population
#' @param lambda Pathogen population linear growth rate
#' @param pi Probability of host being sampled
#' @param sec.p Probability of repeated samples
#' @param sec.t Time between repeated samples
#' @param w.shape Shape parameter of generation time distribution
#' @param w.scale Scale parameter of generation time distribution
#' @param ws.shape Shape parameter of primary sampling time distribution
#' @param ws.scale Scale parameter of primary sampling time distribution
#' @param w.mean Mean parameter of generation time distribution
#' @param w.std Standard deviation of generation time distribution
#' @param ws.mean Mean parameter of primary sampling time distribution
#' @param ws.std Standard deviation of primary sampling time distribution
#' @param outbreak.start Outbreak start date
#' @param obs.start Start time of outbreak sampling
#' @param obs.end Stop time of outbreak sampling
#' @param grid.delta Discrete time step
#'
#' @export
sim_ongoing_lim_t <- function(off.r = 2,
                              off.p = 0.5,
                              kappa = 0.2,
                              lambda = 0.2,
                              pi = 0.8,
                              sec.p = 0.1,
                              sec.t = 14 / 365,
                              w.shape = 2,
                              w.scale = 1,
                              ws.shape = NULL,
                              ws.scale = NULL,
                              w.mean = NULL,
                              w.std = NULL,
                              ws.mean = NULL,
                              ws.std = NULL,
                              outbreak.start = 2000,
                              obs.start = NULL,
                              obs.end = 2010,
                              grid.delta = 1 / 365) {

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

  if (is.null(obs.start)) {

    obs.start <- outbreak.start

  }

  grid <- seq(obs.end, outbreak.start - grid.delta, by = - grid.delta)

  fn_list <- num_approx_disc(grid, off.r, off.p, pi, w.shape, w.scale,
                             ws.shape, ws.scale, obs.start, obs.end)

  omega <- fn_list$omega
  omega_bar <- fn_list$omega_bar
  phi <- fn_list$phi
  pit <- fn_list$pit
  gamma_prob <- fn_list$gamma_prob

  repeat {

    ttree <- matrix(0, 1, 3)

    i <- 1
    n <- 1
    t_inf <- outbreak.start

    ttree[i, 1] <- t_inf

    repeat {

      tidx <- 1 + round((obs.end - t_inf) / grid.delta)

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
    invord <- 1:length(ord)
    invord[ord] <- 1:length(ord)
    ttree <- ttree[ord, , drop = FALSE]
    ttree[ttree[, 3] > 0, 3] <- invord[ttree[ttree[, 3] > 0, 3]]

    obs <- matrix(0, 0, 2)

    for (h in which(ttree[, 2] == 1)) {

      t_inf <- ttree[h, 1]

      repeat {

        t_sam <- t_inf + rgamma(1, shape = ws.shape, scale = ws.scale)

        if (t_sam > obs.start & t_sam < obs.end) {

          break

        }

      }

      obs <- rbind(obs, c(t_sam, h))

      repeat {

        if ((t_sam + sec.t) < obs.end & runif(1) < sec.p) {

          t_sam <- t_sam + sec.t

          obs <- rbind(obs, c(t_sam, h))

        } else {

          break

        }

      }

    }

    phy_ord <- which(!(1:length(ttree[, 1]) %in% ttree[, 3]))

    repeat {

      infectors_inc <- ttree[phy_ord, 3]
      infectors_exc <- ttree[-phy_ord, 3]

      phy_ord_new <- unique(infectors_inc[which(!(infectors_inc %in% infectors_exc))])
      phy_ord_new <- phy_ord_new[which(!(phy_ord_new %in% phy_ord))]

      if (sum(phy_ord_new) == 0) {

        break

      }

      phy_ord <- c(phy_ord, phy_ord_new)

    }


    sams <- length(obs[, 1])

    ctree <- matrix(0, sams, 4)
    ctree[, 1] <- obs[, 1]
    ctree[, 4] <- obs[, 2]

    for (i in 1:length(phy_ord)) {

      host <- phy_ord[i]

      leaves <- which(ctree[, 4] == host)

      if (length(leaves) == 1) {

        ctree <- rbind(ctree, rep(0, 4))

        ctree[length(ctree[, 1]), ] <- c(ttree[host, 1], leaves, 0, ttree[host, 3])

      } else {

        inf_time <- ttree[host, 1]

        leaf_times <- ctree[leaves, 1]
        leaf_order <- order(leaf_times, decreasing = T)

        l <- 1

        lin <- c(leaves[leaf_order[l]])
        t <- leaf_times[leaf_order[l]]

        repeat {

          if (length(lin) == 1) {

            l <- l + 1

            if (l > length(leaves)) {

              ctree <- rbind(ctree, rep(0, 4))

              ctree[length(ctree[, 1]), ] <- c(inf_time, lin, 0, ttree[host, 3])

              break

            } else {

              lin <- c(lin, leaves[leaf_order[l]])
              t <- leaf_times[leaf_order[l]]

            }

          } else {

            if (l == length(leaves)) {

              lim_time <- inf_time

            } else {

              lim_time <- leaf_times[leaf_order[l + 1]]

            }

            if (lambda == 0) {

              coa_time <- t + log(1 - runif(1)) * kappa / choose(length(lin), 2)

            } else {

              coa_time <- inf_time + (1 / lambda) * ((1 - runif(1)) ^ (lambda / choose(length(lin), 2)) * (lambda * (t - inf_time) + kappa) - kappa)

            }

            if (coa_time > lim_time) {

              sam <- sample(1:choose(length(lin), 2), size = 1)

              children_idx <- combn(length(lin), 2)[, sam]
              children <- lin[children_idx]

              ctree <- rbind(ctree, rep(0, 4))

              ctree[length(ctree[, 1]), ] <- c(coa_time, children, host)

              lin <- c(lin[which(!(lin %in% children))], length(ctree[, 1]))
              t <- coa_time

            } else {

              l <- l + 1

              if (l > length(leaves)) {

                for (j in 1:length(lin)) {

                  ctree <- rbind(ctree, rep(0, 4))

                  ctree[length(ctree[, 1]), ] <- c(inf_time, lin[j], 0, ttree[host, 3])

                }

                break

              } else {

                lin <- c(lin, leaves[leaf_order[l]])
                t <- leaf_times[leaf_order[l]]

              }

            }

          }

        }

      }

    }

    if (length(which(ctree[, 4] == 0)) == 1) {

      break

    }

  }

  ord <- c(1:sams, sams + order(ctree[(sams + 1):length(ctree[, 1]), 1], decreasing = T))
  invord = 1:length(ord)
  invord[ord] = 1:length(ord)

  ctree[, ] <- ctree[ord, ]
  ctree[which(ctree[, 2] > 0), 2] <- invord[ctree[which(ctree[, 2] > 0), 2]]
  ctree[which(ctree[, 3] > 0), 3] <- invord[ctree[which(ctree[, 3] > 0), 3]]

  ctree <- order_hosts(ctree)

  nam_host <- ctree[1:sams, 4]
  nam_num <- integer(sams)

  host_count <- integer(length(unique(nam_host)))

  for (i in 1:sams) {

    host_count[nam_host[i]] <- host_count[nam_host[i]] + 1

    nam_num[i] <- host_count[nam_host[i]]

  }

  out <- list(ctree = ctree, nam = paste(nam_host, ".", nam_num, sep = ""))
  class(out) <- 'ctree'

  out <- trim_root(out)

  return(out)

}


#' Simulate an ongoing outbreak with specified number of observations
#'
#' @param off.r Shape parameter for the number of offspring
#' @param off.p Probability parameter for the number of offspring
#' @param kappa Initial pathogen population
#' @param lambda Pathogen population linear growth rate
#' @param pi Probability of host being sampled
#' @param sec.p Probability of repeated samples
#' @param sec.t Time between repeated samples
#' @param w.shape Shape parameter of generation time distribution
#' @param w.scale Scale parameter of generation time distribution
#' @param ws.shape Shape parameter of primary sampling time distribution
#' @param ws.scale Scale parameter of primary sampling time distribution
#' @param w.mean Mean parameter of generation time distribution
#' @param w.std Standard deviation of generation time distribution
#' @param ws.mean Mean parameter of primary sampling time distribution
#' @param ws.std Standard deviation of primary sampling time distribution
#' @param outbreak.start Outbreak start date
#' @param obs.start Start time of outbreak sampling
#' @param grid.delta Discrete time step
#' @param obs.lim Maximum number of samples
#'
#' @export
sim_ongoing_lim_obs <- function(off.r = 2,
                                off.p = 0.5,
                                kappa = 0.2,
                                lambda = 0.2,
                                pi = 0.8,
                                sec.p = 0.1,
                                sec.t = 14 / 365,
                                w.shape = 2,
                                w.scale = 1,
                                ws.shape = NULL,
                                ws.scale = NULL,
                                w.mean = NULL,
                                w.std = NULL,
                                ws.mean = NULL,
                                ws.std = NULL,
                                outbreak.start = 2000,
                                obs.start = NULL,
                                grid.delta = 1 / 365,
                                obs.lim = 100) {


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

  if (is.null(obs.start)) {

    obs.start <- outbreak.start

  }

  attempts <- 0

  repeat {

    attempts <- attempts + 1

    t_lim1 <- Inf
    t_lim2 <- Inf

    ttree <- matrix(0, 1, 3)
    obs <- matrix(0, 0, 2)

    i <- 1
    n <- 1
    t_inf <- outbreak.start

    ttree[i, 1] <- t_inf

    repeat {

      if (runif(1) < pi) {

        t_sam <- t_inf + rgamma(1, shape = ws.shape, scale = ws.scale)

        if (t_sam > obs.start & t_sam < t_lim2) {

          ttree[i, 2] <- 1

          obs <- rbind(obs, c(t_sam, i))

          repeat {

            if ((t_sam + sec.t) < t_lim2 & runif(1) < sec.p) {

              t_sam <- t_sam + sec.t

              obs <- rbind(obs, c(t_sam, i))

            } else {

              break

            }

          }

        } else {

          ttree[i, 2] <- NA

        }

      } else {

        ttree[i, 2] <- NA

      }

      if (dim(obs)[1] > obs.lim) {

        t_lim1 <- sort(obs[, 1])[obs.lim]
        t_lim2 <- sort(obs[, 1])[obs.lim + 1]

      }

      offspring <- rnbinom(1, size = off.r, prob = off.p)

      if (offspring > 0) {

        off_times <- t_inf + rgamma(offspring, shape = w.shape, scale = w.scale)

        off_times <- off_times[which(off_times < t_lim1)]

        offspring_inc <- length(off_times)

        if (offspring_inc > 0) {

          ttree <- rbind(ttree, matrix(0, offspring_inc, 3))

          ttree[n + 1:offspring_inc, 1] <- off_times

          ttree[n + 1:offspring_inc, 3] <- i

          n <- n + offspring_inc

        }

      }

      i <- i + 1

      if (i > nrow(ttree)) {

        break

      } else {

        t_inf <- ttree[i, 1]

      }

    }

    if (length(obs[, 1]) < obs.lim) {

      next

    }

    drop_obs <- which(obs[, 1] > t_lim1)
    keep_obs <- which(obs[, 1] <= t_lim1)

    ttree[!(1:dim(ttree)[1] %in% obs[keep_obs, 2]), 2] <- NA

    obs <- obs[keep_obs, , drop = F]

    repeat {

      torem <- which(is.na(ttree[, 2]) & !(1:dim(ttree)[1] %in% ttree[, 3]))

      if (length(torem) == 0) {

        break

      }

      for (j in 1:length(torem)) {

        h <- torem[j]

        ttree <- ttree[-h, , drop = FALSE]

        ttree[which(ttree[, 3] > h), 3] <- ttree[which(ttree[, 3] > h), 3] - 1

        obs[which(obs[, 2] > h), 2] <- obs[which(obs[, 2] > h), 2] - 1

        torem[which(torem > h)] <- torem[which(torem > h)] - 1

      }

    }

    ord <- c(which(!is.na(ttree[,2])), which(is.na(ttree[,2])))
    invord <- 1:length(ord)
    invord[ord] <- 1:length(ord)
    ttree <- ttree[ord, , drop = FALSE]
    ttree[ttree[, 3] > 0, 3] <- invord[ttree[ttree[, 3] > 0, 3]]
    obs[, 2] <- invord[obs[, 2]]

    if (dim(ttree)[1] == 0) {

      next

    }

    phy_ord <- which(!(1:length(ttree[, 1]) %in% ttree[, 3]))

    repeat {

      infectors_inc <- ttree[phy_ord, 3]
      infectors_exc <- ttree[-phy_ord, 3]

      phy_ord_new <- unique(infectors_inc[which(!(infectors_inc %in% infectors_exc))])
      phy_ord_new <- phy_ord_new[which(!(phy_ord_new %in% phy_ord))]

      if (sum(phy_ord_new) == 0) {

        break

      }

      phy_ord <- c(phy_ord, phy_ord_new)

    }

    sams <- length(obs[, 1])

    ctree <- matrix(0, sams, 4)
    ctree[, 1] <- obs[, 1]
    ctree[, 4] <- obs[, 2]

    for (i in 1:length(phy_ord)) {

      host <- phy_ord[i]

      leaves <- which(ctree[, 4] == host)

      if (length(leaves) == 1) {

        ctree <- rbind(ctree, rep(0, 4))

        ctree[length(ctree[, 1]), ] <- c(ttree[host, 1], leaves, 0, ttree[host, 3])

      } else {

        inf_time <- ttree[host, 1]

        leaf_times <- ctree[leaves, 1]
        leaf_order <- order(leaf_times, decreasing = T)

        l <- 1

        lin <- c(leaves[leaf_order[l]])
        t <- leaf_times[leaf_order[l]]

        repeat {

          if (length(lin) == 1) {

            l <- l + 1

            if (l > length(leaves)) {

              ctree <- rbind(ctree, rep(0, 4))

              ctree[length(ctree[, 1]), ] <- c(inf_time, lin, 0, ttree[host, 3])

              break

            } else {

              lin <- c(lin, leaves[leaf_order[l]])
              t <- leaf_times[leaf_order[l]]

            }

          } else {

            if (l == length(leaves)) {

              lim_time <- inf_time

            } else {

              lim_time <- leaf_times[leaf_order[l + 1]]

            }

            if (lambda == 0) {

              coa_time <- t + log(1 - runif(1)) * kappa / choose(length(lin), 2)

            } else {

              coa_time <- inf_time + (1 / lambda) * ((1 - runif(1)) ^ (lambda / choose(length(lin), 2)) * (lambda * (t - inf_time) + kappa) - kappa)

            }

            if (coa_time > lim_time) {

              sam <- sample(1:choose(length(lin), 2), size = 1)

              children_idx <- combn(length(lin), 2)[, sam]
              children <- lin[children_idx]

              ctree <- rbind(ctree, rep(0, 4))

              ctree[length(ctree[, 1]), ] <- c(coa_time, children, host)

              lin <- c(lin[which(!(lin %in% children))], length(ctree[, 1]))
              t <- coa_time

            } else {

              l <- l + 1

              if (l > length(leaves)) {

                for (j in 1:length(lin)) {

                  ctree <- rbind(ctree, rep(0, 4))

                  ctree[length(ctree[, 1]), ] <- c(inf_time, lin[j], 0, ttree[host, 3])

                }

                break

              } else {

                lin <- c(lin, leaves[leaf_order[l]])
                t <- leaf_times[leaf_order[l]]

              }

            }

          }

        }

      }

    }

    if (length(which(ctree[, 4] == 0)) == 1) {

      break

    }

  }

  ord <- c(1:sams, sams + order(ctree[(sams + 1):length(ctree[, 1]), 1], decreasing = T))
  invord = 1:length(ord)
  invord[ord] = 1:length(ord)

  ctree[, ] <- ctree[ord, ]
  ctree[which(ctree[, 2] > 0), 2] <- invord[ctree[which(ctree[, 2] > 0), 2]]
  ctree[which(ctree[, 3] > 0), 3] <- invord[ctree[which(ctree[, 3] > 0), 3]]

  ctree <- order_hosts(ctree)

  nam_host <- ctree[1:sams, 4]
  nam_num <- integer(sams)

  host_count <- integer(length(unique(nam_host)))

  for (i in 1:sams) {

    host_count[nam_host[i]] <- host_count[nam_host[i]] + 1

    nam_num[i] <- host_count[nam_host[i]]

  }

  l <- list(ctree = ctree, nam = paste(nam_host, ".", nam_num, sep = ""))
  class(l) <- 'ctree'

  l <- trim_root(l)

  return(list(ctree = l, obs.end = 0.5 * (t_lim1 + t_lim2), attempts = attempts))

}



#' Simulate an ongoing outbreak
#'
#' @param off.r Shape parameter for the number of offspring
#' @param off.p Probability parameter for the number of offspring
#' @param kappa Initial pathogen population
#' @param lambda Pathogen population linear growth rate
#' @param pi Probability of host being sampled
#' @param sec.n Number of secondary observations per host
#' @param sec.t Time between repeated observations
#' @param w.shape Shape parameter of generation time distribution
#' @param w.scale Scale parameter of generation time distribution
#' @param ws.shape Shape parameter of primary sampling time distribution
#' @param ws.scale Scale parameter of primary sampling time distribution
#' @param w.mean Mean parameter of generation time distribution
#' @param w.std Standard deviation of generation time distribution
#' @param ws.mean Mean parameter of primary sampling time distribution
#' @param ws.std Standard deviation of primary sampling time distribution
#' @param outbreak.start Outbreak start date
#' @param obs.start Start time of outbreak sampling
#' @param grid.delta Discrete time step
#' @param host.lim Number of sampled hosts
#'
#' @export
sim_ongoing_lim_hosts <- function(off.r = 2,
                                  off.p = 0.5,
                                  kappa = 0.2,
                                  lambda = 0.2,
                                  pi = 0.8,
                                  sec.n = 1,
                                  sec.t = 14 / 365,
                                  w.shape = 2,
                                  w.scale = 1,
                                  ws.shape = NULL,
                                  ws.scale = NULL,
                                  w.mean = NULL,
                                  w.std = NULL,
                                  ws.mean = NULL,
                                  ws.std = NULL,
                                  outbreak.start = 2000,
                                  obs.start = NULL,
                                  grid.delta = 1 / 365,
                                  host.lim = 200) {


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

  if (is.null(obs.start)) {

    obs.start <- outbreak.start

  }

  attempts <- 0

  repeat {

    attempts <- attempts + 1

    t_lim1 <- Inf
    t_lim2 <- Inf

    ttree <- matrix(0,1,3)
    obs <- matrix(0, 0, 2)

    prim_host <- c()
    prim_obs <- c()


    i <- 1
    n <- 1
    t_inf <- outbreak.start

    ttree[i, 1] <- t_inf

    repeat {

      if (runif(1) < pi) {

        t_sam <- t_inf + rgamma(1, shape = ws.shape, scale = ws.scale)

        if (t_sam > obs.start & t_sam < t_lim2) {

          ttree[i, 2] <- 1

          obs <- rbind(obs, c(t_sam, i))

          prim_host <- c(prim_host, i)
          prim_obs <- c(prim_obs, t_sam)

          j <- 1

          while (j <= sec.n) {

            t_sam <- t_sam + sec.t

            obs <- rbind(obs, c(t_sam, i))

            j <- j + 1

          }

        } else {

          ttree[i, 2] <- NA

        }

      } else {

        ttree[i, 2] <- NA

      }

      if (length(prim_obs) > host.lim) {

        t_lim1 <- sort(prim_obs)[host.lim]
        t_lim2 <- sort(prim_obs)[host.lim + 1]

      }

      offspring <- rnbinom(1, size = off.r, prob = off.p)

      if (offspring > 0) {

        off_times <- t_inf + rgamma(offspring, shape = w.shape, scale = w.scale)

        off_times <- off_times[which(off_times < t_lim1)]

        offspring_inc <- length(off_times)

        if (offspring_inc > 0) {

          ttree <- rbind(ttree, matrix(0, offspring_inc, 3))

          ttree[n + 1:offspring_inc, 1] <- off_times

          ttree[n + 1:offspring_inc, 3] <- i

          n <- n + offspring_inc

        }

      }

      i <- i + 1

      if (i > nrow(ttree)) {

        break

      } else {

        t_inf <- ttree[i, 1]

      }

    }

    if (length(prim_obs) < host.lim) {

      next

    }

    keep_hosts <- prim_host[order(prim_obs)[1:host.lim]]

    keep_obs <- which(obs[, 2] %in% keep_hosts)
    drop_obs <- which(!(obs[, 2] %in% keep_hosts))

    ttree[!(1:dim(ttree)[1] %in% obs[keep_obs, 2]), 2] <- NA

    obs <- obs[keep_obs, , drop = F]

    repeat {

      torem <- which(is.na(ttree[, 2]) & !(1:dim(ttree)[1] %in% ttree[, 3]))

      if (length(torem) == 0) {

        break

      }

      for (j in 1:length(torem)) {

        h <- torem[j]

        ttree <- ttree[-h, , drop = FALSE]

        ttree[which(ttree[, 3] > h), 3] <- ttree[which(ttree[, 3] > h), 3] - 1

        obs[which(obs[, 2] > h), 2] <- obs[which(obs[, 2] > h), 2] - 1

        torem[which(torem > h)] <- torem[which(torem > h)] - 1

      }

    }

    ord <- c(which(!is.na(ttree[,2])), which(is.na(ttree[,2])))
    invord <- 1:length(ord)
    invord[ord] <- 1:length(ord)
    ttree <- ttree[ord, , drop=FALSE]
    ttree[ttree[, 3] > 0, 3] <- invord[ttree[ttree[, 3] > 0, 3]]
    obs[, 2] <- invord[obs[, 2]]

    if (dim(ttree)[1] == 0) {

      next

    }

    phy_ord <- which(!(1:length(ttree[, 1]) %in% ttree[, 3]))

    repeat {

      infectors_inc <- ttree[phy_ord, 3]
      infectors_exc <- ttree[-phy_ord, 3]

      phy_ord_new <- unique(infectors_inc[which(!(infectors_inc %in% infectors_exc))])
      phy_ord_new <- phy_ord_new[which(!(phy_ord_new %in% phy_ord))]

      if (sum(phy_ord_new) == 0) {

        break

      }

      phy_ord <- c(phy_ord, phy_ord_new)

    }

    sams <- length(obs[, 1])

    ctree <- matrix(0, sams, 4)
    ctree[, 1] <- obs[, 1]
    ctree[, 4] <- obs[, 2]

    for (i in 1:length(phy_ord)) {

      host <- phy_ord[i]

      leaves <- which(ctree[, 4] == host)

      if (length(leaves) == 1) {

        ctree <- rbind(ctree, rep(0, 4))

        ctree[length(ctree[, 1]), ] <- c(ttree[host, 1], leaves, 0, ttree[host, 3])

      } else {

        inf_time <- ttree[host, 1]

        leaf_times <- ctree[leaves, 1]
        leaf_order <- order(leaf_times, decreasing = T)

        l <- 1

        lin <- c(leaves[leaf_order[l]])
        t <- leaf_times[leaf_order[l]]

        repeat {

          if (length(lin) == 1) {

            l <- l + 1

            if (l > length(leaves)) {

              ctree <- rbind(ctree, rep(0, 4))

              ctree[length(ctree[, 1]), ] <- c(inf_time, lin, 0, ttree[host, 3])

              break

            } else {

              lin <- c(lin, leaves[leaf_order[l]])
              t <- leaf_times[leaf_order[l]]

            }

          } else {

            if (l == length(leaves)) {

              lim_time <- inf_time

            } else {

              lim_time <- leaf_times[leaf_order[l + 1]]

            }

            if (lambda == 0) {

              coa_time <- t + log(1 - runif(1)) * kappa / choose(length(lin), 2)

            } else {

              coa_time <- inf_time + (1 / lambda) * ((1 - runif(1)) ^ (lambda / choose(length(lin), 2)) * (lambda * (t - inf_time) + kappa) - kappa)

            }

            if (coa_time > lim_time) {

              sam <- sample(1:choose(length(lin), 2), size = 1)

              children_idx <- combn(length(lin), 2)[, sam]
              children <- lin[children_idx]

              ctree <- rbind(ctree, rep(0, 4))

              ctree[length(ctree[, 1]), ] <- c(coa_time, children, host)

              lin <- c(lin[which(!(lin %in% children))], length(ctree[, 1]))
              t <- coa_time

            } else {

              l <- l + 1

              if (l > length(leaves)) {

                for (j in 1:length(lin)) {

                  ctree <- rbind(ctree, rep(0, 4))

                  ctree[length(ctree[, 1]), ] <- c(inf_time, lin[j], 0, ttree[host, 3])

                }

                break

              } else {

                lin <- c(lin, leaves[leaf_order[l]])
                t <- leaf_times[leaf_order[l]]

              }

            }

          }

        }

      }

    }

    if (length(which(ctree[, 4] == 0)) == 1) {

      break

    }

  }

  ord <- c(1:sams, sams + order(ctree[(sams + 1):length(ctree[, 1]), 1], decreasing = T))
  invord = 1:length(ord)
  invord[ord] = 1:length(ord)

  ctree[, ] <- ctree[ord, ]
  ctree[which(ctree[, 2] > 0), 2] <- invord[ctree[which(ctree[, 2] > 0), 2]]
  ctree[which(ctree[, 3] > 0), 3] <- invord[ctree[which(ctree[, 3] > 0), 3]]

  ctree <- order_hosts(ctree)

  nam_host <- ctree[1:sams, 4]
  nam_num <- integer(sams)

  host_count <- integer(length(unique(nam_host)))

  for (i in 1:sams) {

    host_count[nam_host[i]] <- host_count[nam_host[i]] + 1

    nam_num[i] <- host_count[nam_host[i]]

  }

  l <- list(ctree = ctree, nam = paste(nam_host, ".", nam_num, sep = ""))
  class(l) <- 'ctree'

  l <- trim_root(l)

  return(list(ctree = l, obs.end = 0.5 * (t_lim1 + t_lim2), attempts = attempts))

}

