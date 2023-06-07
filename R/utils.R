#' Numeric approximations (discretised) of branching process
#'
#' @param grid Discrete grid over which to evaluate functions
#' @param off.r Shape parameter for the number of offspring
#' @param off.p Probability parameter for the number of offspring
#' @param pi Probability of host being observed
#' @param w.shape Shape parameter of generation time distribution
#' @param w.scale Scale parameter of generation time distribution
#' @param ws.shape Shape parameter of primary sampling time distribution
#' @param ws.scale Scale parameter of primary sampling time distribution
#' @param obs.start Start date for observations
#' @param obs.end End date for observations
num_approx_disc <- function(grid,
                            off.r,
                            off.p,
                            pi,
                            w.shape,
                            w.scale,
                            ws.shape,
                            ws.scale,
                            obs.start,
                            obs.end) {

  grid.size <- length(grid)

  omega <- numeric(grid.size)
  omega_bar <- numeric(grid.size)
  phi <- numeric(grid.size)
  pit <- numeric(grid.size)

  omega[1] <- 1
  omega_bar[1] <- 1
  phi[1] <- 1

  pit <- pi * (pgamma(obs.end - grid, shape = ws.shape, scale = ws.scale) -
                 pgamma(obs.start - grid, shape = ws.shape, scale = ws.scale))

  gamma_prob <- pgamma(grid[1] - grid,  shape = w.shape, scale = w.scale) -
    pgamma(grid[2] - grid,  shape = w.shape, scale = w.scale)

  ft <- 1 - cumsum(gamma_prob)

  for (g in 2:grid.size) {

    omega_bar[g] <- ft[g] + sum(gamma_prob[g:2] * omega[1:(g - 1)])

    phi[g] <- (off.p / (1 - (1 - off.p) * omega_bar[g])) ^ off.r

    omega[g] <- (1 - pit[g]) * phi[g]

  }

  return(list(omega = omega,
              omega_bar = omega_bar,
              phi = phi,
              pit = pit,
              gamma_prob = gamma_prob))

}



#' Order hosts in a coloured tree
#'
#' @param ctree Current coloured tree
order_hosts <- function(ctree) {

  unique_hosts <- unique(ctree[, 4])

  n_hosts <- length(unique_hosts) - 1

  ordered_hosts <- order(unique_hosts[1:n_hosts])

  ctree[which(ctree[, 4] > 0), 4] <- ordered_hosts[ctree[which(ctree[, 4] > 0), 4]]

  return(ctree)

}


#' Remove extra root hosts from a coloured tree
#'
#' @param ctree Coloured tree
trim_root <- function(ctree) {

  nam <- ctree$nam
  ctree <- ctree$ctree

  repeat {

    root_row <- which(ctree[, 4] == 0)

    child_row <- ctree[root_row, 2]

    if (ctree[child_row, 2] > 0 & ctree[child_row, 3] == 0) {

      rem_host <- ctree[child_row, 4]

      ctree <- ctree[-root_row, ]

      ctree[child_row, 4] <- 0

      if (rem_host < max(ctree[, 4])) {

        ctree[which(ctree[, 4] > rem_host), 4] <- ctree[which(ctree[, 4] > rem_host), 4] - 1

      }

    } else {

      break

    }

  }

  ctree <- order_hosts(ctree)

  new_ctree <- list(ctree = ctree, nam = nam)
  class(new_ctree) <- 'ctree'

  return(new_ctree)

}



#' Truncate a combined tree
#'
#' @param ctree Coloured tree
#' @param trunc_time Time to truncate tree
#' @param n_obs Number of observations in final tree
trunc_ctree <- function(ctree, trunc_time = NA, n_obs = NA) {

  nam <- ctree$nam
  ctree <- ctree$ctree

  obs_idx <- which(ctree[, 2] == 0 & ctree[, 3] == 0)

  if (!is.na(trunc_time)) {

    if (max(ctree[, 1]) <= trunc_time) {

      stop('Truncation time is larger than the latest observation time.')

    }

    dateT <- trunc_time

  } else if (!is.na(n_obs)) {

    if (length(which(ctree[, 2] == 0 & ctree[, 3] == 0)) <= n_obs) {

      stop('Tree contains fewer observations than n_obs')

    } else {

      obs_ord <- order(ctree[obs_idx, 1], decreasing = F)

      trunc_time <- ctree[obs_ord[n_obs], 1]
      trunc_upper <- ctree[obs_ord[n_obs + 1], 1]

      dateT <- trunc_time + runif(1) * (trunc_upper - trunc_time)

    }

  } else {

    stop('Truncation parameter required.')

  }

  obs_keep <- which(ctree[obs_idx, 1] <= trunc_time)
  obs_rem <- which(ctree[obs_idx, 1] > trunc_time)

  rem_idx <- obs_rem

  while(length(rem_idx) > 0) {

    r_idx <- which(ctree[rem_idx, 1] == max(ctree[rem_idx, 1]))[1]

    r <- rem_idx[r_idx]
    p <- which(ctree[, 2] == r | ctree[, 3] == r)

    rhost <- ctree[r, 4]
    phost <- ctree[p, 4]

    if (ctree[r, 3] == 0) {

      rem_idx <- c(rem_idx, p)

      if (ctree[p, 2] == r) {

        ctree[p, 2] <- -1

      } else {

        ctree[p, 3] <- -1

      }

    } else {

      if (ctree[p, 3] == 0) {

        if (ctree[r, 2] == -1 & ctree[r, 3] == -1) {

          rem_idx <- c(rem_idx, p)

          ctree[p, 2] <- -1

        } else {

          cc <- max(ctree[r, 2:3])

          ctree[p, 2] <- cc

        }

      } else {

        if (ctree[r, 2] == -1 & ctree[r, 3] == -1) {

          rem_idx <- c(rem_idx, p)

          if (ctree[p, 2] == r) {

            ctree[p, 2] <- -1

          } else {

            ctree[p, 3] <- -1

          }

        } else {

          cc <- max(ctree[r, 2:3])

          if (ctree[p, 2] == r) {

            ctree[p, 2] <- cc

          } else {

            ctree[p, 3] <- cc

          }

        }

      }

    }



    ctree <- rbind(ctree[1:(r - 1), ], ctree[(r + 1):length(ctree[, 1]), ])
    ctree[, 2:3][which(ctree[, 2:3] > r)] <- ctree[, 2:3][which(ctree[, 2:3] > r)] - 1
    rem_idx[which(rem_idx > r)] <- rem_idx[which(rem_idx > r)] - 1

    rem_idx <- rem_idx[-r_idx]
    rem_idx <- unique(rem_idx)

  }

  old_hosts <- sort(unique(ctree[, 4]))[-1]

  for (i in 1:length(old_hosts)) {

    ctree[which(ctree[, 4] == old_hosts[i]), 4] <- i

  }

  ctree <- order_hosts(ctree)
  nam <- nam[obs_keep]

  new_ctree <- list(ctree = ctree, nam = nam)
  class(new_ctree) <- 'ctree'

  return(list(ctree = new_ctree, dateT = dateT))

}


#' Calculate primary observation times
#' @param ptree Phylogenetic tree with host data
calc_prim_obs <- function(ptree) {

  host <- ptree$host
  ptree <- ptree$ptree

  prim_times <- numeric(max(host))

  for (i in 1:max(host)) {

    prim_times[i] <- min(ptree[which(host == i), 1])

  }

  return(prim_times)

}
