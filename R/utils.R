#' Order hosts in a coloured tree
#'
#' @param ctree Current coloured tree.
#'
#' @export
order_hosts <- function(ctree) {

  unique_hosts <- unique(ctree[, 4])

  n_hosts <- length(unique_hosts) - 1

  ordered_hosts <- order(unique_hosts[1:n_hosts])

  ctree[which(ctree[, 4] > 0), 4] <- ordered_hosts[ctree[which(ctree[, 4] > 0), 4]]

  return(ctree)

}



#' Truncate a combined tree
#'
#'
#' @export
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
