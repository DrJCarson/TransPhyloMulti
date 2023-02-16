#' Create a transmission tree compatible with the provided phylogenetic tree
#' @param ptree Phylogenetic tree
#'
#' @export
init_ctree <- function(ptree) {

  nam <- ptree$nam
  host <- ptree$host
  ptree <- ptree$ptree

  # Create ctree from ptree
  ctree <- ptree
  ctree <- cbind(ctree, rep(0, nrow(ptree)))
  ctree[1:length(host), 4] <- host


  # Determine parent for each row of ctree
  parents <- rep(NA, nrow(ctree))
  parents[ctree[, 2:3] + 1] <- 1:nrow(ctree)
  parents <- parents[-1]
  parents <- c(parents, NA)

  # For each row of ctree, determine which hosts are downstream
  is_host_downstream <- array(numeric(length(ctree[, 1]) * max(host)),
                              dim = c(length(ctree[, 1]), max(host)))

  for (i in 1:length(host)) {

    h <- host[i]
    b <- i

    repeat {

      is_host_downstream[b, h] <- 1

      b <- parents[b]

      if (is.na(b)) {

        break

      }

    }

  }

  # Determine start and end times for each row / branch of ctree
  branch_low_time <- ctree[parents, 1]
  branch_upp_time <- ctree[1:nrow(ctree), 1]

  # Check that valid tree is possible
  notposs <- 0

  # Track hosts that have been added as infectors
  # Labels may be changed / merged later
  new_host <- max(host) + 1
  host_todo <- c()

  # Loop through sampled hosts
  for (h in 1:max(host)) {

    # Determine rows over which a transmission will occur
    tr_rows <- which(ctree[, 4] == h)

    repeat {

      # Check that all branches for the host overlap
      min_upp_time <- min(branch_upp_time[tr_rows])
      max_low_time <- max(branch_low_time[tr_rows])

      if (min_upp_time > max_low_time) {

        # Calculate infection time
        inf_time <- 0.5 * (min_upp_time + max_low_time)

        # Determine if additional branches need to be added due to downstream hosts
        repeat {

          if (length(tr_rows) == 1) {

            downstream_hosts <- which(is_host_downstream[tr_rows, ] == 1)

          } else {

            downstream_hosts <- which(apply(is_host_downstream[tr_rows, ], 2, max) == 1)

          }

          if (length(downstream_hosts) == 1) {

            inc_branch <- is_host_downstream[, downstream_hosts]

          } else {

            inc_branch <- apply(is_host_downstream[, downstream_hosts], 1, max)

          }

          # Proposed new set of rows
          prop_rows <- which(branch_low_time < inf_time &
                               branch_upp_time > inf_time &
                               inc_branch == 1)

          # If nothing is added, break the loop
          if (length(prop_rows) == length(tr_rows)) {

            break

          } else {

            tr_rows <- prop_rows

          }

        }

        # Add transmission to ctree
        for (i in 1:length(tr_rows)) {

          r <- tr_rows[i]

          anc <- which(ctree[, 2] == r | ctree[, 3] == r)

          new_row <- min(which(ctree[, 1] < inf_time & ctree[, 2] != 0))

          ctree <- rbind(ctree[1:(new_row - 1), ], rep(0, 4), ctree[new_row:nrow(ctree), ])

          # Update downstream hosts
          is_host_downstream <- rbind(is_host_downstream[1:(new_row - 1), ], rep(0, ncol(is_host_downstream)), is_host_downstream[new_row:nrow(is_host_downstream), ])

          is_host_downstream[new_row, ] <- is_host_downstream[r, ]

          ctree[, 2:3][which(ctree[, 2:3] >= new_row)] <- ctree[, 2:3][which(ctree[, 2:3] >= new_row)] + 1

          if (anc >= new_row) {

            anc <- anc + 1

          }

          if (ctree[anc, 2] == r) {

            ctree[anc, 2] <- new_row

          } else {

            ctree[anc, 3] <- new_row

          }

          ctree[new_row, ]  <- c(inf_time, r, 0, new_host)

        }

        # Update list of hosts
        host_todo <- c(host_todo, new_host)
        new_host <- new_host + 1

        # Update parents and branch times
        parents <- rep(NA, nrow(ctree))
        parents[ctree[, 2:3] + 1] <- 1:nrow(ctree)
        parents <- parents[-1]
        parents <- c(parents, NA)

        branch_low_time <- ctree[parents, 1]
        branch_upp_time <- ctree[1:nrow(ctree), 1]

        break

      } else {

        # If some branches do not overlap, iterate backwards in time
        iback <- which(branch_low_time[tr_rows] > min_upp_time)
        rback <- tr_rows[iback]

        prop_rows <- parents[rback]

        # If branch can not be assigned to host, there is no valid ctree
        if (prod(ctree[prop_rows, 4] %in% c(0, h))) {

          ctree[prop_rows, 4] <- h

          tr_rows <- c(tr_rows[-iback], prop_rows)

        } else {

          notposs <- 1

          break

        }

      }

    }

    if (notposs == 1) {

      stop('Not possible to fit a transmission tree to the provided phylogenetic tree.')

    }

  }


  # Loop through infectors and add new transmissions
  repeat {

    if (length(host_todo) == 0) {

      break

    }

    h <- host_todo[1]

    # Determine rows over which transmission occur
    tr_rows <- which(ctree[, 4] == h & ctree[, 3] == 0)

    tr_pars <- parents[tr_rows]

    iback <- which(ctree[tr_pars, 1] == max(ctree[tr_pars, 1]))
    rback <- tr_rows[iback]

    prop_rows <- c(tr_rows[-iback], tr_pars[iback])
    prop_rows <- unique(prop_rows)


    # Check if host is the root host
    if (any(prop_rows == nrow(ctree))) {

      ctree[tr_rows, 4] <- 0

      ctree[which(ctree[, 4] > h), 4] <- ctree[which(ctree[, 4] > h), 4] - 1

      host_todo[which(host_todo > h)] <- host_todo[which(host_todo > h)] - 1

      new_host <- new_host - 1

      host_todo <- host_todo[-1]

      # Check if host needs to be relabeled
    } else if (any(ctree[prop_rows, 4] > 0 & ctree[prop_rows, 4] != h)) {

      prop_host_vec <- unique(ctree[prop_rows, 4])

      ctree[tr_rows, 4] <- prop_host_vec[which(prop_host_vec != 0 & prop_host_vec != h)]

      ctree[which(ctree[, 4] > h), 4] <- ctree[which(ctree[, 4] > h), 4] - 1

      host_todo[which(host_todo > h)] <- host_todo[which(host_todo > h)] - 1

      new_host <- new_host - 1

      host_todo <- host_todo[-1]

    } else {

      anc_rows <- parents[prop_rows]

      min_time <- max(ctree[anc_rows, 1])
      max_time <- min(ctree[prop_rows, 1])

      inf_time <- 0.5 * (min_time + max_time)

      # Determine if additional branches need to be added due to downstream hosts
      repeat {

        if (length(prop_rows) == 1) {

          downstream_hosts <- which(is_host_downstream[prop_rows, ] == 1)

        } else {

          downstream_hosts <- which(apply(is_host_downstream[prop_rows, ], 2, max) == 1)

        }

        if (length(downstream_hosts) == 1) {

          inc_branch <- is_host_downstream[, downstream_hosts]

        } else {

          inc_branch <- apply(is_host_downstream[, downstream_hosts], 1, max)

        }

        # Proposed new set of rows
        prop_rows2 <- which(branch_low_time < inf_time &
                              branch_upp_time > inf_time &
                              inc_branch == 1)

        # If nothing is added, break the loop
        if (length(prop_rows2) == length(prop_rows)) {

          break

        } else {

          prop_rows <- prop_rows2

        }

      }

      ctree[prop_rows, 4] <- h

      for (r in prop_rows) {

        anc <- which(ctree[, 2] == r | ctree[, 3] == r)

        new_row <- min(which(ctree[, 1] < inf_time & ctree[, 2] != 0))

        ctree <- rbind(ctree[1:(new_row - 1), ], rep(0, 4), ctree[new_row:nrow(ctree), ])

        is_host_downstream <- rbind(is_host_downstream[1:(new_row - 1), ], rep(0, ncol(is_host_downstream)), is_host_downstream[new_row:nrow(is_host_downstream), ])

        is_host_downstream[new_row, ] <- is_host_downstream[r, ]

        ctree[, 2:3][which(ctree[, 2:3] >= new_row)] <- ctree[, 2:3][which(ctree[, 2:3] >= new_row)] + 1

        if (anc >= new_row) {

          anc <- anc + 1

        }

        if (ctree[anc, 2] == r) {

          ctree[anc, 2] <- new_row

        } else {

          ctree[anc, 3] <- new_row

        }

        ctree[new_row, ]  <- c(inf_time, r, 0, new_host)

      }

      host_todo <- c(host_todo, new_host)
      new_host <- new_host + 1

      # Update parents and branch times
      parents <- rep(NA, nrow(ctree))
      parents[ctree[, 2:3] + 1] <- 1:nrow(ctree)
      parents <- parents[-1]
      parents <- c(parents, NA)

      branch_low_time <- ctree[parents, 1]
      branch_upp_time <- ctree[1:nrow(ctree), 1]

      host_todo <- host_todo[-1]

    }

  }

  ctree[which(ctree[, 4] == 0), 4] <- new_host

  ctree <- rbind(ctree, c(ctree[nrow(ctree), 1] - 1, nrow(ctree), 0, 0))



  # Attempt bottlenecking
  obs_hosts <- unique(ctree[which(ctree[, 2] == 0 & ctree[, 3] == 0), 4])

  repeat {

    changes <- 0

    tr_rows <- which(ctree[, 2] > 0 & ctree[, 3] == 0)
    host1 <- ctree[tr_rows, 4]
    host2 <- ctree[ctree[tr_rows, 2], 4]

    host_todo <- 1:max(ctree[, 4])
    host_lin <- numeric(length(host_todo))
    host_par <- numeric(length(host_todo))

    for (i in 1:length(host_todo)) {

      h <- host_todo[i]

      host_lin[i] <- length(which(host2 == h))
      host_par[i] <- host1[which(host2 == h)][1]

    }

    ord <- order(host_lin, decreasing = T)

    for (i in 1:length(ord)) {

      h2 <- ord[i]
      h1 <- host_par[h2]

      l2 <- host_lin[h2]
      l1 <- host_lin[h1]

      if (h1 != 0 & !(h1 %in% obs_hosts & h2 %in% obs_hosts)) {

        if (l2 > l1) {

          changes <- changes + 1

          r1 <- which(ctree[, 2] > 0 & ctree[, 3] == 0 & ctree[, 4] == h1)
          rm_rows <- r1[which(ctree[ctree[r1, 2], 4] == h2)]
          rm_rows <- sort(rm_rows, decreasing = T)

          ctree[which(ctree[, 4] == h1), 4] <- h2
          ctree[which(ctree[, 4] > h1), 4] <- ctree[which(ctree[, 4] > h1), 4] - 1

          for (r in rm_rows) {

            nxt_row <- ctree[r, 2]

            prv_col <- 2
            prv_row <- which(ctree[, 2] == r)

            if (length(prv_row) == 0) {

              prv_col <- 3
              prv_row <- which(ctree[, 3] == r)

            }

            ctree[prv_row, prv_col] <- nxt_row

            ctree[, 2:3][which(ctree[, 2:3] > r)] <- ctree[, 2:3][which(ctree[, 2:3] > r)] - 1

            ctree <- ctree[-r, ]

          }

          ord[which(ord == h1)] <- h2
          ord[which(ord > h1)] <- ord[which(ord > h1)] - 1

          host_par[h2] <- host_par[h1]
          host_par[which(host_par == h1)] <- h2
          host_par <- host_par[-h1]
          host_par[which(host_par > h1)] <- host_par[which(host_par > h1)] - 1

          host_lin[h2] <- host_lin[h1]
          host_lin <- host_lin[-h1]

        }

      }

    }

    if (changes == 0) {

      break

    }

  }

  l <- list(ctree = ctree, nam = nam)
  class(l) <- 'ctree'

  return(l)

}
