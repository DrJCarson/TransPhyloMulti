#' Create a transmission tree compatible with the provided phylogenetic tree
#' @param ptree Phylogenetic tree
#'
#' @export
init_ctree <- function(ptree) {

  nam <- ptree$nam
  host <- ptree$host
  ptree <- ptree$ptree

  ctree <- ptree
  ctree <- cbind(ctree, rep(0, nrow(ptree)))
  ctree[1:length(host), 4] <- host


  tr_full <- rep(0, nrow(ctree))

  parents <- rep(NA, nrow(ctree))
  parents[ctree[, 2:3] + 1] <- 1:nrow(ctree)
  parents <- parents[-1]

  branch_low_time <- ctree[parents, 1]
  branch_upp_time <- ctree[1:(nrow(ctree) - 1), 1]

  host_low_time <- rep(0, max(host))
  host_upp_time <- rep(0, max(host))

  notposs <- 0
  for (h in 1:max(host)) {

    tr_rows <- which(ctree[, 4] == h)

    repeat {

      min_upp_time <- min(branch_upp_time[tr_rows])
      max_low_time <- max(branch_low_time[tr_rows])

      if (min_upp_time > max_low_time) {

        host_low_time[h] <- max_low_time
        host_upp_time[h] <- min_upp_time

        tr_full[tr_rows] <- h

        break

      } else {

        iback <- which(branch_low_time[tr_rows] > min_upp_time)
        rback <- tr_rows[iback]

        prop_rows <- parents[rback]

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


  for (h in 1:max(host)) {

    if (length(which(tr_full == h)) == 1) {

      next

    } else {

      tr_rows <- which(tr_full == h)

      repeat {

        iback <- which(branch_low_time[tr_rows] == max(branch_low_time[tr_rows]))
        rback <- tr_rows[iback]

        prop_rows <- parents[rback]
        prop_rows <- unique(prop_rows)

        if (any(is.na(prop_rows))) {

          break

        }

        if (any(ctree[prop_rows, 4] > 0)) {

          break

        } else {

          ctree[prop_rows, 4] <- h

          tr_rows <- c(tr_rows[-iback], prop_rows)

          if (length(tr_rows) == 1) {

            break

          }

        }

      }

      min_upp_time <- min(branch_upp_time[tr_rows])
      max_low_time <- max(branch_low_time[tr_rows])

      host_low_time[h] <- max_low_time
      host_upp_time[h] <- min_upp_time

      tr_full[which(tr_full == h)] <- 0
      tr_full[tr_rows] <- h

    }

  }


  host_inf_time <- 0.5 * (host_low_time + host_upp_time)

  for (h in 1:max(host)) {

    tr_rows <- which(tr_full == h)

    for (r in tr_rows) {

      new_row <- min(which(ctree[, 1] < host_inf_time[h] & ctree[, 2] != 0))

      ctree <- rbind(ctree[1:(new_row - 1), ], rep(0, 4), ctree[new_row:nrow(ctree), ])
      tr_full <- c(tr_full[1:(new_row - 1)], 0, tr_full[new_row:length(tr_full)])

      anc <- which(ctree[, 2] == r | ctree[, 3] == r)

      ctree[, 2:3][which(ctree[, 2:3] >= new_row)] <- ctree[, 2:3][which(ctree[, 2:3] >= new_row)] + 1

      if (ctree[anc, 2] == r) {

        ctree[anc, 2] <- new_row

      } else {

        ctree[anc, 3] <- new_row

      }

      ctree[new_row, ]  <- c(host_inf_time[h], r, 0, max(host) + h)
      tr_full[new_row] <- max(host) + h

    }

  }


  host_todo <- (max(host) + 1):(2 * max(host))

  repeat {

    if (length(host_todo) == 0) {

      break

    }

    h <- host_todo[1]

    tr_rows <- which(tr_full == h)

    if (length(tr_rows) == 1) {

      prop_row <- which(ctree[, 2] == tr_rows | ctree[, 3] == tr_rows)

      if (prop_row == nrow(ctree)) {

        ctree[tr_rows, 4] <- 0

        ctree[which(ctree[, 4] > h), 4] <- ctree[which(ctree[, 4] > h), 4] - 1

        host_todo[which(host_todo > h)] <- host_todo[which(host_todo > h)] - 1

        tr_full[which(tr_full == h)] <- 0
        tr_full[which(tr_full > h)] <- tr_full[which(tr_full > h)] - 1

        host_todo <- host_todo[-1]

      } else {

        if (ctree[prop_row, 4] > 0) {

          ctree[tr_rows, 4] <- ctree[prop_row, 4]

          ctree[which(ctree[, 4] > h), 4] <- ctree[which(ctree[, 4] > h), 4] - 1

          host_todo[which(host_todo > h)] <- host_todo[which(host_todo > h)] - 1

          tr_full[which(tr_full == h)] <- 0
          tr_full[which(tr_full > h)] <- tr_full[which(tr_full > h)] - 1

          host_todo <- host_todo[-1]

        } else {

          r <- prop_row

          tr_full[which(tr_full == h)] <- 0
          tr_full[r] <- h

          ctree[r, 4] <- h

          anc <- which(ctree[, 2] == r | ctree[, 3] == r)

          inf_time <- 0.5 * (ctree[anc, 1] + ctree[r, 1])

          new_row <- min(which(ctree[, 1] < inf_time & ctree[, 2] != 0))

          ctree <- rbind(ctree[1:(new_row - 1), ], rep(0, 4), ctree[new_row:nrow(ctree), ])
          tr_full <- c(tr_full[1:(new_row - 1)], 0, tr_full[new_row:length(tr_full)])

          ctree[, 2:3][which(ctree[, 2:3] >= new_row)] <- ctree[, 2:3][which(ctree[, 2:3] >= new_row)] + 1

          if (anc >= new_row) {

            anc <- anc + 1

          }

          if (ctree[anc, 2] == r) {

            ctree[anc, 2] <- new_row

          } else {

            ctree[anc, 3] <- new_row

          }

          ctree[new_row, ]  <- c(inf_time, r, 0, max(ctree[, 4]) + 1)
          tr_full[new_row] <- max(ctree[, 4])

          host_todo <- c(host_todo, max(ctree[, 4]))
          host_todo <- host_todo[-1]

        }

      }

    } else {

      tr_pars <- numeric(length(tr_rows))

      for (i in 1:length(tr_rows)) {

        tr_pars[i] <- which(ctree[, 2] == tr_rows[i] | ctree[, 3] == tr_rows[i])

      }

      iback <- which(ctree[tr_pars, 1] == max(ctree[tr_pars, 1]))
      rback <- tr_rows[iback]

      prop_rows <- c(tr_rows[-iback], tr_pars[iback])
      prop_rows <- unique(prop_rows)

      if (any(prop_rows == nrow(ctree))) {

        ctree[tr_rows, 4] <- 0

        ctree[which(ctree[, 4] > h), 4] <- ctree[which(ctree[, 4] > h), 4] - 1

        host_todo[which(host_todo > h)] <- host_todo[which(host_todo > h)] - 1

        tr_full[which(tr_full == h)] <- 0
        tr_full[which(tr_full > h)] <- tr_full[which(tr_full > h)] - 1

        host_todo <- host_todo[-1]

      } else if (any(ctree[prop_rows, 4] > 0 & ctree[prop_rows, 4] != h)) {

        ctree[tr_rows, 4] <- max(ctree[prop_rows, 4])

        ctree[which(ctree[, 4] > h), 4] <- ctree[which(ctree[, 4] > h), 4] - 1

        host_todo[which(host_todo > h)] <- host_todo[which(host_todo > h)] - 1

        tr_full[which(tr_full == h)] <- 0
        tr_full[which(tr_full > h)] <- tr_full[which(tr_full > h)] - 1

        host_todo <- host_todo[-1]

      } else {

        anc_rows <- numeric(length(prop_rows))

        for (i in 1:length(prop_rows)) {

          anc_rows[i] <- which(ctree[, 2] == prop_rows[i] | ctree[, 3] == prop_rows[i])

        }

        min_time <- max(ctree[anc_rows, 1])
        max_time <- min(ctree[prop_rows, 1])

        inf_time <- 0.5 * (min_time + max_time)

        tr_full[which(tr_full == h)] <- 0
        tr_full[prop_rows] <- h

        ctree[prop_rows, 4] <- h

        new_host <- max(ctree[, 4]) + 1

        for (r in prop_rows) {

          anc <- which(ctree[, 2] == r | ctree[, 3] == r)

          new_row <- min(which(ctree[, 1] < inf_time & ctree[, 2] != 0))

          ctree <- rbind(ctree[1:(new_row - 1), ], rep(0, 4), ctree[new_row:nrow(ctree), ])
          tr_full <- c(tr_full[1:(new_row - 1)], 0, tr_full[new_row:length(tr_full)])

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
          tr_full[new_row] <- new_host

        }

        host_todo <- c(host_todo, new_host)
        host_todo <- host_todo[-1]

      }

    }

  }

  ctree[which(ctree[, 4] == 0), 4] <- max(ctree[, 4]) + 1

  ctree <- rbind(ctree, c(ctree[nrow(ctree), 1] - 1, nrow(ctree), 0, 0))

  l <- list(ctree = ctree, nam = nam)
  class(l) <- 'ctree'

  return(l)

}
