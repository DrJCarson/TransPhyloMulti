#' Create a transmission tree compatible with the provided phylogenetic tree
#' @param ptree Phylogenetic tree
#'
#' @export
init_ctree <- function(ptree, bottleneck = 1) {

  nam <- ptree$nam
  host <- ptree$host
  ptree <- ptree$ptree

  # Create ctree from ptree
  ctree <- ptree
  ctree <- cbind(ctree, rep(0, nrow(ptree)))
  ctree[1:length(host), 4] <- host


  n_host <- max(host)

  host_count <- numeric(n_host)

  for (i in 1:n_host) {

    host_count[i] <- length(which(host == i))

  }

  n_sam <- length(host)


  new_host <- n_host + 1

  ctree[length(ctree[, 1]), 4] <- new_host
  ctree <- rbind(ctree, c(min(ctree[, 1]) - 1, length(ctree[, 1]), 0, 0))

  # Determine parent for each row of ctree
  parents <- rep(NA, nrow(ctree))
  parents[ctree[, 2:3] + 1] <- 1:nrow(ctree)
  parents <- parents[-1]
  parents <- c(parents, NA)



  # For each row of ctree, determine which hosts are downstream
  host_downstream <- array(numeric(length(ctree[, 1]) * max(host)),
                           dim = c(length(ctree[, 1]), max(host)))

  for (i in 1:n_sam) {

    h <- host[i]
    b <- i

    repeat {

      host_downstream[b, h] <- host_downstream[b, h] + 1

      b <- parents[b]

      if (is.na(b)) {

        break

      }

    }

  }


  branch_upp <- ctree[1:(length(ctree[, 1]) - 1), 1]
  branch_low <- ctree[parents[1:(length(ctree[, 1]) - 1)], 1]


  r <- length(ctree[, 1]) - 1

  repeat {

    if (ctree[r, 2] > 0 & ctree[r, 3] > 0) {

      d1 <- ctree[r, 2]
      d2 <- ctree[r, 3]

      if (ctree[d1, 1] > ctree[d2, 1]) {

        d3 <- d2
        d2 <- d1
        d1 <- d3

      }

      t1 <- 0.5 * (ctree[d1, 1] + ctree[r, 1])

      b1 <- d1
      c1 <- host_downstream[d1, ]
      h1 <- which(c1 > 0)

      ol1 <- which(branch_upp > t1 & branch_low < t1)
      ol1 <- ol1[which(!(ol1 %in% b1))]

      repeat {

        if (length(ol1) == 0) {

          break

        }

        if (length(h1) == 1) {

          to_add <- which(host_downstream[ol1, h1] > 0)

        } else {

          if (length(ol1) == 1) {

            to_add <- which(max(host_downstream[ol1, h1]) > 0)

          } else {

            to_add <- which(apply(host_downstream[ol1, h1], 1, max) > 0)

          }

        }

        if (length(to_add) == 0) {

          break

        }

        b1 <- c(b1, ol1[to_add])
        c1 <- apply(host_downstream[b1, , drop = F], 2, sum)
        h1 <- which(c1 > 0)

        ol1 <- ol1[which(!(ol1 %in% b1))]

      }

      if (all(c1[h1] == host_count[h1])) {

        ctree[which(ctree[, 4] >= new_host), 4] <- ctree[which(ctree[, 4] >= new_host), 4] + 1

        ctree[b1, 4] <- new_host

        for (i in 1:length(b1)) {

          new_row <- min(which(ctree[, 1] < t1 & ctree[, 2] > 0))

          ctree <- rbind(ctree[1:(new_row - 1), ], rep(0, 4), ctree[(new_row):length(ctree[, 1]), ])

          r <- r + 1

          ctree[which(ctree[, 2] >= new_row), 2] <- ctree[which(ctree[, 2] >= new_row), 2] + 1
          ctree[which(ctree[, 3] >= new_row), 3] <- ctree[which(ctree[, 3] >= new_row), 3] + 1

          ctree[which(ctree[, 2] == b1[i]), 2] <- new_row
          ctree[which(ctree[, 3] == b1[i]), 3] <- new_row

          ctree[new_row, ] <- c(t1, b1[i], 0, ctree[r, 4])

          branch_low <- c(branch_low[1:(new_row - 1)], branch_low[b1[i]], branch_low[(new_row):length(ctree[, 1])])
          branch_upp <- c(branch_upp[1:(new_row - 1)], t1, branch_upp[(new_row):length(ctree[, 1])])

          branch_low[b1[i]] <- t1


          host_downstream <- rbind(host_downstream[1:(new_row - 1), ], host_downstream[b1[i], ], host_downstream[(new_row):length(host_downstream[, 1]), ])

        }

      } else {

        ctree[b1, 4] <- ctree[r, 4]

      }

      if (!(d2 %in% b1)) {

        t2 <- 0.5 * (ctree[d2, 1] + ctree[r, 1])

        b2 <- d2
        c2 <- host_downstream[d2, ]
        h2 <- which(c2 > 0)

        ol2 <- which(branch_upp > t2 & branch_low < t2)
        ol2 <- ol2[which(!(ol2 %in% b2))]

        repeat {

          if (length(ol2) == 0) {

            break

          }

          if (length(h2) == 1) {

            to_add <- which(host_downstream[ol2, h2] > 0)

          } else {

            if (length(ol2) == 1) {

              to_add <- which(max(host_downstream[ol2, h2]) > 0)

            } else {

              to_add <- which(apply(host_downstream[ol2, h2], 1, max) > 0)

            }

          }

          if (length(to_add) == 0) {

            break

          }

          b2 <- c(b2, ol2[to_add])
          c2 <- apply(host_downstream[b2, , drop = F], 2, sum)
          h2 <- which(c2 > 0)

          ol2 <- ol2[which(!(ol2 %in% b2))]

        }


        if (all(c2[h2] == host_count[h2])) {

          ctree[which(ctree[, 4] >= new_host), 4] <- ctree[which(ctree[, 4] >= new_host), 4] + 1

          ctree[b2, 4] <- new_host

          for (i in 1:length(b2)) {

            new_row <- min(which(ctree[, 1] < t2 & ctree[, 2] > 0))

            ctree <- rbind(ctree[1:(new_row - 1), ], rep(0, 4), ctree[(new_row):length(ctree[, 1]), ])

            r <- r + 1

            ctree[which(ctree[, 2] >= new_row), 2] <- ctree[which(ctree[, 2] >= new_row), 2] + 1
            ctree[which(ctree[, 3] >= new_row), 3] <- ctree[which(ctree[, 3] >= new_row), 3] + 1

            ctree[which(ctree[, 2] == b2[i]), 2] <- new_row
            ctree[which(ctree[, 3] == b2[i]), 3] <- new_row

            ctree[new_row, ] <- c(t2, b2[i], 0, ctree[r, 4])

            branch_low <- c(branch_low[1:(new_row - 1)], branch_low[b2[i]], branch_low[(new_row):length(ctree[, 1])])
            branch_upp <- c(branch_upp[1:(new_row - 1)], t2, branch_upp[(new_row):length(ctree[, 1])])

            branch_low[b2[i]] <- t2

            host_downstream <- rbind(host_downstream[1:(new_row - 1), ], host_downstream[b2[i], ], host_downstream[(new_row):length(host_downstream[, 1]), ])


          }

        }else {

          ctree[b2, 4] <- ctree[r, 4]

        }

      }

    }


    r <- r - 1

    if (r == 0) {

      break

    }


  }

  #### Order these better

  ctree[which(ctree[, 4] > 0), 4] <- ctree[which(ctree[, 4] > 0), 4] + max(ctree[, 4])

  unique_hosts <- unique(ctree[which(ctree[, 4] > 0), 4])

  for (i in 1:length(unique_hosts)) {

    ctree[which(ctree[, 4] == unique_hosts[i]), 4] <- i

  }

  ctree <- order_hosts(ctree)

  if (!all(ctree[1:n_sam, 4] == host)) {

    stop('Not possible to fit a transmission tree to the provided phylogenetic tree.')

  }


  #################################################################################


  if (bottleneck == 1) {

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

            hmin <- min(c(h1, h2))
            hmax <- max(c(h1, h2))

            ctree[which(ctree[, 4] == hmax), 4] <- hmin
            ctree[which(ctree[, 4] > hmax), 4] <- ctree[which(ctree[, 4] > hmax), 4] - 1

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

            ord[which(ord == hmax)] <- hmin
            ord[which(ord > hmax)] <- ord[which(ord > hmax)] - 1

            host_par[hmin] <- host_par[hmax]
            host_par[which(host_par == hmax)] <- hmin
            host_par <- host_par[-hmax]
            host_par[which(host_par > hmax)] <- host_par[which(host_par > hmax)] - 1

            host_lin[hmin] <- host_lin[hmax]
            host_lin <- host_lin[-hmax]

          }

        }

      }

      if (changes == 0) {

        break

      }

    }

  }


  l <- list(ctree = ctree, nam = nam)
  class(l) <- 'ctree'

  return(l)

}
