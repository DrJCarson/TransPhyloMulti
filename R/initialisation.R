#' Create a transmission tree compatible with the provided phylogenetic tree
#' @param ptree Phylogenetic tree
#'
#' @export
init_ctree <- function(ptree) {

  nam <- ptree$nam
  host <- ptree$host
  ptree <- ptree$ptree

  branch_host <- rep(0, nrow(ptree))
  branch_host[1:length(host)] <- host

  parents <- rep(NA, nrow(ptree))
  parents[ptree[, 2:3] + 1] <- 1:nrow(ptree)
  parents <- parents[-1]

  branch_low_time <- ptree[parents, 1]
  branch_upp_time <- ptree[, 1]

  host_low_time <- rep(0, max(host))
  host_upp_time <- rep(0, max(host))

  notposs <- 0
  for (h in 1:max(host)) {

    repeat {

      r <- which(branch_host == h)

      min_upp_time <- min(branch_upp_time[r])
      max_low_time <- max(branch_low_time[r])

      if (min_upp_time > max_low_time) {

        host_low_time[h] <- max_low_time
        host_upp_time[h] <- min_upp_time

        break

      } else {

        rback <- which(branch_low_time[r] > min_upp_time)

        s <- r[rback]

        s2 <- parents[s]

        branch_host[s] <- 0

        if (prod(branch_host[s2] %in% c(0, h))) {

          branch_host[s2] <- h

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

  ctree <- ptree
  ctree <- cbind(ctree, rep(0, nrow(ptree)))
  ctree[1:length(host), 4] <- host

  new_host <- max(host) + 1

  host_inf_time <- 0.5 * (host_low_time + host_upp_time)

  for (i in 1:length(host)) {

    r <- i

    h <- ctree[i, 4]

    repeat {

      anc <- which(ctree[, 2] == r | ctree[, 3] == r)
      anc_time <- ctree[anc, 4]

      if (ctree[anc, 1] > host_inf_time[h]) {

        r <- anc
        ctree[r, 4] <- h

      } else if (ctree[anc, 3] > 0) {

        j <- min(which(ctree[, 1] < host_inf_time[h] & ctree[, 2] != 0))

        ctree <- rbind(ctree[1:(j - 1), ], rep(0, 4), ctree[j:nrow(ctree), ])

        anc <- anc + 1
        ctree[, 2:3][which(ctree[, 2:3] >= j)] <- ctree[, 2:3][which(ctree[, 2:3] >= j)] + 1
        if (ctree[anc, 2] == r) {

          ctree[anc, 2] <- j

        } else {

          ctree[anc, 3] <- j

        }

        ctree[anc, 4] <- new_host

        ctree[j, 1] <- host_inf_time[h]
        ctree[j, 2] <- r
        ctree[j, 4] <- new_host

        break

      } else {

        break

      }

    }

  }

  ctree[which(ctree[, 4] == 0), 4] <- new_host
  ctree <- rbind(ctree, c(ctree[nrow(ctree), 1] - 1, nrow(ctree), 0, 0))

  l <- list(ctree = ctree, nam = nam)
  class(l) <- 'ctree'

  return(l)

}
