#' Proposal to add a transmission
#'
#' @param ctree Current coloured tree.
#' @param host Proposed infector.
#' @param u1 Random number that determines transmission structure.
#' @param u2 Random number that determines transmission time.
#' @param bn_weight Penalisation term for transmitting multiple lineages.
#'
#' @export
add_transmission <- function(ctree, host, u1, u2, bn_weight = 0.1) {

  nam <- ctree$nam
  ctree <- ctree$ctree

  # Rows for transmission or observations in host
  leaves <- which(ctree[, 3] == 0 & ctree[, 4] == host)
  L <- length(leaves)

  # Event type (observation = 1, transmission = 2) for each leaf
  type <- 1 * (ctree[leaves, 2] == 0 & ctree[leaves, 3] == 0) +
    2 * (ctree[leaves, 2] > 0 & ctree[leaves, 3] == 0)

  # Corresponding host to type (sampled or infected individual)
  host2 <- rep(host, length(leaves))
  host2[which(type == 2)] <- ctree[ctree[leaves[which(type == 2)], 2], 4]

  # All rows assigned to host
  rows_host <- which(ctree[, 4] == host)
  R <- length(rows_host)

  # Active leaves for each row and corresponding time interval
  active <- matrix(0, R, L)
  interval <- matrix(0, R, 2)

  for (l in 1:L) {

    r <- leaves[l]

    repeat {

      active[which(rows_host == r), l] <- 1

      anc <- which(ctree[, 2] == r | ctree[, 3] == r)

      interval[which(rows_host == r), ] <- c(ctree[anc, 1], ctree[r, 1])

      if (ctree[r, 4] != ctree[anc, 4]) {

        break

      } else {

        r <- anc

      }

    }

  }


  # Extend active array to account for non-bottleneck combinations
  active_ex <- active

  # Branch lengths for possible combinations
  len_ex <- interval[, 2] - interval[, 1]

  # Set length to zero for invalid combinations
  for (r in 1:R) {

    l1 <- which(active[r, ] == 1)
    l2 <- which(active[r, ] == 0)

    if (length(host2[l1[which(type[l1] == 1)]]) > 0 & length(host2[l2[which(type[l2] == 1)]]) > 0) {

      if (any(host2[l1[which(type[l1] == 1)]] %in% host2[l2[which(type[l2] == 1)]])) {

        len_ex[r] <- 0

      }

    }


    if (length(host2[l1[which(type[l1] == 2)]]) > 0 & length(host2[l2[which(type[l2] == 2)]]) > 0) {

      if (any(host2[l1[which(type[l1] == 2)]] %in% host2[l2[which(type[l2] == 2)]])) {

        len_ex[r] <- 0

      }

    }

  }

  # Unique identifier for combinations
  bin_ex <- numeric(R)

  for (l in 1:L) {

    bin_ex <- bin_ex + active_ex[, l] * 2^(l - 1)

  }


  # Determine lengths for possible combinations of branches
  if (R > 1) {

    for (r in 2:R) {

      # Possible combinations
      combs <- combn(1:R, r)

      for (c in 1:dim(combs)[2]) {

        inc_branches <- combs[, c]
        inc_active <- apply(active[inc_branches,], 2, max)

        # Do not include invalid combinations
        leaves_1 <- which(inc_active == 1)
        leaves_2 <- which(inc_active == 0)

        if (length(host2[leaves_1[which(type[leaves_1] == 1)]]) > 0 & length(host2[leaves_2[which(type[leaves_2] == 1)]]) > 0) {

          if (any(host2[leaves_1[which(type[leaves_1] == 1)]] %in% host2[leaves_2[which(type[leaves_2] == 1)]])) {

            next

          }

        }

        if (length(host2[leaves_1[which(type[leaves_1] == 2)]]) > 0 & length(host2[leaves_2[which(type[leaves_2] == 2)]]) > 0) {

          if (any(host2[leaves_1[which(type[leaves_1] == 2)]] %in% host2[leaves_2[which(type[leaves_2] == 2)]])) {

            next

          }

        }

        inc_bin <- sum(inc_active * 2^(0:(L - 1)))

        # Define interval over which combinations can occur
        low <- max(interval[inc_branches, 1])
        high <- min(interval[inc_branches, 2])

        # Update lengths
        if (high > low) {

          if (inc_bin %in% bin_ex) {

            k <- which(bin_ex == inc_bin)
            len_ex[k] <- len_ex[k] + (bn_weight ^ (r - 1)) * (high - low)

          } else {

            active_ex <- rbind(active_ex, inc_active)
            len_ex <- c(len_ex, (bn_weight ^ (r - 1)) * (high - low))
            bin_ex <- c(bin_ex, inc_bin)

          }

        }

      }

    }

  }

  sum_len <- sum(len_ex)
  norm_len <- len_ex / sum_len


  # Sample transmissions
  v <- 1
  s <- norm_len[v]

  while (u1 > s) {

    v <- v + 1
    s <- s + norm_len[v]

  }

  tar_bin <- bin_ex[v]
  tar_len <- len_ex[v]

  # Sample transmission time
  u2l <- u2 * tar_len
  sam_time <- NA

  for (r in 1:R) {

    combs <- combn(1:R, r)

    for (c in 1:dim(combs)[2]) {

      inc_branches <- combs[, c]

      if (r > 1) {

        inc_active <- apply(active[inc_branches,], 2, max)

      } else {

        inc_active <- active[inc_branches, ]

      }


      # Skip invalid combinations
      leaves_1 <- which(inc_active == 1)
      leaves_2 <- which(inc_active == 0)

      if (length(host2[leaves_1[which(type[leaves_1] == 1)]]) > 0 & length(host2[leaves_2[which(type[leaves_2] == 1)]]) > 0) {

        if (any(host2[leaves_1[which(type[leaves_1] == 1)]] %in% host2[leaves_2[which(type[leaves_2] == 1)]])) {

          next

        }

      }

      if (length(host2[leaves_1[which(type[leaves_1] == 2)]]) > 0 & length(host2[leaves_2[which(type[leaves_2] == 2)]]) > 0) {

        if (any(host2[leaves_1[which(type[leaves_1] == 2)]] %in% host2[leaves_2[which(type[leaves_2] == 2)]])) {

          next

        }

      }

      # If correct transmissions are formed, sample time
      inc_bin <- sum(inc_active * 2^(0:(L - 1)))

      low <- max(interval[inc_branches, 1])
      high <- min(interval[inc_branches, 2])

      if (inc_bin == tar_bin & high > low) {

        inc_len <- (bn_weight ^ (r - 1)) * (high - low)

        if (inc_len < u2l) {

          u2l <- u2l - inc_len

        } else {

          sam_time <- low + u2l / (bn_weight ^ (r - 1))
          break

        }

      }

    }

    if (!is.na(sam_time)) {

      break

    }

  }

  # Update ctree
  old_host <- host
  new_host <- host + 1

  if (new_host <= max(ctree[, 4])) {

    ctree[which(ctree[, 4] >= new_host), 4] <- ctree[which(ctree[, 4] >= new_host), 4] + 1

  }


  rhl <- rows_host[inc_branches]

  repeat {

    ctree[rhl[1], 4] <- new_host

    if (ctree[rhl[1], 2] != 0) {

      if (ctree[ctree[rhl[1], 2], 4] == old_host) {

        rhl <- c(rhl, ctree[rhl[1], 2])

      }

    }

    if (ctree[rhl[1], 3] != 0) {

      if (ctree[ctree[rhl[1], 3], 4] == old_host) {

        rhl <- c(rhl, ctree[rhl[1], 3])

      }

    }

    rhl <- rhl[-1]

    if (length(rhl) == 0) {

      break

    }

  }


  rhl <- sort(rows_host[inc_branches], decreasing = T)

  tar_row <- min(which(ctree[, 1] < sam_time & ctree[, 2] > 0))

  for (r in rhl) {

    ctree <- rbind(ctree[1:(tar_row - 1), ], rep(0, 4), ctree[(tar_row):length(ctree[, 1]), ])

    ctree[which(ctree[, 2] >= tar_row), 2] <- ctree[which(ctree[, 2] >= tar_row), 2] + 1
    ctree[which(ctree[, 3] >= tar_row), 3] <- ctree[which(ctree[, 3] >= tar_row), 3] + 1

    ctree[which(ctree[, 2] == r), 2] <- tar_row
    ctree[which(ctree[, 3] == r), 3] <- tar_row

    ctree[tar_row, ] <- c(sam_time, r, 0, old_host)

  }

  new_ctree <- list(ctree = ctree, nam = nam)
  class(new_ctree) <- 'ctree'

  return(new_ctree)

}
