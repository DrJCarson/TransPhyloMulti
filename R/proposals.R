#' Proposal to add a transmission
#'
#' @param ctree Current coloured tree.
#' @param host Proposed infector.
#' @param u1 Random number that determines transmission structure.
#' @param u2 Random number that determines transmission time.
#' @param bn_weight Penalisation term for transmitting multiple lineages.
#'
#' @export
add_transmission <- function(ctree, bn_weight = 0.1) {

  # Extract ctree
  nam <- ctree$nam
  ctree <- ctree$ctree

  # Random numbers for sampling
  host <- sample(1:max(ctree[, 4]), size = 1)
  u1 <- runif(1)
  u2 <- runif(1)

  prop_density <- log(1 / max(ctree[, 4]))

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

    # Ensure samples of the same host are assigned to one host
    if (length(host2[l1[which(type[l1] == 1)]]) > 0 & length(host2[l2[which(type[l2] == 1)]]) > 0) {

      if (any(host2[l1[which(type[l1] == 1)]] %in% host2[l2[which(type[l2] == 1)]])) {

        len_ex[r] <- 0

      }

    }

    # Ensure transmitted lineages to the same host have a single infector
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

        # Included branches and leaves in this combination
        inc_branches <- combs[, c]
        inc_active <- apply(active[inc_branches,], 2, max)

        # Do not include invalid combinations
        leaves_1 <- which(inc_active == 1)
        leaves_2 <- which(inc_active == 0)

        # Ensure samples of the same host are assigned to one host
        if (length(host2[leaves_1[which(type[leaves_1] == 1)]]) > 0 & length(host2[leaves_2[which(type[leaves_2] == 1)]]) > 0) {

          if (any(host2[leaves_1[which(type[leaves_1] == 1)]] %in% host2[leaves_2[which(type[leaves_2] == 1)]])) {

            next

          }

        }

        # Ensure transmitted lineages to the same host have a single infector
        if (length(host2[leaves_1[which(type[leaves_1] == 2)]]) > 0 & length(host2[leaves_2[which(type[leaves_2] == 2)]]) > 0) {

          if (any(host2[leaves_1[which(type[leaves_1] == 2)]] %in% host2[leaves_2[which(type[leaves_2] == 2)]])) {

            next

          }

        }

        # Unique identifier for combination
        inc_bin <- sum(inc_active * 2^(0:(L - 1)))

        # Define interval over which combinations can occur
        low <- max(interval[inc_branches, 1])
        high <- min(interval[inc_branches, 2])

        # Update lengths
        if (high > low) {

          # If combination already has weight
          if (inc_bin %in% bin_ex) {

            k <- which(bin_ex == inc_bin)
            len_ex[k] <- len_ex[k] + (bn_weight ^ (r - 1)) * (high - low)

          # If first time combination has occurred
          } else {

            active_ex <- rbind(active_ex, inc_active)
            len_ex <- c(len_ex, (bn_weight ^ (r - 1)) * (high - low))
            bin_ex <- c(bin_ex, inc_bin)

          }

        }

      }

    }

  }

  # Normalise lengths (weights)
  sum_len <- sum(len_ex)
  norm_len <- len_ex / sum_len


  # Sample combination
  v <- 1
  s <- norm_len[v]

  while (u1 > s) {

    v <- v + 1
    s <- s + norm_len[v]

  }

  prop_density <- prop_density + log(norm_len[v])

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

      # If correct transmission tree topology is formed, sample time
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

  prop_density <- prop_density + log((bn_weight ^ (r - 1)) / tar_len)

  tr_lin <- r

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

  # Order hosts
  ctree <- order_hosts(ctree)

  # Log density of reverse move
  # Sample rows and host
  obs_idx <- which(ctree[, 2] == 0 & ctree[, 3] == 0)
  obs_host <- ctree[obs_idx, 4]

  # Transmission rows and hosts
  tr_idx <- which(ctree[, 2] > 0 & ctree[, 3] == 0)
  tr_host1 <- ctree[tr_idx, 4]
  tr_host2 <- ctree[ctree[tr_idx, 2], 4]

  # Exclude transmissions between sampled hosts and root transmission
  tr_idx <- tr_idx[!((tr_host1 %in% obs_host) & (tr_host2 %in% obs_host)) & tr_host1 != 0]

  rev_density <- log(tr_lin / (length(tr_idx)))

  new_ctree <- list(ctree = ctree, nam = nam)
  class(new_ctree) <- 'ctree'

  return(list(ctree = new_ctree, prop_density = prop_density, rev_density = rev_density))

}


#' Proposal to remove a transmission
#'
#' @param ctree Current coloured tree.
#'
#' @export
remove_transmission <- function(ctree, bn_weight = 0.1) {

  # Extract ctree
  nam <- ctree$nam
  ctree <- ctree$ctree

  # Sample rows and host
  obs_idx <- which(ctree[, 2] == 0 & ctree[, 3] == 0)
  obs_host <- ctree[obs_idx, 4]

  # Transmission rows and hosts
  tr_idx <- which(ctree[, 2] > 0 & ctree[, 3] == 0)
  tr_host1 <- ctree[tr_idx, 4]
  tr_host2 <- ctree[ctree[tr_idx, 2], 4]

  # Exclude transmissions between sampled hosts and root transmission
  tr_idx <- tr_idx[!((tr_host1 %in% obs_host) & (tr_host2 %in% obs_host)) & tr_host1 != 0]

  # Sample transmission and corresponding ctree row
  sam_tr <- sample(1:length(tr_idx), size = 1)
  sam_ct <- tr_idx[sam_tr]

  # Sampled transmission information
  sam_time <- ctree[sam_ct, 1]
  sam_host1 <- ctree[sam_ct, 4]
  sam_host2 <- ctree[ctree[sam_ct, 2], 4]

  # All ctree rows associated with transmission
  sam_tr_ex <- which(ctree[tr_idx, 1] == sam_time &
                       ctree[tr_idx, 3] == 0 &
                       ctree[tr_idx, 4] == sam_host1 &
                       ctree[ctree[tr_idx, 2], 4] == sam_host2)

  sam_ct_ex <- tr_idx[sam_tr_ex]

  prop_density <- log(length(sam_tr_ex) / length(tr_idx))

  # Remove infected host
  orig_host <- ctree[, 4]
  ctree[which(ctree[, 4] == sam_host2), 4] <- sam_host1

  for (r in 1:length(sam_ct_ex)) {

    tr_row <- sam_ct_ex[r]
    nxt_row <- ctree[tr_row, 2]

    prv_col <- 2
    prv_row <- which(ctree[, 2] == tr_row)

    if (length(prv_row) == 0) {

      prv_col <- 3
      prv_row <- which(ctree[, 3] == tr_row)

    }

    ctree[prv_row, prv_col] <- nxt_row

  }

  for (r in length(sam_ct_ex):1) {

    ctree[, 2:3][which(ctree[, 2:3] > sam_ct_ex[r])] <-
      ctree[, 2:3][which(ctree[, 2:3] > sam_ct_ex[r])] - 1

  }

  orig_host <- orig_host[-sam_ct_ex]
  ctree <- ctree[-sam_ct_ex, ]

  # Update host numbers
  ctree[, 4][which(ctree[, 4] > sam_host2)] <-
    ctree[, 4][which(ctree[, 4] > sam_host2)] - 1


  rev_density <- log(1 / max(ctree[, 4]))

  host <- sam_host1 - 1

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

    # Ensure samples of the same host are assigned to one host
    if (length(host2[l1[which(type[l1] == 1)]]) > 0 & length(host2[l2[which(type[l2] == 1)]]) > 0) {

      if (any(host2[l1[which(type[l1] == 1)]] %in% host2[l2[which(type[l2] == 1)]])) {

        len_ex[r] <- 0

      }

    }

    # Ensure transmitted lineages to the same host have a single infector
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

  orig_bin <- sum(as.numeric(orig_host[leaves] == sam_host2) * 2 ^ (0:(L - 1)))

  # Determine lengths for possible combinations of branches
  if (R > 1) {

    for (r in 2:R) {

      # Possible combinations
      combs <- combn(1:R, r)

      for (c in 1:dim(combs)[2]) {

        # Included branches and leaves in this combination
        inc_branches <- combs[, c]
        inc_active <- apply(active[inc_branches,], 2, max)

        # Do not include invalid combinations
        leaves_1 <- which(inc_active == 1)
        leaves_2 <- which(inc_active == 0)

        # Ensure samples of the same host are assigned to one host
        if (length(host2[leaves_1[which(type[leaves_1] == 1)]]) > 0 & length(host2[leaves_2[which(type[leaves_2] == 1)]]) > 0) {

          if (any(host2[leaves_1[which(type[leaves_1] == 1)]] %in% host2[leaves_2[which(type[leaves_2] == 1)]])) {

            next

          }

        }

        # Ensure transmitted lineages to the same host have a single infector
        if (length(host2[leaves_1[which(type[leaves_1] == 2)]]) > 0 & length(host2[leaves_2[which(type[leaves_2] == 2)]]) > 0) {

          if (any(host2[leaves_1[which(type[leaves_1] == 2)]] %in% host2[leaves_2[which(type[leaves_2] == 2)]])) {

            next

          }

        }

        # Unique identifier for combination
        inc_bin <- sum(inc_active * 2^(0:(L - 1)))

        # Define interval over which combinations can occur
        low <- max(interval[inc_branches, 1])
        high <- min(interval[inc_branches, 2])

        # Update lengths
        if (high > low) {

          # If combination already has weight
          if (inc_bin %in% bin_ex) {

            k <- which(bin_ex == inc_bin)
            len_ex[k] <- len_ex[k] + (bn_weight ^ (r - 1)) * (high - low)

            # If first time combination has occurred
          } else {

            active_ex <- rbind(active_ex, inc_active)
            len_ex <- c(len_ex, (bn_weight ^ (r - 1)) * (high - low))
            bin_ex <- c(bin_ex, inc_bin)

          }

        }

      }

    }

  }

  # Normalise lengths (weights)
  sum_len <- sum(len_ex)
  norm_len <- len_ex / sum_len

  v <- which(bin_ex == orig_bin)

  rev_density <- rev_density + log(norm_len[v])

  rev_density <- rev_density + log((bn_weight ^ (length(sam_tr_ex) - 1)) / norm_len[v])

  # Order hosts
  ctree <- order_hosts(ctree)

  new_ctree <- list(ctree = ctree, nam = nam)
  class(new_ctree) <- 'ctree'

  return(list(ctree = new_ctree, prop_density = prop_density, rev_density = rev_density))

}


#' Proposal to remove a transmission, and add a new one with the same infector
#'
#' @param ctree Current coloured tree.
#' @param bn_weight Penalisation term for transmitting multiple lineages.
#' @param delta Standard deviation of random walk proposal for root update.
#'
#' @export
remove_add_local <- function(ctree, bn_weight = 0.1, delta = 1) {

  # Random numbers for sampling
  u1 <- runif(1)
  u2 <- runif(1)

  # Extract ctree
  nam <- ctree$nam
  ctree <- ctree$ctree

  # Transmission rows
  tr_idx <- which(ctree[, 2] > 0 & ctree[, 3] == 0)

  # Sample transmission and corresponding row in ctree
  sam_tr <- sample(1:length(tr_idx), size = 1)
  sam_ct <- tr_idx[sam_tr]

  # Transmission information
  sam_time <- ctree[sam_ct, 1]
  sam_host1 <- ctree[sam_ct, 4]
  sam_host2 <- ctree[ctree[sam_ct, 2], 4]

  # Propose new root time
  if (sam_host1 == 0) {

    # Perturb root time
    ctree[sam_ct, 1] <- sam_time + rnorm(1, mean = 0, sd = delta)

    # Reflect if proposed time exceeds next event time
    if (ctree[sam_ct, 1] > ctree[ctree[sam_ct, 2], 1]) {

      rebound <- ctree[sam_ct, 1] - ctree[ctree[sam_ct, 2], 1]

      ctree[sam_ct, 1] <- ctree[ctree[sam_ct, 2], 1] - rebound

    }

    prop_density <- log(1 / length(tr_idx))
    rev_density <- log(1 / length(tr_idx))

  # Remove transmission and replace it
  } else {

    ##################### Remove transmission ##################################

    # Find all ctree rows corresponding to transmission
    sam_tr_ex <- which(ctree[tr_idx, 1] == sam_time &
                         ctree[tr_idx, 3] == 0 &
                         ctree[tr_idx, 4] == sam_host1 &
                         ctree[ctree[tr_idx, 2], 4] == sam_host2)

    sam_ct_ex <- tr_idx[sam_tr_ex]

    prop_density <- log(length(sam_tr_ex) / length(tr_idx))

    # Hosts before removal
    orig_host <- ctree[, 4]

    # Replace removed host in ctree
    ctree[which(ctree[, 4] == sam_host2), 4] <- sam_host1

    # Remove sampled transmission rows from ctree
    for (r in 1:length(sam_ct_ex)) {

      tr_row <- sam_ct_ex[r]
      nxt_row <- ctree[tr_row, 2]

      prv_col <- 2
      prv_row <- which(ctree[, 2] == tr_row)

      if (length(prv_row) == 0) {

        prv_col <- 3
        prv_row <- which(ctree[, 3] == tr_row)

      }

      ctree[prv_row, prv_col] <- nxt_row

    }

    for (r in length(sam_ct_ex):1) {

      ctree[, 2:3][which(ctree[, 2:3] > sam_ct_ex[r])] <-
        ctree[, 2:3][which(ctree[, 2:3] > sam_ct_ex[r])] - 1

    }

    ctree <- ctree[-sam_ct_ex, ]
    orig_host <- orig_host[-sam_ct_ex]

    ##################### Re-add transmission ##################################

    # Rows for transmission or observations in host
    leaves <- which(ctree[, 3] == 0 & ctree[, 4] == sam_host1)
    L <- length(leaves)

    # Event type (observation = 1, transmission = 2) for each leaf
    type <- 1 * (ctree[leaves, 2] == 0 & ctree[leaves, 3] == 0) +
      2 * (ctree[leaves, 2] > 0 & ctree[leaves, 3] == 0)

    # Corresponding host to type (sampled or infected individual)
    host2 <- orig_host[leaves]
    host2[which(type == 2)] <- ctree[ctree[leaves[which(type == 2)], 2], 4]

    # All rows assigned to host
    rows_host <- which(ctree[, 4] == sam_host1)
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

      # Ensure samples of the same host are assigned to one host
      if (length(host2[l1[which(type[l1] == 1)]]) > 0 & length(host2[l2[which(type[l2] == 1)]]) > 0) {

        if (any(host2[l1[which(type[l1] == 1)]] %in% host2[l2[which(type[l2] == 1)]])) {

          len_ex[r] <- 0

        }

      }

      # Ensure transmitted lineages to the same host have a single infector
      if (length(host2[l1[which(type[l1] == 2)]]) > 0 & length(host2[l2[which(type[l2] == 2)]]) > 0) {

        if (any(host2[l1[which(type[l1] == 2)]] %in% host2[l2[which(type[l2] == 2)]])) {

          len_ex[r] <- 0

        }

      }

      # Ensure samples from different hosts are not assigned to a single host
      if (length(host2[l1[which(type[l1] == 1)]]) > 0) {

        if (length(unique(host2[l1[which(type[l1] == 1)]])) > 1) {

          len_ex[r] <- 0

        }

      }

      if (length(host2[l2[which(type[l2] == 1)]]) > 0) {

        if (length(unique(host2[l2[which(type[l2] == 1)]])) > 1) {

          len_ex[r] <- 0

        }

      }

    }

    # Unique identifier for combinations
    bin_ex <- numeric(R)

    for (l in 1:L) {

      bin_ex <- bin_ex + active_ex[, l] * 2^(l - 1)

    }

    orig_bin <- sum(as.numeric(orig_host[leaves] == sam_host2) * 2 ^ (0:(L - 1)))


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

          # Ensure samples of the same host are assigned to one host
          if (length(host2[leaves_1[which(type[leaves_1] == 1)]]) > 0 & length(host2[leaves_2[which(type[leaves_2] == 1)]]) > 0) {

            if (any(host2[leaves_1[which(type[leaves_1] == 1)]] %in% host2[leaves_2[which(type[leaves_2] == 1)]])) {

              next

            }

          }

          # Ensure transmitted lineages to the same host have a single infector
          if (length(host2[leaves_1[which(type[leaves_1] == 2)]]) > 0 & length(host2[leaves_2[which(type[leaves_2] == 2)]]) > 0) {

            if (any(host2[leaves_1[which(type[leaves_1] == 2)]] %in% host2[leaves_2[which(type[leaves_2] == 2)]])) {

              next

            }

          }

          # Ensure samples from different hosts are not assigned to a single host
          if (length(host2[leaves_1[which(type[leaves_1] == 1)]]) > 0) {

            if (length(unique(host2[leaves_1[which(type[leaves_1] == 1)]])) > 1) {

              next

            }

          }

          if (length(host2[leaves_2[which(type[leaves_2] == 1)]]) > 0) {

            if (length(unique(host2[leaves_2[which(type[leaves_2] == 1)]])) > 1) {

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

    prop_density <- prop_density + log(norm_len[v])

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

        # Ensure samples of the same host are assigned to one host
        if (length(host2[leaves_1[which(type[leaves_1] == 1)]]) > 0 & length(host2[leaves_2[which(type[leaves_2] == 1)]]) > 0) {

          if (any(host2[leaves_1[which(type[leaves_1] == 1)]] %in% host2[leaves_2[which(type[leaves_2] == 1)]])) {

            next

          }

        }

        # Ensure transmitted lineages to the same host have a single infector
        if (length(host2[leaves_1[which(type[leaves_1] == 2)]]) > 0 & length(host2[leaves_2[which(type[leaves_2] == 2)]]) > 0) {

          if (any(host2[leaves_1[which(type[leaves_1] == 2)]] %in% host2[leaves_2[which(type[leaves_2] == 2)]])) {

            next

          }

        }

        # Ensure samples from different hosts are not assigned to a single host
        if (length(host2[leaves_1[which(type[leaves_1] == 1)]]) > 0) {

          if (length(unique(host2[leaves_1[which(type[leaves_1] == 1)]])) > 1) {

            next

          }

        }

        if (length(host2[leaves_2[which(type[leaves_2] == 1)]]) > 0) {

          if (length(unique(host2[leaves_2[which(type[leaves_2] == 1)]])) > 1) {

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

    prop_density <- prop_density + log((bn_weight ^ (r - 1)) / tar_len)

    rev_density <- log(r / (length(tr_idx) - length(sam_tr_ex) + r))

    v <- which(bin_ex == orig_bin)

    rev_density <- rev_density + log(norm_len[v])

    rev_density <- rev_density + log((bn_weight ^ (length(sam_tr_ex) - 1)) / norm_len[v])

    # Update ctree
    rhl <- rows_host[inc_branches]

    repeat {

      ctree[rhl[1], 4] <- sam_host2

      if (ctree[rhl[1], 2] != 0) {

        if (ctree[ctree[rhl[1], 2], 4] == sam_host1) {

          rhl <- c(rhl, ctree[rhl[1], 2])

        }

      }

      if (ctree[rhl[1], 3] != 0) {

        if (ctree[ctree[rhl[1], 3], 4] == sam_host1) {

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

      ctree[tar_row, ] <- c(sam_time, r, 0, sam_host1)

    }

    ctree <- order_hosts(ctree)

  }

  new_ctree <- list(ctree = ctree, nam = nam)
  class(new_ctree) <- 'ctree'

  return(list(ctree = new_ctree, prop_density = prop_density, rev_density = rev_density))

}
