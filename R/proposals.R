
#' Proposal to move a transmission whilst maintaining the same host
#'
#' @param ctree Current coloured tree
#' @param epsilon Scale for root time update proposals
remove_add <- function(ctree, pen.prob = 0, pen.size = 1, pen.len = 1, epsilon = 1) {

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

  if (length(tr_host2) > 1) {

    sam_host2 <- sample(unique(tr_host2), size = 1)

  } else {

    sam_host2 <- tr_host2

  }

  sam_host1 <- tr_host1[which(tr_host2 == sam_host2)][1]

  sam_tr_ex <- which(tr_host2 == sam_host2)
  sam_ct_ex <- tr_idx[sam_tr_ex]

  sam_time <- ctree[sam_ct_ex[1], 1]

  # Note if both hosts are observed
  if (sam_host1 %in% obs_host & sam_host2 %in% obs_host) {

    both_obs <- 1

  } else {

    both_obs <- 0

  }

  # Propose new root time
  if (sam_host1 == 0) {

    curr_hosts <- sam_host2
    prop_hosts <- sam_host2

    # Perturb root time
    ctree[sam_ct_ex, 1] <- sam_time + epsilon * (runif(1) - 0.5)

    # Reflect if proposed time exceeds next event time
    if (ctree[sam_ct_ex, 1] > ctree[ctree[sam_ct_ex, 2], 1]) {

      rebound <- ctree[sam_ct_ex, 1] - ctree[ctree[sam_ct_ex, 2], 1]

      ctree[sam_ct_ex, 1] <- ctree[ctree[sam_ct_ex, 2], 1] - rebound

    }

    prop_density <- -log(length(unique(tr_host2)))
    rev_density <- -log(length(unique(tr_host2)))

    # Remove transmission and replace it
  } else {

    curr_hosts <- c(sam_host1, sam_host2)

    ##################### Remove transmission ##################################

    prop_density <- -log(length(unique(tr_host2)))
    rev_density <- -log(length(unique(tr_host2)))

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

    # Vector of hosts at leaves
    unique_hosts <- unique(host2)

    # Number of hosts at leaves
    n_hosts2 <- length(unique_hosts)

    # Number of leaves for each host
    hosts_count <- numeric(n_hosts2)

    # Leaf type for each host
    hosts_type <- numeric(n_hosts2)

    # Loop through possible hosts
    for (i in 1:n_hosts2) {

      # Update host counts and type
      hosts_count[i] <- length(which(host2 == unique_hosts[i]))
      hosts_type[i] <- type[which(host2 == unique_hosts[i])][1]

    }

    # All rows assigned to host
    rows_host <- which(ctree[, 4] == sam_host1)
    R <- length(rows_host)

    # Indicate which leaves are downstream for each branch
    down_leaves <- matrix(0, R, L)

    # Indicate which hosts are downstream for each branch
    # and count the number of relevant leaves
    down_hosts <- matrix(0, R, n_hosts2)

    # Time interval for each branch
    interval <- matrix(0, R, 2)

    # Loop through leaves
    for (l in 1:L) {

      # Initialise at row corresponding to leaf
      r <- leaves[l]

      # Store host corresponding to leaf
      h <- host2[l]

      # Loop backwards until infection time
      repeat {

        # Indicate that current leaf is downstream of this branch
        down_leaves[which(rows_host == r), l] <- 1

        # Update host count for this branch
        down_hosts[which(rows_host == r), which(unique_hosts == h)] <-
          down_hosts[which(rows_host == r), which(unique_hosts == h)] + 1

        # Calculate parent branch
        anc <- which(ctree[, 2] == r | ctree[, 3] == r)

        # Calculate time interval of this branch
        interval[which(rows_host == r), ] <- c(ctree[anc, 1], ctree[r, 1])

        # Determine if we have reached the infection time
        if (ctree[r, 4] != ctree[anc, 4]) {

          # Break loop at infection time
          break

        } else {

          # Iterate backwards if not at infection time
          r <- anc

        }

      }

    }


    # Evaluate break points along mini-tree
    # Leaf times, coalescent times, infection time
    clust_intervals <- sort(interval)
    clust_intervals <- unique(clust_intervals)

    # Determine branch clusters. Branch combinations for which a transmission must
    # occur over all members or none
    clusters <- list()

    # Store possible transmissions according to number of transmitted lineages
    tr_lin <- list()

    for (l in 1:L) {

      tr_lin[[l]] <- list()

    }

    # Store possible transmissions, and length of time over which transmission
    # could occur in each case
    poss_tr <- array(dim = c(0, 2))


    for (i in 1:(length(clust_intervals) - 1)) {

      # Start and end point of interval
      start <- clust_intervals[i]
      end <- clust_intervals[i + 1]

      # Set of clusters for current time interval
      clusters_int <- list()

      # Determine which branches overlap with current interval
      branches_todo <- which(!(interval[, 2] <= clust_intervals[i] | interval[, 1] >= clust_intervals[i + 1]))

      # Number of evaluated clusters for this interval + 1
      j <- 1

      repeat {

        # If all branches have been clustered, break the loops
        if (length(branches_todo) == 0) {

          break

        }

        # Initialise with first un-clustered branch
        # Consider adding a transmission to this branch
        branches_inc <- branches_todo[1]

        # Evaluate number of leaves for each host that would be assigned to the
        # newly added host
        counts_inc <- down_hosts[branches_inc, ]

        # Indicate which hosts would be affected by the newly added transmission
        hosts_inc <- which(counts_inc > 0)

        # Indicate which sampled hosts would be affected by the newly added transmission
        obs_inc <- unique_hosts[hosts_inc[which(hosts_type[hosts_inc] == 1)]]

        # Indicate that branch is now clustered
        branches_todo <- branches_todo[-1]

        repeat {

          # If no branches remain unclustered, break the loop
          if (length(branches_todo) == 0) {

            break

          }

          # If only a single branch remains unclustered
          if (length(branches_todo) == 1) {

            # If only one host would be affected in the current cluster
            if (length(hosts_inc) == 1) {

              # Determine which branches must be added
              to_add <- which(down_hosts[branches_todo, hosts_inc] > 0)

            } else {

              # Determine which branches must be added
              to_add <- which(max(down_hosts[branches_todo, hosts_inc]) > 0)

            }

          } else {

            # If only one host would be affected in the current cluster
            if (length(hosts_inc) == 1) {

              # Determine which branches must be added
              to_add <- which(down_hosts[branches_todo, hosts_inc] > 0)

            } else {

              # Determine which branches must be added
              to_add <- which(apply(down_hosts[branches_todo, hosts_inc], 1, max) > 0)

            }

          }

          # If no branches need to be added to the cluster, break the loop
          if (length(to_add) == 0) {

            break

          }

          # Update branches in cluster
          branches_inc <- c(branches_inc, branches_todo[to_add])

          # Update number of affected leaves for each host
          counts_inc <- apply(down_hosts[branches_inc, , drop = F], 2, sum)

          # Update vector of affected hosts
          hosts_inc <- which(counts_inc > 0)

          # Update vector of affected sampled hosts
          obs_inc <- unique_hosts[hosts_inc[which(hosts_type[hosts_inc] == 1)]]

          # Update vector of unclustered branches
          branches_todo <- branches_todo[-to_add]
          #        branches_todo <- branches_todo[which(!(branches_todo %in% branches_inc))]

        }

        # If all leaves for each host are either included or excluded
        # and at most one observed host is included
        if (all(counts_inc[hosts_inc] == hosts_count[hosts_inc]) &
            length(obs_inc) < 2) {

          # Store cluster information
          # branches: branches included in cluster
          # counts: downstream leaf counts for each host
          # hosts: affected hosts
          # bin: unique numerical identifier for included hosts
          # obs: affected sampled hosts
          clusters_int[[j]] <- list(branches = branches_inc,
                                    counts = counts_inc,
                                    hosts = hosts_inc,
                                    bin = sum(2 ^ (hosts_inc - 1)),
                                    obs = obs_inc)

          # Increment cluster count
          j <- j + 1

        }

      }

      # Store all clusters for current interval
      clusters[[i]] <- clusters_int

      # If there is at least one valid cluster within current interval
      if (length(clusters[[i]]) > 0) {

        # Loop through possible cluster sizes
        for (csize in 1:length(clusters[[i]])) {

          # Possible cluster combinations for given combination size
          combs <- combn(length(clusters[[i]]), csize)

          # Loop through possible cluster combinations for current size
          for (k in 1:dim(combs)[2]) {

            # Resulting branches, numerical ID, and sampled hosts from combination
            temp_branches <- c()
            temp_bin <- 0
            temp_obs <- c()

            # Add relevant clusters
            for (j in combs[, k]) {

              # Update branches, numerical ID, and sampled hosts
              temp_bin <- temp_bin + clusters[[i]][[j]]$bin
              temp_branches <- c(temp_branches, clusters[[i]][[j]]$branches)
              temp_obs <- c(temp_obs, clusters[[i]][[j]]$obs)

            }

            if (both_obs == 0 | length(temp_obs) == 1) {

              # Number of transmitted lineages
              lb <- length(temp_branches)

              tr_lin[[lb]][[length(tr_lin[[lb]]) + 1]] <- list(interval = i,
                                                               clusts = combs[, k],
                                                               branches = temp_branches,
                                                               bin = temp_bin,
                                                               start = start,
                                                               end = end,
                                                               length = end - start)


              # Determine if resulting transmission already has weight
              tr_ind <- which(poss_tr[, 1] == temp_bin)

              if (length(tr_ind) == 0) {

                # If the transmission is novel, add as a possibility
                poss_tr <- rbind(poss_tr, c(temp_bin, (end - start)))

              } else {

                # If transmission currently has weight, update the possible timespan
                poss_tr[tr_ind, 2] <- poss_tr[tr_ind, 2] + (end - start)

              }

            }

          }

        }

      }

    }

    len_lin <- numeric(L)

    for (l in 1:L) {

      temp_len <- sapply(tr_lin[[l]], function(x) x[["length"]])

      if (length(temp_len) > 0) {

        len_lin[l] <- sum(temp_len)

      } else {

        len_lin[l] <- 0

      }

    }

    poss_lin <- which(len_lin > 0)

    if (pen.len == 0) {

      if (pen.prob == 0) {

        prob_lin <- numeric(L)
        prob_lin[poss_lin] <- 1 / length(poss_lin)

      } else {

        prob_lin <- dnbinom(0:(L - 1), prob = pen.prob, size = pen.size) * ((1:L) %in% poss_lin)
        prob_lin <- prob_lin / sum(prob_lin)

      }

    } else {

      if (pen.prob == 0) {

        prob_lin <- len_lin / sum(len_lin)

      } else {

        prob_lin <- dnbinom(0:(L - 1), prob = pen.prob, size = pen.size) * len_lin
        prob_lin <- prob_lin / sum(prob_lin)

      }

    }

    sam_lin_p <- sample(1:L, size = 1, prob = prob_lin)
    sam_lin_r <- length(sam_ct_ex)

    prop_density <- prop_density + log(prob_lin[sam_lin_p])
    rev_density <- rev_density + log(prob_lin[sam_lin_r])

    temp_len_p <- sapply(tr_lin[[sam_lin_p]], function(x) x[["length"]])
    temp_len_r <- sapply(tr_lin[[sam_lin_r]], function(x) x[["length"]])

    prob_clust_p <- temp_len_p / sum(temp_len_p)
    prob_clust_r <- temp_len_r / sum(temp_len_r)

    sam_clust_p <- sample(1:length(prob_clust_p), size = 1, prob = prob_clust_p)

    temp_start_r <- sapply(tr_lin[[sam_lin_r]], function(x) x[["start"]])
    temp_end_r <- sapply(tr_lin[[sam_lin_r]], function(x) x[["end"]])
    temp_bin_r <- sapply(tr_lin[[sam_lin_r]], function(x) x[["bin"]])

    orig_bin <- sum(as.numeric(unique_hosts %in% unique(host2[orig_host[leaves] == sam_host2])) * 2 ^ (0:(length(unique_hosts) - 1)))

    sam_time_r <- sam_time

    sam_clust_r <- which(temp_start_r < sam_time_r & temp_end_r > sam_time_r & temp_bin_r == orig_bin)

    prop_density <- prop_density + log(prob_clust_p[sam_clust_p])
    rev_density <- rev_density + log(prob_clust_r[sam_clust_r])

    # Sample transmission time
    sam_time <- runif(1, min = tr_lin[[sam_lin_p]][[sam_clust_p]]$start,
                      max = tr_lin[[sam_lin_p]][[sam_clust_p]]$end)

    prop_density <- prop_density - log(tr_lin[[sam_lin_p]][[sam_clust_p]]$length)
    rev_density <- rev_density - log(tr_lin[[sam_lin_r]][[sam_clust_r]]$length)

    temp_branches <- tr_lin[[sam_lin_p]][[sam_clust_p]]$branches

    # Update ctree
    rhl <- rows_host[temp_branches]

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


    rhl <- sort(rows_host[temp_branches], decreasing = T)

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

    prop_hosts <- c(ctree[tar_row, 4], ctree[ctree[tar_row, 2], 4])

  }

  new_ctree <- list(ctree = ctree, nam = nam)
  class(new_ctree) <- 'ctree'

  return(list(ctree = new_ctree, prop_density = prop_density, rev_density = rev_density, is_possible = 1,
              curr_hosts = curr_hosts, prop_hosts = prop_hosts))

}


##############################################################################################

#' Proposal to add a transmission
#'
#' @param ctree Current coloured tree
add_transmission <- function(ctree, pen.prob = 0, pen.size = 1, pen.len = 1) {

  # Extract ctree
  nam <- ctree$nam
  ctree <- ctree$ctree

  # Total number of hosts
  n_hosts <- max(ctree[, 4])

  # Probability of sampling each host
  prob_hosts <- rep(1 / n_hosts, n_hosts)

  # Sample a host
  host <- sample(1:n_hosts, size = 1, prob = prob_hosts)

  # Store sampled host for likelihood calculations
  curr_hosts <- host

  # Log density of the proposal
  prop_density <- log(prob_hosts[host])

  # Rows for transmissions or observations in host
  leaves <- which(ctree[, 3] == 0 & ctree[, 4] == host)
  L <- length(leaves)

  # Leaf type (observation = 1, transmission = 2)
  type <- 1 * (ctree[leaves, 2] == 0 & ctree[leaves, 3] == 0) +
    2 * (ctree[leaves, 2] > 0 & ctree[leaves, 3] == 0)

  # Corresponding host to type (sampled or infected individual)
  host2 <- rep(host, L)
  host2[which(type == 2)] <- ctree[ctree[leaves[which(type == 2)], 2], 4]

  # Vector of hosts at leaves
  unique_hosts <- unique(host2)

  # Number of hosts at leaves
  n_hosts2 <- length(unique_hosts)

  # Number of leaves for each host
  hosts_count <- numeric(n_hosts2)

  # Leaf type for each host
  hosts_type <- numeric(n_hosts2)

  # Loop through possible hosts
  for (i in 1:n_hosts2) {

    # Update host counts and type
    hosts_count[i] <- length(which(host2 == unique_hosts[i]))
    hosts_type[i] <- type[which(host2 == unique_hosts[i])][1]

  }


  # All rows / branches assigned to host
  rows_host <- which(ctree[, 4] == host)
  R <- length(rows_host)


  # Indicate which leaves are downstream for each branch
  down_leaves <- matrix(0, R, L)

  # Indicate which hosts are downstream for each branch
  # and count the number of relevant leaves
  down_hosts <- matrix(0, R, n_hosts2)

  # Time interval for each branch
  interval <- matrix(0, R, 2)

  # Loop through leaves
  for (l in 1:L) {

    # Initialise at row corresponding to leaf
    r <- leaves[l]

    # Store host corresponding to leaf
    h <- host2[l]

    # Loop backwards until infection time
    repeat {

      # Indicate that current leaf is downstream of this branch
      down_leaves[which(rows_host == r), l] <- 1

      # Update host count for this branch
      down_hosts[which(rows_host == r), which(unique_hosts == h)] <-
        down_hosts[which(rows_host == r), which(unique_hosts == h)] + 1

      # Calculate parent branch
      anc <- which(ctree[, 2] == r | ctree[, 3] == r)

      # Calculate time interval of this branch
      interval[which(rows_host == r), ] <- c(ctree[anc, 1], ctree[r, 1])

      # Determine if we have reached the infection time
      if (ctree[r, 4] != ctree[anc, 4]) {

        # Break loop at infection time
        break

      } else {

        # Iterate backwards if not at infection time
        r <- anc

      }

    }

  }

  # Evaluate break points along mini-tree
  # Leaf times, coalescent times, infection time
  clust_intervals <- sort(interval)
  clust_intervals <- unique(clust_intervals)

  # Determine branch clusters. Branch combinations for which a transmission must
  # occur over all members or none
  clusters <- list()

  # Store possible transmissions according to number of transmitted lineages
  tr_lin <- list()

  for (l in 1:L) {

    tr_lin[[l]] <- list()

  }

  # Store possible transmissions, and length of time over which transmission
  # could occur in each case
  poss_tr <- array(dim = c(0, 2))

  # Loop through intervals of mini-tree
  for (i in 1:(length(clust_intervals) - 1)) {

    # Start and end point of interval
    start <- clust_intervals[i]
    end <- clust_intervals[i + 1]

    # Set of clusters for current time interval
    clusters_int <- list()

    # Determine which branches overlap with current interval
    branches_todo <- which(!(interval[, 2] <= clust_intervals[i] | interval[, 1] >= clust_intervals[i + 1]))

    # Number of evaluated clusters for this interval + 1
    j <- 1

    # Loop until all branches have been clustered
    repeat {

      # If all branches have been clustered, break the loops
      if (length(branches_todo) == 0) {

        break

      }

      # Initialise with first un-clustered branch
      # Consider adding a transmission to this branch
      branches_inc <- branches_todo[1]

      # Evaluate number of leaves for each host that would be assigned to the
      # newly added host
      counts_inc <- down_hosts[branches_inc, ]

      # Indicate which hosts would be affected by the newly added transmission
      hosts_inc <- which(counts_inc > 0)

      # Indicate which sampled hosts would be affected by the newly added transmission
      obs_inc <- unique_hosts[hosts_inc[which(hosts_type[hosts_inc] == 1)]]

      # Indicate that branch is now clustered
      branches_todo <- branches_todo[-1]

      # Repeatedly test if additional branches must be added to the cluster
      repeat {

        # If no branches remain unclustered, break the loop
        if (length(branches_todo) == 0) {

          break

        }

        # If only a single branch remains unclustered
        if (length(branches_todo) == 1) {

          # If only one host would be affected in the current cluster
          if (length(hosts_inc) == 1) {

            # Determine which branches must be added
            to_add <- which(down_hosts[branches_todo, hosts_inc] > 0)

          } else {

            # Determine which branches must be added
            to_add <- which(max(down_hosts[branches_todo, hosts_inc]) > 0)

          }

        } else {

          # If only one host would be affected in the current cluster
          if (length(hosts_inc) == 1) {

            # Determine which branches must be added
            to_add <- which(down_hosts[branches_todo, hosts_inc] > 0)

          } else {

            # Determine which branches must be added
            to_add <- which(apply(down_hosts[branches_todo, hosts_inc], 1, max) > 0)

          }

        }

        # If no branches need to be added to the cluster, break the loop
        if (length(to_add) == 0) {

          break

        }

        # Update branches in cluster
        branches_inc <- c(branches_inc, branches_todo[to_add])

        # Update number of affected leaves for each host
        counts_inc <- apply(down_hosts[branches_inc, , drop = F], 2, sum)

        # Update vector of affected hosts
        hosts_inc <- which(counts_inc > 0)

        # Update vector of affected sampled hosts
        obs_inc <- unique_hosts[hosts_inc[which(hosts_type[hosts_inc] == 1)]]

        # Update vector of unclustered branches
        branches_todo <- branches_todo[-to_add]
        #        branches_todo <- branches_todo[which(!(branches_todo %in% branches_inc))]

      }

      # If all leaves for each host are either included or excluded
      # and at most one observed host is included
      if (all(counts_inc[hosts_inc] == hosts_count[hosts_inc]) &
          length(obs_inc) < 2) {

        # Store cluster information
        # branches: branches included in cluster
        # counts: downstream leaf counts for each host
        # hosts: affected hosts
        # bin: unique numerical identifier for included hosts
        # obs: affected sampled hosts
        clusters_int[[j]] <- list(branches = branches_inc,
                                  counts = counts_inc,
                                  hosts = hosts_inc,
                                  bin = sum(2 ^ (hosts_inc - 1)),
                                  obs = obs_inc)

        # Increment cluster count
        j <- j + 1

      }

    }

    # Store all clusters for current interval
    clusters[[i]] <- clusters_int

    # If there is at least one valid cluster within current interval
    if (length(clusters[[i]]) > 0) {

      # Loop through possible cluster sizes
      for (csize in 1:length(clusters[[i]])) {

        # Possible cluster combinations for given combination size
        combs <- combn(length(clusters[[i]]), csize)

        # Loop through possible cluster combinations for current size
        for (k in 1:dim(combs)[2]) {

          # Resulting branches, numerical ID, and sampled hosts from combination
          temp_branches <- c()
          temp_bin <- 0
          temp_obs <- c()

          # Add relevant clusters
          for (j in combs[, k]) {

            # Update branches, numerical ID, and sampled hosts
            temp_bin <- temp_bin + clusters[[i]][[j]]$bin
            temp_branches <- c(temp_branches, clusters[[i]][[j]]$branches)
            temp_obs <- c(temp_obs, clusters[[i]][[j]]$obs)

          }

          # Number of transmitted lineages
          lb <- length(temp_branches)

          tr_lin[[lb]][[length(tr_lin[[lb]]) + 1]] <- list(interval = i,
                                                           clusts = combs[, k],
                                                           branches = temp_branches,
                                                           bin = temp_bin,
                                                           start = start,
                                                           end = end,
                                                           length = end - start)


          # Determine if resulting transmission already has weight
          tr_ind <- which(poss_tr[, 1] == temp_bin)

          if (length(tr_ind) == 0) {

            # If the transmission is novel, add as a possibility
            poss_tr <- rbind(poss_tr, c(temp_bin, (end - start)))

          } else {

            # If transmission currently has weight, update the possible timespan
            poss_tr[tr_ind, 2] <- poss_tr[tr_ind, 2] + (end - start)

          }

        }

      }

    }

  }

  len_lin <- numeric(L)

  for (l in 1:L) {

    temp_len <- sapply(tr_lin[[l]], function(x) x[["length"]])

    if (length(temp_len) > 0) {

      len_lin[l] <- sum(temp_len)

    } else {

      len_lin[l] <- 0

    }

  }

  poss_lin <- which(len_lin > 0)

  if (pen.len == 0) {

    if (pen.prob == 0) {

      prob_lin <- numeric(L)
      prob_lin[poss_lin] <- 1 / length(poss_lin)

    } else {

      prob_lin <- dnbinom(0:(L - 1), prob = pen.prob, size = pen.size) * ((1:L) %in% poss_lin)
      prob_lin <- prob_lin / sum(prob_lin)

    }

  } else {

    if (pen.prob == 0) {

      prob_lin <- len_lin / sum(len_lin)

    } else {

      prob_lin <- dnbinom(0:(L - 1), prob = pen.prob, size = pen.size) * len_lin
      prob_lin <- prob_lin / sum(prob_lin)

    }

  }

  sam_lin <- sample(1:L, size = 1, prob = prob_lin)

  prop_density <- prop_density + log(prob_lin[sam_lin])

  temp_len <- sapply(tr_lin[[sam_lin]], function(x) x[["length"]])

  prob_clust <- temp_len / sum(temp_len)

  sam_clust <- sample(1:length(prob_clust), size = 1, prob = prob_clust)

  prop_density <- prop_density + log(prob_clust[sam_clust])

  # Sample transmission time
  sam_time <- runif(1, min = tr_lin[[sam_lin]][[sam_clust]]$start,
                    max = tr_lin[[sam_lin]][[sam_clust]]$end)

  prop_density <- prop_density - log(tr_lin[[sam_lin]][[sam_clust]]$length)

  # Update ctree
  old_host <- host
  new_host <- host + 1

  # Update hosts
  if (new_host <= max(ctree[, 4])) {

    ctree[which(ctree[, 4] >= new_host), 4] <- ctree[which(ctree[, 4] >= new_host), 4] + 1

  }

  temp_branches <- tr_lin[[sam_lin]][[sam_clust]]$branches

  rhl <- rows_host[temp_branches]

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

  # Insert rows
  rhl <- sort(rows_host[temp_branches], decreasing = T)

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

  # Affected hosts in new ctree
  prop_hosts <- c(ctree[tar_row, 4], ctree[ctree[tar_row, 2], 4])

  # Sample rows and host
  obs_idx <- which(ctree[, 2] == 0 & ctree[, 3] == 0)
  obs_host <- ctree[obs_idx, 4]

  # Transmission rows and hosts
  tr_idx <- which(ctree[, 2] > 0 & ctree[, 3] == 0)
  tr_host1 <- ctree[tr_idx, 4]
  tr_host2 <- ctree[ctree[tr_idx, 2], 4]

  # Exclude transmissions between sampled hosts and root transmission
  tr_idx <- tr_idx[!((tr_host1 %in% obs_host) & (tr_host2 %in% obs_host)) & tr_host1 != 0]

  # Hosts that could be removed
  tr_host2 <- ctree[ctree[tr_idx, 2], 4]

  # Log density of reverse move
  rev_density <- -log(length(unique(tr_host2)))

  new_ctree <- list(ctree = ctree, nam = nam)
  class(new_ctree) <- 'ctree'

  return(list(ctree = new_ctree, prop_density = prop_density, rev_density = rev_density, is_possible = 1,
              curr_hosts = curr_hosts, prop_hosts = prop_hosts))

}


#' Proposal to remove a transmission
#'
#' @param ctree Current coloured tree.
remove_transmission <- function(ctree, pen.prob = 0, pen.size = 1, pen.len = 1) {

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

  if (length(tr_idx) == 0) {

    new_ctree <- list(ctree = ctree, nam = nam)
    class(new_ctree) <- 'ctree'

    return(list(ctree = new_ctree, prop_density = 1, rev_density = 1, is_possible = 0,
                curr_hosts = c(), prop_hosts = c()))

  }

  tr_host1 <- ctree[tr_idx, 4]
  tr_host2 <- ctree[ctree[tr_idx, 2], 4]

  if (length(unique(tr_host2)) > 1) {

    sam_host2 <- sample(unique(tr_host2), size = 1)

  } else {

    sam_host2 <- tr_host2[1]

  }

  sam_host1 <- tr_host1[which(tr_host2 == sam_host2)][1]

  sam_tr_ex <- which(tr_host2 == sam_host2)
  sam_ct_ex <- tr_idx[sam_tr_ex]

  sam_time <- ctree[sam_ct_ex[1], 1]

  curr_hosts <- c(sam_host1, sam_host2)

  prop_density <- -log(length(unique(tr_host2)))

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

  if (sam_host1 < sam_host2) {

    host <- sam_host1

  } else {

    host <- sam_host1 - 1

  }


  # Rows for transmission or observations in host
  leaves <- which(ctree[, 3] == 0 & ctree[, 4] == host)
  L <- length(leaves)

  # Event type (observation = 1, transmission = 2) for each leaf
  type <- 1 * (ctree[leaves, 2] == 0 & ctree[leaves, 3] == 0) +
    2 * (ctree[leaves, 2] > 0 & ctree[leaves, 3] == 0)

  # Corresponding host to type (sampled or infected individual)
  host2 <- rep(host, length(leaves))
  host2[which(type == 2)] <- ctree[ctree[leaves[which(type == 2)], 2], 4]

  # Vector of hosts at leaves
  unique_hosts <- unique(host2)

  # Number of hosts at leaves
  n_hosts2 <- length(unique_hosts)

  # Number of leaves for each host
  hosts_count <- numeric(n_hosts2)

  # Leaf type for each host
  hosts_type <- numeric(n_hosts2)

  # Loop through possible hosts
  for (i in 1:n_hosts2) {

    # Update host counts and type
    hosts_count[i] <- length(which(host2 == unique_hosts[i]))
    hosts_type[i] <- type[which(host2 == unique_hosts[i])][1]

  }

  # All rows assigned to host
  rows_host <- which(ctree[, 4] == host)
  R <- length(rows_host)

  # Indicate which leaves are downstream for each branch
  down_leaves <- matrix(0, R, L)

  # Indicate which hosts are downstream for each branch
  # and count the number of relevant leaves
  down_hosts <- matrix(0, R, n_hosts2)

  # Time interval for each branch
  interval <- matrix(0, R, 2)


  # Loop through leaves
  for (l in 1:L) {

    # Initialise at row corresponding to leaf
    r <- leaves[l]

    # Store host corresponding to leaf
    h <- host2[l]

    # Loop backwards until infection time
    repeat {

      # Indicate that current leaf is downstream of this branch
      down_leaves[which(rows_host == r), l] <- 1

      # Update host count for this branch
      down_hosts[which(rows_host == r), which(unique_hosts == h)] <-
        down_hosts[which(rows_host == r), which(unique_hosts == h)] + 1

      # Calculate parent branch
      anc <- which(ctree[, 2] == r | ctree[, 3] == r)

      # Calculate time interval of this branch
      interval[which(rows_host == r), ] <- c(ctree[anc, 1], ctree[r, 1])

      # Determine if we have reached the infection time
      if (ctree[r, 4] != ctree[anc, 4]) {

        # Break loop at infection time
        break

      } else {

        # Iterate backwards if not at infection time
        r <- anc

      }

    }

  }

  # Evaluate break points along mini-tree
  # Leaf times, coalescent times, infection time
  clust_intervals <- sort(interval)
  clust_intervals <- unique(clust_intervals)

  # Determine branch clusters. Branch combinations for which a transmission must
  # occur over all members or none
  clusters <- list()

  # Store possible transmissions according to number of transmitted lineages
  tr_lin <- list()

  for (l in 1:L) {

    tr_lin[[l]] <- list()

  }

  # Store possible transmissions, and length of time over which transmission
  # could occur in each case
  poss_tr <- array(dim = c(0, 2))

  for (i in 1:(length(clust_intervals) - 1)) {

    # Start and end point of interval
    start <- clust_intervals[i]
    end <- clust_intervals[i + 1]

    # Set of clusters for current time interval
    clusters_int <- list()

    # Determine which branches overlap with current interval
    branches_todo <- which(!(interval[, 2] <= clust_intervals[i] | interval[, 1] >= clust_intervals[i + 1]))

    # Number of evaluated clusters for this interval + 1
    j <- 1

    repeat {

      # If all branches have been clustered, break the loops
      if (length(branches_todo) == 0) {

        break

      }

      # Initialise with first un-clustered branch
      # Consider adding a transmission to this branch
      branches_inc <- branches_todo[1]

      # Evaluate number of leaves for each host that would be assigned to the
      # newly added host
      counts_inc <- down_hosts[branches_inc, ]

      # Indicate which hosts would be affected by the newly added transmission
      hosts_inc <- which(counts_inc > 0)

      # Indicate which sampled hosts would be affected by the newly added transmission
      obs_inc <- unique_hosts[hosts_inc[which(hosts_type[hosts_inc] == 1)]]

      # Indicate that branch is now clustered
      branches_todo <- branches_todo[-1]

      repeat {

        # If no branches remain unclustered, break the loop
        if (length(branches_todo) == 0) {

          break

        }

        # If only a single branch remains unclustered
        if (length(branches_todo) == 1) {

          # If only one host would be affected in the current cluster
          if (length(hosts_inc) == 1) {

            # Determine which branches must be added
            to_add <- which(down_hosts[branches_todo, hosts_inc] > 0)

          } else {

            # Determine which branches must be added
            to_add <- which(max(down_hosts[branches_todo, hosts_inc]) > 0)

          }

        } else {

          # If only one host would be affected in the current cluster
          if (length(hosts_inc) == 1) {

            # Determine which branches must be added
            to_add <- which(down_hosts[branches_todo, hosts_inc] > 0)

          } else {

            # Determine which branches must be added
            to_add <- which(apply(down_hosts[branches_todo, hosts_inc], 1, max) > 0)

          }

        }

        # If no branches need to be added to the cluster, break the loop
        if (length(to_add) == 0) {

          break

        }

        # Update branches in cluster
        branches_inc <- c(branches_inc, branches_todo[to_add])

        # Update number of affected leaves for each host
        counts_inc <- apply(down_hosts[branches_inc, , drop = F], 2, sum)

        # Update vector of affected hosts
        hosts_inc <- which(counts_inc > 0)

        # Update vector of affected sampled hosts
        obs_inc <- unique_hosts[hosts_inc[which(hosts_type[hosts_inc] == 1)]]

        # Update vector of unclustered branches
        branches_todo <- branches_todo[-to_add]
        #        branches_todo <- branches_todo[which(!(branches_todo %in% branches_inc))]


      }

      # If all leaves for each host are either included or excluded
      # and at most one observed host is included
      if (all(counts_inc[hosts_inc] == hosts_count[hosts_inc]) &
          length(obs_inc) < 2) {

        # Store cluster information
        # branches: branches included in cluster
        # counts: downstream leaf counts for each host
        # hosts: affected hosts
        # bin: unique numerical identifier for included hosts
        # obs: affected sampled hosts
        clusters_int[[j]] <- list(branches = branches_inc,
                                  counts = counts_inc,
                                  hosts = hosts_inc,
                                  bin = sum(2 ^ (hosts_inc - 1)),
                                  obs = obs_inc)

        # Increment cluster count
        j <- j + 1

      }

    }

    # Store all clusters for current interval
    clusters[[i]] <- clusters_int

    # If there is at least one valid cluster within current interval
    if (length(clusters[[i]]) > 0) {

      # Loop through possible cluster sizes
      for (csize in 1:length(clusters[[i]])) {

        # Possible cluster combinations for given combination size
        combs <- combn(length(clusters[[i]]), csize)

        # Loop through possible cluster combinations for current size
        for (k in 1:dim(combs)[2]) {

          # Resulting branches, numerical ID, and sampled hosts from combination
          temp_branches <- c()
          temp_bin <- 0
          temp_obs <- c()

          # Add relevant clusters
          for (j in combs[, k]) {

            # Update branches, numerical ID, and sampled hosts
            temp_bin <- temp_bin + clusters[[i]][[j]]$bin
            temp_branches <- c(temp_branches, clusters[[i]][[j]]$branches)
            temp_obs <- c(temp_obs, clusters[[i]][[j]]$obs)

          }

          # Number of transmitted lineages
          lb <- length(temp_branches)

          tr_lin[[lb]][[length(tr_lin[[lb]]) + 1]] <- list(interval = i,
                                                           clusts = combs[, k],
                                                           branches = temp_branches,
                                                           bin = temp_bin,
                                                           start = start,
                                                           end = end,
                                                           length = end - start)


          # Determine if resulting transmission already has weight
          tr_ind <- which(poss_tr[, 1] == temp_bin)

          if (length(tr_ind) == 0) {

            # If the transmission is novel, add as a possibility
            poss_tr <- rbind(poss_tr, c(temp_bin, (end - start)))

          } else {

            # If transmission currently has weight, update the possible timespan
            poss_tr[tr_ind, 2] <- poss_tr[tr_ind, 2] + (end - start)

          }

        }

      }

    }

  }

  len_lin <- numeric(L)

  for (l in 1:L) {

    temp_len <- sapply(tr_lin[[l]], function(x) x[["length"]])

    if (length(temp_len) > 0) {

      len_lin[l] <- sum(temp_len)

    } else {

      len_lin[l] <- 0

    }

  }

  poss_lin <- which(len_lin > 0)

  if (pen.len == 0) {

    if (pen.prob == 0) {

      prob_lin <- numeric(L)
      prob_lin[poss_lin] <- 1 / length(poss_lin)

    } else {

      prob_lin <- dnbinom(0:(L - 1), prob = pen.prob, size = pen.size) * ((1:L) %in% poss_lin)
      prob_lin <- prob_lin / sum(prob_lin)

    }

  } else {

    if (pen.prob == 0) {

      prob_lin <- len_lin / sum(len_lin)

    } else {

      prob_lin <- dnbinom(0:(L - 1), prob = pen.prob, size = pen.size) * len_lin
      prob_lin <- prob_lin / sum(prob_lin)

    }

  }

  sam_lin <- length(sam_ct_ex)

  rev_density <- rev_density + log(prob_lin[sam_lin])

  temp_len <- sapply(tr_lin[[sam_lin]], function(x) x[["length"]])

  prob_clust <- temp_len / sum(temp_len)

  temp_start <- sapply(tr_lin[[sam_lin]], function(x) x[["start"]])
  temp_end <- sapply(tr_lin[[sam_lin]], function(x) x[["end"]])
  temp_bin <- sapply(tr_lin[[sam_lin]], function(x) x[["bin"]])

  orig_bin <- sum(as.numeric(unique_hosts %in% unique(host2[orig_host[leaves] == sam_host2])) * 2 ^ (0:(length(unique_hosts) - 1)))

  sam_clust <- which(temp_start < sam_time & temp_end > sam_time & temp_bin == orig_bin)

  rev_density <- rev_density + log(prob_clust[sam_clust])

  rev_density <- rev_density - log(tr_lin[[sam_lin]][[sam_clust]]$length)

  # Order hosts
  ctree <- order_hosts(ctree)

  prop_hosts <- ctree[rows_host[1], 4]

  new_ctree <- list(ctree = ctree, nam = nam)
  class(new_ctree) <- 'ctree'

  return(list(ctree = new_ctree, prop_density = prop_density, rev_density = rev_density, is_possible = 1,
              curr_hosts = curr_hosts, prop_hosts = prop_hosts))

}

