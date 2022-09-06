#' Extracts phylogenetic tree from a combined phylogenetic/transmission tree
#' @param ctree Combined tree
#'
#' @export
extractPTreeM <- function(ctree)  {

  tree <- ctree$ctree
  nam <- ctree$nam

  n <- sum(tree[ ,2] + tree[ ,3] == 0)

  tra <- n + 1
  while (tra <= nrow(tree))  { #Removing transmission events from the tree one by one

    if (tree[tra,3] != 0)  {

      tra <- tra + 1
      next

    }

    t <- tree[ ,2:3]
    f <- which(t == tra)
    t[f] <- tree[tra,2]
    f <- which( t > tra )
    t[f] <- t[f]-1
    tree[ ,2:3] <- t
    tree <- tree[-tra, , drop=FALSE]
    tra <- n + 1

  }

  ptree <- tree[, 1:(ncol(tree) - 1), drop=FALSE]

  l <- list(ptree = ptree, host = tree[1:n, 4], nam = nam)
  class(l) <- 'ptreem'

  return(l)

}


#' Extracts transmission tree from a combined phylogenetic/transmission tree
#' @param ctree Combined tree
#' @export
extractTTreeM <- function(ctree)  {

  nam <- ctree$nam
  ctree <- ctree$ctree

  host <- ctree[ ,4]

  tr_rows <- which(ctree[, 2] != 0 & ctree[, 3] == 0)
  tr_times <- ctree[tr_rows, 1]
  infectors <- ctree[tr_rows, 4]
  infecteds <- ctree[ctree[tr_rows, 2], 4]

  obs_rows <- which(ctree[, 2] == 0 & ctree[, 3] == 0)
  obs_host <- ctree[obs_rows, 4]
  nsam <- length(obs_rows)

  ttree <- array(max(host) * 3, dim = c(max(host), 3))
  obs <- array(nsam * 2, dim = c(nsam, 2))

  for (i in 1:max(host)) {

    h_ind <- which(infecteds == i)[1]

    ttree[i, ] <- c(tr_times[h_ind],
                    length(which(obs_host == i)),
                    infectors[h_ind])

  }

  for (i in 1:nsam) {

    obs[i, ] <- c(ctree[i, 1],
                  ctree[i, 4])

  }

  l <- list(ttree = ttree, obs = obs, nam = nam)
  class(l)<-'ttreem'

  return(l)

}
