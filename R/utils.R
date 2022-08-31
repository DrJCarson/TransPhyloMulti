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
