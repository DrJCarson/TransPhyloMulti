#' Print function for ptreem objects
#' @param x Object of class ptreem, ie a phylogenetic tree
#' @param ... Additional parameters are passed on
#' @export
print.ptreem <- function(x, ...) {

  stopifnot(inherits(x, "ptreem"))

  cat( 'Phylogenetic tree\n')
  cat(sprintf('Number of sampled hosts = %d\n', max(x$host)))
  cat(sprintf('Number of samples = %d\n', length(x$nam)))

  invisible(x)

}

#' Print function for ttreem objects
#' @param x Object of class ttreem, ie a transmission tree
#' @param ... Additional parameters are passed on
#' @export
print.ttreem <- function(x, ...) {

  stopifnot(inherits(x, "ttreem"))

  cat('Transmission tree\n')
  cat(sprintf('Number of hosts = %d\n', nrow(x$ttree)))
  cat(sprintf('Number of sampled hosts = %d\n', length(unique(x$obs[, 2]))))
  cat(sprintf('Number of samples = %d\n', length(x$nam)))

  invisible(x)
}

#' Plotting for ptreem
#' @param x Object of class ptreem, ie  a phylogenetic tree
#' @param ... Additional parameters are passed on to ape::plot.phylo
#' @return Plot of ptreem
#' @export
plot.ptreem <- function(x, ...) {

  stopifnot(inherits(x, "ptreem"))

  phy <- phyloFromPTreeM(x)

  ape::plot.phylo(phy, ...)
  ape::axisPhylo(backward = F)

}

#' Plotting for ttreem
#' @param x Object of class ttreem, ie  a transmission tree
#' @param type Type of plot to display, can be 'detailed' or 'summarised' (default)
#' @param w.shape Shape parameter of the generation time, needed for detailed plot only
#' @param w.scale Scale parameter of the generation time, needed for detailed plot only
#' @param ... Additional parameters are passed on
#' @export
plot.ttreem <- function(x, type='summarised', w.shape=NA, w.scale=NA, ...) {

  stopifnot(inherits(x, "ttreem"))

  if (type == 'summarised') {

    plotTTree_summ(x,...)

  } else if (type == 'detailed') {

    if (is.na(w.shape) || is.na(w.scale)) {

      stop('You need to specify w.shape and w.scale to display this plot.')

    } else {

      plotTTree_det(x, w.shape = w.shape, w.scale = w.scale, ...)

    }

  }

}

#' Plotting for resTransPhylo
#' @param x Output from inferTTree
#' @param ... Additional parameters are passed on
#' @return Plot of TransPhylo results
#' @export
plot.resTransPhyloM <- function(x, ...) {

  stopifnot(inherits(x, "resTransPhyloM"))

  plotTracesM(x)

}

