#' Plot a transmission tree in an economic format
#' @param ttree Transmission tree
#' @param showLabels Boolean for whether or not to show the labels
#' @param maxTime Maximum value of time to show on x axis
#' @param cex Expansion factor
#' @export
plotTTree_summ = function(ttree, showLabels = T, maxTime = NA, cex = 1) {

  nam <- ttree$nam
  ttree <- ttree$ttree

  #Determine ys
  n <- nrow(ttree)

  ys <- rep(0, n)
  scale <- rep(1, n)

  todo <- c(which(ttree[, 3] == 0))
  while (length(todo) > 0) {

    f <- which(ttree[, 3] == todo[1])
    o <- rank(-ttree[f, 1], ties.method = "first")
    f[o] <- f

    for (i in f) {

      ys[i] <- ys[todo[1]] + scale[todo[1]] * which(f == i) / (length(f) + 1)
      scale[i] <- scale[todo[1]] / (length(f) + 1)
      todo <- c(todo, i)

    }

    todo <- todo[-1]

  }

  ys <- rank(ys)


  #Do the plot
  oldpar <- par('yaxt', 'bty')
  on.exit(par(oldpar))
  par(yaxt='n', bty='n')

  mi <- min(ttree[which(!is.na(ttree[, 1])), 1])
  ma <- max(ttree[which(!is.na(ttree[, 1])), 1])

  if (!is.na(maxTime)) {

    ma <- maxTime

  }

  plot(c(), c(), xlim = c(mi, ma), ylim=c(0, n + 1), xlab = '', ylab = '')

  for (i in 1:n) {

    if (ttree[i, 3] != 0) {

      arrows(ttree[ttree[i, 3], 1], ys[ttree[i, 3]], ttree[i, 1], ys[i], length = 0)

    }

    if (showLabels && ttree[i, 2] > 0) {

      text(ttree[i, 1], ys[i], i, pos = 4, cex = cex)

    }

  }

  for (i in 1:n) {

    points(ttree[i, 1], ys[i], pch = 21, bg = ifelse(ttree[i, 2] == 0, 'white', 'black'), cex = cex)

  }

  return(invisible(ttree))

}


#' Plot a transmission tree in a detailed format
#'
#' @param ttree Transmission tree
#' @param w.shape Shape parameter of the Gamma probability density function representing the generation time
#' @param w.scale Scale parameter of the Gamma probability density function representing the generation time
#' @param showLabels Whether or not to show the labels
#' @param maxTime Maximum value of time to show on x axis
#' @param cex Expansion factor
#' @export
plotTTree_det = function(ttree, w.shape, w.scale, showLabels = TRUE, maxTime = NA, cex=1) {

  nam <- ttree$nam
  obs <- ttree$obs
  ttree <- ttree$ttree

  n <- nrow(ttree)

  #Determine ys
  ys <- rep(0, n)
  scale <- rep(1, n)

  todo <- c(which(ttree[, 3] == 0))
  while (length(todo) > 0) {

    f <- which(ttree[, 3] == todo[1])
    o <- rank(-ttree[f, 1])
    f[o] <- f

    for (i in f) {

      ys[i] <- ys[todo[1]] + scale[todo[1]] * which(f == i) / (length(f) + 1)
      scale[i] <- scale[todo[1]] / (length(f) + 1)
      todo <- c(todo, i)

    }

    todo <- todo[-1]

  }

  ys <- rank(ys)

  oldpar <- par('yaxt', 'bty')
  on.exit(par(oldpar))
  par(yaxt='n', bty='n')

  mi <- min(ttree[, 1])

  if (is.na(maxTime)) {

    ma <- max(obs[, 1])

  }  else {

    ma <- maxTime

  }

  xstep <- (ma - mi) / 2000

  plot(c(), c(), xlim = c(mi - (ma - mi) * 0.05, ma + (ma - mi) * 0.05), ylim=c(0, n + 1), xlab = '', ylab = '')

  maxcol <- max(dgamma(seq(0, ma - mi, xstep), shape = w.shape, scale = w.scale))

  for (i in 1:n) {

    as <- seq(ttree[i, 1], ma, xstep)
    bs <- rep(ys[i], length(as))
    cs <- abs((maxcol - dgamma(as - ttree[i, 1], shape = w.shape, scale = w.scale)) / maxcol)
    cs <- grDevices::gray(cs)

    segments(as, bs, x1 = as + xstep, col = cs)

    host_obs <- which(obs[, 2] == i)

    if (length(host_obs) > 0) {

      obs_times <- obs[host_obs, 1]

      points(obs_times, rep(ys[i], length(obs_times)), col = 'red')

    }

    if (showLabels) {

      text(ma + (ma - mi) * 0.05, ys[i], i, cex = cex)

    }

    if (ttree[i, 3] == 0) {

      next

    }

    arrows(ttree[i, 1], ys[ttree[i, 3]], ttree[i, 1], ys[i], length = 0.1)

  }

  return(invisible(ttree))

}


#' Plot MCMC traces
#'
#' @param record Output from inferTTree function
#' @param burnin Proportion of the MCMC output to be discarded as burnin
#' @return Returns invisibly the first parameter
#' @export
plotTracesM <- function(record, burnin = 0) {

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  par(mfrow = c(3, 2))

  record <- record[max(1, round(length(record) * burnin)):length(record)]

  plot(sapply(record, function(x) x$pTTree + x$pPTree), ylab = 'Posterior probability',
       xlab = 'MCMC iterations', type = 'l')

  plot(sapply(record, function(x) x$pi), ylab = 'Sampling proportion pi',
       xlab = 'MCMC iterations', type = 'l')

  plot(sapply(record, function(x) x$off.r), ylab = 'off.r',
       xlab = 'MCMC iterations', type = 'l')

  plot(sapply(record, function(x) x$off.p), ylab = 'off.p',
       xlab = 'MCMC iterations', type = 'l')

  plot(sapply(record, function(x) x$kappa), ylab = 'Within-host initial population kappa',
       xlab = 'MCMC iterations', type = 'l')

  plot(sapply(record, function(x) x$lambda), ylab = 'Within-host population growth rate lambda',
       xlab = 'MCMC iterations', type = 'l')

  return(invisible(record))

}
