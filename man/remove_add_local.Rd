% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/proposals.R
\name{remove_add_local}
\alias{remove_add_local}
\title{Proposal to remove a transmission, and add a new one with the same infector}
\usage{
remove_add_local(ctree, bn_weight = 0.1, delta = 1)
}
\arguments{
\item{ctree}{Current coloured tree.}

\item{bn_weight}{Penalisation term for transmitting multiple lineages.}

\item{delta}{Standard deviation of random walk proposal for root update.}
}
\description{
Proposal to remove a transmission, and add a new one with the same infector
}
