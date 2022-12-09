% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/likelihoods.R
\name{log_likelihood_coalescence_linear}
\alias{log_likelihood_coalescence_linear}
\title{Likelihood evaluation for the linear growth model}
\usage{
log_likelihood_coalescence_linear(
  infected_time,
  final_time,
  start_time,
  lm_const,
  lm_rate,
  branch_combs,
  coalescence
)
}
\arguments{
\item{infected_time}{Time at which individual was infected.}

\item{final_time}{Lower bound of time period.}

\item{start_time}{Upper bound of time period.}

\item{lm_const}{Initial population vale.}

\item{lm_rate}{Growth rate.}

\item{branch_combs}{Number of possible coalescence possibilities.}

\item{coalescence}{Whether or not a coalescence occurs at final_time.}
}
\description{
Likelihood evaluation for the linear growth model
}