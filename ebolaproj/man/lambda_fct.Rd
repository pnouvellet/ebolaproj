% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lambda_fct.R
\name{lambda_fct}
\alias{lambda_fct}
\title{'force of infection'}
\usage{
lambda_fct(param, I, N_l, ws, SItrunc)
}
\arguments{
\item{I}{matrix of observed incidence, same dimension as lambda}

\item{N_l}{integer of  numbers of locations}

\item{ws}{vector reversed serial interval distribution (output from SI_gamma_dist_EpiEstim reverse in the MCMC function) )}

\item{SItrunc}{integer, threshold of serial interval distribution (see SI_gamma_dist_EpiEstim)}

\item{param:}{'force of infection' matrix (incidence weighted by serial interval),
column number of days, row: number of locations}
}
\description{
return incidence weighted by serial interval for the time window of interest 
internal to the MCMC
}
\details{
lambda incidence weighted by serial interval
}
