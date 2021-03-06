% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mctnetwork.r
\name{mctnetwork_propogate_from_t}
\alias{mctnetwork_propogate_from_t}
\title{Compute backward and forward flow propgation of metacell probability from time t}
\usage{
mctnetwork_propogate_from_t(mct, t, mc_p)
}
\arguments{
\item{mct}{mct network object}

\item{t}{flows will be computed for the (time,time+1) interval}

\item{mc_p}{probabilities at time t}
}
\value{
a list with two elements: probs is a matrix of probabilities over metacells (rows)  and time (columns). step_m is a list of sparse matrices inferred flows between metacells per time.
}
\description{
Compute backward and forward flow propgation of metacell probability from time t
}
