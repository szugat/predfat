% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/gradH.R
\name{gradH}
\alias{gradH}
\title{Derivation of quantile from hypoexponential distribution with respect to theta}
\usage{
gradH(x, y, theta, lambda, gradient, type)
}
\arguments{
\item{x}{values of influental variable for the link function}

\item{y}{values of depedent variable}

\item{theta}{numeric vector of length four with link function's parameters}

\item{lambda}{[\code{function(theta, x)}]\cr
link function for exponential distribution}

\item{gradient}{[\code{function(x, theta, ...)}]\cr
gradient of link function}

\item{type}{[\code{integer}]\cr
if link function is not given a collection of given link function is available, see \code{\link{linkfun}}}
}
\value{
derivation of quantile from hypoexponential distribution with respect to theta
}
\description{
Derivation of quantile from hypoexponential distribution with respect to theta
}

