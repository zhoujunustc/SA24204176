% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Mirror.R
\name{GM}
\alias{GM}
\title{Gauss Mirrors}
\usage{
GM(y, X, q)
}
\arguments{
\item{y}{vector of the response variable}

\item{X}{the design matrix}

\item{q}{the designated FDR level to be controlled at}
}
\value{
a list with the selected features and the mirror statistic of all features.
}
\description{
Use Gaussian Mirrors to select features while controlling FDR for linear regression model.
}
\examples{
\dontrun{
n <- 100; p <- 5
beta <- c(0,1,1,2,0)
X <- matrix(rnorm(n*p),n,p)
y <- X \%*\% beta + rnorm(n)
GM(y,X,0.1)$feature
}
}
