\name{lrf}
\alias{lrf}
\alias{lrf.default}
\alias{lrf.ssa}
\alias{lrf.1d-ssa}
\title{Calculate the Linear Recurrent Formula}

\description{
  Calculates the Linear Recurrent Formula given the 'ssa' object.
}

\usage{
\method{lrf}{ssa}(x, group)
}

\arguments{
  \item{x}{SSA object holding the decomposition}
  \item{group}{indices of eigenvectors used to derive the LRF}
}

\details{
  Produces the linear recurrent formula from the series. Check
  \code{\link{basis2lrf}} for calculation details.
}

\value{
  Object of class 'lrf'
}

\seealso{
  \code{\link[Rssa:basis2lrf]{basis2lrf}}
  \code{\link[Rssa:plot.lrf]{plot.lrf}}
  \code{\link[Rssa:roots.lrf]{roots.lrf}}
}

\examples{
# Decompose 'co2' series with default parameters
s <- new.ssa(co2, L = 24)
# Calculate the LRF out of first 3 eigentriples
lrf(s, group = c(1:3))
}
