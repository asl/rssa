\name{lrr}
\alias{lrr}
\alias{lrr.default}
\alias{lrr.ssa}
\alias{lrr.1d-ssa}
\title{Calculate the Linear Recurrence Relation}

\description{
  Calculates the Linear Recurrence Relation given the 'ssa' object.
}

\usage{
\method{lrr}{ssa}(x, group)
}

\arguments{
  \item{x}{SSA object holding the decomposition}
  \item{group}{indices of eigenvectors used to derive the LRR}
}

\details{
  Produces the linear recurrence relation from the series. Check
  \code{\link{basis2lrr}} for calculation details.
}

\value{
  Object of class 'lrr'
}

\seealso{
  \code{\link[Rssa:basis2lrr]{basis2lrr}}
  \code{\link[Rssa:plot.lrr]{plot.lrr}}
  \code{\link[Rssa:roots.lrr]{roots.lrr}}
}

\examples{
# Decompose 'co2' series with default parameters
s <- new.ssa(co2, L = 24)
# Calculate the LRR out of first 3 eigentriples
lrr(s, group = c(1:3))
}
