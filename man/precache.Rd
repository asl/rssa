\name{precache}
\alias{precache}
\alias{precache.ssa}
\title{Calculates and caches elementary series inside SSA object}

\description{
  Calculates all the elementary series and saves inside SSA
  object. After this the grouping procedure can be performed much
  faster.
}
\usage{
\method{precache}{ssa}(this, n, \dots)
}

\arguments{
  \item{this}{SSA object.}
  \item{n}{integer, number of series to calculate and save.}
  \item{\dots}{further arguments passed to the reconstruction routines.}
}

\note{
  In most cases it is not necessary to call this routine directly. By
  default functions from the package collect all elementary
  subseries they encounter during the calculations.
}

\seealso{
  \code{\link[Rssa:reconstruct.ssa]{reconstruct}}
}

\examples{
# Decompose 'co2' series with default parameters
s <- new.ssa(co2)
summary(s)
# Precache the stuff
precache(s)
summary(s)
} 
