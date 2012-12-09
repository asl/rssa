\name{reconstruct}
\alias{reconstruct}
\alias{reconstruct.ssa}
\title{Perform a series reconstruction}

\description{
  Reconstruct the series given the SSA decomposition and the desired
  grouping of the elementary series.
}

\usage{
\method{reconstruct}{ssa}(x, groups, \dots, drop = FALSE, cache = TRUE)
}

\arguments{
  \item{x}{SSA object}
  \item{groups}{list of numeric vectors, indices of elementary series
    used for reconstruction.}
  \item{\dots}{further arguments passed to routines (e.g. to
  \code{decompose} routine if the continuation is desired).}
  \item{drop}{logical, if 'TRUE' then the attributes of the input series
  are not copied to the reconstructed ones.}
  \item{cache}{logical, if 'TRUE' then intermediate results will be
    cached in the SSA object.}
}

\value{
  List of vectors of reconstructed series. Elements of the list are
  named 'F1', 'F2', and so on.
}

\note{
  By default (argument \code{drop}) the routine copies all the
  attributes of the input series to the reconstructed series. Thus, the
  reconstruction result of 'ts' object is the 'ts' object with the same
  time scale as well.
}

\seealso{
  \code{\link{Rssa}} for an overview of the package, as well as,
  \code{\link[Rssa:ssa]{ssa}},
  \code{\link[Rssa:plot.reconstruction]{plot.reconstruction}},
}

\examples{
# Decompose 'co2' series with default parameters
s <- ssa(co2)
# Reconstruct the series, grouping elementary series.
r <- reconstruct(s, groups = list(c(1, 4), c(2,3), c(5, 6)))
plot(r)
# 'groups' argument might contain duplicate entries as well
r <- reconstruct(s, groups = list(1, 1:4, 1:6))
plot(r)
}
