\name{reconstruct}
\alias{reconstruct}
\alias{reconstruct.ssa}
\title{Perform a series reconstruction}

\description{
  Reconstruct the series given the SSA decomposition and the desired
  grouping of the elementary subseries.
}

\usage{
\method{reconstruct}{ssa}(this, groups, \dots, cache = TRUE)
}

\arguments{
  \item{this}{SSA object}
  \item{groups}{list of numeric vectors, indices of elementary subseries
    used for reconstruction.}
  \item{\dots}{further arguments passed to routines (e.g. to
  \code{decompose} routine if the continuation is desired).}
  \item{cache}{logical, if 'TRUE' then intermediate results will be
    cached in the SSA object.}
}

\value{
  List of vectors of reconstructed series. Elements of the list are
  named 'F1'', 'F2', and so on.
}

\examples{
# Decompose 'co2' series with default parameters
s <- new.ssa(co2)
# Reconstruct the series, grouping elementary subseries 2, 3 and 4, 5.
r <- reconstruct(s, groups = list(1, c(2,3), c(4,5)))
\dontrun{plot(r$F1, col = "black")}
\dontrun{lines(r$F1 + r$F2, col = "red")}
\dontrun{lines(r$F1 + r$F2 + r$F3, col = "blue")}
# 'groups' argument might contain duplicate entires as well
r <- reconstruct(s, groups = list(1, 1:3, 1:5))
\dontrun{plot(r$F1, col = "black")}
\dontrun{lines(r$F2, col = "red")}
\dontrun{lines(r$F3, col = "blue")}
}
