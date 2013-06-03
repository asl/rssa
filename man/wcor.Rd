\name{wcor}
\alias{wcor}
\alias{wcor.ssa}
\alias{wcor.default}
\alias{plot.wcor.matrix}

\title{Calculate the W-correlation matrix}
\description{
  Function calculates the W-correlation matrix for the series.
}
\usage{
\method{wcor}{ssa}(x, groups, \dots, cache = TRUE)
\method{wcor}{default}(x, L = (N + 1) \%/\% 2, \dots)
\method{plot}{wcor.matrix}(x,
     grid = c(),
     \dots,
     cuts = 20,
     zlim = range(abs(x), 0, 1))
}

\arguments{
  \item{x}{the input object. This might be ssa object for \emph{ssa}
    method, or just a matrix with elementary series in columns for
  \emph{default} implementation.}
  \item{L}{window length.}
  \item{groups}{list of numeric vectors, indices of elementary series
    used for reconstruction.}
  \item{\dots}{further arguments passed to \code{reconstruct} routine
  for \code{wcor} or to \code{plot} for \code{plot.wcor.matrix}}
  \item{cache}{logical, if 'TRUE' then intermediate results will be
    cached in 'ssa' object.}
  \item{grid}{numeric vector, indices of matrix blocks (groups)
    which will be separated by grid line.
    Lines will be drawn on the left of and under noted blocks.
    Also this argument can be list of two numeric vectors
    with names 'x' and 'y', for control vertical and horizontal
    grid lines separately.}
  \item{cuts}{integer, the number of levels the range of W-cor
    values will be divided into.}
  \item{zlim}{range of displaying W-cor values.
    Values outside the range will be absolutely not shown.}
}
\details{
  W-correlation matrix is a standard way of checking for weak
  separability between the elementary time series. In particular, the
  strongly correlated elementary series should be placed into the same
  group.The function calculates such a matrix either directly from 'ssa'
  object or from the matrix of elementary series.

  Additional (non-standard) graphical parameters which can be transfered via \dots:
  \describe{
    \item{useRaster}{logical. For wcor-plot, indicating whether raster representations
      should be used. 'TRUE' by default.}
  }
}

\value{
  Object of type 'wcor.matrix'.
}

\references{
  Golyandina, N., Nekrutkin, V. and Zhigljavsky, A. (2001): \emph{Analysis of
    Time Series Structure: SSA and related techniques.} Chapman and
  Hall/CRC. ISBN 1584881941
}
\seealso{
  \code{\link[Rssa:reconstruct.ssa]{reconstruct}}
}

\examples{
# Decompose co2 series with default parameters
s <- ssa(co2)
# Calculate the w-correlation matrix between first 20 series
# for a guess for grouping
w <- wcor(s, groups = 1:20)
plot(w, grid = c(2,4, 5,7))
# Calculate the w-correlation matrix for the chosen groups
# to check separability
w <- wcor(s, groups = list(c(1,4), c(2,3), c(5,6)))
}

