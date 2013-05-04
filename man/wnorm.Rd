\name{wnorm}
\alias{wnorm}
\alias{wnorm.1d.ssa}
\alias{wnorm.toeplitz.ssa}
\alias{wnorm.default}

\title{Calculate weighted norm of series}
\description{
  Function calculates the W-norm for input series or for series stored in input ssa obect.
}
\usage{
\method{wnorm}{1d.ssa}(x, ...)
\method{wnorm}{toeplitz.ssa}(x, ...)
\method{wnorm}{default}(x, L = (N + 1) \%/\% 2, ...)
}

\arguments{
  \item{x}{the input object. This might be ssa object for \emph{ssa}
    method, or just a series.}
  \item{L}{window length.}
  \item{\dots}{arguments to be passed to methods.}
}

\details{
  \code{L}-weighted norm of series is Frobenius norm of its \code{L}-trajectory matrix. So,
  if \code{x} is vector (series), the result of \code{wnorm(x, L)} is equal to \code{sqrt(sum(hankel(x, L)^2)}, but
  would be evaluated much more efficiently.
  For 1d SSA and Toeplitz SSA \code{wnorm(x)} calculates weighted norm for stored original input series and stored window length.

  For 2d SSA this method is not implemented yet.
}

\references{
  Golyandina, N., Nekrutkin, V. and Zhigljavsky, A. (2001): \emph{Analysis of
    Time Series Structure: SSA and related techniques.} Chapman and
  Hall/CRC. ISBN 1584881941
}

\seealso{
  \code{\link[Rssa:hankel]{hankel}}
  \code{\link[Rssa:wcor.ssa]{wcor}}
}

\examples{
wnorm(co2, 20)
ss <- ssa(co2)
wnorm(ss)
}
