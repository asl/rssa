\name{wnorm}
\alias{wnorm}
\alias{wnorm.1d.ssa}
\alias{wnorm.nd.ssa}
\alias{wnorm.toeplitz.ssa}
\alias{wnorm.mssa}
\alias{wnorm.default}
\alias{wnorm.complex}

\title{Calculate Weighted Norm of series}
\description{
  Function calculates the W-norm for input objects or for objects stored in input ssa obect.
}
\usage{
\method{wnorm}{1d.ssa}(x, ...)
\method{wnorm}{nd.ssa}(x, ...)
\method{wnorm}{toeplitz.ssa}(x, ...)
\method{wnorm}{mssa}(x, ...)
\method{wnorm}{default}(x, L = (N + 1) \%/\% 2, ...)
\method{wnorm}{complex}(x, L = (N + 1) \%/\% 2, ...)
}

\arguments{
  \item{x}{the input object. This might be ssa object for \emph{ssa}
    method, or just a series.}
  \item{L}{window length.}
  \item{\dots}{arguments to be passed to methods.}
}

\details{
  \code{L}-weighted norm of series is Frobenius norm of its
  \code{L}-trajectory matrix. So, if \code{x} is vector (series), the
  result of \code{wnorm(x, L)} is equal to \code{sqrt(sum(hankel(x,
  L)^2)}, but in fact is calculated much more efficiently. For 1d SSA and
  Toeplitz SSA \code{wnorm(x)} calculates weighted norm for stored
  original input series and stored window length.

  \code{L}-weighted norm of 2d array is Frobenius norm of its \code{L[1]
  * L[2]}-trajectory hankel-block-hankel matrix.  For 2d SSA this method
  calculates weighted norm for stored original input array and stored
  2d-window lengths.
}

\references{
  Golyandina, N., Nekrutkin, V. and Zhigljavsky, A. (2001): \emph{Analysis of
    Time Series Structure: SSA and related techniques.} Chapman and
  Hall/CRC. ISBN 1584881941
}

\seealso{
  \code{\link[Rssa:ssa-input]{ssa-input}},
  \code{\link[Rssa:hankel]{hankel}},
  \code{\link[Rssa:wcor]{wcor}}
}

\examples{
wnorm(co2, 20)
# Construct ssa-object for 'co2' with default parameters but don't decompose
ss <- ssa(co2, force.decompose = FALSE)
wnorm(ss)

# Artificial image for 2D SSA
\donttest{
mx <- outer(1:50, 1:50,
            function(i, j) sin(2*pi * i/17) * cos(2*pi * j/7) + exp(i/25 - j/20)) +
      rnorm(50^2, sd = 0.1)
# Construct ssa-object for 'mx' with default parameters but don't decompose
s <- ssa(mx, kind = "2d-ssa", force.decompose = FALSE)
wnorm(s)
}
}
