\name{wnorm}
\alias{wnorm}
\alias{wnorm.ssa}
\alias{wnorm.default}

\title{Calculate weighted norm of series}
\description{
  Function calculates the W-norm for input series or for series stored in input ssa obect.
}
\usage{
\method{wnorm}{ssa}(x, ...)
\method{wnorm}{default}(x, L = (N - 1) \%/\% 2, ...)
}

\arguments{
  \item{x}{the input object. This might be ssa object for \emph{ssa}
    method, or just a series.}
  \item{L}{window length.}
  \item{\dots}{arguments to be passed to methods.}
}

\details{
  If \code{x} is vector (series), the result of \code{wnorm(x, L)} is equal to \code{sqrt(sum(hankel(x, L)^2)}, but
  would be evaluated much more efficiently.
  For 1d SSA and Toeplitz SSA \code{wnorm(x)} just calls \emph{default} method for stored series and window length.

  For 2d SSA this method is not implemented yet.
}

\examples{
wnorm(co2, 20)
ss <- ssa(co2)
wnorm(ss)
}
