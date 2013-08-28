\name{calc.v}
\alias{calc.v}
\alias{calc.v.ssa}
\alias{calc.v.1d.ssa}
\alias{calc.v.toeplitz.ssa}
\alias{calc.v.2d.ssa}
\alias{calc.v.mssa}

\title{Calculate Factor Vector(s)}

\description{
  Generic function for the factor vector calculation given the SSA
  decomposition.
}

\usage{
\method{calc.v}{1d.ssa}(x, idx, \dots)
\method{calc.v}{toeplitz.ssa}(x, idx, \dots)
\method{calc.v}{mssa}(x, idx, \dots)
\method{calc.v}{2d.ssa}(x, idx, \dots)
}

\arguments{
  \item{x}{SSA object holding the decomposition.}
  \item{idx}{indices of the factor vectors to compute.}
  \item{\dots}{additional arguments to 'calc.v'.}
}

\details{
  Factor vector is a column of the factor matrix V, which is calculated
  as follows:
  \deqn{%
    V = \Sigma^{-1} X^{T} U,
  }{%
    V = \Sigma^(-1) X^T U,
  }
  where X is a Hankel trajectory matrix, U is the matrix of eigenvectors
  and Sigma is a matrix of singular values.
}

\value{
  A numeric vector of suitable length (usually depends on SSA method
  and window length).
}

\examples{
# Decompose 'co2' series with default parameters
s <- ssa(co2)
# Calculate the 5th factor vector
v <- calc.v(s, 5)
}

\keyword{algebra}

