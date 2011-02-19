\name{calc.v}
\alias{calc.v}
\alias{calc.v.ssa}
\alias{calc.v.1d-ssa}
\alias{calc.v.2d-ssa}

\title{Calculate Factor Vector(s)}

\description{
  Generic function for the factor vector calculation given the SSA
  decomposition.
}

\usage{
\method{calc.v}{ssa}(this, idx, env = .GlobalEnv, \dots)
}

\arguments{
  \item{this}{SSA object holding the decomposition.}
  \item{idx}{indices of the factor vectors to compute.}
  \item{env}{temporary environment used for precalculation.}
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

  Note that for several SVD methods (`eigen' and `svd' currently) the
  matrix X is not stored inside the `ssa' object and must be calculated
  out of series which might be pretty time-consuming, especially if
  multiple calls to the routine are desired. User might provide a
  temporary garbage-collected environment, which the routine will use to
  store different values which will be shared between the calls. It is
  the user responsibility to ensure that the environment is valid
  (e.g. it was populated via the call of calc.v with the same `ssa'
  object as an argument).
  
  Temporary environment is not required for fast FFT-based SVD method
  (`propack' and `nutrlan').
}

\value{
  A numeric vector of suitable length (usually depends on SSA method
  and window length).
}

\examples{
# Decompose 'co2' series with default parameters
s <- new.ssa(co2)
# Calculate the 5th factor vector
v <- calc.v(s, 5)
}

\keyword{algebra}

