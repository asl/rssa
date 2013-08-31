\name{cadzow}
\alias{cadzow}
\alias{cadzow.ssa}
\alias{cadzow.1d.ssa}

\title{
  Cadzow Iterations
}
\description{
  Perform the finite rank approximation of the series via Cadzow iterations
}

\usage{
\method{cadzow}{ssa}(x, rank, eps = 1e-6, numiter = 0, \dots, cache = TRUE)
\method{cadzow}{1d.ssa}(x, rank, correct = TRUE, eps = 1e-6, numiter = 0, \dots, cache = TRUE)
}

\arguments{
  \item{x}{input SSA object}
  \item{rank}{desired rank of approximation}
  \item{correct}{logical, if 'TRUE' then additional correction per
  (Gillard et al, 2013) is performed}
  \item{eps}{tolerance value used for convergence criteria}
  \item{numiter}{number of iterations to perform, if zero then
    iterations are performed until the convergence}
  \item{\dots}{further arguments passed to \code{reconstruct}}
  \item{cache}{logical, if 'TRUE' then intermediate results will be
    cached in the SSA object.}
}

\details{
  Cadzow iterations aim to solve the problem of the approximation of the
  input series by a series of finite rank. The idea of the algorithm is
  quite simple: the sequence of projections onto the space of Hankel
  matricies of finite rank is performed until the convergence.

  Unfortunately, being simple, the method often yields the solution which is
  far away from the optimum.
}

\seealso{
  \code{\link{Rssa}} for an overview of the package, as well as,
  \code{\link[Rssa:reconstruct]{reconstruct}}
}

\references{
  Cadzow J. A. (1988) Signal enhancement a composite property mapping algorithm,
  IEEE Transactions on Acoustics, Speech, and Signal Processing, 36,
  49-62.

  Gillard, J. and Zhigljavsky, A. (2013) Stochastic optimization
  algorithms for Hankel structured low-rank approximation. Unpublished
  Manuscript. Cardiff School of Mathematics. Cardiff.
}
  
\examples{
# Decompose 'co2' series with default parameters
s <- ssa(co2)
# Now make rank 3 approximation using the Cadzow iterations
F <- cadzow(s, rank = 3)
# Compare with SSA reconstruction
Fr <- reconstruct(s, groups = list(1:3))$F1
plot(F - Fr, type = "l")
}
