\name{basis2lrf}
\alias{basis2lrf}
\title{Derive LRF From Set of Basis Vectors}

\description{
  Function to calculate the Linear Recurrence Formula (LRF) given a set
  of basis vectors.
}

\usage{
basis2lrf(U, eps = sqrt(.Machine$double.eps))
}

\arguments{
  \item{U}{Matrix containing the basis vectors in columns}
  \item{eps}{Tolerance for verticality checking}
}

\details{
  Denote by \eqn{U_i} the columns of matrix \eqn{U}. Denote by
  \eqn{\tilde{U}_{i}} the same vector \eqn{U_i} but without the
  last coordinate. Denote the last coordinate of \eqn{U_i} by
  \eqn{\pi_i}. The returned value is
  \deqn{
    \mathcal{R} = \frac{1}{1-\nu^2}\sum_{i=1}^{d}{\pi_i \tilde{U}_{i}},
  }
  where
  \deqn{
    \nu^2 = \pi_1^2 + \dots + \pi_d^2.
  }
}

\value{
  a 'numeric' vector of length equals to number of rows in U minus 1.
}

\note{
  In most cases you should not use this function directly. Use
  \code{\link{lrf}} routine instead.
}

\seealso{
  \code{\link{lrf}}
}

\references{
  Golyandina, N., Nekrutkin, V. and Zhigljavsky, A. (2001): \emph{Analysis of
    Time Series Structure: SSA and related techniques.} Chapman and Hall/CRC. ISBN 1584881941
}

\examples{
# Decompose 'co2' series with default parameters
s <- new.ssa(co2, neig = 20)
# Derive the LRF which governs the whole series
print(basis2lrf(s$U))
}

\keyword{algebra}


