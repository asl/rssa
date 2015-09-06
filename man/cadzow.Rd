\name{cadzow}
\alias{cadzow}
\alias{cadzow.ssa}

\title{
  Cadzow Iterations
}
\description{
  Perform the finite rank approximation of the series via Cadzow iterations
}

\usage{
\method{cadzow}{ssa}(x, rank, correct = TRUE, tol = 1e-6, maxiter = 0,
         norm = function(x) max(abs(x)),
         trace = FALSE, \dots, cache = TRUE)
}

\arguments{
  \item{x}{input SSA object}
  \item{rank}{desired rank of approximation}
  \item{correct}{logical, if 'TRUE' then additional correction as in
  (Gillard et al, 2013) is performed}
  \item{tol}{tolerance value used for convergence criteria}
  \item{maxiter}{number of iterations to perform, if zero then
    iterations are performed until the convergence}
  \item{norm}{distance function used for covergence criterion}
  \item{trace}{logical, indicates whether the convergence process should be traced}
  \item{\dots}{further arguments passed to \code{reconstruct}}
  \item{cache}{logical, if 'TRUE' then intermediate results will be
    cached in the SSA object.}
}

\details{
  Cadzow iterations aim to solve the problem of the approximation of the
  input series by a series of finite rank. The idea of the algorithm is
  quite simple: alternating projections of the trajectory matrix to
  Hankel and low-rank matrices are performed which hopefully converge to
  a Hankel low-rank matrix.

  Note that the results of one Cadzow iteration with no correction
  coincides with the result of reconstruction by the leading \code{rank}
  components.
  
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
# Decompose co2 series with default parameters
s <- ssa(co2)
# Now make rank 3 approximation using the Cadzow iterations
F <- cadzow(s, rank = 3, tol = 1e-10)
library(lattice)
xyplot(cbind(Original = co2, Cadzow = F), superpose = TRUE)
# All but the first 3 eigenvalues are close to 0
plot(ssa(F))

# Compare with SSA reconstruction
F <- cadzow(s, rank = 3, maxiter = 1, correct = FALSE)
Fr <- reconstruct(s, groups = list(1:3))$F1
print(max(abs(F - Fr)))

# Cadzow with and without weights
set.seed(3)
N <- 60
L <- 30
K <- N - L + 1
alpha <- 0.1

sigma <- 0.1
signal <- cos(2*pi * seq_len(N) / 10)
x <- signal + rnorm(N, sd = sigma)

weights <- rep(alpha, K)
weights[seq(1, K, L)] <- 1
salpha <- ssa(x, L = L,
              column.oblique = "identity",
              row.oblique = weights)
calpha <- cadzow(salpha, rank = 2)

cz <- cadzow(ssa(x, L = L), rank = 2)

print(mean((cz - signal)^2))
print(mean((calpha - signal)^2))
}
