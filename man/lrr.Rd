\name{lrr}
\alias{lrr}
\alias{lrr.default}
\alias{lrr.ssa}
\alias{lrr.1d.ssa}
\alias{lrr.toeplitz.ssa}
\alias{roots}
\alias{roots.lrr}
\alias{plot.lrr}

\title{Calculate the Linear Recurrence Relation}

\description{
  Calculates the Linear Recurrence Relation given the one-dimensional 'ssa' object.
}

\usage{
\method{lrr}{1d.ssa}(x, groups, \dots, drop = TRUE)
\method{lrr}{toeplitz.ssa}(x, groups, \dots, drop = TRUE)
\method{lrr}{default}(U, eps = sqrt(.Machine$double.eps), \dots)
\method{roots}{lrr}(x, ..., method = c("companion", "polyroot"))
\method{plot}{lrr}(x, ..., raw = FALSE)
}

\arguments{
  \item{x}{SSA object holding the decomposition for \code{lrr} call or 'lrr'
    object itself for other function calls}
  \item{groups}{list, the grouping of eigentriples used to derive the LRR}
  \item{\dots}{further arguments to be passed to \code{decompose} or
    \code{plot} call, if necessary}
  \item{drop}{logical, if 'TRUE' then the result is coerced to lrr object
    itself, when possible (length of 'groups' is one)}
  \item{U}{Matrix containing the orthonormal basis vectors in columns}
  \item{eps}{Tolerance for verticality checking}
  \item{method}{methods used for calculation of the polynomial roots: via eigenvalues
    of companion matrix or R's standard \code{polyroot} routine}
  \item{raw}{logical, if 'TRUE' then \code{plot} routine will not add any
    additional plot components (e.g. unit circle)}
}

\details{
  Produces the linear recurrence relation from the series. The default
  implementation works as follows.

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

  For \code{lrr.ssa} case the matrix \eqn{U} used is the matrix of basis
  vector corresponding to the selected elementary series.
}

\value{
  Named list of object of class 'lrr' for \code{lrr} function call,
  where elements have the same names as elements of \code{groups}
  (if group is unnamed, corresponding component gets name `Fn',
  where `n' is its index in \code{groups} list).
  Or the object itself if 'drop = TRUE' and groups has length one.

  Vector with the roots of the of the characteristic
  polynomial of the LRR for \code{roots} function call. Roots are
  ordered by moduli decreasing.
}

\examples{
# Decompose 'co2' series with default parameters
s <- ssa(co2, L = 24)
# Calculate the LRR out of first 3 eigentriples
l <- lrr(s, groups = list(1:3))
# Calculate the roots of the LRR
r <- roots(l)
# Moduli of the roots
Mod(r)
# Periods of three roots with maximal moduli
2*pi/Arg(r)[1:3]
# Plot the roots
plot(l)
}
