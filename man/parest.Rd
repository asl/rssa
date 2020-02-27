\name{parestimate}
\alias{parestimate}
\alias{parestimate.1d.ssa}
\alias{parestimate.toeplitz.ssa}
\alias{parestimate.mssa}
\alias{parestimate.cssa}
\alias{parestimate.nd.ssa}
\title{Estimate periods from (set of) eigenvectors}

\description{
  Function to estimate the parameters (frequencies and rates) given a set of SSA eigenvectors.
}

\usage{
\method{parestimate}{1d.ssa}(x, groups, method = c("esprit", "pairs"),
            subspace = c("column", "row"),
            normalize.roots = NULL,
            dimensions = NULL,
            solve.method = c("ls", "tls"),
            \dots,
            drop = TRUE)
\method{parestimate}{toeplitz.ssa}(x, groups, method = c("esprit", "pairs"),
            subspace = c("column", "row"),
            normalize.roots = NULL,
            dimensions = NULL,
            solve.method = c("ls", "tls"),
            \dots,
            drop = TRUE)
\method{parestimate}{mssa}(x, groups, method = c("esprit", "pairs"),
            subspace = c("column", "row"),
            normalize.roots = NULL,
            dimensions = NULL,
            solve.method = c("ls", "tls"),
            \dots,
            drop = TRUE)
\method{parestimate}{cssa}(x, groups, method = c("esprit", "pairs"),
            subspace = c("column", "row"),
            normalize.roots = NULL,
            dimensions = NULL,
            solve.method = c("ls", "tls"),
            \dots,
            drop = TRUE)
\method{parestimate}{nd.ssa}(x, groups,
            method = c("esprit"),
            subspace = c("column", "row"),
            normalize.roots = NULL,
            dimensions = NULL,
            solve.method = c("ls", "tls"),
            pairing.method = c("diag", "memp"),
            beta = 8,
            \dots,
            drop = TRUE)
}

\arguments{
  \item{x}{SSA object}
  \item{groups}{list of indices of eigenvectors to estimate from}
  \item{\dots}{further arguments passed to 'decompose' routine, if
    necessary}
  \item{drop}{logical, if 'TRUE' then the result is coerced to lowest
    dimension, when possible (length of \code{groups} is one)}
  \item{dimensions}{a vector of dimension indices to perform ESPRIT along. 'NULL' means all dimensions.}
  \item{method}{For 1D-SSA, Toeplitz SSA, and MSSA:
      parameter estimation method,
      'esprit' for 1D-ESPRIT (Algorithm 3.3 in Golyandina et al (2018)),
      'pairs' for rough estimation based on pair of eigenvectors (Algorithm 3.4 in Golyandina et al (2018)).
      For nD-SSA: parameter estimation method. 
      For now only 'esprit' is supported (Algorithm 5.6 in Golyandina et al (2018)).}
  \item{solve.method}{approximate matrix equation solving method, 'ls' for least-squares, 'tls' for total-least-squares.}
  \item{pairing.method}{method for esprit roots pairing, 'diag' for `2D-ESPRIT diagonalization', 'memp' for ``MEMP with an
    improved pairing step'}
  \item{subspace}{which subspace will be used for parameter estimation}
  \item{normalize.roots}{logical vector or 'NULL', force signal roots to lie on unit circle.
    'NULL' means automatic selection: normalize iff circular topology OR Toeplitz SSA used}
  \item{beta}{In nD-ESPRIT, coefficient(s) in convex linear combination of
    shifted matrices. The length of \code{beta} should be \code{ndim - 1}, where \code{ndim} is the number of independent dimensions.
    If only one value is passed, it is expanded to a geometric progression.}
}

\value{
  For 1D-SSA (and Toeplitz), a list of objects of S3-class `fdimpars.1d'. Each object is a list with 5 components:
  \describe{
    \item{roots}{complex roots of minimal LRR characteristic polynomial}
    \item{periods}{periods of dumped sinusoids}
    \item{frequencies}{frequencies of dumped sinusoids}
    \item{moduli}{moduli of roots}
    \item{rates}{rates of exponential trend (\code{rates == log(moduli)})}
  }

  For 'method' = 'pairs' all moduli are set equal to 1 and all rates equal to 0.

  For nD-SSA, a list of objects of S3-class `fdimpars.nd'. Each object
  is named list of \code{n} `fdimpars.1d' objects, each for corresponding
  spatial coordinate.

  In all cases elements of the list have the same names as elements of
  \code{groups}. If group is unnamed, corresponding component gets name
  `Fn', where `n' is its index in \code{groups} list.

  If 'drop = TRUE' and length of 'groups' is one, then corresponding
  list of estimated parameters is returned.
}

\details{
  See Sections 3.1 and 5.3 in Golyandina et al (2018) for full details.

  Briefly, the time series is assumed to satisfy the model
  \deqn{
    x_n = \sum_k{C_k\mu_k^n}
  }
  for complex \eqn{\mu_k} or, alternatively,
  \deqn{
    x_n = \sum_k{A_k \rho_k^n \sin(2\pi\omega_k n + \phi_k)}.
  }
  The return value are the estimated moduli and arguments of complex
  \eqn{\mu_k}, more precisely, \eqn{\rho_k} ('moduli') and \eqn{T_k =
  1/\omega_k} ('periods').

  For images, the model
  \deqn{
    x_{ij}=\sum_k C_k \lambda_k^i \mu_k^j
  }
  is considered.

  Also `print' and `plot' methods are implemented for classes
  `fdimpars.1d' and `fdimpars.nd'.
}


\references{
  Golyandina N., Korobeynikov A., Zhigljavsky A. (2018):
  \emph{Singular Spectrum Analysis with R.}  Use R!.
    Springer, Berlin, Heidelberg.

  Roy, R., Kailath, T., (1989): \emph{ESPRIT: estimation of signal parameters via
  rotational invariance techniques}. IEEE Trans. Acoust. 37, 984--995.

  Rouquette, S., Najim, M. (2001): \emph{Estimation of frequencies and damping factors by two-
   dimensional esprit type methods}. IEEE Transactions on Signal Processing 49(1), 237--245.

  Wang, Y., Chan, J-W., Liu, Zh. (2005): \emph{Comments on ``estimation of frequencies and
    damping factors by two-dimensional esprit type methods''}.
    IEEE Transactions on Signal Processing 53(8), 3348--3349.

  Shlemov A, Golyandina N (2014) Shaped extensions of Singular Spectrum Analysis. In: 21st
  international symposium on mathematical theory of networks and systems, July 7--11, 2014.
  Groningen, The Netherlands, pp 1813--1820.
}

\seealso{
  \code{\link{Rssa}} for an overview of the package, as well as,
  \code{\link[Rssa:ssa]{ssa}},
  \code{\link[Rssa:lrr]{lrr}},
}

\examples{
# Decompose 'co2' series with default parameters
s <- ssa(co2, neig = 20)
# Estimate the periods from 2nd and 3rd eigenvectors using 'pairs' method
print(parestimate(s, groups = list(c(2, 3)), method = "pairs"))
# Estimate the peroids from 2nd, 3rd, 5th and 6th eigenvectors using ESPRIT
pe <- parestimate(s, groups = list(c(2, 3, 5, 6)), method = "esprit")
print(pe)
plot(pe)

\donttest{
# Artificial image for 2D SSA
mx <- outer(1:50, 1:50,
            function(i, j) sin(2*pi * i/17) * cos(2*pi * j/7) + exp(i/25 - j/20)) +
      rnorm(50^2, sd = 0.1)
# Decompose 'mx' with default parameters
s <- ssa(mx, kind = "2d-ssa")
# Estimate parameters
pe <- parestimate(s, groups = list(1:5))
print(pe)
plot(pe, col = c("green", "red", "blue"))

# Real example: Mars photo
data(Mars)
# Decompose only Mars image (without background)
s <- ssa(Mars, mask = Mars != 0, wmask = circle(50), kind = "2d-ssa")
# Reconstruct and plot texture pattern
plot(reconstruct(s, groups = list(c(13,14, 17, 18))))
# Estimate pattern parameters
pe <- parestimate(s, groups = list(c(13,14, 17, 18)))
print(pe)
plot(pe, col = c("green", "red", "blue", "black"))
}
}
