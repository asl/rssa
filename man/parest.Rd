\name{parestimate}
\alias{parestimate}
\alias{parestimate.1d.ssa}
\alias{parestimate.toeplitz.ssa}
\alias{parestimate.2d.ssa}
\alias{parestimate.shaped2d.ssa}
\title{Estimate periods from (set of) eigenvectors}

\description{
  Function to estimate sine periods given a set of SSA eigenvectors.
}

\usage{
\method{parestimate}{1d.ssa}(x, groups, method = c("pairs", "esprit-ls", "esprit-tls"), \dots, drop = TRUE)
\method{parestimate}{toeplitz.ssa}(x, groups, method = c("pairs", "esprit-ls", "esprit-tls"), \dots, drop = TRUE)
\method{parestimate}{2d.ssa}(x, groups,
            method = c("esprit-diag-ls", "esprit-diag-tls",
                       "esprit-memp-ls", "esprit-memp-tls"),
            \dots,
            beta = 8,
            drop = TRUE)
\method{parestimate}{shaped2d.ssa}(x, groups,
            method = c("esprit-diag-ls", "esprit-diag-tls",
                       "esprit-memp-ls", "esprit-memp-tls"),
            \dots,
            beta = 8,
            drop = TRUE)
}

\arguments{
  \item{x}{SSA object}
  \item{groups}{list of indices of eigenvectors to estimate from}
  \item{\dots}{further arguments passed to 'decompose' routine, if
    necessary}
  \item{drop}{logical, if 'TRUE' then the result is coerced to lowest
    dimension, when possible (length of 'groups' is one)}
  \item{method}{estimation method. For 1dSSA: 'pairs' for rough estimation based on
    pair of eigenvectors, 'esprit-ls' for least-squares 1D-ESPRIT and
    'esprit-tls' for total-least-squares 1D-ESPRIT.
    For 2dSSA: 'esprit-diag-ls' or 'esprit-diag-tls' for `Two-Dimensional ESPRIT Method'
    (2D-ESPRIT diagonalization) and
    'esprit-memp-ls' or `esprit-memp-tls' for `MEMP with an improved pairing step'.
    Suffix '-ls' or '-tls' means approximate matrix equation solve method, as in 1D-ESPRIT case.}
  \item{beta}{In 2D-ESPRIT, coefficient in considered shift matrices convex linear combination.}
}

\value{
  For 1D-SSA (and Toeplitz), a list of objects of S3-class `fdimpars.1d'. Each object is a list with 5 components:
  'roots' - complex roots of minimal LRR characteristic polynom,
  'periods' - periods of dumped sinusoids,
  'frequencies' - frequencies of dumped sinusoids,
  'moduli' - moduli of roots and
  'rates' - rates of exponential trend (\code{rates == log(moduli)}).

  For 'method' = 'pairs' all moduli equal to 1 and all rates equal to 0.

  For 2d SSA and Shaped 2d SSA, a list of objects of S3-class `fdimpars.2d'. Each object is unnamed list of two
  `fdimpars.1d' objects, each for corresponding spatial coordinate.

  If 'drop = TRUE' and length of 'groups' is one, then the list of
  periods itself.
}

\details{
  The time series is assumed to satisfy the model
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

  Also `print' and `plot' methods are implemented for classes `fdimpars.1d' and `fdimpars.2d'.
}

\references{
  Golyandina, N., Zhigljavsky, A. (2013): \emph{Singular Spectrum
    Analysis for time series}. Springer Briefs in Statistics. Springer.

  Roy, R., Kailath, T., (1989): \emph{ESPRIT: estimation of signal parameters via
  rotational invariance techniques}. IEEE Trans. Acoust. 37, 984--995.

  Rouquette, S., Najim, M. (2001): \emph{Estimation of frequencies and damping factors by two-
   dimensional esprit type methods}. IEEE Transactions on Signal Processing 49(1), 237--245.

  Wang, Y., Chan, J-W., Liu, Zh. (2005): \emph{Comments on ``estimation of frequencies and
    damping factors by two-dimensional esprit type methods''}.
    IEEE Transactions on Signal Processing 53(8), 3348--3349.
}

\examples{
# Decompose 'co2' series with default parameters
s <- ssa(co2, neig = 20)
# Estimate the periods from 2nd and 3rd eigenvectors using default
# 'pairs' method
print(parestimate(s, groups = list(c(2, 3))))
# Estimate the peroids from 2nd, 3rd, 5th and 6th eigenvectors using ESPRIT
pe <- parestimate(s, groups = list(c(2, 3, 5, 6)), method = "esprit-ls")
print(pe)
plot(pe)

# Artificial field for 2dSSA
mx <- outer(1:50, 1:50,
            function(i, j) sin(2*pi * i/17) * cos(2*pi * j/7) + exp(i/25 - j/20)) +
      rnorm(50^2, sd = 0.1)
# Decompose 'mx' with default parameters
s <- ssa(mx, kind = "2d-ssa")
# Estimate parameters
pe <- parestimate(s, groups = list(1:5))
print(pe)
plot(pe)

# Decompose 'mx' with circular window
circle <- 5
s <- ssa(mx, kind = "shaped2d-ssa", umask = circle(I(circle)), neig = 10)
# Estimate parameters
pe <- parestimate(s, groups = list(1:5), method = "esprit-memp-tls")
print(pe)
plot(pe)
}

