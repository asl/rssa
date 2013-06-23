\name{parestimate}
\alias{parestimate}
\alias{parestimate.1d.ssa}
\alias{parestimate.toeplitz.ssa}
\title{Estimate periods from (set of) eigenvectors}

\description{
  Function to estimate sine periods given a set of SSA eigenvectors.
}

\usage{
\method{parestimate}{1d.ssa}(x, groups, method = c("pairs", "esprit-ls"), \dots, drop = TRUE)
\method{parestimate}{toeplitz.ssa}(x, groups, method = c("pairs", "esprit-ls"), \dots, drop = TRUE)
}

\arguments{
  \item{x}{SSA object}
  \item{groups}{list of indices of eigenvectors to estimate from}
  \item{\dots}{further arguments passed to 'decompose' routine, if
    necessary}
  \item{drop}{logical, if 'TRUE' then the result is coerced to lowest
    dimension, when possible (length of 'groups' is one)}
  \item{method}{estimation method, 'pairs' for rough estimation based on
    pair of eigenvectors, 'esprit-ls' for least-squares ESPRIT}
}

\value{
  a list of objects of S3-class `fdimpars.1d'. Each object is a list with 5 components:
  'roots' - complex roots of minimal LRR characteristic polynom,
  'periods' - periods of dumped sinusoids,
  'frequencies' - frequencies of dumped sinusoids,
  'moduli' - moduli of roots and
  'rates' - rates of exponential trend (\code{rates == log(moduli)}).

  For 'method' = 'pairs' all moduli equal to 1 and all rates equal to 0.

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

  Also `print' and `plot' methods are implemented for class `fdimpars.1d'.
}

\references{
  Golyandina, N., Zhigljavsky, A. (2013): \emph{Singular Spectrum
    Analysis for time series}. Springer Briefs in Statistics. Springer.

  Roy, R., Kailath, T., (1989): \emph{ESPRIT: estimation of signal parameters via
  rotational invariance techniques}. IEEE Trans. Acoust. 37, 984--995.
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
}

