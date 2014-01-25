\name{iossa.result}
\alias{iossa.result}
\alias{print.iossa.result}
\alias{summary.iossa.result}
\title{Summary of I-OSSA results}

\description{
  Various routines to print Iterative Oblique SSA results
}

\usage{
  \method{print}{iossa.result}(x, ...)
  \method{summary}{iossa.result}(object, ...)
}

\arguments{
  \item{x, object}{object of class `iossa.result'}
  \item{\dots}{further arguments passed to method}
}

\details{
  An object of class `iossa.result' is a list with the following fields:
  \describe{
    \item{F}{initial series of whole decomposition}
    \item{Fs}{list of series, the final decomposition}
    \item{Cond}{numeric vector with two elements, condition numbers of the final column and row inner products}
    \item{gwcor}{generalized W-correlations matrix of the final decomposition}
    \item{wcor}{W-correlations matrix of the final decomposition}
    \item{initial.wcor}{W-correlation matrix of the initial decomposition}
    \item{hrr}{numeric vector, proportions of high rank components contribution for each of final series
      (denotes how well the series is approximated by a series of finite rank)}
    \item{initial.hrr}{proportions of high rank components contribution for each of initial series}
    \item{initial.rec}{list of initial series}
    \item{iter}{the number of OSSA iterations}
    \item{converged}{logical, whether algorithm has been converged}
    \item{kappa, maxiter, tol}{I-OSSA procedure parameters}
  }
}

\seealso{
  \code{\link[Rssa:iossa]{iossa}},
  \code{\link[Rssa:summary.ssa]{summary.ssa}}.
}

\examples{
# Separate three non-separable sines with different amplitudes
N <- 150
L <- 70

omega1 <- 0.05
omega2 <- 0.06
omega3 <- 0.07

F <- 4*sin(2*pi*omega1 * (1:N)) + 2*sin(2*pi*omega2 * (1:N)) + sin(2*pi*omega3 * (1:N))
s <- ssa(F, L)
ios <- iossa(s, nested.groups = list(1:2, 3:4, 5:6), kappa = NULL, maxIter = 200, tol = 1e-8)

print(ios)
print(ios$iossa.result)
plot(ios$iossa.result)
}
