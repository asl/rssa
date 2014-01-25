\name{iossa.result}
\alias{iossa.result}
\alias{print.iossa.result}
\alias{summary.iossa.result}
\alias{plot.iossa.result}
\title{S3 class for I-OSSA algorithm results storaging}

\description{
  Class for storaging I-OSSA procedure results and methods
  for displaying them.
}

\usage{
  \method{print}{iossa.result}(x, ...)
  \method{summary}{iossa.result}(object, ...)
  \method{plot}{iossa.result}(x, ...)
}

\arguments{
  \item{x, object}{object of `iossa.result' class}
  \item{\dots}{further arguments passed to method}
}

\details{
  `iossa.result' class object is a list of following fields:
  \describe{
    \item{F}{initial series of whole decomposition}
    \item{Fs}{list of series, the final decomposition}
    \item{Cond}{numeric vector of two elements, condition numbers of the final column and row inner products}
    \item{gwcor}{Generalized W-correlations matrix of the final decomposition}
    \item{wcor}{Standard W-correlations matrix of the final decomposition}
    \item{initial.wcor}{Standard W-correlation matrix of the initial decomposition}
    \item{hrr}{numeric vector, proportions of high rank components contribution for each of final series
      (measure of finitrankness)}
    \item{initial.hrr}{proportions of high rank components contribution for each of initial series}
    \item{initial.rec}{list of initial series}
    \item{iter}{the number of iterations}
    \item{converged}{logical, whether algorithm has been converged}
    \item{kappa, maxiter, tol}{corresponding I-OSSA procedure parameters}
  }

  For ssa-objects returned from \code{iossa} routine, method \code{summary} also prints I-OSSA result
  after common SSA info.
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
