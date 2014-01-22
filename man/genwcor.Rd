\name{genwcor}
\alias{genwcor}
\title{Calculate generalized (oblique) W-correlation matrix}

\description{
  Function calculates oblique W-correlation matrix for the series.
}

\usage{
  genwcor(x, groups, basis, ..., cache = TRUE)
}

\arguments{
  \item{x}{the input object of `ossa' class}
  \item{groups}{list of numeric vectors, indices of elementary components
    used for reconstruction.}
  \item{basis}{vector of elementary component indices which are used for
    inner products construction. If missing, I-OSSA components will be used}
  \item{\dots}{further arguments passed to \code{reconstruct} routine}
  \item{cache}{logical, if 'TRUE' then intermediate results will be
    cached in 'ssa' object.}
}

\value{
  Object of class `wcor.matrix'
}

\details{
  Matrix of generalized (oblique) weighted correlations will be computed.
  For two series, generalized W-covariation defained as follows:
  \deqn{%
    \mathrm{gwcov}(F_1, F_2) = \mathrm{trace}(L^\dagger X_1 (R^\dagger)^\mathrm{T}
      R^\dagger X_2^\mathrm{T} (L^\dagger)^\mathrm{T}),
  }{%
    gwcov(F_1, F_2) = trace(L* X_1 (R*)^T R* X_2^T (L*)^T),
  }
  where \eqn{X_1, X_2} --- trajectory matrices of series \eqn{F_1, F_2}
  correspondingly, \eqn{L = [U_{b_1} : ... : U_{b_r}], R = [V_{b_1}: ... V_{b_r}]},
  where \eqn{b_1, \dots, b_r} are elementary component indices denoted by \code{basis} argument and
  `\eqn{\dagger}{*}' denotes Moore-Penrose pseudo-inverse matrix.

  And generalized W-corraltaion defined the following way:
  \deqn{%
    \mathrm{gwcor}(F_1, F_2) = \frac{\mathrm{gwcov}(F_1, F_2)}
      {\sqrt{\mathrm{gwcov}(F_1, F_1) \cdot \mathrm{gwcov(F_2, F_2)}}}
  }{%
    gwcor(F_1, F_2) = gwcov(F_1, F_2) / sqrt(gwcov(F_1, F_1) gwcov(F_2, F_2))
  }
}

\examples{
# Separate two non-separable sines
N <- 150
L <- 70

omega1 <- 0.06
omega2 <- 0.065

F <- 4*sin(2*pi*omega1 * (1:N)) + sin(2*pi*omega2 * (1:N))
s <- ssa(F, L)
ios <- iossa(s, nested.groups = list(1:2, 3:4), kappa = NULL, maxIter = 200, tol = 1e-8)

plot(wcor(ios, groups = list(1:2, 3:4)))
plot(genwcor(ios, groups = list(1:2, 3:4)))
}
