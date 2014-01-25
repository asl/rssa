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
  For two series, generalized W-covariation is defined as follows:
  \deqn{%
    \mathrm{gwcov}(F_1, F_2) =
      \langle L^\dagger X_1 (R^\dagger)^\mathrm{T},
              L^\dagger X_2 (R^\dagger)^\mathrm{T} \rangle_\mathrm{F},
  }{%
    gwcov(F_1, F_2) = <L* X_1 (R*)^T, L* X_2 (R*)^T>_F,
  }
  where
  \eqn{X_1, X_2} --- trajectory matrices of series \eqn{F_1, F_2}
  correspondingly, \eqn{L = [U_{b_1} : ... : U_{b_r}], R = [V_{b_1}: ... V_{b_r}]},
  where \eqn{b_1, \dots, b_r} are elementary component indices denoted by \code{basis} argument,
  `\eqn{\langle \cdot, \cdot \rangle_\mathrm{F}}{<., .>_F}' denotes Frobenius matrix inner product and
  `\eqn{\dagger}{*}' denotes Moore-Penrose pseudo-inverse matrix.

  And generalized W-correltaion is defined the following way:
  \deqn{%
    \mathrm{gwcor}(F_1, F_2) = \frac{\mathrm{gwcov}(F_1, F_2)}
      {\sqrt{\mathrm{gwcov}(F_1, F_1) \cdot \mathrm{gwcov(F_2, F_2)}}}
  }{%
    gwcor(F_1, F_2) = gwcov(F_1, F_2) / sqrt(gwcov(F_1, F_1) gwcov(F_2, F_2))
  }

  Generalized W-correlation is OSSA analogue of W-correlation --- measure of series separability.
  If I-OSSA procedure separates series exactly, their generalized W-correlation is almost zero.
}

\references{
  Golyandina N. and Shlemov A. (2013): \emph{Variations of Singular Spectrum Analysis
  for separability improvement: non-orthogonal decompositions of time series},
  \url{http://arxiv.org/abs/1308.4022}
}

\seealso{
  \code{\link[Rssa:wcor]{wcor}},
  \code{\link[Rssa:iossa]{iossa}}.
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

p.wcor <- plot(wcor(ios, groups = list(1:2, 3:4)))
p.gwcor <- plot(genwcor(ios, groups = list(1:2, 3:4)), main = "GW-correlation matrix")
print(p.wcor, split = c(1, 1, 2, 1), more = TRUE)
print(p.gwcor, split = c(2, 1, 2, 1))
}
