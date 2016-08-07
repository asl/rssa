\name{owcor}
\alias{owcor}
\title{Calculate generalized (oblique) W-correlation matrix}

\description{
  Function calculates oblique W-correlation matrix for the series.
}

\usage{
  owcor(x, groups, ..., cache = TRUE)
}

\arguments{
  \item{x}{the input object of `ossa' class}
  \item{groups}{list of numeric vectors, indices of elementary components
    used for reconstruction. The elementary components must belong to
    the current OSSA component set}
  \item{\dots}{further arguments passed to \code{reconstruct} routine}
  \item{cache}{logical, if 'TRUE' then intermediate results will be
    cached in 'ssa' object.}
}

\value{
  Object of class `wcor.matrix'
}

\details{
  Matrix of oblique weighted correlations will be computed.
  For two series, oblique W-covariation is defined as follows:
  \deqn{%
    \mathrm{owcov}(F_1, F_2) =
      \langle L^\dagger X_1 (R^\dagger)^\mathrm{T},
              L^\dagger X_2 (R^\dagger)^\mathrm{T} \rangle_\mathrm{F},
  }{%
    owcov(F_1, F_2) = <L* X_1 (R*)^T, L* X_2 (R*)^T>_F,
  }
  where
  \eqn{X_1, X_2} denotes the trajectory matrices of series \eqn{F_1, F_2}
  correspondingly, \eqn{L = [U_{b_1} : ... : U_{b_r}], R = [V_{b_1}: ... V_{b_r}]},
  where \eqn{\\\{b_1, \dots, b_r\\\}}{\{b_1, \dots, \b_r\}} is current OSSA component set
  (see description of `ossa.set' field of `ossa' object),
  `\eqn{\langle \cdot, \cdot
  \rangle_\mathrm{F}}{<., .>_F}' denotes Frobenius matrix inner product
  and `\eqn{\dagger}{*}' denotes Moore-Penrose pseudo-inverse matrix.

  And oblique W-correlation is defined the following way:
  \deqn{%
    \mathrm{owcor}(F_1, F_2) = \frac{\mathrm{owcov}(F_1, F_2)}
      {\sqrt{\mathrm{owcov}(F_1, F_1) \cdot \mathrm{owcov(F_2, F_2)}}}
  }{%
    owcor(F_1, F_2) = owcov(F_1, F_2) / sqrt(owcov(F_1, F_1) owcov(F_2, F_2))
  }

  Oblique W-correlation is OSSA analogue of W-correlation, that is, a
  measure of series separability. If I-OSSA procedure separates series
  exactly, their oblique W-correlation will be equal to zero.
}

\references{
  Golyandina N. and Shlemov A. (2015): \emph{Variations of Singular Spectrum Analysis
  for separability improvement: non-orthogonal decompositions of time series}, 
  Statistics and Its Interface. Vol.8, No 3, P.277-294.
  \url{http://arxiv.org/abs/1308.4022}
}

\seealso{
  \code{\link{Rssa}} for an overview of the package, as well as,
  \code{\link[Rssa:wcor]{wcor}},
  \code{\link[Rssa:iossa]{iossa}},
  \code{\link[Rssa:iossa]{fossa}}.
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
p.owcor <- plot(owcor(ios, groups = list(1:2, 3:4)), main = "OW-correlation matrix")
print(p.wcor, split = c(1, 1, 2, 1), more = TRUE)
print(p.owcor, split = c(2, 1, 2, 1))
}
