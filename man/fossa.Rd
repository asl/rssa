\name{fossa}
\alias{fossa}
\alias{fossa.ssa}
\title{Nested Filter-adjusted O-SSA decomposition}

\description{
  Perform nested decomposition by Filter-adjusted O-SSA (FOSSA).
}

\usage{
  \method{fossa}{ssa}(x, nested.groups, filter = c(-1, 1), gamma = Inf, normalize = TRUE, \dots)
}

\arguments{
  \item{x}{SSA object holding SSA decomposition}
  \item{nested.groups}{vector of numbers of eigentriples from full decomposition
    for nested decomposition. The argument is coerced to a vector, if
    necessary}
  \item{filter}{numeric vector or array of reversed impulse response (IR)
    coefficients for filter adjustment or list of such vectors or arrays}
  \item{gamma}{weight of filter adjustment. See `Details' and `References'}
  \item{normalize}{logical, whether to normalize left decomposition vectors before filtering}
  \item{\dots}{additional arguments passed to \code{\link{decompose}}
    routines}
}

\value{
  Object of class `ossa'. The field `ossa.set' contains the vector of indices
  of elementary components used in Filter-adjusted O-SSA (that is, used in \code{nested.groups}).
}

\details{
  FOSSA serves for decomposition of series components that are
  mixed due to equal contributions of their elementary components, e.g. of
  sinusoids with equal amplitudes or of complex-form trend and
  periodics. FOSSA performs a new decomposition of a part of the
  ssa-object, which is given by a set of eigentriples. Note that
  eigentriples that do not belong to the chosen set are not changed.

  In particular, Filter-adjusted O-SSA performs a nested decomposition
  specified by a number of eigentriples via Oblique SSA with a specific
  inner product in the row space:
  \deqn{%
    \langle x, y \rangle = (x, y) + \gamma^2(\Phi(x), \Phi(y)),
  }{<x, y> = (x, y) + gamma^2(Phi(x), Phi(y))}
  where \eqn{(\cdot, \cdot)}{(., .)} denotes conventional inner product and
  `\eqn{\Phi}{Phi}' is linear filtration which is specified by \code{filter} argument.

  The default value of \eqn{\Phi}{Phi} corresponds to sequential
  differences, that is, to derivation. Such version of Filter-adjusted
  O-SSA is called `DerivSSA'. See `References' for more details.

  \subsection{\code{filter} argument}{
    For 1D-SSA, Toeplitz-SSA and MSSA:
    Filter can be given by a vector or a list of vectors. Each vector corresponds to reversed IR
    for a filter, these filters are applied independently and their results are stacked such that the matrix
    \eqn{[X:\Phi_1(X):\Phi_2(X)]}{[X:Phi_1(X):Phi_2(X)]} is decomposed.

    For 2D-SSA: the following variants are possible: (1) a list of vectors.
    Each vector corresponds to reversed IR for a filter. Each filter is
    applied to different dimensions, the first to columns, the second to rows, and the results are stacked.
    (2) single vector. Given vector corresponds to one-dimensional filter applied to both dimensions, the same as list of two equal vectors.
    (3) a list of matrices, where each matrix provides 2d filter coefficients and the results are stacked.
    (4) single matrix. Given matrix corresponds to two-dimensional filter applied once, the same as list of one matrix.

    For nD-SSA: the same as for 2D-SSA, a list of vectors for filters by directions, single vector,
    a list of arrays (matroids) for nD filters or single array.
  }


  \subsection{Normalization}{
  Let us explain for the 1D case.
  Let \eqn{X}{X} be the reconstructed matrix, corresponding to the selected eigentriples
  \eqn{\{(\sigma_i,U_i,V_i)\}}{\{sigma_i,U_i,V_i\}}, \eqn{\Psi(X)}{Psi(X)} is the matrix, where the filter is applied to
  each row of \eqn{X}{X}.

  Then \code{normalize = FALSE} corresponds to finding the basis in the column space of
  \eqn{X}{X} by means of the SVD of \eqn{[X, \Psi(X)]}{[X, Phi(X)]}, while
  \code{normalize = TRUE} (by default)
  corresponds to finding the basis by the SVD of \eqn{[V, \Phi(V)]}{[V, Phi(V)]},
  where the rows of
  matrix \eqn{V}{V} are \eqn{V_i}{V_i}.
  The value by default \code{TRUE} guaranties that the contributions
  of sine waves will be ordered by decreasing of frequencies, although can
  slightly worsen the weak separability}
}

\references{
  Golyandina N. and Shlemov A. (2015): \emph{Variations of Singular Spectrum Analysis
  for separability improvement: non-orthogonal decompositions of time series}, 
  Statistics and Its Interface. Vol.8, No 3, P.277-294.
  \url{http://arxiv.org/abs/1308.4022}
}

\seealso{
  \code{\link{Rssa}} for an overview of the package, as well as,
  \code{\link[Rssa:iossa]{iossa}}.
}

\examples{
# Separation of two mixed sine-waves with equal amplitudes
N <- 150
L <- 70
omega1 <- 1/15
omega2 <- 1/10

v <- sin(2*pi*omega1 * (1:N)) + sin(2*pi*omega2 * (1:N))
s <- ssa(v, L)
fs <- fossa(s, nested.groups = 1:4, gamma = 100)

# Rssa does most of the plots via lattice
ws <- plot(wcor(s, groups = 1:4))
wfs <- plot(wcor(fs, groups = 1:4))
plot(ws, split = c(1, 1, 2, 1), more = TRUE)
plot(wfs, split = c(2, 1, 2, 1), more = FALSE)

opar <- par(mfrow = c(2, 1))
plot(reconstruct(s, groups = list(1:2, 3:4)))
plot(reconstruct(fs, groups = list(1:2, 3:4)))
par(opar)

# Real-life example: Australian Wine Sales
\donttest{
data(AustralianWine)
s <- ssa(AustralianWine[1:120, "Fortified"], L = 60)
fs <- fossa(s, nested.groups = list(6:7, 8:9, 10:11), gamma = 10)

plot(reconstruct(fs, groups = list(6:7, 8:9, 10:11)))
plot(wcor(s, groups = 6:11))
plot(wcor(fs, groups = 6:11))

# Real life example: improving of strong separability
data(USUnemployment)
unempl.male <- USUnemployment[, "MALE"]
s <- ssa(unempl.male)
fs <- fossa(s, nested.groups = 1:13, gamma = 1000)

# Comparison of reconstructions
rec <- reconstruct(s, groups = list(c(1:4, 7:11), c(5:6, 12:13)))
frec <- reconstruct(fs, groups <- list(5:13, 1:4))
# Trends
matplot(data.frame(frec$F1, rec$F1, unempl.male), type= 'l',
        col=c("red","blue","black"), lty=c(1,1,2))
# Seasonalities
matplot(data.frame(frec$F2, rec$F2), type = 'l', col=c("red","blue"), lty=c(1,1))

# W-cor matrices before and after FOSSA
ws <- plot(wcor(s, groups = 1:30), grid = 14)
wfs <- plot(wcor(fs, groups = 1:30), grid = 14)
plot(ws, split = c(1, 1, 2, 1), more = TRUE)
plot(wfs, split = c(2, 1, 2, 1), more = FALSE)

# Eigenvectors before and after FOSSA
plot(s, type = "vectors", idx = 1:13)
plot(fs, type = "vectors", idx = 1:13)

# 2D plots of periodic eigenvectors before and after FOSSA
plot(s, type = "paired", idx = c(5, 12))
plot(fs, type = "paired", idx = c(1, 3))

# Compare FOSSA with and without normalize
N <- 150
L <- 70
omega1 <- 1/15
omega2 <- 1/10

v <- 3*sin(2*pi*omega1 * (1:N)) + 2*sin(2*pi*omega2 * (1:N))
s <- ssa(v, L)
fs <- fossa(s, nested.groups = 1:4, gamma = 100)
fs.norm <- fossa(s, nested.groups = 1:4, gamma = 100, normalize = TRUE)
opar <- par(mfrow = c(2, 1))
plot(reconstruct(fs, groups = list(1:2, 3:4)))
plot(reconstruct(fs.norm, groups = list(1:2, 3:4)))
par(opar)

# 2D example
data(Mars)
s <- ssa(Mars)
plot(s, "vectors", idx = 1:50)
plot(s, "series", idx = 1:50)
fs <- fossa(s, nested.groups = 1:50, gamma = Inf)
plot(fs, "vectors", idx = 1:14)
plot(fs, "series", groups = 1:13)

# Filters example, extracting horizontal and vertical stripes
data(Mars)
s <- ssa(Mars)
fs.hor <- fossa(s, nested.groups = 1:50, gamma = Inf,
                filter = list(c(-1, 1), c(1)))
plot(fs.hor, "vectors", idx = 1:14)
plot(fs.hor, "series", groups = 1:13)
fs.ver <- fossa(s, nested.groups = 1:50, gamma = Inf,
                filter = list(c(1), c(-1, 1)))
plot(fs.ver, "vectors", idx = 1:14)
plot(fs.ver, "series", groups = 1:13)
}
}
