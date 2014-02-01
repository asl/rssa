\name{fossa}
\alias{fossa}
\title{Nested Filter-adjusted O-SSA decomposition}

\description{
  Perform nested decomposition by Filter-adjusted O-SSA.
}

\usage{
  fossa(x, nested.groups, FILTER = diff, gamma = 1, \dots)
}

\arguments{
  \item{x}{SSA object holding the one-dimensional SSA (Basic 1D-SSA, Toeplitz SSA) decomposition}
  \item{nested.groups}{vector of elementary component indices for
    nested decomposition. The argument is coerced to a vector, if necessary}
  \item{FILTER}{function or numeric vector of impulse response
    coefficients for filter adjustment}
  \item{gamma}{weight of filter adjustment. See `Details' and `References'}
  \item{\dots}{additional arguments passed to \code{\link{decompose}}
    routines}
}

\value{
  Object of class `ossa'. The field `ossa.set' contains the vector of indices
  of elementary components used in Filter-adjusted O-SSA (that is, used in nested.groups).
}

\details{
  Filter-adjusted O-SSA performs a nested decomposition specified by a
  number eigentriples via Oblique SSA with a specific inner product in
  the row space:
  \deqn{%
    \langle x, y \rangle = (x, y) + \gamma^2(\Phi(x), \Phi(y)),
  }{<x, y> = (x, y) + gamma^2(Phi(x), Phi(y))}
  where \eqn{(\cdot, \cdot)}{(., .)} denotes conventional inner product and
  `\eqn{\Phi}{Phi}' is function (supposed to be linear operator, but
  right now no additional checks are performed)
  or linear filtration which is specified by \code{FILTER} argument.

  The default value of \eqn{\Phi}{Phi} corresponds to sequential
  differences, that is, to derivation. Such version of Filter-adjusted
  O-SSA is called `DerivSSA'. See `References' for more details.
}

\references{
  Golyandina N. and Shlemov A. (2013): \emph{Variations of Singular Spectrum Analysis
  for separability improvement: non-orthogonal decompositions of time series},
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
omega1 <- 1/5
omega2 <- 1/10

v <- sin(2*pi*omega1 * (1:N)) + sin(2*pi*omega2 * (1:N))
s <- ssa(v, L)
fs <- fossa(s, nested.groups = 1:4, gamma = 100)

plot(wcor(s, groups = 1:4))
plot(wcor(fs, groups = 1:4))


# Real-life example: Australian Wine Sales
data(AustralianWine)
s <- ssa(AustralianWine[1:120, "Fortified"], L = 60)
fs <- fossa(s, nested.groups = list(6:7, 8:9, 10:11), gamma = 10)

plot(reconstruct(fs, groups = list(6:7, 8:9, 10:11)))
plot(wcor(s, groups = 6:11))
plot(wcor(fs, groups = 6:11))


# Real life example: improving of strong separability
data(USUnemployment)
s <- ssa(USUnemployment[, "MALE"])
fs <- fossa(s, nested.groups = 1:13, gamma = 1000)

# w-cor matrices before and after FOSSA
plot(wcor(s, groups = 1:13))
plot(wcor(fs, groups = 1:13))

# eigenvectors before and after FOSSA
plot(s, type = "vectors", idx = 1:13)
plot(fs, type = "vectors", idx = 1:13)

# 2D plots of periodic eigenvectors before and after FOSSA
plot(s, type = "paired", idx = c(5, 12))
plot(fs, type = "paired", idx = c(1, 3))
}
