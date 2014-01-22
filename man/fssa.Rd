\name{fssa}
\alias{fssa}
\title{Nested FSSA redecomposition}

\description{
  Perform Filtered SSA nested redecomposition of the ssa-object.
}


\usage{
  fssa(x, nested.groups, FILTER = diff, gamma = 1, \dots)
}

\arguments{
  \item{x}{SSA object holding the one-dimensional SSA (classic 1d-SSA, Toeplitz SSA or PSSA) decomposition}
  \item{nested.groups}{vector of elementary components indices for nested redecomposition. If list passed,
    \code{unlist} and then \code{unique} will be applyed}
  \item{FILTER}{function or numeric vector of impulse response coefficients for filtration}
  \item{gamma}{weight coefficient of filtration addition. See `Details' and `References'}
  \item{\dots}{additional arguments passed to \code{\link{decompose}}
    routines}
}

\value{
  Object of `ossa' class --- modified input object.
}

\details{
  Filtered SSA performs nested (low-rank) Quotient SVD redecomposition of specified eigentriples.
  Inner product is conventional in column space and in row space it has form:
  \deqn{%
    \langle x, y \rangle = (x, y) + \gamma^2(\Phi(x), \Phi(y)),
  }{<x, y> = (x, y) + gamma^2(Phi(x), Phi(y))}
  where \eqn{(\cdot, \cdot)}{(., .)} denotes conventional inner product and
  `\eqn{\Phi}{Phi}' is function (supposed to be linear operator, but not any checkings performed)
  or linear filtration, specified by \code{FILTER} argument.

  If \eqn{\Phi}{Phi} is difference (by default), such Filtered SSA called `DerivSSA'
  (this very case is dicussed in referenced article).
}

\references{
  Golyandina N. and Shlemov A. (2013): \emph{Variations of Singular Spectrum Analysis
  for separability improvement: non-orthogonal decompositions of time series},
  \url{http://arxiv.org/abs/1308.4022}
}

\seealso{
  \code{\link[Rssa:iossa]{iossa}}.
}

\examples{
# Separate two non-strongly separated series
N <- 150
L <- 70
omega1 <- 1/5
omega2 <- 1/10

v <- sin(2*pi*omega1 * (1:N)) + sin(2*pi*omega2 * (1:N))
s <- ssa(v, L)
fs <- fssa(s, nested.groups = 1:4, gamma = 100)

plot(wcor(s, groups = 1:4))
plot(wcor(fs, groups = 1:4))

# Real-life example: Australian Wine Sales
data(AustralianWine)
s <- ssa(AustralianWine[1:120, "Fortified"], L = 60)
fs <- fssa(s, nested.groups = list(6:7, 8:9, 10:11), gamma = 10)

plot(reconstruct(fs, groups = list(6:7, 8:9, 10:11)))
plot(wcor(s, groups = 6:11))
plot(wcor(fs, groups = 6:11))
}
