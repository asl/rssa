\name{iossa}
\alias{iossa}
\title{Iterative O-SSA nested decomposition}

\description{
  Perform Iterative O-SSA (IOSSA) algorithm.
}

\usage{
iossa(x, nested.groups, \dots, tol = 1e-5, kappa = 2,
      maxiter = 100,
      norm = function(x) sqrt(mean(x^2)),
      trace = FALSE,
      kappa.balance = 0.5,
      svd.method = "auto")
}

\arguments{
  \item{x}{SSA object holding the one-dimensional SSA (classic 1D-SSA, Toeplitz SSA) decomposition}
  \item{nested.groups}{list or named list, describes initial grouping for IOSSA iterations}
  \item{tol}{tolerance for IOSSA iterations}
  \item{kappa}{`kappa' parameter for sigma-correction (see `Details' and `References') procedure.
    If 'NULL', sigma-correction will not be performed}
  \item{maxiter}{upper bound for the number of iterations}
  \item{norm}{function, series norm for convergence detection}
  \item{trace}{logical, indicates whether the convergence process should be traced}
  \item{kappa.balance}{proportion of the distribution of sigma-correction multiplier
    between colums and rows inner products}
  \item{svd.method}{singular value decomposition method, passed to \code{ssa} routine (IOSSA method
    uses Basic 1D SSA for rank and signal subspace estimations)}
  \item{\dots}{additional arguments passed to \code{\link{decompose}}
    routines}
}

\value{
  Object of `ossa' class. In additional to usual `ssa' class
  fields, it also contains the following fields:
  \describe{
    \item{iossa.result}{object of `iossa.result' class, a list which contains algorithm parameters,
      final inner products condition numbers, separability
      measures, the number of iterations and convergence status (see
      \code{\link{iossa.result}})}
    \item{iossa.groups}{list, describes final (after sigma-correction)
      grouping of series components}
    \item{iossa.groups.all}{list, describes cumulative grouping after
      several sequentional IOSSA applications}
    \item{ossa.set}{vector of elementary component indices involved in
      procedure}
  }
}

\details{
  Iterative Oblique SSA (IOSSA) is an iterative (EM-like) method for
  improving separability in SSA.  At each iteration signal subspaces of
  series components are estimated by series decomposition from previous
  iteration, then oblique reduced inner products are constructed from
  the estimated signal subspaces (and optionally sigma-correction
  procedure is applied). Finally nested (low-rank) generalized SVD
  decomposition with new inner products is performed for specified
  eigentriples. Such a procedure is performed until the stationary
  point.

  Sigma-correction procudure does the renormalization of new inner
  products. This prevents the mixing of the components during the next
  iteration. Such approach makes the whole procedure more stable and can
  solve lack of strong separability problem (see References).
}

\references{
  Golyandina N. and Shlemov A. (2013): \emph{Variations of Singular Spectrum Analysis
  for separability improvement: non-orthogonal decompositions of time series},
  \url{http://arxiv.org/abs/1308.4022}
}

\seealso{
  \code{\link[Rssa:fossa]{fossa}},
  \code{\link[Rssa:owcor]{owcor}},
  \code{\link[Rssa:iossa.result]{iossa.result}}.
}

\examples{
# Separate three non-separable sine series with different amplitudes
N <- 150
L <- 70

omega1 <- 0.05
omega2 <- 0.06
omega3 <- 0.07

F <- 4*sin(2*pi*omega1 * (1:N)) + 2*sin(2*pi*omega2 * (1:N)) + sin(2*pi*omega3 * (1:N))
s <- ssa(F, L)
ios <- iossa(s, nested.groups = list(1:2, 3:4, 5:6), kappa = NULL, maxiter = 100, tol = 1e-3)

plot(reconstruct(ios, groups = ios$iossa.groups))
summary(ios)

\donttest{
# Separate two non-separable sines with equal amplitudes
N <- 200
L <- 100
omega1 <- 0.07
omega2 <- 0.06

F <- sin(2*pi*omega1 * (1:N)) + sin(2*pi*omega2 * (1:N))
s <- ssa(F, L)

# Apply FOSSA and then IOSSA
fs <- fossa(s, nested.groups = 1:4)
ios <- iossa(fs, nested.groups = list(1:2, 3:4), maxiter = 100)

opar <- par(mfrow = c(3, 1))
plot(reconstruct(s, groups = list(1:2, 3:4)))
plot(reconstruct(fs, groups = list(1:2, 3:4)))
plot(reconstruct(ios, groups = ios$iossa.groups))
par(opar)
summary(ios)

plot(wcor(ios, groups = 1:4))
plot(owcor(ios, groups = 1:4))
}


\donttest{
# Real life example: improving of strong separability
data(USUnemployment)
s <- ssa(USUnemployment[, "MALE"])
ios <- iossa(s, nested.groups = list(c(1:4, 7:11), c(5:6, 12:13)))

print(ios$iossa.result)

# w-cor matrix before and after IOSSA
plot(wcor(s, groups = 1:13))
plot(wcor(ios, groups = 1:13))

# eigenvectors before and after Iterative O-SSA
plot(s, type = "vectors", idx = 1:13)
plot(ios, type = "vectors", idx = 1:13)

# 2D plots of periodic eigenvectors before and after Iterative O-SSA
plot(s, type = "paired", idx = c(5, 12))
plot(ios, type = "paired", idx = c(10, 12), plot.contrib = FALSE)
}
}
