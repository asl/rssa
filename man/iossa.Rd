\name{iossa}
\alias{iossa}
\title{Iterative OSSA nested redecomposition}

\description{
  Perform Iterative OSSA (I-OSSA) algorithm to the ssa-object.
}

\usage{
iossa(x, nested.groups, \dots, tol = 1e-5, kappa = 1.2,
      maxiter = 1000,
      initial.approx = reconstruct(x, nested.groups),
      norm = function(x) sqrt(mean(x^2)),
      trace = FALSE,
      kappa.balance = 0.5)
}

\arguments{
  \item{x}{SSA object holding the one-dimensional SSA (classic 1d-SSA, Toeplitz SSA or PSSA) decomposition}
  \item{nested.groups}{list or named list, describes initial grouping for I-OSSA iterations}
  \item{tol}{tolerance for I-OSSA iterations}
  \item{kappa}{`kappa' parameter for sigma-correction (see `Details' and `References') procedure.
    If 'NULL', sigma-correction will not be performed}
  \item{maxiter}{upper bound for the number of iterations}
  \item{initial.approx}{list of series the same length as \code{nested.groups}, using for initial signal
    subspaces' approximations}
  \item{norm}{function, norm series norm for convergence detection}
  \item{trace}{logical, indicating whether convergence process should be traced}
  \item{kappa.balance}{proportion of distribution sigma-correction multiplier
    between colums and rows inner probucts}
  \item{\dots}{additional arguments passed to \code{\link{decompose}}
    routines}
}

\value{
  Object of `ossa' class --- modificated input object. Also is contains following fiels:
  \describe{
    \item{ossa.result}{object of `ossa.result' class, list which contains algorithm traces,
      final inner products condition numbers and so on. This object can be printed or plotted}
    \item{iossa.groups}{list, describes final (after sigma-correction) component grouping}
    \item{iossa.groups.all}{list, describes cumulative grouping after several sequentional I-OSSA applications}
  }
}

\details{
  TODO Write some short description of algorithm.
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

plot(reconstruct(ios, groups = ios$iossa.groups))
print(ios$ossa.result)
plot(ios$ossa.result)

# Separate two non-separable sines with equal amplitudes
N <- 200
L <- 100
omega1 <- 0.07
omega2 <- 0.06

F <- sin(2*pi*omega1 * (1:N)) + sin(2*pi*omega2 * (1:N))
s <- ssa(F, L)
# Apply FSSA and then I-OSSA
fs <- fssa(s, nested.groups = 1:4)
ios <- iossa(fs, nested.groups = list(1:2, 3:4), maxIter = 1000)

plot(reconstruct(s, groups = list(1:2, 3:4)))
plot(reconstruct(ios, groups = ios$iossa.groups))
plot(ios$ossa.result)
print(ios$ossa.result)

plot(wcor(ios, groups = 1:4))
plot(genwcor(ios, groups = 1:4))
}
