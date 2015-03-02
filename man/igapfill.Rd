\name{igapfill}
\alias{igapfill}
\alias{igapfill.ssa}
\alias{igapfill.1d.ssa}
\alias{igapfill.cssa}
\alias{igapfill.toeplitz.ssa}
\alias{igapfill.nd.ssa}
\title{Perform SSA gapfilling via iterative reconstruction}

\description{
  Perform iterative gapfilling of the series.
}

\usage{
\method{igapfill}{1d.ssa}(x, groups, fill = NULL, tol = 1e-6, maxiter = 0,
          norm = function(x) sqrt(max(x^2)),
          base = c("original", "reconstructed"), \dots, trace = FALSE,
          drop = TRUE, drop.attributes = FALSE, cache = TRUE)
\method{igapfill}{cssa}(x, groups, fill = NULL, tol = 1e-6, maxiter = 0,
          norm = function(x) sqrt(max(x^2)),
          base = c("original", "reconstructed"), \dots, trace = FALSE,
          drop = TRUE, drop.attributes = FALSE, cache = TRUE)
\method{igapfill}{toeplitz.ssa}(x, groups, fill = NULL, tol = 1e-6, maxiter = 0,
          norm = function(x) sqrt(max(x^2)),
          base = c("original", "reconstructed"), \dots, trace = FALSE,
          drop = TRUE, drop.attributes = FALSE, cache = TRUE)
\method{igapfill}{nd.ssa}(x, groups, fill = NULL, tol = 1e-6, maxiter = 0,
          norm = function(x) sqrt(max(x^2)),
          base = c("original", "reconstructed"), \dots, trace = FALSE,
          drop = TRUE, drop.attributes = FALSE, cache = TRUE)
}

\arguments{
  \item{x}{Shaped SSA object holding the decomposition}
  \item{groups}{list, the grouping of eigentriples to be used in the forecast}
  \item{fill}{initial values for missed entries, recycled if necessary;
      if missed, then average of the series will be used}
  \item{tol}{tolerance for reconstruction iterations}
  \item{maxiter}{upper bound for the number of iterations}
  \item{norm}{distance function used for covergence criterion}
  \item{base}{series used as a 'seed' for gapfilling: original or
    reconstructed according to the value of \code{groups} argument}
  \item{\dots}{additional arguments passed to \code{\link{reconstruct}}
    routines}
  \item{trace}{logical, indicates whether the convergence process should be traced}
  \item{drop}{logical, if 'TRUE' then the result is coerced to series
    itself, when possible (length of 'groups' is one)}
  \item{drop.attributes}{logical, if 'TRUE' then the attributes of the input series
    are not copied to the reconstructed ones.}
  \item{cache}{logical, if 'TRUE' then intermediate results will be
    cached in the SSA object.}
}

\details{
  Iterative gapfilling starts from filling missed entries with initial
  values, then the missed values are imputed from the successive
  reconstructions. This process continues until convergence up to a
  stationary point (e.g. filling / reconstruction does not change missed
  values at all).
}

\note{
  The method is very sensitive to the initial value of missed entries
  ('fill' argument). If the series are not stationary (e.g. contains
  some trend) than the method may be prohibitely slow, or even fail to
  converge or produce bogus results.
}

\value{
  List of objects with gaps filled in. Elements of the list have the
  same names as elements of \code{groups}. If group is unnamed,
  corresponding component gets name `Fn', where `n' is its index in
  \code{groups} list.

  Or, the forecasted object itself, if length of groups is one and 'drop = TRUE'.
}

\references{
  Kondrashov, D. & Ghil, M. (2006) \emph{Spatio-temporal filling of
    missing points in geophysical data sets}. Nonlinear Processes In
  Geophysics, Vol. 13(2), pp. 151-159.
}

\seealso{
  \code{\link{Rssa}} for an overview of the package, as well as,
  \code{\link[Rssa:igapfill]{gapfill}},
}

\examples{
# Produce series with gaps
F <- co2; F[100:200] <- NA
# Perform shaped SSA
s <- ssa(F, L = 72)
# Fill in gaps using the trend and 2 periodicty components
# Due to trend, provide a linear filler to speedup the process
fill <- F; fill[100:200] <- F[99] + (1:101)/101*(F[201] - F[99])
g <- igapfill(s, groups = list(1:6), fill = fill, maxit = 50)
# Compare the result
plot(g)
lines(co2, col = "red")
}
