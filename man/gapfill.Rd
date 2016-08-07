\name{gapfill}
\alias{gapfill}
\alias{gapfill.1d.ssa}
\alias{gapfill.toeplitz.ssa}
\alias{gapfill.mssa}
\alias{gapfill.cssa}

\title{Perform SSA gapfilling via forecast}

\description{
  Perform SSA gapfilling of the series.
}

\usage{
\method{gapfill}{1d.ssa}(x, groups, base = c("original", "reconstructed"),
        method = c("sequential", "simultaneous"),
        alpha = function(len) seq.int(0, 1, length.out = len), \dots,
        drop = TRUE, drop.attributes = FALSE, cache = TRUE)
\method{gapfill}{mssa}(x, groups, base = c("original", "reconstructed"),
        alpha = function(len) seq.int(0, 1, length.out = len), \dots,
        drop = TRUE, drop.attributes = FALSE, cache = TRUE)
\method{gapfill}{cssa}(x, groups, base = c("original", "reconstructed"),
        method = c("sequential", "simultaneous"),
        alpha = function(len) seq.int(0, 1, length.out = len), \dots,
        drop = TRUE, drop.attributes = FALSE, cache = TRUE)
\method{gapfill}{toeplitz.ssa}(x, groups, base = c("original", "reconstructed"),
        method = c("sequential", "simultaneous"),
        alpha = function(len) seq.int(0, 1, length.out = len), \dots,
        drop = TRUE, drop.attributes = FALSE, cache = TRUE)
}

\arguments{
  \item{x}{Shaped SSA object holding the decomposition}
  \item{groups}{list, the grouping of eigentriples to be used in the forecast}
  \item{base}{series used as a 'seed' for gapfilling: original or
    reconstructed according to the value of \code{groups} argument}
  \item{method}{method used for gapfilling, "sequential" means to
    filling by a recurrent forecast from complete parts; "simultaneous"
    tries to build a projections onto the signal subspace. See
    'References' for more info.}
  \item{alpha}{weight used for combining forecasts from left and right
    when method = "sequential"; 0.5 means that the forecasts are averaged,
    0 (1) means that only forecast from the left (right correspondingly)
    is used, arbitrary function could be specified; by default linear
    weights are used.}
  \item{\dots}{additional arguments passed to \code{\link{reconstruct}}
    routines}
  \item{drop}{logical, if 'TRUE' then the result is coerced to series
    itself, when possible (length of 'groups' is one)}
  \item{drop.attributes}{logical, if 'TRUE' then the attributes of the input series
    are not copied to the reconstructed ones.}
  \item{cache}{logical, if 'TRUE' then intermediate results will be
    cached in the SSA object.}
}

\details{
  The function fills in the missed entries in the series. Both methods
  described in (Golyandina, Osipov, 2007) are implemented:
  \itemize{
    \item method = "sequential" performs forecast from complete chunks
    onto incomplete. For internal gaps forecast is performed from both
    sides of the gap and average is taken in order to reduce the
    forecast error. For gaps in the beginning or end of the series the
    method coincides with ordinary recurrent forecast;
    \item method = "simultaneous" performs gap filling via projections
    onto signal subspace. The method may fail if insufficient complete
    observations are provided.
  }
}

\value{
  List of objects with gaps filled in. Elements of the list have the
  same names as elements of \code{groups}. If group is unnamed,
  corresponding component gets name `Fn', where `n' is its index in
  \code{groups} list.

  Or, the forecasted object itself, if length of groups is one and 'drop = TRUE'.
}

\references{
  N. Golyandina, E. Osipov (2007): \emph{The "Caterpillar"-SSA method
  for analysis of time series with missing values}. Journal of
  Statistical Planning and Inference, Vol. 137, No. 8, Pp 2642--2653
  \url{http://www.gistatgroup.com/cat/mvssa1en.pdf}
}

\seealso{
  \code{\link{Rssa}} for an overview of the package, as well as,
  \code{\link[Rssa:rforecast]{rforecast}},
  \code{\link[Rssa:igapfill]{igapfill}},
  \code{\link[Rssa:clplot]{clplot}},
  \code{\link[Rssa:summarize.gaps]{summarize.gaps}},
}

\examples{
# Produce series with gaps
F <- co2; F[100:200] <- NA
# Perform shaped SSA
s <- ssa(F, L = 72)
# Fill in gaps using the trend and 2 periodicty components
g <- gapfill(s, groups = list(1:6))
# Compare the result
plot(g)
lines(co2, col = "red")
}
