\name{bforecast}
\alias{bforecast}
\alias{bforecast.default}
\alias{bforecast.ssa}
\alias{bforecast.1d.ssa}
\alias{bforecast.toeplitz.ssa}
\title{Perform bootstrap SSA forecasting of the series}

\description{
  Perform bootstrap SSA forecasting of the one-dimensional series.
}

\usage{
\method{bforecast}{1d.ssa}(x, groups, len = 1, R = 100, level = 0.95,
          type = c("recurrent", "vector"), \dots,
          drop = TRUE, drop.attributes = FALSE, cache = TRUE)
\method{bforecast}{toeplitz.ssa}(x, groups, len = 1, R = 100, level = 0.95,
          type = c("recurrent", "vector"), \dots,
          drop = TRUE, drop.attributes = FALSE, cache = TRUE)
}

\arguments{
  \item{x}{SSA object holding the decomposition}
  \item{groups}{list, the grouping of eigentriples to be used in the forecast}
  \item{len}{the desired length of the forecasted series}
  \item{R}{number of bootstrap replications}
  \item{level}{vector of confidence levels for bounds}
  \item{type}{the type of forecast method to be used during bootstrapping}
  \item{\dots}{additional arguments passed to forecasting routines}
  \item{drop}{logical, if 'TRUE' then the result is coerced to series
    itself, when possible (length of 'groups' is one)}
  \item{drop.attributes}{logical, if 'TRUE' then the attributes of the input series
    are not copied to the reconstructed ones}
  \item{cache}{logical, if 'TRUE' then intermediate results will be
    cached in the SSA object}
}

\details{
  The routine uses the reconstruction residuals in order to calculate
  their empirical distribution (the residuals are assumed to be
  stationary). Empirical distribution of the residuals is used to
  perform bootstrap series simulation. Such bootsrapped series are then
  extended via selected forecast method. Finally, the distribution of
  forecasted values is used to calculate bootstrap estimate of series
  forecast and confidence bounds.
}

\value{
  List of matricies. Each matrix has 1 + 2*length(level) columns and
  'len' rows. First column contains the forecasted values, remaining
  columns --- low and upper bootstrap confidence bounds for average
  forecasted values.

  The matrix itself, if length of groups is one and 'drop = TRUE'.
}

\seealso{
  \code{\link{Rssa}} for an overview of the package, as well as,
  \code{\link[Rssa:rforecast]{rforecast}},
  \code{\link[Rssa:bforecast]{bforecast}},
  \code{\link[Rssa:forecast]{forecast}}.
}

\examples{
# Decompose 'co2' series with default parameters
s <- ssa(co2)
# Produce 24 forecasted values and confidence bounds of the series using
# the first 3 eigentriples as a base space for the forecast.
\donttest{
f <- bforecast(s, groups = list(1:3), len = 24, R = 50)
matplot(f, col = c("black", "red", "red"), type='l')}
}
