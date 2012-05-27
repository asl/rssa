\name{bforecast}
\alias{bforecast}
\alias{bforecast.default}
\alias{bforecast.ssa}
\alias{bforecast.1d-ssa}
\title{Perform bootstrap SSA forecasting of the series}

\description{
  Perform bootstrap SSA forecasting of the series.
}

\usage{
\method{bforecast}{ssa}(x, group, len = 1, R = 100, level = 0.95, \dots, cache = TRUE)
}

\arguments{
  \item{x}{SSA object holding the decomposition}
  \item{group}{the grouping of eigentriples to be used in the forecast}
  \item{len}{the desired length of the forecasted series}
  \item{R}{number of bootstrap replications}
  \item{level}{confidence level for bounds}
  \item{\dots}{additional arguments passed to \code{\link{reconstruct}} routines}
  \item{cache}{logical, if 'TRUE' then intermediate results will be
    cached in the SSA object.}
}

\details{
  The routine uses the reconstruction residuals in order to calculate
  their empirical distribution (the residuals are assumed to be
  stationary). Empirical distribution of the residuals is used to
  perform bootstrap series simulation. Such bootsrapped series are then
  extended via recurrent forecast. Finally, the distribution of
  forecasted values is used to calculate bootstrap estimate of series
  forecast and confidence bounds.
}

\value{
  Matrix with 3 columns and 'len' rows. First column contains the
  forecasted values, second and third columns - low and upper bootstrap
  confidence bounds for forecasted values.
}

\seealso{
  \code{\link[svd:svd]{svd}},
  \code{\link[Rssa:new.ssa]{new.ssa}},
  \code{\link[Rssa:reconstruct]{reconstruct}},
  \code{\link[Rssa:rforecast]{rforecast}}
}

\examples{
# Decompose 'co2' series with default parameters
s <- new.ssa(co2)
# Produce 5 forecasted values and confidence bounds of the series using
# the first 3 eigentriples as a base space for the forecast.
bforecast(s, group = 1:3, len = 5)
}
