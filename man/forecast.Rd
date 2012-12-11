\name{forecast}
\alias{forecast.ssa}
\alias{forecast.1d.ssa}
\alias{forecast.toeplitz.ssa}
\alias{predict.ssa}
\alias{predict.1d.ssa}
\alias{predict.toeplitz.ssa}
\title{Perform SSA forecasting of series}

\description{
  All-in-one function to perform SSA forecasting of series.
}

\usage{
\method{forecast}{1d.ssa}(object,
         group, len = 1,
         method = c("recurrent", "vector",
                    "bootstrap-recurrent", "bootstrap-vector"),
         \dots,
         drop = FALSE, cache = TRUE)
\method{forecast}{toeplitz.ssa}(object,
         group, len = 1,
         method = c("recurrent", "vector",
                    "bootstrap-recurrent", "bootstrap-vector"),
         \dots,
         drop = FALSE, cache = TRUE)
\method{predict}{1d.ssa}(object,
        group, len = 1,
        method = c("recurrent", "vector",
                   "bootstrap-recurrent", "bootstrap-vector"),
        \dots,
        drop = FALSE, cache = TRUE)
\method{predict}{toeplitz.ssa}(object,
        group, len = 1,
        method = c("recurrent", "vector",
                   "bootstrap-recurrent", "bootstrap-vector"),
        \dots,
        drop = FALSE, cache = TRUE)
}

\arguments{
  \item{object}{SSA object holding the decomposition}
  \item{group}{vector of integers, numbers of eigentriples to be used in the forecast}
  \item{len}{the desired length of the forecasted series}
  \item{method}{method of forecasting to be used. The confidence bounds
    are available only for bootstrap-based methods}
  \item{\dots}{further arguments passed for forecast routines
  (e.g. \code{level} argument to \code{bforecast})}
  \item{drop}{logical, if 'TRUE' then the forecast routines do not try
  to infer the time index arguments for the forecasted series.}
  \item{cache}{logical, if 'TRUE' then intermediate results will be
    cached in the SSA object.}
}

\details{
  This function is a convenient wrapper over other forecast routines (see
  'See Also') turning their value into object of type 'forecast' which
  can be used with the routines from \pkg{forecast} package.
}

\value{
  object of class 'forecast' for \code{forecast} function call,
  predicted series for \code{predict} call.
}

\seealso{
  \code{\link{Rssa}} for an overview of the package, as well as,
  \code{\link[Rssa:rforecast]{rforecast}},
  \code{\link[Rssa:vforecast]{vforecast}},
  \code{\link[Rssa:bforecast]{bforecast}},
  \code{\link[forecast:forecast]{forecast (package)}}
}

\examples{
s <- ssa(co2)
# Calculate 24-point forecast using first 6 components as a base
f <- forecast(s, group = c(1:6), method = "bootstrap-recurrent", len = 24, R = 50)
# Plot the result including the last 24 points of the series
plot(f, include = 24, shadecols = "green", type = "l")
# Use of predict() for prediction
p <- predict(s, group = c(1:6), method = "recurrent", len = 24)
# Simple plotting
plot(p)
}
