\name{forecast}
\alias{forecast.ssa}
\alias{forecast.1d.ssa}
\alias{forecast.toeplitz.ssa}
\alias{predict.ssa}
\alias{predict.1d.ssa}
\alias{predict.mssa}
\alias{predict.toeplitz.ssa}
\title{Perform SSA forecasting of series}

\description{
  All-in-one function to perform SSA forecasting of one-dimensional series.
}

\usage{
\method{forecast}{1d.ssa}(object,
         groups, h = 1,
         method = c("recurrent", "vector"),
         interval = c("none", "confidence", "prediction"),
         only.intervals = TRUE,
         \dots,
         drop = TRUE, drop.attributes = FALSE, cache = TRUE)
\method{forecast}{toeplitz.ssa}(object,
         groups, h = 1,
         method = c("recurrent", "vector"),
         interval = c("none", "confidence", "prediction"),
         only.intervals = TRUE,
         \dots,
         drop = TRUE, drop.attributes = FALSE, cache = TRUE)
\method{predict}{1d.ssa}(object,
        groups, len = 1,
        method = c("recurrent", "vector"),
        interval = c("none", "confidence", "prediction"),
        only.intervals = TRUE,
        \dots,
        drop = TRUE, drop.attributes = FALSE, cache = TRUE)
\method{predict}{toeplitz.ssa}(object,
        groups, len = 1,
        method = c("recurrent", "vector"),
        interval = c("none", "confidence", "prediction"),
        only.intervals = TRUE,
        \dots,
        drop = TRUE, drop.attributes = FALSE, cache = TRUE)
\method{predict}{mssa}(object,
        groups, len = 1,
        method = c("recurrent", "vector"),
        direction = c("column", "row"),
        \dots,
        drop = TRUE, drop.attributes = FALSE, cache = TRUE)
}

\arguments{
  \item{object}{SSA object holding the decomposition}
  \item{groups}{list, the grouping of eigentriples to be used in the forecast}
  \item{h,len}{the desired length of the forecasted series}
  \item{method}{method of forecasting to be used}
  \item{interval}{type of interval calculation}
  \item{only.intervals}{logical, if 'TRUE' then bootstrap method is used
    for confidence bounds only, otherwise --- mean bootstrap forecast is
    returned as well}
  \item{direction}{direction of forecast in multichannel SSA case, "column"
    stands for so-called L-forecast and "row" stands for K-forecast}
  \item{\dots}{further arguments passed for forecast routines
    (e.g. \code{level} argument to \code{bforecast})}
  \item{drop}{logical, if 'TRUE' then the result is coerced to series
    itself, when possible (length of 'groups' is one)}
  \item{drop.attributes}{logical, if 'TRUE' then the forecast routines do not try
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
f <- forecast(s, groups = list(1:6), method = "recurrent", bootstrap = TRUE, len = 24, R = 10)
\donttest{
# Plot the result including the last 24 points of the series
plot(f, include = 24, shadecols = "green", type = "l")
# Use of predict() for prediction
p <- predict(s, groups = list(1:6), method = "recurrent", len = 24)
# Simple plotting
plot(p, ylab = "Forecasteed Values")
}
}
