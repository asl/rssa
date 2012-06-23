\name{sforecast}
\alias{sforecast}
\alias{sforecast.default}
\alias{sforecast.ssa}
\alias{sforecast.1d-ssa}
\title{Perform sliding SSA forecasting of the series}

\description{
  Perform sliding SSA forecasting of the series.
}

\usage{
\method{sforecast}{ssa}(x, group, len = 1, \dots, cache = TRUE)
}

\arguments{
  \item{x}{SSA object holding the decomposition}
  \item{group}{the grouping of eigentriples to be used in the forecast}
  \item{len}{integer, the desired length of the forecasted series}
  \item{\dots}{additional arguments passed to \code{\link{decompose}} routines}
  \item{cache}{logical, if 'TRUE' then intermediate results will be
    cached in the SSA object.}
}

\details{
  The routing performs the so-called multistart sliding forecast
  applying the linear recurrence relation to sliding segments of the
  original series. The resulting series can be used to predict the
  estimation quality of the LRR on the series.
}

\value{
  Vector of forecasted series of length \eqn{N - L - len}, where N is
  the length of the series and L is the window length. The series
  denotes the LRR forecast of last \eqn{N - L - len} values of the
  series itself.
}

\references{
  Golyandina, N., Nekrutkin, V. and Zhigljavsky, A. (2001): \emph{Analysis of
    Time Series Structure: SSA and related techniques.} Chapman and Hall/CRC. ISBN 1584881941
}
\seealso{
  \code{\link[Rssa:rforecast]{rforecast}}
}

\examples{
# Decompose 'co2' series with default parameters
s <- new.ssa(co2)
# Produce the series forecast from 5th forecasted values derived from
# sliding segement of the series using the first 3 eigentriples as a
# base space for the forecast.
sforecast(s, group = 1:3, len = 5)
}
