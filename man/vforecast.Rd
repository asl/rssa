\name{vforecast}
\alias{vforecast}
\alias{vforecast.default}
\alias{vforecast.ssa}
\alias{vforecast.1d-ssa}
\title{Perform vector SSA forecasting of the series}

\description{
  Perform vector SSA forecasting of the series.
}

\usage{
\method{vforecast}{ssa}(x, groups, len = 1, \dots)
}

\arguments{
  \item{x}{SSA object holding the decomposition}
  \item{groups}{the grouping of eigentriples to be used in the forecast}
  \item{len}{the desired length of the forecasted series}
  \item{\dots}{additional arguments passed to \code{\link{decompose}} routines}
}

\details{
  The routines applies the vectors SSA forecasting algorithm to
  produce the new series which are expected to 'continue' the current
  series on the basis of the decomposition given. Vector forecast
  differs from recurrent forecast in such way that it tries to extend
  the embedded subspace and then - derive the series out of this
  extended subspace.
}

\value{
  List of vectors of forecasted series. Elements of the list are
  named 'F1', 'F2', and so on.
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
# Produce 5 forecasted values of the series using the first 3
# eigentriples as a base space for the forecast.
vforecast(s, groups = list(1:3), len = 5)
}
