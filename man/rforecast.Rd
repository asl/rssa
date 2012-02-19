\name{rforecast}
\alias{rforecast}
\alias{rforecast.default}
\alias{rforecast.ssa}
\alias{rforecast.1d-ssa}
\title{Perform recurrent SSA forecasting of the series}

\description{
  Perform recurrent SSA forecasting of the series.
}

\usage{
\method{rforecast}{ssa}(x, groups, len = 1, base = c("reconstructed", "original"), \dots, cache = TRUE)
}

\arguments{
  \item{x}{SSA object holding the decomposition}
  \item{groups}{the grouping of eigentriples to be used in the forecast}
  \item{len}{the desired length of the forecasted series}
  \item{base}{series used as a 'seed' of forecast: original or
    reconstructed according to the value of \code{groups} argument}
  \item{\dots}{additional arguments passed to \code{\link{reconstruct}} routines}
  \item{cache}{logical, if 'TRUE' then intermediate results will be
    cached in the SSA object.}
}

\details{
  The routines applies the recurrent SSA forecasting algorithm to
  produce the new series which are expected to 'continue' the current
  series on the basis of the decomposition given. The algorithm
  sequentialy projects the incomplete embedding vectors (either original
  or from reconstructed series) onto the subspace spanned by the
  selected eigentriples of the decomposition to derive the missed
  (ending) values of the such vectors.

  In such a way the forecasted elements of the series are produced on
  one-by-one basis.
}

\note{
  In fact the current implementation of the algorithm instead of using
  the direct projections of the series calculates the so-called Linear
  Recurrent Formula which governs the series and uses it for the
  forecast.
}

\value{
  List of vectors of forecasted series. Elements of the list are
  named 'F1', 'F2', and so on.
}

\seealso{
  \code{\link[svd:svd]{svd}},
  \code{\link[Rssa:new.ssa]{new.ssa}},
  \code{\link[Rssa:reconstruct]{reconstruct}}
}

\examples{
# Decompose 'co2' series with default parameters
s <- new.ssa(co2)
# Produce 5 forecasted values of the series using the first 3
# eigentriples as a base space for the forecast.
rforecast(s, groups = list(1:3), len = 5)
}
