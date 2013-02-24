\name{vforecast}
\alias{vforecast}
\alias{vforecast.default}
\alias{vforecast.ssa}
\alias{vforecast.1d.ssa}
\alias{vforecast.toeplitz.ssa}
\title{Perform vector SSA forecasting of the series}

\description{
  Perform vector SSA forecasting of the series.
}

\usage{
\method{vforecast}{1d.ssa}(x, groups, len = 1, only.new = TRUE, \dots,
          drop = TRUE, drop.attributes = FALSE)
\method{vforecast}{toeplitz.ssa}(x, groups, len = 1, only.new = TRUE, \dots,
          drop = TRUE, drop.attributes = FALSE)
}

\arguments{
  \item{x}{SSA object holding the decomposition}
  \item{groups}{list, the grouping of eigentriples to be used in the forecast}
  \item{len}{integer, the desired length of the forecasted series}
  \item{only.new}{logical, if 'TRUE' then only forecasted values are returned,
    whole series otherwise}
  \item{\dots}{additional arguments passed to \code{\link{decompose}}
    routines}
  \item{drop}{logical, if 'TRUE' then the result is coerced to series
    itself, when possible (length of 'groups' is one)}
  \item{drop.attributes}{logical, if 'TRUE' then the attributes of the input series
    are not copied to the reconstructed ones.}
}

\details{
  The routines applies the vectors SSA forecasting algorithm to
  produce the new series which is expected to 'continue' the current
  series on the basis of the decomposition given. Vector forecast
  differs from recurrent forecast in such way that it continues the set
  of vectors in the subspace spanning the chosen eigenvectors (the same
  formula as described in \code{\link{lrr}} is used for constructing of
  the last components of the new vectors) and then derive the series out
  of this extended set of vectors.
}

\value{
  List of vectors of forecasted series. Elements of the list are
  named 'F1', 'F2', and so on.

  The series itself, if length of groups is one and 'drop = TRUE'.
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
# Produce 24 forecasted values of the series using different sets of eigentriples
# as a base space for the forecast.
vfor <- vforecast(s, groups = list(c(1,4), 1:4), len = 24, only.new=FALSE)
matplot(data.frame(c(co2, rep(NA, 24)), vfor), type="l")
}
