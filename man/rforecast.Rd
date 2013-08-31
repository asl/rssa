\name{rforecast}
\alias{rforecast}
\alias{rforecast.default}
\alias{rforecast.ssa}
\alias{rforecast.1d.ssa}
\alias{rforecast.toeplitz.ssa}
\alias{rforecast.mssa}
\title{Perform recurrent SSA forecasting of the series}

\description{
  Perform recurrent SSA forecasting of the series.
}

\usage{
\method{rforecast}{1d.ssa}(x, groups, len = 1, base = c("reconstructed", "original"),
          only.new = TRUE, \dots, drop = TRUE, drop.attributes = FALSE, cache = TRUE)
\method{rforecast}{toeplitz.ssa}(x, groups, len = 1, base = c("reconstructed", "original"),
          only.new = TRUE, \dots, drop = TRUE, drop.attributes = FALSE, cache = TRUE)
\method{rforecast}{mssa}(x, groups, len = 1, base = c("reconstructed", "original"),
          direction = c("row", "column"), only.new = TRUE, \dots, drop = TRUE,
          drop.attributes = FALSE, cache = TRUE)
}

\arguments{
  \item{x}{SSA object holding the decomposition}
  \item{groups}{list, the grouping of eigentriples to be used in the forecast}
  \item{len}{integer, the desired length of the forecasted series}
  \item{base}{series used as a 'seed' of forecast: original or
    reconstructed according to the value of \code{groups} argument}
  \item{direction}{direction of forecast in multichannel SSA case, "column"
    stands for so-called L-forecast and "row" stands for K-forecast}
  \item{only.new}{logical, if 'TRUE' then only forecasted values are returned,
    whole series otherwise}
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
  The routines applies the recurrent SSA forecasting algorithm to
  produce the new series which is expected to 'continue' the current
  series on the basis of the decomposition given. The algorithm
  sequentialy projects the incomplete embedding vectors (either original
  or from reconstructed series) onto the subspace spanned by the
  selected eigentriples of the decomposition to derive the missed
  (ending) values of the such vectors.

  In such a way the forecasted elements of the series are produced on
  one-by-one basis.

  In particular, the \eqn{m}-th step of the forecast is calculated by
  means of linear recurrence relation (see \code{\link[Rssa:lrr]{lrr}}) as
  \eqn{y_{n+m}=\sum_{k=1}^{L-1} a_k y_{n+m-k}} where the starting points
  \eqn{y_{n-(L-2)}}, \ldots, \eqn{y_{n}} are taken from the reconstructed
  time series (\code{base="reconstructed"}) or from the initial
  (\code{base="initial"}) time series.

  For multichannel SSA the column forecast is obtained via applying the
  LRR to each series separately. Forecast uses the formulae from
  (Golyandina and Stepanov, 2005).
}

\value{
  List of forecasted objects. Elements of the list have the same names
  as elements of \code{groups}. If group is unnamed, corresponding
  component gets name `Fn', where `n' is its index in \code{groups}
  list.

  Or, the forecasted object itself, if length of groups is one and 'drop = TRUE'.
}

\references{
  Golyandina, N., Nekrutkin, V. and Zhigljavsky, A. (2001): \emph{Analysis of
    Time Series Structure: SSA and related techniques.} Chapman and
  Hall/CRC. ISBN 1584881941

  Golyandina, N. and Stepanov, D. (2005): \emph{SSA-based approaches to
    analysis and forecast of multidimensional time series}. In
    Proceedings of the 5th St.Petersburg Workshop on Simulation, June
    26-July 2, 2005, St. Petersburg State University, St. Petersburg,
    293--298. \url{http://www.gistatgroup.com/gus/mssa2.pdf}
}
\seealso{
  \code{\link{Rssa}} for an overview of the package, as well as,
  \code{\link[Rssa:forecast]{forecast}},
  \code{\link[Rssa:vforecast]{vforecast}},
  \code{\link[Rssa:bforecast]{bforecast}}.
}

\examples{
# Decompose 'co2' series with default parameters
s <- ssa(co2)
# Produce 24 forecasted values of the series using different sets of eigentriples
# as a base space for the forecast.
rfor <- rforecast(s, groups = list(c(1,4), 1:4), len = 24, only.new=FALSE)
matplot(data.frame(c(co2, rep(NA, 24)), rfor), type="l")
}
