\name{rforecast}
\alias{rforecast}
\alias{rforecast.default}
\alias{rforecast.ssa}
\alias{rforecast.1d.ssa}
\alias{rforecast.toeplitz.ssa}
\alias{rforecast.mssa}
\alias{rforecast.cssa}
\alias{rforecast.pssa.1d.ssa}
\title{Perform recurrent SSA forecasting of the series}

\description{
  Perform recurrent SSA forecasting of the series.
}

\usage{
\method{rforecast}{1d.ssa}(x, groups, len = 1, base = c("reconstructed", "original"),
          only.new = TRUE, reverse = FALSE, \dots,
          drop = TRUE, drop.attributes = FALSE, cache = TRUE)
\method{rforecast}{toeplitz.ssa}(x, groups, len = 1, base = c("reconstructed", "original"),
          only.new = TRUE, reverse = FALSE, \dots,
          drop = TRUE, drop.attributes = FALSE, cache = TRUE)
\method{rforecast}{mssa}(x, groups, len = 1, base = c("reconstructed", "original"),
          direction = c("row", "column"), only.new = TRUE, \dots, drop = TRUE,
          drop.attributes = FALSE, cache = TRUE)
\method{rforecast}{cssa}(x, groups, len = 1, base = c("reconstructed", "original"),
          only.new = TRUE, reverse = FALSE, \dots,
          drop = TRUE, drop.attributes = FALSE, cache = TRUE)
\method{rforecast}{pssa.1d.ssa}(x, groups, len = 1, base = c("reconstructed", "original"),
          only.new = TRUE, reverse = FALSE, \dots,
          drop = TRUE, drop.attributes = FALSE, cache = TRUE)
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
  \item{reverse}{logical, direction of forecast in 1D SSA case, 'FALSE'
    (default) means that the forecast moves forward in the time and
    'TRUE' means the opposite}
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
  series on the basis of a given decomposition. The algorithm
  sequentialy projects the incomplete embedding vectors (from either the original
  or the reconstructed series) onto the subspace spanned by the
  selected eigentriples of the decomposition to derive the missed
  (last) values of the such vectors. Then the filled value 
  In such a way the forecasting elements are produced on
  one-by-one basis.

  It is shown in Golyandina et al (2001) that this approach corresponds
  to application of a linear recurrence formula (the same formula as 
  described in \code{\link{lrr}}) to initial data taken from either the original
  or the reconstructed series.

  In particular, the \eqn{m}-th step of the forecast is calculated by
  means of linear recurrence relation (see \code{\link[Rssa:lrr]{lrr}}) as
  \eqn{y_{n+m}=\sum_{k=1}^{L-1} a_k y_{n+m-k}} where the starting points
  \eqn{y_{n-(L-2)}}, \ldots, \eqn{y_{n}} are taken from the reconstructed
  time series (\code{base="reconstructed"}) or from the initial
  (\code{base="initial"}) time series.

  For multichannel SSA the column forecast is obtained via applying the
  LRR to each series separately. The row forecast is more complicated
  and is based on a multivariate LRR. Forecast uses the formulae from
  Golyandina and Stepanov (2005) and Golyandina et.al (2015).

  For details of 1D-SSA recurrent forecasting, see Section 3.2.1.2 and 
  Algorithm 3.5 in Golyandina et al (2018).
  For details of MSSA recurrent forecasting, see Section 4.3.1.2  and 
  Algorithm 4.4 (column forecasting).
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
  
  Golyandina, N., Korobeynikov, A., Shlemov, A. and Usevich, K. (2015):
  \emph{Multivariate and 2D Extensions of Singular Spectrum Analysis
    with the Rssa Package}. Journal of Statistical Software, Vol. 67, Issue 2.
  \doi{10.18637/jss.v067.i02}

  Golyandina, N. and Stepanov, D. (2005): \emph{SSA-based approaches to
    analysis and forecast of multidimensional time series}. In
    Proceedings of the 5th St.Petersburg Workshop on Simulation, June
    26-July 2, 2005, St. Petersburg State University, St. Petersburg,
    293--298. \url{https://www.gistatgroup.com/gus/mssa2.pdf}
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
matplot(data.frame(c(co2, rep(NA, 24)), rfor), type = "l")

# Forecast `co2' trend by SSA with projections
s <- ssa(co2, column.projector = 2, row.projector = 2)
len <- 100
rfor <- rforecast(s, groups = list(trend = seq_len(nspecial(s))), len = len, only.new = FALSE)
matplot(data.frame(c(co2, rep(NA, len)), rfor), type = "l")

# Forecast finite rank series with polynomial component by SSA with projections
v <- 5000 * sin(2*pi / 13 * (1:100)) +  (1:100)^2 + 10000
s <- ssa(v, row.projector = 2, column.projector = 2)
plot(rforecast(s, groups = list(all = 1:6), len = 100, only.new = FALSE), type = "l")
}
