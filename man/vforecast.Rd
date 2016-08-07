\name{vforecast}
\alias{vforecast}
\alias{vforecast.default}
\alias{vforecast.ssa}
\alias{vforecast.1d.ssa}
\alias{vforecast.toeplitz.ssa}
\alias{vforecast.mssa}
\alias{vforecast.cssa}
\alias{vforecast.pssa.1d.ssa}
\title{Perform vector SSA forecasting of the series}

\description{
  Perform vector SSA forecasting of the series.
}

\usage{
\method{vforecast}{1d.ssa}(x, groups, len = 1, only.new = TRUE, \dots,
          drop = TRUE, drop.attributes = FALSE)
\method{vforecast}{toeplitz.ssa}(x, groups, len = 1, only.new = TRUE, \dots,
          drop = TRUE, drop.attributes = FALSE)
\method{vforecast}{toeplitz.ssa}(x, groups, len = 1, only.new = TRUE, \dots,
          drop = TRUE, drop.attributes = FALSE)
\method{vforecast}{mssa}(x, groups, len = 1,
          direction = c("row", "column"),
          only.new = TRUE, \dots,
          drop = TRUE, drop.attributes = FALSE)
\method{vforecast}{cssa}(x, groups, len = 1, only.new = TRUE, \dots,
          drop = TRUE, drop.attributes = FALSE)
\method{vforecast}{pssa.1d.ssa}(x, groups, len = 1, only.new = TRUE, \dots,
          drop = TRUE, drop.attributes = FALSE)
}

\arguments{
  \item{x}{SSA object holding the decomposition}
  \item{groups}{list, the grouping of eigentriples to be used in the forecast}
  \item{len}{integer, the desired length of the forecasted series}
  \item{direction}{direction of forecast in multichannel SSA case, "column"
    stands for so-called L-forecast and "row" stands for K-forecast}
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
  
  For multichannel SSA, forecast uses the formulae from
  (Golyandina et al, 2015).
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
    \url{https://www.jstatsoft.org/article/view/v067i02}

  Golyandina, N. and Stepanov, D. (2005): \emph{SSA-based approaches to
    analysis and forecast of multidimensional time series}. In
    Proceedings of the 5th St.Petersburg Workshop on Simulation, June
    26-July 2, 2005, St. Petersburg State University, St. Petersburg,
    293--298. \url{http://www.gistatgroup.com/gus/mssa2.pdf}
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

# Forecast `co2' trend by SSA with projections
s <- ssa(co2, column.projector = 2, row.projector = 2)
len <- 100
vfor <- vforecast(s, groups = list(trend = seq_len(nspecial(s))), len = len, only.new = FALSE)
matplot(data.frame(c(co2, rep(NA, len)), vfor), type = "l")

# Forecast finite rank series with polynomial component by SSA with projections
v <- 5000 * sin(2*pi / 13 * (1:100)) +  (1:100)^2 + 10000
s <- ssa(v, row.projector = 2, column.projector = 2)
plot(vforecast(s, groups = list(all = 1:6), len = 100, only.new = FALSE), type = "l")
}
