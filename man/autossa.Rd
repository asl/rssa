\name{grouping.auto.pgram}
\alias{grouping.auto.pgram}
\alias{grouping.auto.pgram.ssa}
\alias{grouping.auto.pgram.toeplitz.ssa}
\alias{grouping.auto.pgram.1d.ssa}

\title{
  Group elementary series using periodogram
}

\description{
  Group elementary components automatically using their frequency contributions
}

\usage{
  \method{grouping.auto.pgram}{1d.ssa}(x, groups,
             base = c("series", "eigen", "factor"),
             freq.bins = 2,
             threshold = 0,
             method = c("constant", "linear"),
             \dots,
             drop = TRUE)
}

\arguments{
  \item{x}{SSA object}
  \item{groups}{indices of elementary components for grouping}
  \item{base}{input for periodogram: elementary reconstructed series, eigenvectors or factor vectors}
  \item{freq.bins}{single integer number > 1 (the number of intervals),
    vector of frequency breaks (of length >=2) or list of frequency ranges.
    For each range, if only one element provided it will
    be used as the upper bound and the lower bound will be zero}
  \item{threshold}{contribution threshold. If zero then dependent grouping approach will be used}
  \item{method}{method of periodogram interpolation}
  \item{\dots}{additional arguments passed to \code{\link{reconstruct}}
    routines}
  \item{drop}{logical, whether to exclude empty groups from resulted list}
}

\value{
  list of groups (indices of elementary components)
}

\details{
  Elementary components are grouped using their frequency contribution (periodogram).
  Optionally (see argument 'base') periodogram of eigen or factor vectors may be used.

  For each elementary component and for each frequency interval
  (which are specified by 'freq.bins' argument)
  relative (from 0 till 1) contribution
  is computed using one of two methods:
  'constant' (periodogram is considered as a sequence of separate bars)
  or 'linear' (periodogram is linearly interpolated).

  Two approaches of grouping is implemented:
  \describe{
    \item{'independent' or 'threshold'}{Each group includes components
      with frequency contribution in correspondent interval is greater than specified threshold;
      resulted groups can intersect.
      If 'threshold' is a vector, correspondent value of threshold will be using
      for each interval.}
    \item{'dependent'}{Elementary components are separated to disjoint subsets;
      for each component interval with the highest contribution is selected.}
  }

  If 'freq.bins' is named, result groups will take the same names.

  If drop = 'TRUE' (by default), empty groups will be excluded from result.

  See the paper in 'References' for the details of the algorithm.
}

\references{
  Alexandrov, Th., Golyandina, N. (2005):
  \emph{Automatic extraction and forecast of time series cyclic components
    within the framework of SSA.}
    In Proceedings of the 5th St.Petersburg Workshop on Simulation,
    June 26 -- July 2, 2005, St.Petersburg State University, St.Petersburg, Pp. 45--50
    \url{http://www.gistatgroup.com/gus/autossa2.pdf}
}


\seealso{
  \code{\link{Rssa}} for an overview of the package, as well as,
  \code{\link[Rssa:reconstruct]{reconstruct}},
  \code{\link[Rssa:rforecast]{rforecast}},
  \code{\link[Rssa:vforecast]{vforecast}},
  \code{\link[Rssa:parestimate]{parestimate}}
}

\examples{
  ss <- ssa(co2)
  plot(ss, type = "vectors", idx = 1:12)
  plot(ss, type = "vectors", vectors = "factor", idx = 1:12)
  plot(ss, type = "series", groups = 1:12)

  grouping.auto(ss, base = "series", freq.bins = list(0.005), threshold = 0.95)
  grouping.auto(ss, base = "eigen", freq.bins = 2, threshold = 0)
  grouping.auto(ss, base = "factor", freq.bins = list(c(0.1), c(0.1, 0.2)), threshold = 0, method = "linear")
  grouping.auto(ss, freq.bins = c(0.1, 0.2), threshold = 0)

  g <- grouping.auto(ss, freq.bins = 8, threshold = 0)
  plot(reconstruct(ss, groups = g))

  g <- grouping.auto(ss, freq.bins = list(0.1, 0.2, 0.3, 0.4, 0.5), threshold = 0.95)
  plot(reconstruct(ss, groups = g))
}
