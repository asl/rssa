\name{grouping.auto}
\alias{grouping.auto}
\alias{grouping.auto.ssa}
\alias{grouping.auto.toeplitz.ssa}
\alias{grouping.auto.1d.ssa}

\title{
  Group decomposition components automatically using periodograms
}

\description{
  Group elementary components automatically using their frequency contributions
}

\usage{
  \method{grouping.auto}{1d.ssa}(x, groups,
             base = c("series", "vectors"),
             vectors = c("eigen", "factor"),
             freq.bins = 2,
             threshold = 0,
             drop = TRUE,
             \dots)
}

\arguments{
  \item{x}{SSA object}
  \item{groups}{indices of elementary components for grouping}
  \item{base}{what use for periodogram computation: elementary reconstructed series or decomposition vectors}
  \item{vectors}{for \code{base = 'vectors'} choose vectors to use}
  \item{freq.bins}{single integer number > 1 (the number of intervals),
    vector of frequency breaks (of length >=2) or list of frequency ranges.
    For each range, if only one element provided it will
    be used as the upper bound and the lower bound will be zero}
  \item{threshold}{contribution threshold. If zero then dependent grouping approach will be used}
  \item{drop}{logical, whether to exclude empty groups from resulted list}
  \item{\dots}{additional arguments passed to \code{\link{reconstruct}}
    routines}
}

\value{
  list of groups (indices of elementary components)
}

\details{
  TODO
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
  grouping.auto(ss, base = "series", freq.bins = 2, threshold = 0)
  grouping.auto(ss, base = "series", freq.bins = list(c(0.1), c(0.1, 0.2)), threshold = 0)
  grouping.auto(ss, base = "series", freq.bins = c(0.1, 0.2), threshold = 0)

  g <- grouping.auto(ss, base = "series", freq.bins = 8, threshold = 0)
  plot(reconstruct(ss, groups = g))

  g <- grouping.auto(ss, base = "series", freq.bins = list(0.1, 0.2, 0.3, 0.4, 0.5), threshold = 0.95)
  plot(reconstruct(ss, groups = g))
}
