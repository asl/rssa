\name{lrf-routines}
\alias{roots}
\alias{roots.lrf}
\alias{plot.lrf}

\title{Miscellaneous LRF rouintes}

\description{
  Small helper functions to operate on LRFs.
}

\usage{
\method{roots}{lrf}(x)
\method{plot}{lrf}(x, ..., raw = FALSE)
}

\arguments{
  \item{x}{'LRF' object}
  \item{\dots}{other arguments forwarded to \code{plot}}
  \item{raw}{logical, if 'TRUE' then \code{plot} routine will not add any
    additional plot components (e.g. unit circle)}
}

\details{
  \code{roots} function calculates the roots of the characteristic
  polynomial of the LRF.

  \code{plot} just plots these roots adding unit circle if \code{raw}
  agrument is not 'TRUE'.
}

\value{
  vector with the roots of the of the characteristic
  polynomial of the LRF for \code{roots} function call.
}

\seealso{
  \code{\link[Rssa:lrf]{lrf}}
}

\examples{
# Decompose 'co2' series with some parameters
s <- new.ssa(co2, L = 24)
# Evaluate the LRF
l <- lrf(s, group = c(1:3))
# Calculate the roots of the LRF
roots(l)
# Plot the roots
\dontrun{plot(l)}
}
