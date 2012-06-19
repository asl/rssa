\name{lrr-routines}
\alias{roots}
\alias{roots.lrr}
\alias{plot.lrr}

\title{Miscellaneous LRR rouintes}

\description{
  Small helper functions to operate on LRRs.
}

\usage{
\method{roots}{lrr}(x, ..., method = c("companion", "polyroot"))
\method{plot}{lrr}(x, ..., raw = FALSE)
}

\arguments{
  \item{x}{'LRR' object}
  \item{method}{methods used for poly roots calculation: via eigenvalues
    of companion matrix or R's standard polyroot routine}
  \item{\dots}{other arguments forwarded to \code{plot}}
  \item{raw}{logical, if 'TRUE' then \code{plot} routine will not add any
    additional plot components (e.g. unit circle)}
}

\details{
  \code{roots} function calculates the roots of the characteristic
  polynomial of the LRR.

  \code{plot} just plots these roots adding unit circle if \code{raw}
  agrument is not 'TRUE'.
}

\value{
  vector with the roots of the of the characteristic
  polynomial of the LRR for \code{roots} function call.
}

\seealso{
  \code{\link[Rssa:lrr]{lrr}}
}

\examples{
# Decompose 'co2' series with some parameters
s <- new.ssa(co2, L = 24)
# Evaluate the LRR
l <- lrr(s, group = c(1:3))
# Calculate the roots of the LRR
roots(l)
# Plot the roots
\dontrun{plot(l)}
}
