\name{summarize.gaps}
\alias{summarize.gaps}
\alias{summarize.gaps.ssa}
\alias{summarize.gaps.1d.ssa}
\alias{summarize.gaps.toeplitz.ssa}
\alias{summarize.gaps.cssa}
\alias{summarize.gaps.default}

\title{Summarize Gaps in a Series}

\description{Provide a summary about the gaps in a series given desired
  window length, namely whether the gap is internal or not, whether it
  is sparce or dense, etc.}

\usage{
\method{summarize.gaps}{1d.ssa}(x, L = NULL)
\method{summarize.gaps}{toeplitz.ssa}(x, L = NULL)
\method{summarize.gaps}{cssa}(x, L = NULL)
\method{summarize.gaps}{default}(x, L)
}

\arguments{
  \item{x}{SSA object}
  \item{L}{vector of window lengths, if missing or NULL, then all viable
    window lengths are considered}
}

\value{
  Object of type 'ssa.gaps': a list with entries which correspond to
  every window length. For each window length, entry is a list of gaps
  with their descriptions.
}

\seealso{
  \code{\link{Rssa}} for an overview of the package, as well as,
  \code{\link[Rssa:gapfill]{gapfill}},
  \code{\link[Rssa:igapfill]{igapfill}},
  \code{\link[Rssa:clplot]{clplot}},
}

\examples{
# Produce series with gaps
F <- co2; F[c(12, 100:200, 250)] <- NA
# Summarize the gaps
s <- ssa(F, L = 72)
g <- summarize.gaps(s, L = c(36, 72, 144))
# Print the results
print(g)
# Plot the proportion of complete lag-vectors
plot(g)
}

