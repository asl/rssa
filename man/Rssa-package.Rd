\name{Rssa-package}
\alias{Rssa-package}
\alias{Rssa}
\docType{package}
\title{A collection of methods for singular spectrum analysis}
\description{
  Singular Spectrum Analysis (SSA, in short) is a modern
  non-parametric method for the timeseries analysis. This
  package provides a set of fast and reliable implementations of various
  routines to perform the decomposition and reconstruction of the
  timeseries.
}

\details{
  Typically the use of the package starts with the \emph{decomposition} of the
  time series using \code{\link{new.ssa}}. After this a suitable
  \emph{grouping} of the elementary time series is required. This can be done
  heuristically, for example, via looking at the eigenplots of the
  decomposition (\code{\link[Rssa:plot.ssa]{plot}}). Alternatively, one can
  examine the so-called w-correlation matrix (\code{\link{wcor}}). Next step
  includes the \emph{reconstruction} of the time-series using the selected
  grouping (\code{\link[Rssa:reconstruct.ssa]{reconstruct}}).
}

\author{Anton Korobeynikov \email{anton@korobeynikov.info}}
\references{
  Golyandina, N., Nekrutkin, V. and Zhigljavsky, A. (2001): \emph{Analysis of
    Time Series Structure: SSA and related techniques.} Chapman and Hall/CRC. ISBN 1584881941
 
  Korobeynikov, A. (2010) \emph{Computation- and space-efficient implementation of
    SSA.} Statistics and Its Interface, Vol. 3, No. 3, Pp. 257-268
}
\keyword{package}
\seealso{
    \code{\link{new.ssa}}, \code{\link[Rssa:decompose.ssa]{decompose}},
    \code{\link[Rssa:reconstruct.ssa]{reconstruct}},
    \code{\link{wcor}}, \code{\link[Rssa:plot.ssa]{plot}}
}
\examples{
s <- new.ssa(co2) # Perform the decomposition using the default window length
summary(s)        # Show various information about the decomposition
plot(s)           # Show the plot of the eigenvalues
f <- reconstruct(s, groups = list(1, c(2, 3), 4)) # Reconstruct into 3 series
plot(co2)         # Plot the original series
lines(f$F1, col = "blue")            # Extract the trend
lines(f$F1+f$F2, col = "red")        # Add the periodicity
lines(f$F1+f$F2+f$F3, col = "green") # Add slow-varying component
}
