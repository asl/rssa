\name{Rssa-package}
\alias{Rssa-package}
\alias{Rssa}
\docType{package}
\title{A collection of methods for singular spectrum analysis}
\description{
  Singular Spectrum Analysis (SSA, in short) is a modern non-parametric
  method for the analysis of time series and digital images. This
  package provides a set of fast and reliable implementations of various
  routines to perform decomposition, reconstruction and forecasting.
}

\details{
  Typically the use of the package starts with the \emph{decomposition}
  of the time series using \code{\link{ssa}}. After this a suitable
  \emph{grouping} of the elementary time series is required. This can be
  done heuristically, for example, via looking at the plots of the
  decomposition (\code{\link[Rssa:plot.ssa]{plot}}). Alternatively, one
  can examine the so-called w-correlation matrix
  (\code{\link{wcor}}). Next step includes the \emph{reconstruction} of
  the time-series using the selected grouping
  (\code{\link[Rssa:reconstruct.ssa]{reconstruct}}). One ends with
  frequency estimation (\code{\link[Rssa:parestimate]{parestimate}}) and
  series forecasting (\code{\link[Rssa:forecast.ssa]{forecast}}).
}

\author{Anton Korobeynikov \email{anton@korobeynikov.info}}
\references{
  Golyandina, N., Nekrutkin, V. and Zhigljavsky, A. (2001):
  \emph{Analysis of Time Series Structure: SSA and related techniques.}
  Chapman and Hall/CRC. ISBN 1584881941
  
  Golyandina, N. and Usevich, K. (2009). \emph{2D-extensions of singular
  spectrum analysis: algorithm and elements of theory.} In Matrix
  Methods: Theory, Algorithms, Applications. World Scientific
  Publishing, 450-474.
  
  Korobeynikov, A. (2010) \emph{Computation- and space-efficient
    implementation of SSA.} Statistics and Its Interface, Vol. 3, No. 3,
  Pp. 257-268

  Golyandina, N., Zhigljavsky, A. (2013): \emph{Singular Spectrum
  Analysis for time series}. Springer Briefs in Statistics. Springer.
}
\keyword{package}
\seealso{
    \code{\link{ssa}}, \code{\link[Rssa:decompose.ssa]{decompose}},
    \code{\link[Rssa:reconstruct.ssa]{reconstruct}},
    \code{\link{wcor}}, \code{\link[Rssa:plot.ssa]{plot}},
    \code{\link[Rssa:parestimate]{parestimate}},
    \code{\link[Rssa:forecast.ssa]{forecast}}
}
\examples{
s <- ssa(co2) # Perform the decomposition using the default window length
summary(s)        # Show various information about the decomposition
plot(s)           # Show the plot of the eigenvalues
f <- reconstruct(s, groups = list(c(1, 4), c(2, 3))) # Reconstruct into 2 series
plot(co2)         # Plot the original series
plot(f)           # Plot the reconstruction
}
