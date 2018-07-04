\name{reconstruct}
\alias{reconstruct}
\alias{reconstruct.ssa}
\title{Perform a series reconstruction}

\description{
  Reconstruct the data given the SSA decomposition and the desired
  grouping of the elementary components.
}

\usage{
\method{reconstruct}{ssa}(x, groups, \dots, drop.attributes = FALSE, cache = TRUE)
}

\arguments{
  \item{x}{SSA object}
  \item{groups}{list of numeric vectors, indices of elementary components
    used for reconstruction, the entries of the list can be named, see
    'Value' for more information}
  \item{\dots}{further arguments passed to routines (e.g. to
    \code{decompose} routine if the continuation is desired).}
  \item{drop.attributes}{logical, if 'TRUE' then the attributes of the
    input objects are not copied to the reconstructed ones.}
  \item{cache}{logical, if 'TRUE' then intermediate results will be
    cached in the SSA object.}
}

\details{
  Reconstruction is performed in a common form for different types of input objects.
  See Section 1.1.2.6 in Golyandina et al (2018) for the explanation.
  Formal algorithms are described in this book in Algorithm 2.2 for 1D-SSA, Algorithm 4.3 for MSSA,
  Algorithm 5.2 for 2D-SSA and Algorithm 5.6 for Shaped 2D-SSA.

  Fast implementation of reconstruction with the help of FFT is described in Korobeynikov (2010)
  for the 1D case and in Section 6.2 (Rank-one quasi-hankelization) of Golyandina et al (2015) 
  for the general case.
}

\value{
  List of reconstructed objects. Elements of the list have the same
  names as elements of \code{groups}. If the group is unnamed, then
  corresponding component will obtain name `Fn', where `n' is its index
  in \code{groups} list.
}

\note{
  By default (argument \code{drop.attributes}) the routine tries to
  preserve all the attributes of the input object. This way, for
  example, the reconstruction result of 'ts' object is the 'ts' object
  with the same time scale.
}

\references{
  Golyandina N., Korobeynikov A., Zhigljavsky A. (2018):
  \emph{Singular Spectrum Analysis with R.}  Use R!.
    Springer, Berlin, Heidelberg.

  Korobeynikov, A. (2010): \emph{Computation- and space-efficient
    implementation of SSA.} Statistics and Its Interface, Vol. 3, No. 3,
  Pp. 257-268

  Golyandina, N., Korobeynikov, A., Shlemov, A. and Usevich, K. (2015):
  \emph{Multivariate and 2D Extensions of Singular Spectrum Analysis
    with the Rssa Package}. Journal of Statistical Software, Vol. 67, Issue 2.
    \url{https://www.jstatsoft.org/article/view/v067i02}
}

\seealso{
  \code{\link{Rssa}} for an overview of the package, as well as,
  \code{\link[Rssa:ssa-input]{ssa-input}},  
  \code{\link[Rssa:ssa]{ssa}},
  \code{\link[Rssa:plot.reconstruction]{plot.reconstruction}},
}

\examples{
# Decompose 'co2' series with default parameters
s <- ssa(co2)
# Reconstruct the series, grouping elementary series.
r <- reconstruct(s, groups = list(Trend = c(1, 4), Season1 = c(2,3), Season2 = c(5, 6)))
plot(r)
# 'groups' argument might contain duplicate entries as well
r <- reconstruct(s, groups = list(1, 1:4, 1:6))
plot(r)

\donttest{
# Real example: Mars photo
data(Mars)
# Decompose only Mars image (without backgroud)
s <- ssa(Mars, mask = Mars != 0, wmask = circle(50), kind = "2d-ssa")
# Reconstruct and plot trend
plot(reconstruct(s, 1), fill.uncovered = "original")
# Reconstruct and plot texture pattern
plot(reconstruct(s, groups = list(c(13, 14, 17, 18))))

# Decompose 'EuStockMarkets' series with default parameters
s <- ssa(EuStockMarkets, kind = "mssa")
r <- reconstruct(s, groups = list(Trend = 1:2))
# Plot original series, trend and residuals superimposed
plot(r, plot.method = "xyplot", superpose = TRUE,
     auto.key = list(columns = 3),
     col = c("blue", "green", "red", "violet"),
     lty = c(rep(1, 4), rep(2, 4), rep(3, 4)))
}
}

