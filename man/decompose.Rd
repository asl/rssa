\name{decompose}
\alias{decompose}
\alias{decompose.default}
\alias{decompose.ssa}
\alias{decompose.1d.ssa}
\alias{decompose.2d.ssa}
\alias{decompose.toeplitz.ssa}
\title{Calculate SSA decomposition of the series}

\description{
  Performs the SSA decomposition of the series.
}

\usage{
\method{decompose}{1d.ssa}(x, neig = min(50, L, K), \dots, force.continue = FALSE)
\method{decompose}{toeplitz.ssa}(x, neig = min(50, L, K), \dots, force.continue = FALSE)
\method{decompose}{2d.ssa}(x, neig = min(50, prod(L), prod(K)), \dots)
}

\arguments{
  \item{x}{SSA object holding the decomposition.}
  \item{neig}{number of desired eigentriples.}
  \item{\dots}{additional arguments passed to SVD routines.}
  \item{force.continue}{logical, if TRUE then continuation of the decomposition is explicitly requested}
}

\details{
  This is the main function which does the decomposition of the series
  trajectory matrix. Depending on the SVD method selected in the
  \code{\link{ssa}} different SVD implementations are called. This
  might be the ordinary full SVD routines or fast methods which exploit
  the Hankel / Toeplitz / Hankel with Hankel blocks matrix structure and
  allow the calculation of first few eigentriples.

  Some SVD methods support continuation of the decomposition: if the
  'ssa' object already holds some decomposition and more eigentriples
  are requested, then the decomposition continues using the current
  values as a starting point reducing the computation time dramatically.
}

\value{
  The SSA object.
}

\seealso{
  \code{\link{Rssa}} for an overview of the package, as well as,
  \code{\link[svd:svd]{svd}}, \code{\link[Rssa:ssa]{ssa}}.
}

\examples{
# Decompose 'co2' series with default parameters and decomposition turned off.
s <- ssa(co2, force.decompose = FALSE)
# Perform the decomposition
decompose(s, neig = 50)
# Continue the decomposition
decompose(s, neig = 100)
}
