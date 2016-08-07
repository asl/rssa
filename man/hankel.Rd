\name{hmat}
\alias{new.hmat}
\alias{is.hmat}
\alias{hcols}
\alias{hrows}
\alias{hankel}
\alias{hmatmul}

\title{Hankel matrices operations.}

\description{
  A set of routines to operate on Hankel matrices stored in
  compact FFT-based form.
}

\usage{
new.hmat(F, L = (N + 1)\%/\%2, circular = FALSE, wmask = NULL,
         fmask = NULL, weights = NULL, fft.plan = NULL)
is.hmat(h)
hcols(h)
hrows(h)
hmatmul(hmat, v, transposed = FALSE)
hankel(X, L)
}

\arguments{
  \item{F}{series to construct the trajectory matrix for.}
  \item{fft.plan}{internal hint argument, should be NULL in most cases}
  \item{wmask, fmask, weights}{special parameters for shaped SSA case (see \code{\link[Rssa:ssa]{ssa}}).
    \code{wmask} and \code{fmask} are logical vectors, window and factor masks respectively.
    \code{weights} is integer vector which denotes hankel weights for array elements. If 'NULL',
    parameters for simple 1D SSA case are used.}
  \item{circular}{logical vector of one element, describes series topology.
    'TRUE' means circularity by time.}
  \item{L}{the window length.}
  \item{h, hmat}{matrix to operate on.}
  \item{transposed}{logical, if 'TRUE' the multiplication is performed
    with the transposed matrix.}
  \item{v}{vector to multiply with.}
  \item{X}{series to construct the trajectory matrix for or matrix for hankelization}
}

\details{
  Fast Fourier Transform provides a very efficient matrix-vector
  multiplication routine for Hankel matrices. See the paper in
  'References' for the details of the algorithm.
}

\references{
  Korobeynikov, A. (2010) \emph{Computation- and space-efficient implementation of
    SSA.} Statistics and Its Interface, Vol. 3, No. 3, Pp. 257-268
}

\seealso{
  \code{\link{Rssa}} for an overview of the package, as well as,
  \code{\link[Rssa:ssa]{ssa}},
  \code{\link[Rssa:decompose.ssa]{decompose}},
}

\examples{
# Construct the Hankel trajectory matrix for 'co2' series
h <- new.hmat(co2, L = 10)
# Print number of columns and rows
print(hrows(h))
print(hcols(h))
}
