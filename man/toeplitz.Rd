\name{tmat}
\alias{new.tmat}
\alias{is.tmat}
\alias{tcols}
\alias{trows}
\alias{tmatmul}

\title{Toeplitz matrices operations.}

\description{
  A set of routines to operate on Toeplitz matrices stored in compact
  FFT-based form.
}

\usage{
new.tmat(F, L = (N + 1) \%/\% 2, circular = FALSE, fft.plan = NULL)
is.tmat(t)
tcols(t)
trows(t)
tmatmul(tmat, v, transposed = FALSE)
}

\arguments{
  \item{F}{series to construct the Toeplitz version of L x L autocovariance matrix.}
  \item{fft.plan}{internal hint argument, should be NULL in most cases}
  \item{L}{the window length.}
  \item{circular}{logical vector of one element, describes series topology.
    'TRUE' means series circularity}
  \item{t, tmat}{matrix to operate on.}
  \item{transposed}{logical, if 'TRUE' the multiplication is performed
    with the transposed matrix.}
  \item{v}{vector to multiply with.}
}

\details{
  Fast Fourier Transform provides a very efficient matrix-vector
  multiplication routine for Toeplitz matrices. See the paper in
  'References' for the details of the algorithm.
}

\references{
  Korobeynikov, A. (2010) \emph{Computation- and space-efficient implementation of
    SSA.} Statistics and Its Interface, Vol. 3, No. 3, Pp. 257-268
}

\seealso{
  \code{\link{Rssa}} for an overview of the package, as well as,
  \code{\link[Rssa:ssa]{ssa}},
}

\examples{
# Construct the Toeplitz version of the autocovariance matrix for 'co2' series
h <- new.tmat(co2, L = 10)
# Print the number of columns and rows
print(trows(h)); print(tcols(h))
}
