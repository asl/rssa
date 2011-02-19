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
new.hmat(F, L = (N - 1) \%/\% 2)
is.hmat(h)
hcols(h)
hrows(h)
hmatmul(hmat, v, transposed = FALSE)
}

\arguments{
  \item{F}{series to construct the trajectory matrix for.}
  \item{L}{the window length.}
  \item{h, hmat}{matrix to operate on.}
  \item{transposed}{logical, if 'TRUE' the multiplication is performed
    with the transposed matrix.}
  \item{v}{vector to multiply with.}
}

\details{
  Fast Fourier Transforms provides a very efficient matrix-vector
  multiplication routine for Hankel matrices. See the paper in
  'References' for the details of the algorithm.
}

\references{
  Korobeynikov, A. (2010) \emph{Computation- and space-efficient implementation of
    SSA.} Statistics and Its Interface, Vol. 3, No. 3, Pp. 257-268
}

\examples{
# Construct the Hankel trajectory matrix for 'co2' series
h <- new.hmat(co2, 10)
# Print number of columns and rows
print(hrows(h))
print(hcols(h))
}
