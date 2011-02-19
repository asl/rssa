\name{hbhmat}
\alias{new.hbhmat}
\alias{is.hbhmat}
\alias{hbhcols}
\alias{hbhrows}
\alias{hbhankel}
\alias{hbhmatmul}
\title{Hankel with Hankel block matrices operations.}

\description{
  A set of routines to operate on Hankel with Hankel block matrices
  stored in compact FFT-based form.
}

\usage{
new.hbhmat(F, L = (N - 1) \%/\% 2)
is.hbhmat(h)
hbhcols(h)
hbhrows(h)
hbhmatmul(hmat, v, transposed = FALSE)
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
  multiplication routine for Hankel with Hankel blocks matrices. See the
  paper in 'References' for the details of the algorithm.
}

\references{
  Korobeynikov, A. (2010) \emph{Computation- and space-efficient implementation of
    SSA.} Statistics and Its Interface, Vol. 3, No. 3, Pp. 257-268
}

\author{Konstantin Usevich}
