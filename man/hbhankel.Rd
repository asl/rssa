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
new.hbhmat(F, L = (N + 1) \%/\% 2,
           wmask = NULL, fmask = NULL, weights = NULL,
           circular = FALSE)
is.hbhmat(h)
hbhcols(h)
hbhrows(h)
hbhmatmul(hmat, v, transposed = FALSE)
}

\arguments{
  \item{F}{array to construct the trajectory matrix for.}
  \item{L}{the window length.}
  \item{wmask, fmask, weights}{special parameters for shaped SSA case (see \code{\link[Rssa:ssa]{ssa}}).
    \code{wmask} and \code{fmask} are logical matrices, window and factor masks respectively.
    \code{weights} is integer matrix which denotes hankel weights for array elements. If 'NULL',
    parameters for simple rectangular 2D SSA case are used.}
  \item{circular}{logical vector of one or two elements, describes field topology.
    'TRUE' means circularity by a corresponding coordinate. If vector has only one element,
    this element will be used twice.}
  \item{h, hmat}{matrix to operate on.}
  \item{transposed}{logical, if 'TRUE' the multiplication is performed
    with the transposed matrix.}
  \item{v}{vector to multiply with.}
}

\details{
  Fast Fourier Transform provides a very efficient matrix-vector
  multiplication routine for Hankel with Hankel blocks matrices. See the
  paper in 'References' for the details of the algorithm.
}

\references{
  Korobeynikov, A. (2010) \emph{Computation- and space-efficient implementation of
    SSA.} Statistics and Its Interface, Vol. 3, No. 3, Pp. 257-268
}

\author{Konstantin Usevich}
