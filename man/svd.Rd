\name{svd}
\alias{trlan.svd}
\alias{propack.svd}
\alias{trlan.eigen}

\title{
  Generic Singular Value Decomposition of a Matrix
}

\description{
  Compute the singular-value decomposition of a real rectangular matrix.
}

\usage{
propack.svd(X, neig = min(m, n), opts = list())
trlan.svd(X, neig = min(m, n), opts = list(), lambda = NULL, U = NULL)
}

\arguments{
  \item{X}{the matrix to be decomposed. This can be either normal matrix
  or 'external matrix' object (e.g. one, created via 'new.hmat' function).}
  \item{neig}{number of desired eigentriples}
  \item{opts}{different options for eigensolver. See 'Details' section
  for more information}
  \item{lambda}{set of already computed singular values (used for
  continuation of the decomposition).}
  \item{U}{matrix of already computed eigenvectors (used for
  continuation of the decomposition).}
}

\details{
  These routines provides an interface to two state-of-art
  implementations of truncated SVD.

  PROPACK does this via the implicitly restarted Lanczos
  bidiagonalization with partial reorthogonalization. nu-TRLAN does the
  thick-restart Lanczos eigendecomposition of cross-product matrix.
  
  'opts' is a list of different options which can be passed to the
  routines. Note that by default more or less suitable values for these
  options are set by the routines automatically.

  The options for PROPACK are:
  \describe{
    \item{kmax}{integer, maximum number of iterations.}
    \item{dim}{integer, dimension of Krylov subspace.}
    \item{p}{integer, number of shits per restart.}
    \item{maxiter}{integer. maximum number of restarts.}
    \item{tol}{numeric, tolerance level.}
    \item{verbose}{logical, if 'TRUE', provide verbose output.}
  }
  
  The options for nu-TRLAN are:
  \describe{
    \item{kmax}{integer, maximum number of iterations.}
    \item{maxiter}{integer. maximum number of matrix-vector products.}
    \item{tol}{numeric, tolerance level.}
    \item{verbose}{integer, verboseness level.}
  }
}

\value{
  The returned value is a list with components
  \item{d}{a vector containing the singular values of 'x'}
  \item{u}{a matrix whose columns contain the left singular vectors of
    'X'}
  \item{v}{a matrix whose columns contain the right singular vectors of
    'X' (only for 'propack.svd')}
}

\references{
  Wu, K. and Simon, H. (2000). \emph{Thick-restart Lanczos method for
  large symmetric eigenvalue problems.} SIAM J. Matrix Anal. Appl. 22, 2, 602-616.

  Yamazaki, I., Bai, Z., Simon, H., Wang, L.-W., and Wu,
  K. (2008). \emph{Adaptive projection subspace dimension for the thick
  restart Lanczos method.} Tech. rep., Lawrence Berkeley National
  Laboratory, University of California, One Cyclotron road, Berkeley, California 94720.  
  
  Larsen, R. M. (1998). \emph{Efficient algorithms for helioseismic inversion.}
  Ph.D. thesis, University of Aarhus, Denmark.

  Korobeynikov, A. (2010) \emph{Computation- and space-efficient implementation of
  SSA.} Statistics and Its Interface, Vol. 3, No. 3, Pp. 257-268
}

\keyword{array}
\keyword{algebra}
