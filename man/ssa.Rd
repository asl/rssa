\name{ssa}
\alias{nlambda}
\alias{nu}
\alias{nv}
\title{Properties of SSA object}

\description{
  Functions to query for number of singular values, eigenvectors and
  factor vectors stored in the SSA object.
}

\usage{
nlambda(this, ...)
nu(this, ...)
nv(this, ...)
}

\arguments{
  \item{this}{SSA object to query.}
  \item{\dots}{further arguments}
}

\details{
  The numbers returned solely depend on used singular value
  decomposition method and parameters of this method (e.g. 'neig'
  argument for 'propack' and 'nutrlan' SVD methods).
}

\value{
  an 'integer' of length 1.
}

\examples{
# Decompose 'co2' series with default parameters
s <- new.ssa(co2, neig = 20)
# Show the number of eigentriples saved in the 's'
print(nlambda(s))
}

\keyword{array}
\keyword{algebra}


