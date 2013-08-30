\name{ssa-object}
\alias{ssa-object}
\alias{nlambda}
\alias{nu}
\alias{nv}
\alias{summary.ssa}
\alias{$.ssa}

\title{Properties of SSA object}

\description{
  Functions to access various fields of SSA object, query
  for number of singular values, eigenvectors and factor vectors stored
  in the SSA object and other miscellaneous info.
}

\usage{
nlambda(x)
nu(x)
nv(x)
\method{summary}{ssa}(object, digits = max(3, getOption("digits") - 3), \dots)
\method{$}{ssa}(x, name)
}

\arguments{
  \item{x}{SSA object to query}
  \item{object}{an object for which a summary is desired}
  \item{digits}{integer, used for number formatting}
  \item{\dots}{additional arguments affecting the summary produced}
  \item{name}{field of SSA object to extract. See 'Details' for list of
  the fields}
}

\details{
  The internals of SSA object is inheritely opaque, because depending on
  the selected SVD method and other conditions it might contains
  different fields.

  However, it is possible to extract some fields out of it using the
  operator \code{$}, in particular, the following values of argument
  \code{name} are supported:
  \describe{
    \item{lambda}{Vector of singular values}
    \item{U}{The matrix of eigenvectors}
    \item{V}{The matrix of factor vectors. Might not exist depending on the selected SVD method}
  }
}

\value{
  an 'integer' of length 1 for \code{nu}, \code{nv}, \code{nlambda}
  routines, matrix or vector for \code{$} operator.
}

\note{
  For \code{nlambda}, \code{nu}, \code{nv}, \code{$}  routines, the
  values returned solely depend on used singular value decomposition
  method and parameters of this method (e.g. 'neig' argument for
  'propack' and 'nutrlan' SVD methods).
}

\seealso{
  \code{\link{Rssa}} for an overview of the package, as well as,
  \code{\link[Rssa:ssa]{ssa}},
  \code{\link[Rssa:calc.v]{calc.v}}
}

\examples{
# Decompose 'co2' series with default parameters
s <- ssa(co2, neig = 20)
# Show the number of eigentriples saved in the 's'
print(nlambda(s))
# Show the summary 
summary(s)
# Show the set of singular values
print(s$lambda)
# Show the first eigenvector
print(s$U[, 1])
}
