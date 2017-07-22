\name{ssa-object}
\alias{ssa-object}
\alias{nlambda}
\alias{nsigma}
\alias{nu}
\alias{nv}
\alias{nspecial}
\alias{nspecial.ssa}
\alias{summary.ssa}
\alias{$.ssa}
\alias{contributions}

\title{Properties of SSA object}

\description{
  Functions to access various fields of SSA object, query
  for number of singular values, eigenvectors, factor vectors and `special' decomposition
  triples (now, ProjectionSSA triples) in the SSA object and other miscellaneous info.
}

\usage{
nsigma(x)
nu(x)
nv(x)
\method{nspecial}{ssa}(x)
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
    \item{sigma}{Vector of singular values}
    \item{U}{The matrix of eigenvectors}
    \item{V}{The matrix of factor vectors. Might not exist depending on the selected SVD method}
  }

  If SSA with projections is being performed, then the eigentriples are ordered as follows:
  at first, row projection triples, then column projection triples and then
  SVD-triples. Non-SVD triples (like projection triples) are called `special triples'.
  The number of special triples can be obtained by calling \code{nspecial} method.
  Also, one can use the following fields of the SSA object:
  \describe{
    \item{nPR}{the number of row projection triples, may be NULL}
    \item{nPL}{the number of column projection triples, may be NULL}
  }
}

\section{About decompositions}{
  The result of Decomposition step of SSA and its modifications can be written down in the following form:
  \deqn{%
    {(*)} \quad \mathbf{X} = \sum_i \mathbf{X}_i, \qquad \mathbf{X}_i = \sigma_i U_i V_i^\mathrm{T},}{%
    (*) X = sum_i X_i, X_i = sigma_i U_i (V_i)^T,}
  where \eqn{\mathbf{X}}{X} is the trajectory matrix,
  \eqn{U_i \in R^L}{U_i in R^L}, \eqn{V_i \in R^K}{V_i in R^K}, \eqn{\sigma_i} are non-negative numbers.
  Also, we assume that \eqn{\|U_i\| = 1}{||U_i|| = 1}, \eqn{\|V_i\| = 1}{||V_i|| = 1}.

  The Singular Value Decomposition is a particular case of \eqn{(*)} and corresponds to orthonormal systems of
  \eqn{\{U_i\}}{{U_i}} and \eqn{\{V_i\}}{{V_i}}.
  We call \eqn{(\sigma_i, U_i, V_i)} eigentriple, \eqn{\sigma_i}
  are singular values, \eqn{U_i} are left singular values or eigenvectors,
  \eqn{V_i} are right singular vectors or factor vectors, by analogy with the SVD.

  For the most of SSA decompositions, \eqn{U_i} belongs to the column space of \eqn{\mathbf{X}}{X},
  while \eqn{V_i} belongs to the row space of \eqn{\mathbf{X}}{X}.
  Therefore, let us consider such decompositions called \emph{consistent.}

  Note that \eqn{(*)} is a decomposition of \eqn{\mathbf{X}}{X} into a sum of rank-one matrices.
  If the systems \eqn{\{U_i\}}{{U_i}} and \eqn{\{V_i\}}{{V_i}} are linearly-independent,
  then the decomposition \eqn{(*)} is minimal (has minimal possible number of addends).

  If at least one of the systems is not linear independent,
  the decomposition \eqn{(*)} is not minimal.
  If both \eqn{\{U_i\}}{{U_i}} and \eqn{\{V_i\}}{{V_i}} are orthonormal,
  then the decomposition \eqn{(*)} is called bi-orthogonal.
  If \eqn{\{U_i\}}{{U_i}} is orthonormal, the decomposition is called left-orthogonal;
  If \eqn{\{V_i\}}{{V_i}} is orthonormal, the decomposition is called right-orthogonal.

  Let \eqn{r} be rank of \eqn{\mathbf{X}}{X}. Minimal decomposition has exactly \eqn{r} addends.
  Introduce the Frobenius-inner product as
  \eqn{\langle \mathbf{Z}, \mathbf{Y} \rangle _\mathrm{F} = \sum_{i, j} z_{i, j} \cdot y_{i, j}}{<Z, Y>_F = sum_{i, j} (z_{i, j} y_{i, j})}.
  Thus, we can say about F-orthogonality and F-orthogonal decompositions if \eqn{\mathbf{X}_i}{X_i}
  are F-orthogonal. For F-orthogonality, left or right orthogonality is sufficient.

  Generally, \eqn{\|\mathbf{X}\|^2}{||X||^2} can be not equal to \eqn{\sum_i \|\mathbf{X}_i\|^2}{sum_i ||X_i||^2}.
  For F-orthogonal decompositions, \eqn{\|\mathbf{X}\|^2 = \sum_i \|\mathbf{X}_i\|^2}{||X||^2 = sum_i ||X_i||^2}.

  The contribution of k-th matrix component is defined as
  \eqn{{\|\mathbf{X}_k\|^2} / {\|\mathbf{X}\|^2} = {\sigma_k^2} / (\sum_i \sigma_i^2)}{||X_k||^2 / ||X||^2 = sigma_k^2 / (sum_i sigma_i^2)}.

  For F-orthogonal decompositions, the sum of component contributions is equal to 1.
  Otherwise, this sum can considerably differ from 1
  (e.g., the sum of component contributions can be 90\% or 146\%).

  \strong{Remark.} If the system \eqn{\{U_i\}}{{U_i}} (or \eqn{\{V_i\}}{{V_i}}) has vectors
  that do not belong to the column (or row) spaces, then the decomposition can be not minimal even
  if \eqn{\{U_i\}}{{U_i}} (or \eqn{\{V_i\}}{{V_i}}) are linearly independent,
  since these projections on the column (or row) space can be dependent.
}

\section{Decompositions for different SSA modifications}{
  \describe{
    \item{Basic SSA}{the SVD, consistent, minimal, bi-orthogonal and
      therefore F-orthogonal decomposition. Implemented in
      \code{\link[Rssa:ssa]{ssa}} with \code{kind='1d-ssa'}}
    \item{FOSSA}{consistent, minimal F-orthogonal
      decomposition. Implemented in \code{\link[Rssa:fossa]{fossa}}}
    \item{IOSSA}{consistent, minimal oblique decomposition. Implemented
      in \code{\link[Rssa:iossa]{iossa}}}
    \item{SSA with projections}{non-consistent if at least one basis
      vector used for the projection does not belong to the column (row)
      trajectory space, F-orthogonal decomposition. The components,
      which are obtained by projections, are located at the beginning of
      the decomposition and have numbers \eqn{1, \dots,
        n_\mathrm{special}}{1, \dots,
        n_special}. Implemented in \code{\link[Rssa:ssa]{ssa}} with
      \code{kind='1d-ssa'} and non-\code{NULL} \code{row.projector} or
      \code{column.projector} arguments}
    \item{Toeplitz SSA}{generally, non-consistent, non-minimal
      F-orthogonal decomposition. Implemented in
      \code{\link[Rssa:ssa]{ssa}} with \code{kind='toeplitz-ssa'}}
  }
}

\value{
  an 'integer' of length 1 for \code{nu}, \code{nv}, \code{nsigma},
  \code{nspecial} routines, matrix or vector for \code{$} operator.
}

\note{
  For \code{nsigma}, \code{nu}, \code{nv}, \code{$}  routines, the
  values returned solely depend on used singular value decomposition
  method and parameters of this method (e.g. 'neig' argument for
  'propack' and 'nutrlan' SVD methods).
}

\seealso{
  \code{\link{Rssa}} for an overview of the package, as well as,
  \code{\link[Rssa:ssa]{ssa}},
  \code{\link[Rssa:calc.v]{calc.v}},
  \code{\link[Rssa:iossa]{iossa}},
  \code{\link[Rssa:fossa]{fossa}},
}

\examples{
# Decompose 'co2' series with default parameters
s <- ssa(co2, neig = 20)
# Show the number of eigentriples saved in the 's'
print(nsigma(s))
# Show the summary 
summary(s)
# Show the set of singular values
print(s$sigma)
# Show the first eigenvector
print(s$U[, 1])

# Decompose `co2' series with polynomial projections
s <- ssa(co2, row.projector = 1, column.projector = 2)
print(nspecial(s))
print(c(s$nPL, s$nPR))
# Reconstruct a polynomial trend
plot(reconstruct(s, groups = list(trend = seq_len(nspecial(s)))))
}
