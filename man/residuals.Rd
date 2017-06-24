\name{residuals}
\alias{residuals.ssa}
\alias{residuals.ssa.reconstruction}
\title{Obtain the residuals from SSA reconstruction}
\description{
  Obtain the residuals from SSA reconstruction
}

\usage{
\method{residuals}{ssa}(object, groups, \dots, cache = TRUE)
\method{residuals}{ssa.reconstruction}(object, \dots)
}

\arguments{
  \item{object}{input object}
  \item{groups}{list of numeric vectors, indices of elementary components
    used for reconstruction, the entries of the list can be named.}
  \item{\dots}{further arguments passed to \code{reconstruct} routine}
  \item{cache}{logical, if 'TRUE' then intermediate results will be
    cached in the SSA object.}
}

\details{
  This function calculates the residuals either from SSA object
  corresponding to reconstruction using \code{groups} arguments, or just
  extracts the residuals from reconstruction object.
}

\value{
  residuals object
}

\seealso{
  \code{\link{Rssa}} for an overview of the package, as well as,
  \code{\link[Rssa:reconstruct]{reconstruct}}.
}

\examples{
# Decompose 'co2' series with default parameters
s <- ssa(co2)
# Reconstruct the series, grouping elementary series.
r <- reconstruct(s, groups = list(c(1, 4), c(2,3), c(5, 6)))
print(residuals(r))

# If there are several groups, then the residuals are calculated as
# residuals for the model corresponding to the combined model.
r <- reconstruct(s, groups = list(c(6, 7), c(6,7), c(8, 9)))
r1 <- reconstruct(s, groups = list(6:9))
max(abs(residuals(r) - residuals(r1))) # 0 
max(abs(co2 - (r1$F1 + residuals(r1)))) # 0
}
