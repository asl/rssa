\name{wcor}
\alias{wcor}
\alias{wcor.ssa}
\alias{wcor.default}

\title{Calculate the W-correlation matrix}
\description{
  Function calculates the W-correlation matrix for the series.
}
\usage{
\method{wcor}{ssa}(X, groups, \dots, cache = TRUE)
\method{wcor}{default}(X, L = (N + 1) \%/\% 2, \dots)
}

\arguments{
  \item{X}{the input object. This might be ssa object for \emph{ssa}
    method, or just a matrix with elementary series in columns for
  \emph{default} implementation.}
  \item{L}{window length.}
  \item{groups}{list of numeric vectors, indices of elementary subseries
    used for reconstruction.}
  \item{\dots}{further arguments passed to \code{reconstruct} routine.}
  \item{cache}{logical, if 'TRUE' then intermediate results will be
    cached in 'ssa' object.}
}
\details{
  W-correlation matrix is a standard way of checking for weak
  separability between the elementary time series (see Golyandina et
  all, 2001 for more information). The function calculates such a matrix
  either directly from 'ssa' object or from the matrix of elementary series.
}

\value{
  Object of type 'wcor.matrix'.
}

\references{
  Golyandina, N., Nekrutkin, V. and Zhigljavsky, A. (2001): \emph{Analysis of
    Time Series Structure: SSA and related techniques.} Chapman and Hall/CRC. ISBN 1584881941
}
\seealso{
  \code{\link[Rssa:reconstruct.ssa]{reconstruct}}
}

\examples{
# Decompose co2 series with default parameters
s <- new.ssa(co2)
# Calculate the w-correlation matrix between first 10 subseries
w <- wcor(s, groups = 1:10)
print(w)
plot(w)
}

