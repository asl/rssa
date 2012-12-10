\name{clusterify}
\alias{clusterify}
\alias{clusterify.ssa}
\title{Automatic grouping of the elementary series}

\description{Try to perform the grouping step of the SSA via clusterization.}

\usage{
\method{clusterify}{ssa}(x, group, nclust = length(group) / 2,
           \dots,
           type = c("wcor"), cache = TRUE)
}

\arguments{
  \item{x}{SSA object}
  \item{group}{vector of integers, initial set of elementary series numbers}
  \item{nclust}{integer, desired number of output series}
  \item{\dots}{further arguments passed to clusterization methods}
  \item{type}{algorithm of clusterization, see 'Details' for more info}
  \item{cache}{logical, if 'TRUE' then the results of the intermediate
    calculations are saved and reused as much as possible}
}

\value{
  List of integer vectors holding the indices of the elementary series
  forming each grouped series
}

\details{
  Right now the only clusterization method is supported: use of
  w-correlation matrix as the dissimilarity matrix. Standard
  \code{hclust} routine is used to perform the grouping of the
  elementary series.
}

\seealso{
  \code{\link{hclust}}, \code{\link{wcor}}
}

\examples{
# Decompose 'co2' series with default parameters
s <- ssa(co2)
# Form 3 series from the initial 6 ones:
lst <- clusterify(s, group = 1:6, nclust=3)
# Automatic grouping:
print(lst)
# Check separability
w <- wcor(s, groups = lst)
plot(w)
}

