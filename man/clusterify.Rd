\name{grouping.auto.wcor}
\alias{grouping.auto.wcor}
\alias{grouping.auto.wcor.ssa}
\title{Group Elementary Series Using W-correlation Matrix}

\description{Group elemenatry series automatically via
  the hierarchical clustering with w-correlation matrix as a proximity
  matrix}

\usage{
\method{grouping.auto.wcor}{ssa}(x, groups, nclust = length(groups) / 2, \dots)
}

\arguments{
  \item{x}{SSA object}
  \item{groups}{list of numeric vectors, indices of elementary components
    used for reconstruction}
  \item{nclust}{integer, desired number of output series}
  \item{\dots}{further arguments passed to \code{hclust}}
}

\value{
  List of integer vectors holding the indices of the elementary components
  forming each grouped objects
}

\details{
  Standard \code{hclust} routine is used to perform the grouping
  of the elementary components.
}

\seealso{
  \code{\link{hclust}}, \code{\link{wcor}}
}

\examples{
# Decompose 'co2' series with default parameters
s <- ssa(co2)
# Form 3 series from the initial 6 ones:
lst <- grouping.auto(s, grouping.method = "wcor",
                     groups = 1:6, nclust=3)
# Automatic grouping:
print(lst)
plot(lst)
# Check separability
w <- wcor(s, groups = lst)
plot(w)
}

