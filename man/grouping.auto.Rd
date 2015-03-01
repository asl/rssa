\name{grouping.auto}
\alias{grouping.auto}
\title{Group Elementary Series}

\description{The `grouping.auto' function performs the Grouping Step of
  SSA using different approaches.}

\usage{
grouping.auto(x, \dots, grouping.method = c("pgram", "wcor"))
}

\arguments{
  \item{x}{SSA object}
  \item{grouping.method}{String specifying the method used to perform
    the grouping. Allowed methods are `"pgram"' (the default) and
    `"wcor"'}
  \item{\dots}{Further arguments to specific methods}
}

\value{
  List of integer vectors holding the indices of the elementary components
  forming each grouped objects.
}

\details{
 `grouping.auto' is a wrapper function which calls the methods
 `grouping.auto.pgram' and `grouping.auto.wcor'.
}

\seealso{
  \code{\link{grouping.auto.pgram}}, \code{\link{grouping.auto.wcor}}
}
