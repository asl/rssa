\name{ssa.capabilities}
\alias{ssa.capabilities}
\title{
  SSA methods and capabilities check
}
\description{
  Not all SSA algorithms and methods could be applied to SSA objects of
  any kind (e.g. gapfilling requires shaped SSA object, one cannot
  forecast for 3D-SSA and so on). This function allows one to determine
  a set of methods allowed to be applied to a particular SSA object
}
\usage{
ssa.capabilities(x)
}

\arguments{
  \item{x}{SSA object holding the decomposition}
}

\value{
  Logical vector, indicating which methods are allowed
}

\examples{
# Decompose 'co2' series with default parameters
s <- ssa(co2)
# Since this is 1d SSA object, everything should be supported except
# gapfilling
print(ssa.capabilities(s))
}
