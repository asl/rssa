\name{cleanup}
\alias{cleanup}
\alias{cleanup.ssa}

\title{
  Cleanup of all cached data from SSA objects
}
\description{
  Function to copy SSA objects
}

\usage{
cleanup(x)
}

\arguments{
  \item{x}{object to be cleaned}
}

\details{
  For the sake of memory efficiency SSA objects hold references to the
  data, not the data itself. That is why they can hold huge amount of
  data and passing them by value is still cheap.

  Also, SSA routines tend to save some intermediate information which
  can be used later inside SSA object. This includes (but not limited
  to) elementary series, etc.

  \code{cleanup} call deletes all pre-cached stuff freeing memory necessary
  for calculations.
}
