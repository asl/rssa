\name{clone}
\alias{clone}
\alias{clone.ssa}

\title{
  Cloning of SSA objects
}
\description{
  Function to copy SSA objects
}

\usage{
\method{clone}{ssa}(this, copy.cache = TRUE, \dots)
}

\arguments{
  \item{this}{object to be cloned}
  \item{copy.cache}{enable/disable copying of the set of pre-cached subseries.}
  \item{\dots}{additional arguments to 'clone'.}
}
\details{
  For the sake of memory efficiency SSA objects hold references to the
  data, not the data itself. That is why they can hold huge amount of
  data and passing them by value is still cheap.

  However, this means that one cannot safely copy the object using
  normal assignment operator, since freeing of references in one object
  would yield stale references in another. The `clone' method provides
  safe `deep copy' of SSA objects.
}

\examples{
# Decompose 'co2' series with default parameters
s <- new.ssa(co2);
# Perform 'normal copy' of SSA object
s1 <- s;
# Perform 'deep copy' of SSA object
s2 <- clone(s);
# Add some data to 's'
reconstruct(s);
# Now 's1' also contains this data, but 's2' - not
summary(s1);
summary(s2);
}
