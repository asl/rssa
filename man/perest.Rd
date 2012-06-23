\name{perestimate}
\alias{perestimate}
\alias{perestimate.ssa}
\title{Estimate periods from (set of) eigenvectors}

\description{
  Function to estimate sine periods given a set of SSA eigenvectors.
}

\usage{
\method{perestimate}{ssa}(this, group, \dots, method = c("pairs", "esprit-ls"))
}

\arguments{
  \item{this}{SSA object}
  \item{group}{indices of eigenvectors to estimate from}
  \item{\dots}{further arguments passed to 'decompose' routine, if
    necessary}
  \item{method}{estimation method, 'pairs' for rough estimation based on
    pair of eigenvectors, 'esprit-ls' for least-squares ESPRIT}
}

\value{
  a 'list' with 2 components named 'periods' containing the periods and
  'modules' containing the modules of eigen roots (may not be present).
}

\references{
  Golyandina, N., Nekrutkin, V. and Zhigljavsky, A. (2001): \emph{Analysis of
    Time Series Structure: SSA and related techniques.} Chapman and Hall/CRC. ISBN 1584881941
}

\examples{
# Decompose 'co2' series with default parameters
s <- new.ssa(co2, neig = 20)
# Estimate the periods from 2nd and 3rd eigenvectors using default
# 'pairs' method
print(perestimate(s, group = c(2, 3)))
# Estimate the peroids from 2nd and 3rd eigenvectors using ESPRIT
print(perestimate(s, group = c(2, 3), method = "esprit-ls"))
}

