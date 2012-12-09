\name{hmatr}
\alias{hmatr}
\alias{plot.hmatr}

\title{Calculate the heterogeneity matrix.}
\description{
  Function calculates the heterogeneity matrix for the series.
}
\usage{hmatr(F, \dots,
      B = N \%/\% 4, T = N \%/\% 4, L = B \%/\% 2,
      neig = 10)

\method{plot}{hmatr}(x,
     col = rev(heat.colors(256)),
     main = "Heterogeneity Matrix", xlab = "", ylab = "", \dots)
}

\arguments{
  \item{F}{the series to be checked for structural changes}
  \item{\dots}{further arguments passed to \code{ssa} routine for
  \code{hmatr} call or \code{image} for \code{plot.hmatr} call}
  \item{B}{integer, length of base series}
  \item{T}{integer, length of tested series}
  \item{L}{integer, window length for the decomposition of the base
    series}
  \item{neig}{integer, number of eigentriples to consider for
    calculating projections}
  \item{x}{'hmatr' object}
  \item{col}{color palette to use}
  \item{main}{plot title}
  \item{xlab,ylab}{labels for 'x' and 'y' axis}
}
\details{
  The \emph{heterogeneity matrix} (H-matrix) provides a
  consistent view on the structural discrepancy between different parts of the
  series. Denote by \eqn{F_{i,j}} the subseries of F of the form: \eqn{F_{i,j} =
  \left(f_{i},\dots,f_{j}\right)}. Fix two integers \eqn{B > L} and \eqn{T \geq L}. Let
  these integers denote the lengths of \emph{base} and \emph{test} subseries,
  respectively. Introduce the H-matrix \eqn{G_{B,T}} with the elements \eqn{g_{ij}} as
  follows:
  \deqn{
    g_{ij} = g(F_{i,i+B}, F_{j,j+T}),
  }
  for \eqn{i=1,\dots,N-B+1} and \eqn{j=1,\dots,N-T+1}, that is we split the series
  F into subseries of lengths B and T and calculate the heterogeneity index
  between all possible pairs of the subseries.

  The heterogeneity index \eqn{g(F^{(1)}, F^{(2)})} between the series
  \eqn{F^{(1)}} and \eqn{F^{(2)}} can be calculated as follows: let
  \eqn{U_{j}^{(1)}}, \eqn{j=1,\dots,L} denote the eigenvectors of the
  SVD of the trajectory matrix of the series \eqn{F^{(1)}}. Fix I to be a
  subset of \eqn{\left\{1,\dots,L\right\}} and denote \eqn{\mathcal{L}^{(1)} =
  \mathrm{span}\,\left(U_{i},\, i \in I\right)}. Denote by
  \eqn{X^{(2)}_{1},\dots,X^{(2)}_{K_{2}}} (\eqn{K_{2} = N_{2} - L + 1}) the
  L-lagged vectors of the series \eqn{F^{(2)}}. Now define
  \deqn{
    g(F^{(1)},F^{(2)})
    = \frac{\sum_{j=1}^{K_{2}}{\mathrm{dist}\,^{2}\left(X^{(2)}_{j},
        \mathcal{L}^{(1)}\right)}}
    {\sum_{j=1}^{K_{2}}{\left\|X^{(2)}_{j}\right\|^{2}}}, } where
  \eqn{\mathrm{dist}\,(X,\mathcal{L})} denotes the Euclidean distance between the
  vector X and the subspace \eqn{\mathcal{L}}. One can easily see that
  \eqn{0 \leq g \leq 1}.
}

\value{
  object of type 'hmatr'
}

\references{
  Golyandina, N., Nekrutkin, V. and Zhigljavsky, A. (2001): \emph{Analysis of
    Time Series Structure: SSA and related techniques.} Chapman and Hall/CRC. ISBN 1584881941
}
\seealso{
  \code{\link[Rssa:ssa]{ssa}}
}

\examples{
# Calculate H-matrix for co2 series
h <- hmatr(co2, L = 24)
# Plot the matrix
plot(h)
}

