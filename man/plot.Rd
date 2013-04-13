\name{plot}
\alias{plot.ssa}

\title{Plot SSA object}
\description{
  This function plots various sorts of figures related to the SSA method.
}
\usage{
\method{plot}{ssa}(x,
     type = c("values", "vectors", "paired", "series", "wcor"),
     \dots,
     plot.contrib = TRUE,
     numvalues = nlambda(x),
     numvectors = min(nlambda(x), 10),
     idx = 1:numvectors,
     idy,
     groups)
}

\arguments{
  \item{x}{SSA object holding the decomposition}
  \item{type}{Type of the plot (see 'Details' for more information)}
  \item{\dots}{Arguments to be passed to methods, such as graphical parameters}
  \item{plot.contrib}{logical. If 'TRUE' (the default), the contribution
    of the component to the total variance is plotted}
  \item{numvalues}{Number of eigenvalues to plot (for type = 'values')}
  \item{numvectors}{Total number of eigenvectors to plot (for type = 'vectors')}
  \item{idx}{Indices of eigenvectors to plot (for type = 'vectors')}
  \item{idy}{Second set of indices of eigenvectors to plot (for type = 'paired')}
  \item{groups}{Grouping used for the decomposition (see \code{\link[Rssa:reconstruct.ssa]{reconstruct}})}
}

\details{
  This function is the single entry to various plots of SSA objects. Right
  now this includes:
  \describe{
    \item{values}{plot the graph of eigenvalues.}
    \item{vectors}{plot the eigenvectors.}
    \item{paired}{plot the pairs of eigenvectors (useful for the
      detection of periodic components).}
    \item{series}{plot the reconstructed series.}
    \item{wcor}{plot the W-correlation matrix for the reconstructed series.}
  }

  Additional (non-standard) graphical parameters can be transfered via \dots:
  \describe{
    \item{plot.type}{just lattice plot type. This parameter is introduced in order to
      solve conflict with argument \code{type}.}
  }
}

\examples{
# Decompose 'co2' series with default parameters
s <- ssa(co2)
# Plot the eigenvalues
plot(s, type = "values")
# Plot W-cor matrix for first 10 reconstructed components
plot(s, type = "wcor", groups = 1:10)
# Plot the paired plot for first 6 eigenvectors
plot(s, type = "paired", idx = 1:6)
# Plot the first 4 reconstructed components
plot(s, type = "series", groups = list(1:4))
}
