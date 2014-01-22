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
     vectors = c("eigen", "factor"),
     plot.contrib = TRUE,
     numvalues = nsigma(x),
     numvectors = min(nsigma(x), 10),
     idx = 1:numvectors,
     idy,
     groups)
}

\arguments{
  \item{x}{SSA object holding the decomposition}
  \item{type}{Type of the plot (see 'Details' for more information)}
  \item{\dots}{Arguments to be passed to methods, such as graphical
    parameters}
  \item{vectors}{For type = 'vectors', choose the vectors to plot}
  \item{plot.contrib}{logical. If 'TRUE' (the default), the contribution
    of the component to the total variance is plotted.
    For `ossa' class, Frobenius orthogonality checking of elementary matrices is performed.
    If not all matrices are orthogonal, corresponding warning is risen}
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
    \item{values}{plot the graph of the singular values.}
    \item{vectors}{plot the eigenvectors.}
    \item{paired}{plot the pairs of eigenvectors (useful for the
      detection of periodic components).}
    \item{series}{plot the reconstructed series.}
    \item{wcor}{plot the W-correlation matrix for the reconstructed objects.}
  }

  Additional (non-standard) graphical parameters which can be transfered via \dots:
  \describe{
    \item{plot.type}{lattice plot type. This argument will be transfered as \code{type}
    argument to function \code{panel.xyplot}.}
    \item{ref}{logical. Whether to plot zero-level lines in series-plot, eigenvectors-plot and paired-plot.
      Zero-level isolines will be plotted for 2d-eigenvectors-plot.}
    \item{symmetric}{logical. Whether to use symmetric scales in series-plot, eigenvectors-plot and paired-plot.}
    \item{useRaster}{logical. For 2d-eigenvector-plot and wcor-plot, indicating whether raster representations
      should be used. 'TRUE' by default.}
    \item{at}{for 2d-eigenvectors-plot, a numeric vector giving breakpoints along the range of \code{z}.
      Also can be a character string. Values 'free' and 'same' are allowed. 'free' means
      special breakpoints' vectors (will be evaluated automatically, see description of \code{cuts}
      argument in 'Details') for each component. 'same' means one breakpoints' vector for all
      component (will be evaluated automatically too).}
    \item{cuts}{for 2d-reconstruction-plot, the number of levels the range of \code{z} would be divided into.}
    \item{fill.color}{color or 'NULL'. Defines background color for shaped 2d-eigenvectors plot. If 'NULL', standard white
      background will be used.}
  }
}

\examples{
\donttest{
# Decompose 'co2' series with default parameters
s <- ssa(co2)
# Plot the eigenvalues
plot(s, type = "values")
# Plot W-cor matrix for first 10 reconstructed components
plot(s, type = "wcor", groups = 1:10)
# Plot the paired plot for first 6 eigenvectors
plot(s, type = "paired", idx = 1:6)
# Plot eigenvectors for first 6 components
plot(s, type = "vectors", idx = 1:6)
# Plot the first 4 reconstructed components
plot(s, type = "series", groups = list(1:4))
# Plot the eigenvalues by points only
plot(s, type = "values", plot.type = "p")

# Artificial image for 2dSSA
mx <- outer(1:50, 1:50,
            function(i, j) sin(2*pi * i/17) * cos(2*pi * j/7) + exp(i/25 - j/20)) +
      rnorm(50^2, sd = 0.1)
# Decompose 'mx' with default parameters
s <- ssa(mx, kind = "2d-ssa")
# Plot the eigenvalues
plot(s, type = "values")
# Plot eigenvectors for first 6 components
plot(s, type = "vectors", idx = 1:6,
     ref = TRUE, at = "same", cuts = 50,
     plot.contrib = TRUE, symmetric = TRUE)
# Plot factor vectors for first 6 components
plot(s, type = "vectors", vectors = "factor", idx = 1:6,
     ref = TRUE, at = "same", cuts = 50,
     plot.contrib = TRUE, symmetric = TRUE)
# Plot wcor for first 12 components
plot(s, type = "wcor", groups = 1:12, grid = c(2, 6))
}
}
