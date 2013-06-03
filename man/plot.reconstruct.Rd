\name{plot.reconstruction}
\alias{plot.reconstruction}
\alias{plot.ssa.reconstruction}
\alias{plot.1d.ssa.reconstruction}
\alias{plot.toeplitz.ssa.reconstruction}
\alias{plot.2d.ssa.reconstruction}

\title{Plot the results of SSA reconstruction}
\description{
  Plot the result of SSA Reconstruction step
}
\usage{
\method{plot}{1d.ssa.reconstruction}(x, \dots,
    type = c("raw", "cumsum"),
    plot.method = c("native", "matplot"),
    base.series = NULL,
    add.original = TRUE,
    add.residuals = TRUE)
\method{plot}{toeplitz.ssa.reconstruction}(x, \dots,
    type = c("raw", "cumsum"),
    plot.method = c("native", "matplot"),
    base.series = NULL,
    add.original = TRUE,
    add.residuals = TRUE)
\method{plot}{2d.ssa.reconstruction}(x, \dots,
    type = c("raw", "cumsum"),
    base.series = NULL,
    add.original = TRUE,
    add.residuals = TRUE,
    add.ranges,
    at)
}

\arguments{
  \item{x}{SSA object holding the decomposition}
  \item{type}{Type of the plot (see 'Details' for more information)}
  \item{\dots}{Arguments to be passed to methods, such as graphical
    parameters}
  \item{plot.method}{Plotting method to use: either ordinary all-in-one
    via matplot, or native plotting method of the input time series}
  \item{base.series}{another SSA reconstruction object, the series of
    which should be considered as an original. Useful for plotting the
    results of sequential SSA}
  \item{add.original}{logical, if 'TRUE' then the original series are added
    to the plot}
  \item{add.residuals}{logical, if 'TRUE' then the residuals are added
  to the plot}
  \item{at}{A numeric vector giving breakpoints along the range of \code{z}.
    Also can be a character string. Values 'free' and 'same' are allowed. 'free' means
    special breakpoints' vectors (will be evaluated automatically, see description of \code{cuts}
    argument in 'Details') for each component. 'same' means one breakpoints' vector for all
    component (will be evaluated automatically too)}
  \item{add.ranges}{logical, if 'TRUE', the range of the components values
    will be printed in panels captions}
}

\details{
  Additional (non-standard) graphical parameters which can be transfered via \dots:
  \describe{
    \item{cuts}{the number of levels the range of \code{z} would be divided into.}
    \item{ref}{logical. Whether to plot zero-level isolines in 2d-reconstruction-plots.}
    \item{symmetric}{logical. Whether to use symmetric \code{z}-scale in 2d-reconstruction-plots.}
    \item{useRaster}{logical. For 2d-reconstruction-plot, indicating whether raster representations
      should be used. 'TRUE' by default.}
  }
}

\examples{
# Decompose 'co2' series with default parameters
s <- ssa(co2)
r <- reconstruct(s, groups = list(c(1, 4), c(2, 3), c(5, 6)))
# Plot full 'co2' reconstruction into trend, periodic components and noise
plot(r)

# Artificial field for 2dSSA
mx <- outer(1:50, 1:50,
            function(i, j) sin(2*pi * i/17) * cos(2*pi * j/7) + exp(i/25 - j/20)) +
      rnorm(50^2, sd = 0.1)
# Decompose 'mx' with default parameters
s <- ssa(mx, kind = "2d-ssa")
# Reconstruct
r <- reconstruct(s, groups = list(1, 2:5, 6))
# Plot components, original field and residuals
plot(r)
# Plot cumulative sum of components only
plot(r, type = "cumsum", add.residuals = FALSE, add.original = FALSE)
}
