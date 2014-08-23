\name{plot.reconstruction}
\alias{plot.reconstruction}
\alias{plot.ssa.reconstruction}
\alias{plot.1d.ssa.reconstruction}
\alias{plot.toeplitz.ssa.reconstruction}
\alias{plot.mssa.reconstruction}
\alias{plot.2d.ssa.reconstruction}
\alias{plot.nd.ssa.reconstruction}

\title{Plot the results of SSA reconstruction}
\description{
  Plot the result of SSA Reconstruction step
}
\usage{
\method{plot}{1d.ssa.reconstruction}(x, \dots,
     type = c("raw", "cumsum"),
     plot.method = c("native", "matplot", "xyplot"),
     base.series = NULL,
     add.original = TRUE,
     add.residuals = TRUE)
\method{plot}{toeplitz.ssa.reconstruction}(x, \dots,
     type = c("raw", "cumsum"),
     plot.method = c("native", "matplot", "xyplot"),
     base.series = NULL,
     add.original = TRUE,
     add.residuals = TRUE)
\method{plot}{mssa.reconstruction}(x,
     slice = list(),
     \dots,
     type = c("raw", "cumsum"),
     plot.method = c("native", "matplot", "xyplot"),
     na.pad = c("left", "right"),
     base.series = NULL,
     add.original = TRUE,
     add.residuals = TRUE)
\method{plot}{2d.ssa.reconstruction}(x, \dots,
     type = c("raw", "cumsum"),
     base.series = NULL,
     add.original = TRUE,
     add.residuals = TRUE,
     add.ranges,
     col = grey(c(0, 1)),
     zlim,
     at)
\method{plot}{nd.ssa.reconstruction}(x, slice, \dots)
}

\arguments{
  \item{x}{SSA object holding the decomposition}
  \item{slice}{for `mssa': list with elements named 'series' and 'components';
    for `nd.ssa': list with elements named 'i', 'j', 'k' or 'x', 'y', 'z', 't'
    or 'd1', 'd2', ... or `1`, `2`, ...; works like '['-operator,
    allows one to select which components from the reconstruction of
    multivariate time series or which subarray from reconstruction of
    multidimentional array to draw.}
  \item{type}{Type of the plot (see 'Details' for more information)}
  \item{\dots}{Arguments to be passed to methods, such as graphical
    parameters}
  \item{plot.method}{Plotting method to use: either ordinary all-in-one
    via matplot or xyplot, or native plotting method of the input time
    series}
  \item{na.pad}{select how to pad the series of unequal length with NA's}
  \item{base.series}{another SSA reconstruction object, the series of
    which should be considered as an original. Useful for plotting the
    results of sequential SSA}
  \item{add.original}{logical, if 'TRUE' then the original series are added
    to the plot}
  \item{add.residuals}{logical, if 'TRUE' then the residuals are added
  to the plot}
  \item{col}{color vector for colorscale,
    given by two or more colors,
    the first color corresponds to the minimal value,
    while the last one corresponds to the maximal value (will be interpolated by \code{colorRamp})}
  \item{zlim}{for 2d-plot, range of displayed values}
  \item{at}{for 2d-eigenvectors-plot, a numeric vector giving breakpoints along the range of \code{z},
    a list of such vectors
    or a character string.
    If a list is given, corresponding list element (with recycling) will be used for each
    plot panel.
    For character strings, values 'free' and 'same' are allowed: 'free' means
    special breakpoints' vectors (will be evaluated automatically, see description of \code{cuts}
    argument in 'Details') for each component. 'same' means one breakpoints' vector for all
    component (will be evaluated automatically too)}
  \item{add.ranges}{logical, if 'TRUE', the range of the components values
    will be printed in panels captions}
}

\details{
  Additional (non-standard) graphical parameters applicable to 2D SSA plots can be transfered via \dots:
  \describe{
    \item{cuts}{the number of levels the range of image would be divided into.}
    \item{ref}{logical, whether to plot zero-level isolines}
    \item{symmetric}{logical, whether to use symmetric image range scale}
    \item{useRaster}{logical, indicates whether raster representations should be used. 'TRUE' by default.}
    \item{fill.uncovered}{single number, matrix, one of the following strings: 'mean', 'original', 'void'
      or a list of such objects.
      For shaped 2d-reconstruction-plot this argument defines filling
      method for uncovered by window array elements on components and residuals plots.
      If number, all uncovered elements will be replaced by it. If matrix, all uncovered elements will be replaced by
      corresponding matrix elements. If 'mean', they will be replaced by mean value of current component. If 'original',
      they will be replaced by corresponding elements of original array. 'void' (by default) means no filling.
    If list is given, corresponding list element (with recycling) will be used for each plot panel.}
    \item{fill.color}{color or 'NULL'. Defines background color for shaped 2d-reconstruction plot. If 'NULL', standard white
      background will be used.}
  }
}

\examples{
\donttest{
# Decompose 'co2' series with default parameters
s <- ssa(co2)
r <- reconstruct(s, groups = list(c(1, 4), c(2, 3), c(5, 6)))
# Plot full 'co2' reconstruction into trend, periodic components and noise
plot(r)

# Artificial image for 2dSSA
mx <- outer(1:50, 1:50,
            function(i, j) sin(2*pi * i/17) * cos(2*pi * j/7) + exp(i/25 - j/20)) +
      rnorm(50^2, sd = 0.1)
# Decompose 'mx' with default parameters
s <- ssa(mx, kind = "2d-ssa")
# Reconstruct
r <- reconstruct(s, groups = list(1, 2:5))
# Plot components, original image and residuals
plot(r)
# Plot cumulative sum of components only
plot(r, type = "cumsum", add.residuals = FALSE, add.original = FALSE)

# Real example: Mars photo
data(Mars)
# Decompose only Mars image (without backgroud)
s <- ssa(Mars, mask = Mars != 0, wmask = circle(50), kind = "2d-ssa")
# Reconstruct and plot trend
plot(reconstruct(s, 1), fill.uncovered = "original")
# Reconstruct and plot texture pattern
plot(reconstruct(s, groups = list(c(13, 14, 17, 18))))

# Decompose 'EuStockMarkets' series with default parameters
s <- ssa(EuStockMarkets, kind = "mssa")
r <- reconstruct(s, groups = list(Trend = 1:2))
# Plot original series, trend and residuals superimposed
plot(r, plot.method = "xyplot", superpose = TRUE,
     auto.key = list(columns = 3),
     col = c("blue", "green", "red", "violet"),
     lty = c(rep(1, 4), rep(2, 4), rep(3, 4)))
# Plot the series separately
plot(r, plot.method = "xyplot", add.residuals = FALSE,
     screens = list(colnames(EuStockMarkets)),
     col = c("blue", "green", "red", "violet"),
     lty = c(rep(1, 4), rep(2, 4), rep(3, 4)))

# 3D-SSA example (2D-MSSA)
data(Barbara)
ss <- ssa(Barbara, L = c(50, 50, 1))
plot(reconstruct(ss, groups = 1), slice = list(k = 1))
}
}
