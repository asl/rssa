\name{plot.reconstruction}
\alias{plot.reconstruction}
\alias{plot.ssa.reconstruction}
\alias{plot.1d.ssa.reconstruction}
\alias{plot.toeplitz.ssa.reconstruction}

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
}

\examples{
# Decompose 'co2' series with default parameters
s <- ssa(co2)
r <- reconstruct(s, groups = list(c(1, 4), c(2, 3), c(5, 6)))
# Plot full 'co2' reconstruction into trend, periodic components and noise
plot(r) 
}
