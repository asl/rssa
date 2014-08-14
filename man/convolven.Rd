\name{convolven}
\alias{convolven}

\title{Convolution of arrays via mdFFT}

\description{
  Compute different kinds of convolutions of two numeric arrays using multidimensional
  Fast Fourier Transform
}

\usage{
convolven(x, y, conj = TRUE, type = "circular")
}

\arguments{
  \item{x,y}{numeric (real) vectors, matrices or arrays of the same rank (the number of dimensions)}
  \item{conj}{logical, whether perform complex conjugation to the Fourier transform
    of the second argument before element-wise multiplication. 'TRUE' (default) for
    conventional convolution}
  \item{type}{character vector, convolution type for each dimension. Possible values are
    'circular' (default), 'open' and 'filter'. Circular convolution is allowed even for arrays with
     different corresponding dimensions. If argument hasn't got sufficient length, recycling will be used}
}

\value{
  Vector, matrix or array of the same rank as input arrays.
}

\details{
  This function is the direct generalization of standard \code{\link[stats:convolve]{convolve}} function
  for multidimensional case. FFTW library (\url{http://www.fftw.org}) is used for multidimensional FFT.
}

\examples{
  data(Mars)
  # Smooth using ordinary square window
  Mars.smoothed <- convolven(Mars, matrix(1/400, 20, 20), type = "filter")

  library(lattice)
  p.original <- levelplot(Mars, main = "original")
  p.smoothed <- levelplot(Mars.smoothed, main = "smoothed")
  plot(p.original, split = c(1, 1, 2, 1), more = TRUE)
  plot(p.smoothed, split = c(2, 1, 2, 1), more = FALSE)
}
