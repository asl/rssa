\name{ssa-input}
\alias{ssa-input}
\title{Input Data Formats Used by SSA Routines}

\description{
  The inputs of SSA can be quite different depending on the kind of SSA
  used. However, there is a common of all the variants of SSA and all
  the routines. The package tries hard to preserve the specifics of
  input object as much as possible. This means, that all the attributes,
  etc. are copied back to the reconstructed objects. This way, the
  result of the SSA decomposition of a 'ts' object is a 'ts' object as
  well.

  For forecasting, it is not possible in general to preserve the
  attributes of the input objects. However, \code{Rssa} knows about some
  common time series classes (e.g. 'ts') and tries to infer the time
  scales for forecasted objects as well.

  The input formats are as follows:
  \subsection{1d SSA and Toeplitz SSA}{
    Input is assumed to be a simple vector, or vector-like object
    (e.g. univariare 'ts' or 'zooreg' object). Everything else is
    coerced to vector.
  }

  \subsection{2d SSA}{
    Input assumed to be a matrix. If there are any \code{NA}'s then the
    shaped variant of 2d SSA will be used. All non-\code{NA} elements
    will be used as a mask.
  }

  \subsection{nd SSA}{
    Input assumed to be an array of arbitrary dimension. If there are any
    \code{NA}'s then the shaped variant will be used.
  }

  \subsection{MSSA}{
    While the representation of a one dimensional time series in
    R is pretty obvious, there are multiple possible ways of
    defining the multivariate time series. Let us outline some common
    choices.

    \itemize{
      \item Matrix with separate series in the columns. Optionally,
      additional time structure like in 'mts' objects, can be embedded.

      \item Matrix-like (e.g. a 'data.frame') object with series in the
      columns. In particular, 'data.frame' would be a result of reading
      the series from the file via 'read.table' function.

      \item List of separate time series objects (e.g. a 'list' of 'ts'
      or 'zoo' objects).
    }
    
    Also, the time scales of the individual time series can be
    normalized via head or tail padding with \code{NA} (for example, as
    a result of the \code{ts.union} call), or specified via time series
    attributes. Or, everything can be mixed all together.

    The \code{ssa} routine with 'kind = mssa' allows one to provide any
    of the outlined multivariate series formats. As usual, all the
    attributes, names of the series, NA padding, etc. is carefully
    preserved.

  }

  \subsection{CSSA}{
    Complex vectors are assumed at the input.
  }
}

\seealso{
  \code{\link[Rssa:ssa]{ssa}}
}

\examples{
s <- ssa(co2) # Perform the decomposition using the default window length
r <- reconstruct(s, groups = list(Trend = c(1, 4),
                                  Seasonality = c(2, 3))) # Reconstruct into 2 series
class(r$Trend) # Result is 'ts' object

# Simultaneous trend extraction using MSSA
s <- ssa(EuStockMarkets, kind = "mssa")
r <- reconstruct(s, groups = list(Trend = c(1,2)))
class(r$Trend) # Result is 'mts' object

# Trend forecast
f <- rforecast(s, groups = list(Trend = c(1, 2)), len = 50, only.new = FALSE)
class(f) # For 'ts' objects the time scales are inferred automatically

# Artificial image for 2dSSA
mx <- outer(1:50, 1:50,
            function(i, j) sin(2*pi * i/17) * cos(2*pi * j/7) + exp(i/25 - j/20)) +
      rnorm(50^2, sd = 0.1)
# Decompose 'mx' with circular window
s <- ssa(mx, kind = "2d-ssa", wmask = circle(5), neig = 10)
# Reconstruct
r <- reconstruct(s, groups = list(1, 2:5))
# Plot components, original image and residuals
plot(r)

# 3D-SSA example (2D-MSSA)
data(Barbara)
\donttest{
ss <- ssa(Barbara, L = c(50, 50, 1))
plot(ss)
}
}
