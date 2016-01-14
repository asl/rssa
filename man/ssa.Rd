\name{ssa}
\alias{ssa}
\alias{new.ssa}
\title{Create a new SSA object}

\description{
  Set up the SSA object and perform the decomposition, if necessary.
}

\usage{
ssa(x,
    L = (N + 1) \%/\% 2,
    neig = NULL,
    mask = NULL,
    wmask = NULL,
    column.projector = "none",
    row.projector = "none",
    column.oblique = "identical",
    row.oblique = "identical",
    \dots,
    kind = c("1d-ssa", "2d-ssa", "nd-ssa", "toeplitz-ssa", "mssa", "cssa"),
    circular = FALSE,
    svd.method = c("auto", "nutrlan", "propack", "svd", "eigen"),
    force.decompose = TRUE)
}

\arguments{
  \item{x}{object to be decomposed. See
  \code{\link[Rssa:ssa-input]{ssa-input}} for more information}
  \item{L}{integer, window length. Fixed to half of the series length by
    default. Should be vector of length 2 for 2d SSA}
  \item{neig}{integer, number of desired eigentriples. If 'NULL', then
    sane default value will be used, see 'Details'}
  \item{mask}{for shaped 2d SSA case only. Logical matrix with same dimension as
    \code{x}. Specifies form of decomposed array. If 'NULL', then all
    non-NA elements will be used}
  \item{wmask}{for shaped 2d SSA case only. Logical matrix which specifies window form.
    See `Details' for more information about the window shape selection}
  \item{\dots}{further arguments passed to \code{\link[Rssa:decompose.ssa]{decompose}} }
  \item{kind}{SSA method. This includes ordinary 1d SSA, 2d SSA,
    Toeplitz variant of 1d SSA, multichannel variant of SSA and complex
    SSA}
  \item{circular}{logical vector of one or two elements, describes series topology for
    1d SSA and Toeplitz SSA or field topology for 2d SSA. 'TRUE' means
    series circularity for 1d case or circularity by a corresponding
    coordinate for 2d case. See (Shlemov, 2014) for more info}
  \item{svd.method}{singular value decomposition method. See 'Details'
    for more info}
  \item{column.projector, row.projector}{column and row signal subspaces projectors for
    SSA with projection. See `Details' for information about methods of
    projectors specification}
  \item{column.oblique, row.oblique}{column and row matrix weights for 
  Weighted Oblique SSA. See `Details' for information about how to use this feature}
  \item{force.decompose}{logical, if 'TRUE' then the decomposition is
    performed before return.}
}

\value{
  Object of class `ssa'. The precise layout of the object is mostly
  meant opaque and subject to change in different version of the
  package. See \code{\link[Rssa:ssa-object]{ssa-object}} for details.
}

\details{
  This is the main entry point to the package. This routine
  constructs the SSA object filling all necessary internal structures
  and performing the decomposition if necessary.

  \subsection{Variants of SSA}{
    The following implementations of the SSA method are supported
    (corresponds to different values of \code{kind} argument):
    \describe{
      \item{1d-ssa}{Basic 1d SSA as described in Chapter 1 of (Golyandina et al,
        2001). This is also known as Broomhead-King variant of SSA or
        BK-SSA, see (Broomhead and King, 1986).}
      \item{toeplitz-ssa}{Toeplitz variant of 1d SSA. See Section 1.7.2
        in (Golyandina et al, 2001). This is also knows as Vatuard-Gill
        variant of SSA or VG-SSA for analysis of stationary time series,
        see (Vautard and Ghil, 1989).}
      \item{mssa}{Multichannel SSA for simultaneous decomposition of
        several time series (possible of unequal length). See
        (Golyandina and Stepanov, 2005).}
      \item{cssa}{Complex variant of 1d SSA.}
      \item{2d-ssa}{2d SSA for decomposition of images and arrays. See
        (Golyandina and Usevich, 2009) for more information.}
      \item{nd-ssa}{Multidimensional SSA decomposition for arrays (tensors).}
    }
  }

  \subsection{Window shape selection (for shaped 2d SSA)}{
    Window shape may be specified by argument \code{wmask}. If \code{wmask} is 'NULL',
    then standard rectangular window (specified by \code{L}) will be used.

    Also in \code{wmask} one may use following functions:
    \describe{
      \item{circle(R)}{circular mask of radius \code{R}}
      \item{triangle(side)}{mask in form of isosceles right-angled triangle with
        cathetus \code{side}. Right angle lay on topleft corner of container square
        matrix}
    }
    These functions are not exported, they defined only for \code{wmask} expression.
    If one has objects with the same names and wants to use them rather than these functions,
    one should use special wrapper function \code{I()} (see 'Examples'). 
  }

  \subsection{Projectors specification for SSA with projection}{
    Projectors are specified by means of \code{column.projector} and
    \code{row.projector} arguments. Each may be a matrix of orthonormal
    (otherwise QR orthonormalization process will be perfomed) basis of
    projection subspace, or single integer, which will be interpreted as
    dimension of orthogonal polynomial basis (note that the dimension
    equals to degree plus 1, e.g. quadratic basis has dimension 3), or
    one of following character strings (or unique prefix): 'none',
    'constant' (or 'centering'), 'linear', 'quadratic' or 'qubic' for
    orthonormal bases of the corresponding functions.

    Here is the the list of the most used options
    \describe{
      \item{both projectors are 'none'}{corresponds to ordinary 1D SSA,}
      \item{column.projector='centering'}{corresponds to 1D SSA
        with centering,}
      \item{column.projector='centering' and
        row.projector='centering'}{corresponds to 1D SSA with double
        centering.}
    }

    SSA with centering and double centering may improve the separation
    of linear trend (see (Golyandina, 2001) for more information).
  }

  \subsection{Weighted Oblique SSA}{
    Corresponding matrix norm weights may be specified for ordinary 1D SSA case 
    by means of \code{column.oblique} and \code{row.oblique} arguments. These 
    arguments should be either 'identical' or positive numeric vectors of length 
    \code{L} and \code{N - L + 1} for \code{column.oblique} and 
    \code{row.oblique} respectively.


    Weighted Oblique SSA inside \code{\link[Rssa:cadzow]{Cadzow}} iterations 
    may improve finite-rank estimation of signal (see (Zvonarev, 2016) for more 
    information).
  }

  \subsection{SVD methods}{
    The main step of the SSA method is the singular decomposition of the
    so-called series trajectory matrix. Package provides several
    implementations of this procedure (corresponds to different values of
    \code{svd.method}) argument:
    \describe{
      \item{auto}{Automatic method selection depending on the series length,
        window length, SSA kind and number of eigenvalues
        requested.}
      \item{nutrlan}{Thick-restart Lanczos eigensolver which operates on
        cross-product matrix. This methods exploits the Hankel structure of
        the trajectory matrix efficiently and is really fast. The method
        allows the truncated SVD (only specifid amount of eigentriples to be
        computed) and the continuation of the decomposition. See
        (Korobeynikov, 2010) for more information.}
      \item{propack}{SVD via implicitly restarted Lanczos
        bidiagonalization with partial reothogonalization. This methods
        exploits the Hankel structure of the trajectory matrix efficiently
        and is really fast. This is the 'proper' SVD implementation (the
        matrix of factor vectors are calculated), thus the memory
        requirements of the methods are higher than for nu-TRLAN. Usually
        the method is slightly faster that nu-TRLAN and more numerically
        stable. The method allows the truncated SVD (only specifid amount
        of eigentriples to be computed). See (Korobeynikov, 2010) for more
        information.}
      \item{svd}{Full SVD as provided by LAPACK DGESDD
        routine. Neither continuation of the decomposition nor the
        truncated SVD is supported. The method does not assume anything
        special about the trajectory matrix and thus is slow.}
      \item{eigen}{Full SVD via eigendecompsition of the cross-product
        matrix. In many cases faster than previous method, but still
        really slow for more or less non-trivial matrix sizes.}
    }

    Usually the \code{ssa} function tries to provide the best SVD
    implementation for given series length and the window size. In
    particular, for small series and window sizes it is better to use
    generic black-box routines (as provided by 'svd' and 'eigen'
    methods). For long series special-purpose routines are to be used.
  }
}

\references{
  Broomhead, D.S., and King, G.P. (1986a): \emph{Extracting qualitative
  dynamics from experimental data}, Physica D, 20, 217--236.

  Vautard, R., and Ghil, M. (1989): \emph{Singular spectrum analysis in
  nonlinear dynamics, with applications to paleoclimatic time series},
  Physica D, 35, 395--424.

  Golyandina, N., Nekrutkin, V. and Zhigljavsky, A. (2001): \emph{Analysis of
    Time Series Structure: SSA and related techniques.} Chapman and
  Hall/CRC. ISBN 1584881941

  Golyandina, N. and Stepanov, D. (2005): \emph{SSA-based approaches to
    analysis and forecast of multidimensional time series}. In
    Proceedings of the 5th St.Petersburg Workshop on Simulation, June
    26-July 2, 2005, St. Petersburg State University, St. Petersburg,
    293--298. \url{http://www.gistatgroup.com/gus/mssa2.pdf}
  
  Golyandina, N. and Usevich, K. (2009). \emph{2D-extensions of singular
  spectrum analysis: algorithm and elements of theory.} In Matrix
  Methods: Theory, Algorithms, Applications. World Scientific
  Publishing, 450-474.

  Korobeynikov, A. (2010) \emph{Computation- and space-efficient
    implementation of SSA.} Statistics and Its Interface, Vol. 3, No. 3,
  Pp. 257-268

  Golyandina, N., Korobeynikov, A. (2012, 2014) \emph{Basic Singular Spectrum
    Analysis and Forecasting with R.} Computational Statistics and Data
  Analysis, Vol. 71, Pp. 934-954. \url{http://arxiv.org/abs/1206.6910}

  Golyandina, N., Zhigljavsky, A. (2013): \emph{Singular Spectrum
    Analysis for time series}. Springer Briefs in Statistics. Springer.

  Golyandina, N., Korobeynikov, A., Shlemov, A. and Usevich, K.
  \emph{Multivariate and 2D Extensions of Singular Spectrum Analysis
    with the Rssa Package}. \url{http://arxiv.org/abs/1309.5050}

  Shlemov, A. and Golyandina, N. \emph{Shaped extensions of singular
  spectrum analysis}. \url{http://arxiv.org/abs/1401.4980}

  Zvonarev, N. and Golyandina, N. \emph{Iterative algorithms for weighted and 
  unweighted finite-rank time-series approximations}. 
  \url{http://arxiv.org/abs/1507.02751}
}
\seealso{
  \code{\link[svd:svd]{svd}},
  \code{\link[Rssa:ssa-object]{ssa-object}},
  \code{\link[Rssa:ssa-input]{ssa-input}},
  \code{\link[Rssa:decompose.ssa]{decompose}},
  \code{\link[Rssa:reconstruct]{reconstruct}},
  \code{\link[Rssa:plot]{plot}},
  \code{\link[Rssa:forecast]{forecast}},
}

\examples{
# Decompose 'co2' series with default parameters
s <- ssa(co2)
# Show the summary
summary(s)
# Reconstruct the series, with suitable grouping
r <- reconstruct(s, groups = list(c(1, 4), c(2, 3), c(5, 6)))
\donttest{
plot(r)
}

# Decompose 'EuStockMarkets' series with default parameters
s <- ssa(EuStockMarkets, kind = "mssa")
r <- reconstruct(s, groups = list(Trend = 1:2))
# Plot original series, trend and residuals superimposed
\donttest{
plot(r, plot.method = "xyplot", superpose = TRUE,
     auto.key = list(columns = 3),
     col = c("blue", "green", "red", "violet"),
     lty = c(rep(1, 4), rep(2, 4), rep(3, 4)))
}

# Artificial image for 2dSSA
mx <- outer(1:50, 1:50,
            function(i, j) sin(2*pi * i/17) * cos(2*pi * j/7) + exp(i/25 - j/20)) +
      rnorm(50^2, sd = 0.1)
# Decompose 'mx' with circular window
s <- ssa(mx, kind = "2d-ssa", wmask = circle(5), neig = 10)
# Reconstruct
r <- reconstruct(s, groups = list(1, 2:5))
\donttest{
# Plot components, original image and residuals
plot(r)
}

# Real example: Mars photo
data(Mars)
# Decompose only Mars image (without backgroud)
s <- ssa(Mars, mask = Mars != 0, wmask = circle(50), kind = "2d-ssa")
\donttest{
# Plot eigenarrays
plot(s, type = "vectors", idx = 1:25)
# Plot factor arrays
plot(s, type = "vectors", vectors = "factor", idx = 1:25)
# Reconstruct and plot trend
plot(reconstruct(s, 1), fill.uncovered = "original")
# Reconstruct and plot texture pattern
plot(reconstruct(s, groups = list(c(13,14, 17, 18))))

# I()-wrapper demo
circle <- 50
s <- ssa(Mars, wmask = circle(R = I(circle)))

# CSSA-based trend extraction
s <- ssa(EuStockMarkets[, 1] + 1.0i*EuStockMarkets[, 2], kind = "cssa")
r <- reconstruct(s, groups = list(Trend = 1:2))
plot(r)

# `co2' decomposition with double projection to linear functions
s <- ssa(co2, column.projector = "centering", row.projector = "centering")
plot(reconstruct(s, groups = list(trend = seq_len(nspecial(s)))))
}

# Artificial 2d example with double projection
ii <- matrix(1:100, 100, 100); jj <- t(ii)
x <- ii + 2 * jj
s <- ssa(x, column.projector = "centering", row.projector = "centering")
plot(s)
plot(reconstruct(s, groups = list(trend = seq_len(nspecial(s)))))

# 3D-SSA example (2D-MSSA)
data(Barbara)
Barbara.noised <- Barbara

# Corrupt image by regular noise
noise <- outer(seq_len(dim(Barbara)[1]),
               seq_len(dim(Barbara)[2]),
               function(i, j) sin(2*pi * (i/13 + j/23)))
Barbara.noised[,, 1] <- Barbara.noised[,, 1] + 10 * noise
Barbara.noised[,, 2] <- Barbara.noised[,, 2] + 30 * noise
Barbara.noised[,, 3] <- Barbara.noised[,, 3] + 5 * noise

# Normalize image for plotting
Barbara.noised <- (Barbara.noised - min(Barbara.noised)) / diff(range(Barbara.noised))

plot(c(0, 1), c(0, 1), type = "n", xlab = "", ylab = "")
rasterImage(Barbara.noised, 0, 0, 1, 1, interpolate = FALSE)

\donttest{
# Multichannel 2D-SSA
ss <- ssa(Barbara.noised, L = c(50, 50, 1))
plot(ss, type = "series", groups = 1:18, slice = list(k = 1))
plot(ss, type = "vectors", idx = 1:12, slice = list(k = 1))
plot(ss, type = "vectors", vectors = "factor", idx = 1:12, slice = list(k = 3))
# Denoise image
Barbara.rec <- residuals(ss, groups = 5:6)
plot(c(0, 1), c(0, 1), type = "n", xlab = "", ylab = "")
rasterImage(Barbara.rec, 0, 0, 1, 1, interpolate = FALSE)
}
}
