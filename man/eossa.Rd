\name{eossa}
\alias{eossa}
\alias{eossa.ssa}
\title{ESPRIT-based O-SSA nested decomposition}

\description{
  Perform ESPRIT-based O-SSA (EOSSA) algorithm.
}

\usage{
\method{eossa}{ssa}(x, nested.groups, k = 2,
      subspace = c("column", "row"),
      dimensions = NULL,
      solve.method = c("ls", "tls"),
      beta = 8,
      \dots)
}

\arguments{
  \item{x}{SSA object holding SSA decomposition}
  \item{nested.groups}{list or named list of numbers of eigentriples
    from full decomposition, describes elementary components for EOSSA nested redecomposition}
  \item{k}{the number of components in desired resultant decomposition}
  \item{subspace}{which subspace will be used for oblique matrix construction}
  \item{dimensions}{a vector of dimension indices to construct shift matrices along. 'NULL' means all dimensions}
  \item{solve.method}{approximate matrix equation solving method, 'ls' for least-squares, 'tls' for total-least-squares.}
  \item{beta}{In multidimensional (nD) case, coefficient(s) in convex linear combination of
    shifted matrices. The length of \code{beta} should be \code{ndim - 1},
    where \code{ndim} is the number of independent dimensions.
    If only one value is passed, it is expanded to a geometric progression.}
  \item{\dots}{additional arguments passed to \code{\link{decompose}}
    routines}
}

\value{
  Object of `ossa' class.
}

\details{
  EOSSA is an experimental signal separation method working in Nested Oblique SSA setting.
  As opposed to \code{\link[Rssa:iossa]{iossa}},
  this method does not require initial approximate decomposition.
  Moreover, it can be used for initial decomposition construction for IOSSA.

  EOSSA is motivated by parametric model of finite-dimensional signal,
  however it does not exploit this model directly and does not estimate the parameters.
  Therefore, it works for wider class of time series.
  According to the experiments, it works for series that could be locally
  approximated by a series of finite dimension, but at this moment there is no any theoretical results for this.

  EOSSA constructs shift matrix estimation by the same way is in ESPRIT
  (see \code{\link[Rssa:parestimate]{parestimate}}) method and
  uses its eigenspace to build separating scalar products
  (see \code{\link[Rssa:iossa]{iossa}} for more information about Oblique SSA decompositions).
  Consequently, the method ideally separates signals of finite dimension with absence of noise.
  With presence of noise it provides approximate results due to continuity.
  The method performs eigenvectors clustering inside (for now \code{\link[stats:hclust]{hclust}} is used),
  the number of components (argument \code{k}) should be passed.
}

\references{
  Shlemov A. (2017): \emph{The method of signal separation using the eigenspaces of the shift matrices (in Russian)},
  In Proceedings of the SPISOK-2017 conference, April 26--28, Saint Petersburg, Russia.
}

\seealso{
  \code{\link{Rssa}} for an overview of the package, as well as,
  \code{\link[Rssa:ssa-object]{ssa-object}},
  \code{\link[Rssa:parestimate]{ESPRIT}},
  \code{\link[Rssa:iossa]{iossa}},
  \code{\link[Rssa:fossa]{fossa}},
  \code{\link[Rssa:owcor]{owcor}},
  \code{\link[Rssa:iossa.result]{iossa.result}}.
}

\examples{
# Separability of three finite-dimensional series, EOSSA vs Basic SSA
N <- 150
L <- 70
omega1 <- 0.065
omega2 <- 0.07
omega3 <- 0.02
sigma <- 0.5

F1.real <- 2*sin(2*pi*omega1*(1:N))
F2.real <- 4*sin(2*pi*omega2*(1:N))
F3.real <- sin(2*pi*omega3*(1:N))

noise <- rnorm(N, sd = sigma)
F <- F1.real + F2.real + F3.real + noise

ss <- ssa(F, L)
eoss <- eossa(ss, nested.groups = list(1:2, 3:4, 5:6), k = 3)

print(eoss)

plot(ss, type = "series", groups = list(1:2, 3:4, 5:6))
plot(eoss, type = "series", groups = eoss$iossa.groups)

plot(reconstruct(ss,
                 groups = list(1:2, 3:4, 5:6)),
     add.residuals = TRUE, plot.method = "xyplot", main = "",
     xlab = "")

plot(reconstruct(eoss, groups = list(1:2, 3:4, 5:6)),
     add.residuals = TRUE, plot.method = "xyplot", main = "",
     xlab = "")

plot(reconstruct(ss,
                 groups = list(Reconstructed = 1:6, F1 = 1:2, F2 = 3:4, F3 = 5:6)),
     add.residuals = TRUE, plot.method = "xyplot", main = "",
     xlab = "")

plot(reconstruct(eoss,
                 groups = list(Reconstructed = 1:6, F1 = 1:2, F2 = 3:4, F3 = 5:6)),
     add.residuals = TRUE, plot.method = "xyplot", main = "",
     xlab = "")

rec.ideal <- reconstruct(ss,
                         groups = list(Signal = 1:6, F1 = 1:2, F2 = 3:4, F3 = 5:6))
rec.ideal$Signal <- F1.real + F2.real + F3.real
rec.ideal$F1 <- F2.real
rec.ideal$F2 <- F1.real
rec.ideal$F3 <- F3.real

plot(rec.ideal,
     add.residuals = TRUE, plot.method = "xyplot", main = "",
     xlab = "")

# Real-life example (co2), EOSSA vs Basic SSA
sigma <- 0.05 
ss <- ssa(co2)
plot(ss, type = "vector")
eoss <- eossa(ss, 1:6, k = 4)
eoss$iossa.groups

plot(eoss)
rec <- reconstruct(eoss, groups = eoss$iossa.groups)
plot(rec)

plot(reconstruct(ss,
                 groups = list(ET1 = 1,ET2 = 2,ET3 = 3,ET4 = 4,ET5 = 5,ET6 = 6)),
     add.residuals = TRUE, plot.method = "xyplot", main = "",
     xlab = "")

plot(reconstruct(eoss,
                 groups = eoss$iossa.groups),
     add.residuals = TRUE, plot.method = "xyplot", main = "",
     xlab = "")

# Sine wave with phase shift, EOSSA vs Basic SSA
omega1 <- 0.06
omega2 <- 0.07
sigma <- 0.25

F1.real <- sin(2*pi*omega1*(1:N))
F2.real <- sin(2*pi*omega2*(1:N))
v <- c(F1.real,  F2.real)
v <- v + rnorm(v, sd = sigma)
# v <- c(F1.real,  F2.real)

ss <- ssa(v, L = 35)

eoss <- eossa(ss, 1:4, 2)
ioss <- iossa(ss, list(1:2, 3:4))

plot(reconstruct(eoss, groups = eoss$iossa.groups))

plot(reconstruct(eoss,
     groups = eoss$iossa.groups), plot.method = "xyplot", main = "",
     xlab = "")

plot(reconstruct(ss, groups = list(1:2, 3:4)),
     plot.method = "xyplot",
     main = "", xlab = "")
plot(reconstruct(ss, groups = list(1,2, 3,4)),
     plot.method = "xyplot",
     main = "", xlab = "")
}
