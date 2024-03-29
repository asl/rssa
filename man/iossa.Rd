\name{iossa}
\alias{iossa}
\alias{iossa.ssa}
\title{Iterative O-SSA nested decomposition}

\description{
  Perform Iterative O-SSA (IOSSA) algorithm.
}

\usage{
\method{iossa}{ssa}(x, nested.groups, \dots, tol = 1e-5, kappa = 2,
      maxiter = 100,
      norm = function(x) sqrt(mean(x^2)),
      trace = FALSE,
      kappa.balance = 0.5)
}

\arguments{
  \item{x}{SSA object holding SSA decomposition}
  \item{nested.groups}{list or named list of numbers of eigentriples
    from full decomposition, describes initial grouping for IOSSA iterations}
  \item{tol}{tolerance for IOSSA iterations}
  \item{kappa}{`kappa' parameter for sigma-correction (see `Details' and `References') procedure.
    If 'NULL', sigma-correction will not be performed}
  \item{maxiter}{upper bound for the number of iterations}
  \item{norm}{function, calculates a norm of a vector; this norm is
    applied to the difference between the reconstructed series at
    sequential iterations and is used for convergence detection}
  \item{trace}{logical, indicates whether the convergence process should be traced}
  \item{kappa.balance}{sharing proportion of sigma-correction multiplier
    between column and row inner products}
  \item{\dots}{additional arguments passed to \code{\link{decompose}}
    routines}
}

\value{
  Object of `ossa' class. In addition to usual `ssa' class
  fields, it also contains the following fields:
  \describe{
    \item{iossa.result}{object of `iossa.result' class, a list which
      contains algorithm parameters, condition numbers, separability
      measures, the number of iterations and convergence status (see
      \code{\link{iossa.result}})}
    \item{iossa.groups}{list of groups within the nested decomposition;
      numbers of components correspond to their numbers in the full
      decomposition}
    \item{iossa.groups.all}{list, describes cumulative grouping after
      after sequential Iterative O-SSA decompositions in the case of
      non-intersecting \code{nested.groups}. Otherwise,
      \code{iossa.groups.all} coincides with \code{iossa.groups}}
    \item{ossa.set}{vector of the indices of elementary components used
      in Iterative O-SSA (that is, used in \code{nested.groups})}
  }
}

\details{
  See Golyandina N. and Shlemov A. (2015) and Section 2.4 in Golyanina et al (2018) 
  for full details in the 1D case and 
  p.250-252 from the same book for an example in the 2D case.

  Briefly, Iterative Oblique SSA (IOSSA) is an iterative (EM-like) method for
  improving separability in SSA. In particular, it serves for separation
  of mixed components, which are not  orthogonal, e.g., of sinusoids
  with close frequencies or for trend separation for short series. IOSSA
  performs a new decomposition of a part of the ssa-object, which is
  given by a set of eigentriples. Note that eigentriples that do not
  belong to the chosen set are not changed.

  Oblique SSA can make many series orthogonal by the choice of inner
  product. Iterative O-SSA find the separating inner products by
  iterations that are hopefully converges to a stationary point. See
  References for more details.

  Sigma-correction procedure does the renormalization of new inner
  products. This prevents the mixing of the components during the next
  iteration. Such approach makes the whole procedure more stable and can
  solve the problem of lack of strong separability.

  Details of the used algorithms can be found in Golyandina et al (2018), 
  Algorithms 2.7 and 2.8.
}

\references{
  Golyandina N., Korobeynikov A., Zhigljavsky A. (2018):
  \emph{Singular Spectrum Analysis with R.}  Use R!.
    Springer, Berlin, Heidelberg.

  Golyandina N. and Shlemov A. (2015): \emph{Variations of Singular Spectrum Analysis
  for separability improvement: non-orthogonal decompositions of time series}, 
  Statistics and Its Interface. Vol.8, No 3, P.277-294.
  \url{https://arxiv.org/abs/1308.4022}
}

\seealso{
  \code{\link{Rssa}} for an overview of the package, as well as,
  \code{\link[Rssa:ssa-object]{ssa-object}},
  \code{\link[Rssa:fossa]{fossa}},
  \code{\link[Rssa:owcor]{owcor}},
  \code{\link[Rssa:iossa.result]{iossa.result}}.
}

\examples{
# Separate three non-separable sine series with different amplitudes
N <- 150
L <- 70

omega1 <- 0.05
omega2 <- 0.06
omega3 <- 0.07

F <- 4*sin(2*pi*omega1 * (1:N)) + 2*sin(2*pi*omega2 * (1:N)) + sin(2*pi*omega3 * (1:N))
s <- ssa(F, L)
ios <- iossa(s, nested.groups = list(1:2, 3:4, 5:6), kappa = NULL, maxiter = 100, tol = 1e-3)

plot(reconstruct(ios, groups = ios$iossa.groups))
summary(ios)

\donttest{
# Separate two non-separable sines with equal amplitudes
N <- 200
L <- 100
omega1 <- 0.07
omega2 <- 0.06

F <- sin(2*pi*omega1 * (1:N)) + sin(2*pi*omega2 * (1:N))
s <- ssa(F, L)

# Apply FOSSA and then IOSSA
fs <- fossa(s, nested.groups = 1:4)
ios <- iossa(fs, nested.groups = list(1:2, 3:4), maxiter = 100)
summary(ios)

opar <- par(mfrow = c(3, 1))
plot(reconstruct(s, groups = list(1:2, 3:4)))
plot(reconstruct(fs, groups = list(1:2, 3:4)))
plot(reconstruct(ios, groups = ios$iossa.groups))
par(opar)

wo <- plot(wcor(ios, groups = 1:4))
gwo <- plot(owcor(ios, groups = 1:4))
plot(wo, split = c(1, 1, 2, 1), more = TRUE)
plot(gwo, split = c(2, 1, 2, 1), more = FALSE)
}

\donttest{
data(USUnemployment)
unempl.male <- USUnemployment[, "MALE"]

s <- ssa(unempl.male)
ios <- iossa(s, nested.groups = list(c(1:4, 7:11), c(5:6, 12:13)))
summary(ios)

# Comparison of reconstructions
rec <- reconstruct(s, groups = list(c(1:4, 7:11), c(5:6, 12:13)))
iorec <- reconstruct(ios, groups <- ios$iossa.groups)
# Trends
matplot(data.frame(iorec$F1, rec$F1, unempl.male), type='l',
        col=c("red","blue","black"), lty=c(1,1,2))
# Seasonalities
matplot(data.frame(iorec$F2, rec$F2), type='l', col=c("red","blue"),lty=c(1,1))

# W-cor matrix before IOSSA and w-cor matrix after it
ws <- plot(wcor(s, groups = 1:30), grid = 14)
wios <- plot(wcor(ios, groups = 1:30), grid = 14)
plot(ws, split = c(1, 1, 2, 1), more = TRUE)
plot(wios, split = c(2, 1, 2, 1), more = FALSE)

# Eigenvectors before and after Iterative O-SSA
plot(s, type = "vectors", idx = 1:13)
plot(ios, type = "vectors", idx = 1:13)

# 2D plots of periodic eigenvectors before and after Iterative O-SSA
plot(s, type = "paired", idx = c(5, 12))
plot(ios, type = "paired", idx = c(10, 12), plot.contrib = FALSE)

data(AustralianWine)
Fortified <- AustralianWine[, "Fortified"]
s <- ssa(window(Fortified, start = 1982 + 5/12, end = 1986 + 5/12), L = 18)
ios <- iossa(s, nested.groups = list(trend = 1, 2:7),
             kappa = NULL,
             maxIter = 1)
fs <- fossa(s, nested.groups = 1:7, gamma = 1000)

rec.ssa <- reconstruct(s, groups = list(trend = 1, 2:7))
rec.iossa <- reconstruct(ios, groups = ios$iossa.groups);
rec.fossa <- reconstruct(fs, groups = list(trend = 7, 1:6))

Fort <- cbind(`Basic SSA trend` = rec.ssa$trend,
              `Iterative O-SSA trend` = rec.iossa$trend,
              `DerivSSA trend` = rec.fossa$trend,
              `Full series` = Fortified)

library(lattice)
xyplot(Fort, superpose = TRUE, col = c("red", "blue", "green4", "black"))
}


# Shaped 2D I. O-SSA separates finite rank fields exactly
mx1 <- outer(1:50, 1:50,
             function(i, j) exp(i/25 - j/20))
mx2 <- outer(1:50, 1:50,
             function(i, j) sin(2*pi * i/17) * cos(2*pi * j/7))

mask <- matrix(TRUE, 50, 50)
mask[23:25, 23:27] <- FALSE
mask[1:2, 1] <- FALSE
mask[50:49, 1] <- FALSE
mask[1:2, 50] <- FALSE

mx1[!mask] <- mx2[!mask] <- NA

s <- ssa(mx1 + mx2, kind = "2d-ssa", L = c(10, 10))
plot(reconstruct(s, groups = list(1, 2:5)))

ios <- iossa(s, nested.groups = list(1, 2:5), kappa = NULL)
plot(reconstruct(ios, groups = ios$iossa.groups))

\donttest{
# I. O-SSA for MSSA
N.A <- 150
N.B <- 120
L <- 40

omega1 <- 0.05
omega2 <- 0.055

tt.A <- 1:N.A
tt.B <- 1:N.B
F1 <- list(A = 2 * sin(2*pi * omega1 * tt.A), B = cos(2*pi * omega1 * tt.B))
F2 <- list(A = 1 * sin(2*pi * omega2 * tt.A), B = cos(2*pi * omega2 * tt.B))

F <- list(A = F1$A + F2$A, B = F1$B + F2$B)

s <- ssa(F, kind = "mssa")
plot(reconstruct(s, groups = list(1:2, 3:4)), plot.method = "xyplot")

ios <- iossa(s, nested.groups = list(1:2, 3:4), kappa = NULL)
plot(reconstruct(ios, groups = ios$iossa.groups), plot.method = "xyplot")
}
}
