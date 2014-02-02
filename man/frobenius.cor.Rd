\name{frobenius.cor}
\alias{frobenius.cor}

\title{Calculate Frobenius correlations of the component matrices}

\description{Function calculates Frobenius correlations between grouped matrices
  from the SSA matrix decomposition}

\usage{
frobenius.cor(x, groups, ...)
}

\arguments{
  \item{x}{input SSA object, supposed to be of class `ossa'}
  \item{groups}{list of numeric vectors, indices of elementary matrix components
    in the SSA matrix decomposition}
  \item{\dots}{further arguments passed to \code{decompose}}
}

\details{
  Function computes matrix of Frobenius correlations between grouped matrices from the SSA matrix decomposition.
  For group \eqn{\mathcal{I} = \{i_1, \dots, i_s\}}{I = {i_1, \dots, i_s}} the group matrix is defined as
  \eqn{\mathbf{X}_\mathcal{I} = \sum_{i \in \mathcal{I}} \sigma_i U_i V_i^\mathrm{T}}{X_i = sum_{i in I} (sigma_i U_i (V_i)^T)}.

  Frobenius correlation of two matrices is defined as follows:
  \deqn{%
    \mathrm{fcor}(\mathbf{Z}, \mathbf{Y}) =
      \frac{\langle \mathbf{Z}, \mathbf{Y} \rangle_\mathrm{F}}
           {\|\mathbf{Z}\|_\mathrm{F} \cdot \|\mathbf{Y}\|_\mathrm{F}}.
  }{%
    fcor (Z, Y) = <Z, Y>_F / (||Z||_F ||Y||_F)).
  }

  Frobenius correlation is a measure of Frobenius orthogonality of the components.
  If grouped matrices are correlated then the w-correlations of the corresponding reconstructed series is not
  relevant measure of separability (and one should use \code{\link[Rssa:owcor]{owcor}} instead).
  Also, if the elementary matrices \eqn{\mathbf{X}_i = \sigma_i U_i V_i^\mathrm{T}}{X_i = sigma_i U_i (V_i)^T}
  of the decomposition are not F-orthogonal,
  then \eqn{\sigma_i} do not reflect their true contributions into the matrix decomposition.

  This function normally should be used only for object of class `ossa'.
  Otherwise it always returns identical matrix (for disjoint groups).
}

\value{
  Object of type 'wcor.matrix'.
}

\seealso{
  \code{\link[Rssa:wcor]{wcor}},
  \code{\link[Rssa:owcor]{owcor}},
  \code{\link[Rssa:iossa]{iossa}}.
}

\examples{
# Separation of two mixed sine-waves with equal amplitudes
N <- 150
L <- 70
omega1 <- 1/5
omega2 <- 1/10

v <- sin(2*pi*omega1 * (1:N)) + sin(2*pi*omega2 * (1:N))
s <- ssa(v, L)
fs <- fossa(s, nested.groups = 1:4, gamma = 100)

# Decomposition is F-orthogonal
plot(frobenius.cor(fs, groups = 1:4), main = "F-correlation matrix")

plot(wcor(s, groups = 1:4))
plot(wcor(fs, groups = 1:4))


# Separate two non-separable sine series with different amplitudes
\donttest{
N <- 150
L <- 70

omega1 <- 0.07
omega2 <- 0.0675

F <- 2*sin(2*pi*omega1 * (1:N)) + 2*sin(2*pi*omega2 * (1:N))
s <- ssa(F, L)
ios <- iossa(s, nested.groups = list(1:2, 3:4),
             kappa = NULL, maxiter = 1000, tol = 1e-5)

plot(reconstruct(ios, groups = ios$iossa.groups))
summary(ios)

# Decomposition is really oblique
plot(frobenius.cor(ios, groups = 1:4), main = "F-correlation matrix")

plot(wcor(ios, groups = 1:4))
plot(owcor(ios, groups = list(1:2, 3:4)), main = "Oblique W-correlation matrix")
}


\donttest{
data(USUnemployment)
unempl.male <- USUnemployment[, "MALE"]

s <- ssa(unempl.male)
ios <- iossa(s, nested.groups = list(c(1:4, 7:11), c(5:6, 12:13)))
summary(ios)

# W-cor matrix before IOSSA and w-cor matrix after it
plot(wcor(s, groups = 1:30))
plot(wcor(ios, groups = 1:30))

# Confirmation of the indicated max value in the above warning
plot(frobenius.cor(ios, groups = 1:30), main = "F-correlation matrix")
}
}

