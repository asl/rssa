\name{clplot}
\alias{clplot}

\title{
  Ratio of complete lag vectors in dependence on window length
}

\description{
  Function to plot the dependence of ratios of complete lagged vectors on window lengths.
}

\usage{
clplot(x, \dots)
}

\arguments{
  \item{x}{input series}
  \item{\dots}{further arguments passed to plotting functions}
}

\details{
  The function plots the dependence of ratios of complete lagged vectors on window lengths.
  This information can be used for the choice of window length, since only complete lagged vectors
  are used for construction of the SVD expansion in SSA. See page 89 (Chapter 2) in Golyandina et al (2018).
}

\references{
  Golyandina N., Korobeynikov A., Zhigljavsky A. (2018):
  \emph{Singular Spectrum Analysis with R.}  Use R!.
    Springer, Berlin, Heidelberg.
}

\seealso{
  \code{\link{Rssa}} for an overview of the package, as well as,
  \code{\link[Rssa:igapfill]{igapfill}},
  \code{\link[Rssa:gapfill]{gapfill}}
  \code{\link[Rssa:summarize.gaps]{summarize.gaps}},
}
