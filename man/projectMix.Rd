\name{projectMix}
\alias{projectMix}
\title{Cell Mixture Projection (reference-based)}
\description{Constrained linear projection for estimating cell mixture or related coefficients.
}
\usage{projectMix(Y, Xmat, nonnegative=TRUE, sumLessThanOne=TRUE, lessThanOne=!sumLessThanOne)}
\arguments{
  \item{Y}{Matrix (m CpGs x n Subjects) of DNA methylation beta values}
  \item{Xmat}{Matrix (m CpGs x K cell types) of cell-type specific methylomes}
  \item{nonnegative}{All coefficients >=0?}
  \item{sumLessThanOne}{Coefficient rows should sum to less than one?}
  \item{lessThanOne}{Every value should be less than one (but possibly sum to value greater than one)?}
}
\details{Function for projecting methylation values (Y) onto space of methyomes (Xmat), with various constraints.
This is the reference-based
method described in Houseman et al. (2012) and also appearing in the minfi package.}
\value{
  Projection coefficients resulting from constrained projection
}
\references{Houseman EA, Accomando WP et al. DNA methylation arrays as surrogate measures of cell mixture distribution,
BMC Bioinformatics, 2012.}
\author{E. Andres Houseman}

