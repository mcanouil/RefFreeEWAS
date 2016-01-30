\name{RefFreeCellMixArrayDevianceBoot}
\alias{RefFreeCellMixArrayDevianceBoot}
\title{RefFreeCellMixArrayDevianceBoot}
\description{
Vector of bootstrapped deviances corresponding to an array of reference-free cell-mixture decompositions
}
\usage{
RefFreeCellMixArrayDevianceBoot(rfArray, Y, EPSILON=1E-9, bootstrapIterations=5)
}
\arguments{
  \item{rfArray}{list of RefFreeCellMix objects (e.g. from RefFreeCellMixArray)}
  \item{Y}{Methylation matrix on which x was based}
  \item{EPSILON}{Minimum value of variance (zero variances will be reset to this value)}
  \item{bootstrapIterations}{Number of RefFreeCellMix iterations to use in bootstrap}
}
\details{
  Vector of bootstrapped deviances corresponding to an array of reference-free cell-mixture decompositions,
  used to determine optimal number of cell types. This function returns one bootstrapped vector.  
  See \code{\link{RefFreeCellMixArrayDevianceBoots}} for more than one bootstrapped vector.
  The bootstrapped deviance is based on normal distribution applied to errors of Y after accounting for cell mixture effect, Mu   Omega^T. 
  See \code{\link{RefFreeCellMixArray}} for example.
}

