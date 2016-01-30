\name{deviance.RefFreeCellMix}
\alias{deviance.RefFreeCellMix}
\title{deviance.RefFreeCellMix}
\description{
Deviance method for objects of type RefFreeCellMix.
}
\usage{
\method{deviance}{RefFreeCellMix}(object, Y, Y.oob=NULL, EPSILON=1E-9,
  bootstrapIterations=0, bootstrapIndices=NULL, ...)
}
\arguments{
  \item{object}{RefFreeCellMix object to summarize}
  \item{Y}{Methylation matrix on which x was based}
  \item{Y.oob}{Alternate ("out-of-box") methylation matrix for which to calculate deviance, based on x}
  \item{EPSILON}{Minimum value of variance (zero variances will be reset to this value)}
  \item{bootstrapIterations}{Number of RefFreeCellMix iterations to use in bootstrap (see details)}
  \item{bootstrapIndices}{Bootstrap indices (see details)}
  \item{\dots}{(Unused).}
}
\details{
  Deviance based on normal distribution applied to errors of Y after accounting for cell mixture effect,
  Mu Omega^T. Since RefFreeCellMix does not save the original data Y in the resulting object x, Y must
  be supplied here.  However, deviance may be calculated for an alternative "out-of-bag" methylation matrix,
  Y.oob.  If bootstrapIterations=0, this is what is done.  If bootstrapIterations>0, then x$Mu is used to
  initialize a new value of x via RefFreeCellMix executed on a bootstrap sample of Y 
  with the number of indicated iterations.  If bootstrapIndices is provided, the bootstrap will be based on
  these indices, otherwise the indices will be sampled randomly with replacement from 1:ncol(Y).
  See \code{\link{RefFreeCellMix}} for example.
}

