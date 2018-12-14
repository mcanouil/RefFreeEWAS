\name{SVDwithMissing}
\alias{SVDwithMissing}
\title{SVD with missing values}
\description{Compute singular value decomposition on a matrix with missing values, using a naive/simple method for imputing missing values by row-mean}
\usage{SVDwithMissing(Y)}
\arguments{
  \item{Y}{Matrix for which to compute SVD.}
}
\details{Computes singular value decomposition on a matrix with missing values, using a naive/simple method for imputing missing values by row-mean.  Not recommended for matrices with very large numbers of missing values.}
\value{singular value decomposition (as returned by svd function)}

