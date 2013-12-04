\name{svdSafe}
\alias{svdSafe}
\title{Safe SVD-like matrix decomposition}
\description{SVD that traps errors and switches to QR when necessary
}
\usage{svdSafe(X)}
\arguments{
  \item{X}{Matrix to decompose}
}
\details{
This function traps errors in the svd function due to numerically zero singular values, and replaces the operation with a QR decomposition.  Technically, the R component of the decomposition fails the orthogonality constraint required for the SVD decomposition, but this function exists to save bootstraps from rudely failing; since the critical component of the SVD (in this application) is the left orthogonal matrix, this is a reasonable approximation for bootstrap purposes.  If there are too many svd failures (which will will be reported by the function) then it is worth looking into the design matrix.
}
\value{
 A list as in what svd produces: U and V matrices as well as the d vector of singular values. 
}
\author{E. Andres Houseman}
\seealso{\code{\link{svd}}}
\keyword{svd}

