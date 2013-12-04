\name{RefFreeEwasModel}
\alias{RefFreeEwasModel}
\title{Reference-Free EWAS Model}
\description{
 Reference-free method for conducting EWAS while deconvoluting DNA methylation arising as mixtures of cell types.
}
\usage{RefFreeEwasModel(Y, X, K, smallOutput=FALSE)}
\arguments{
  \item{Y}{Matrix of DNA methylation beta values (CpGs x subjects).  Missing values *are* supported.}
  \item{X}{Design matrix (subjects x covariates).}
  \item{K}{Latent variable dimension (d in Houseman et al., 2013, technical report)}
  \item{smallOutput}{Smaller output?  (Should be FALSE if you intend to run bootstraps.)}
}
\details{
Reference-free method for conducting EWAS while deconvoluting DNA methylation arising as mixtures of cell types.  This method is similar to surrogate variable analysis (SVA and ISVA), except that it makes additional use of a biological mixture assumption.  Returns mixture-adjusted Beta and unadjusted Bstar, as well as estimates of various latent quantities.
}
\value{
  A list object of class \dQuote{RefFreeEwasModel}. The most important elements are Beta and Bstar.
}
\references{Houseman EA, Molitor J, and Marsit CJ (2013), Reference-Free Cell Mixture Adjustments in Analysis of DNA Methylation Data.  Currently a tech report, in revision for publication.}
\author{E. Andres Houseman}
\seealso{\code{\link{BootRefFreeEwasModel}}}
\examples{

data(RefFreeEWAS)

\dontrun{
  tmpDesign <- cbind(1, rfEwasExampleCovariate)
  tmpBstar <- rfEwasExampleBetaValues %*% tmpDesign %*% solve(t(tmpDesign)%*%tmpDesign)
  
  EstDimRMT(rfEwasExampleBetaValues-tmpBstar %*% t(tmpDesign))$dim  
}

test <- RefFreeEwasModel(
  rfEwasExampleBetaValues, 
  cbind(1,rfEwasExampleCovariate),
  4)

testBoot <- BootRefFreeEwasModel(test,10)
summary(testBoot)

}
\keyword{deconvolution}
\keyword{DNA methylation}
\keyword{EWAS}
\keyword{surrogate variable}
\keyword{cell mixture}
\keyword{svd}


