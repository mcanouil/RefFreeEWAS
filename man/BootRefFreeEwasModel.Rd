\name{BootRefFreeEwasModel}
\alias{BootRefFreeEwasModel}
\title{Bootstrap for Reference-Free EWAS Model}
\description{
 Bootstrap procedure for reference-free method for conducting EWAS while deconvoluting DNA methylation arising as mixtures of cell types.
}
\usage{BootRefFreeEwasModel(mod, nboot)}
\arguments{
  \item{mod}{model object of class RefFreeEwasModel (generated with smallOutput=FALSE).}
  \item{nboot}{Number of bootstrap samples to generate}
}
\details{
Generates the bootstrap samples for the reference-free method for conducting EWAS while deconvoluting DNA methylation arising as mixtures of cell types.
}
\value{
  An array object of class \dQuote{BootRefFreeEwasModel}. Bootstraps are generated for both Beta and Bstar.
}
\references{Houseman EA, Molitor J, and Marsit CJ (2013), Reference-Free Cell Mixture Adjustments in Analysis of DNA Methylation Data.  Currently a tech report, in revision for publication.}
\author{E. Andres Houseman}
\seealso{\code{\link{RefFreeEwasModel}}}
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
\keyword{bootstrap}
