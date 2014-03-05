\name{omnibusBoot}
\alias{omnibusBoot}
\title{Bootstrap-based omnibus test of significance across all features}
\description{
 Support for bootstrap-based omnibus test of significance accounting for correlation.
}
\usage{omnibusBoot(est, boots, denDegFree)}
\arguments{
  \item{est}{Vector of m estimates, one for each of m features.}
  \item{boots}{Matrix (m x R) of bootstrap samples corresponding to the estimates}
  \item{denDegFree}{Single number representing the denominator degrees-of-freedom for computing p-values}
}
\details{
Returns one omnibus p-value based on Kolmogorov-Smirnov distance from a uniform distribution
}
\value{
A single number representing the p-value for the omnibus test over all features.
}
\references{Houseman EA, Molitor J, and Marsit CJ (2014), Reference-Free Cell Mixture Adjustments in Analysis of DNA Methylation Data.  Bioinformatics, doi: 10.1093/bioinformatics/btu029.}
\author{E. Andres Houseman}
\seealso{\code{\link{RefFreeEwasModel}}}
\examples{

data(RefFreeEWAS)

test <- RefFreeEwasModel(
  rfEwasExampleBetaValues,
  cbind(1,rfEwasExampleCovariate),
  4)

testBoot <- BootRefFreeEwasModel(test,10)
summary(testBoot)
omnibusBoot(test$Beta[,2], testBoot[,2,"B",],-diff(dim(test$X))) 
omnibusBoot(test$Bstar[,2], testBoot[,2,"B*",],-diff(dim(test$X)))
}
\keyword{bootstrap}
\keyword{omnibus}
\keyword{kolmogorov-smirnov}

