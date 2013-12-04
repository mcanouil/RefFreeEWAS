\name{EstDimIC}
\alias{EstDimIC}
\title{Dimension estimation by AIC and BIC}
\description{Method for estimating latent dimension by AIC and BIC.
}
\usage{EstDimIC(Rmat,Krange=0:25)}
\arguments{
  \item{Rmat}{Residual matrix for which to estimate latent dimension.}
  \item{Krange}{Vector of integers representing candidate dimensions to consider}
}
\details{Method for estimating latent dimension by AIC and BIC.  Inferior to the RMT method in the isva package, but it appears here because it's mentioned in our paper.
}
\value{
  A list containing AIC and BIC for candidate dimensions, as well as the best dimension for each. 
}
\references{Houseman EA, Molitor J, and Marsit CJ (2013), Reference-Free Cell Mixture Adjustments in Analysis of DNA Methylation Data.  Currently a tech report, in revision for publication.}
\author{E. Andres Houseman}
\seealso{\code{\link{EstDimRMT}}}
\examples{
data(RefFreeEWAS)

\dontrun{
  tmpDesign <- cbind(1, rfEwasExampleCovariate)
  tmpBstar <- rfEwasExampleBetaValues %*% tmpDesign %*% solve(t(tmpDesign)%*%tmpDesign)

  EstDimIC(rfEwasExampleBetaValues-tmpBstar %*% t(tmpDesign))
}
}

