\name{BootOneRefFreeEwasModel}
\alias{BootOneRefFreeEwasModel}
\title{One Bootstrap sample for Reference-Free EWAS Model}
\description{
 Bootstrap generation procedure for reference-free method for conducting EWAS while deconvoluting DNA methylation arising as mixtures of cell types.
}
\usage{BootOneRefFreeEwasModel(mod)}
\arguments{
  \item{mod}{model object of class RefFreeEwasModel (generated with smallOutput=FALSE).}
}
\details{
Generates one bootstrapped data set for the reference-free method for conducting EWAS while deconvoluting DNA methylation arising as mixtures of cell types.  Typically not run by user.
}
\value{
  A matrix representing a bootstrap sample of an DNA methylation assay matrix.
}
\references{Houseman EA, Molitor J, and Marsit CJ (2013), Reference-Free Cell Mixture Adjustments in Analysis of DNA Methylation Data.  Currently a tech report, in revision for publication.}
\author{E. Andres Houseman}
\seealso{\code{\link{BootRefFreeEwasModel}}}
\keyword{deconvolution}
\keyword{DNA methylation}
\keyword{EWAS}
\keyword{surrogate variable}
\keyword{cell mixture}

