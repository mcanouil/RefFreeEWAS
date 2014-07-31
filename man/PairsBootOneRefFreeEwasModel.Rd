\name{PairsBootOneRefFreeEwasModel}
\alias{PairsBootOneRefFreeEwasModel}
\title{One Bootstrap Sample for Reference-Free EWAS Model, Accounting for Paired Data}
\description{
 Bootstrap generation procedure for reference-free method for conducting EWAS while deconvoluting DNA methylation arising as mixtures of cell types.
 This version accounts for paired data (e.g. twin data)
}
\usage{PairsBootOneRefFreeEwasModel(mod, pairID)}
\arguments{
  \item{mod}{model object of class RefFreeEwasModel (generated with smallOutput=FALSE).}
  \item{pairID}{Pair IDs (one unique value per pair).}
}
\details{
Generates one bootstrapped data set for the reference-free method for conducting EWAS while deconvoluting DNA methylation arising as mixtures of cell types.  This version facilitates the estimation of robust standard errors to account for paired data (e.g. twin data) using a strategy similar to that employed by Generalized Estimating Equations (GEEs).  Specifically, in bootstrapping the errors, the pairs are sampled rather than individual arrays.  Typically not run by user. 
}
\value{
  A matrix representing a bootstrap sample of an DNA methylation assay matrix.
}
\references{Houseman EA, Molitor J, and Marsit CJ (Bioinformatics,2014), Reference-Free Cell Mixture Adjustments in Analysis of DNA Methylation Data.  Bioinformatics, doi: 10.1093/bioinformatics/btu029.}
\author{E. Andres Houseman}
\seealso{\code{\link{BootRefFreeEwasModel}},\code{\link{BootOneRefFreeEwasModel}}}
\keyword{deconvolution}
\keyword{DNA methylation}
\keyword{EWAS}
\keyword{surrogate variable}
\keyword{cell mixture}

