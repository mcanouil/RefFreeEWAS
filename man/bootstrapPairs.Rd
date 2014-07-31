\name{bootstrapPairs}
\alias{bootstrapPairs}
\title{One Bootstrap Sample for Pairs}
\description{
 Bootstrap generation procedure for sampling paired data (e.g. twin data)
}
\usage{bootstrapPairs(obs, pairID)}
\arguments{
  \item{obs}{Observation ids (numeric vector).}
  \item{pairID}{Pair IDs (one unique value per pair).}
}
\details{
Generates one bootstrapped set of ids corresponding to pairs for the method for conducting EWAS while deconvoluting DNA methylation arising as mixtures of cell types.  Typically not run by user. 
}
\value{
  A vector of IDs corresponding to bootstrapped pairs
}
\references{Houseman EA, Molitor J, and Marsit CJ (Bioinformatics, 2014).}
\author{E. Andres Houseman}
\seealso{\code{\link{BootRefFreeEwasModel}},\code{\link{PairsBootRefFreeEwasModel}}}
\keyword{deconvolution}
\keyword{DNA methylation}
\keyword{EWAS}
\keyword{surrogate variable}
\keyword{cell mixture}

