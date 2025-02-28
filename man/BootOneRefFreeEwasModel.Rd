% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BootRefFreeEwasModel.R
\name{BootOneRefFreeEwasModel}
\alias{BootOneRefFreeEwasModel}
\title{One Bootstrap sample for Reference-Free EWAS Model}
\usage{
BootOneRefFreeEwasModel(mod)
}
\arguments{
\item{mod}{model object of class RefFreeEwasModel
(generated with \code{smallOutput = FALSE}).}
}
\value{
A matrix representing a bootstrap sample of an DNA methylation assay matrix.
}
\description{
Bootstrap generation procedure for reference-free method for conducting EWAS
while deconvoluting DNA methylation arising as mixtures of cell types.
}
\details{
Generates one bootstrapped data set for the reference-free method
for conducting EWAS while deconvoluting DNA methylation arising as
mixtures of cell types.  Typically not run by user.
}
\seealso{
BootRefFreeEwasModel
}
\keyword{DNA}
\keyword{EWAS,}
\keyword{cell}
\keyword{deconvolution,}
\keyword{methylation,}
\keyword{mixture}
\keyword{surrogate}
\keyword{variable,}
