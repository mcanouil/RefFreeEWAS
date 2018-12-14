\name{RefFreeCellMixArrayWithCustomStart}
\alias{RefFreeCellMixArrayWithCustomStart}
\title{Reference-Free Cell Mixture Projection - Custom Initialization}
\description{Array of reference-free cell-mixture decompositions of a DNA methylation data set, with custom initialization}
\usage{RefFreeCellMixArrayWithCustomStart(Y,mu.start,Klist=1:5,iters=10,
   Yfinal=NULL,verbose=FALSE)}
\arguments{
  \item{Y}{Matrix (m CpGs x n Subjects) of DNA methylation beta values}
  \item{mu.start}{matrix of starting values for Mu: number of columns must be at least the maximum in Klist}
  \item{Klist}{List of K values (each K = assumed number of cell types)}
  \item{iters}{Number of iterations to execute for each value of K}
  \item{Yfinal}{Matrix (m* CpGs x n Subjects) of DNA methylation beta values on which to base final methylomes}
  \item{verbose}{Report summary of errors after each iteration for each fit?}
}
\details{List of Reference-free decompositions for a range of K values.  For each value of K, the decomposition is initialized by using the first K columns of mu.start. Note that for each K, the decomposition will be based on Y, but Yfinal (=Y by default) will be used to determine the final value of Mu based on the last iterated value of Omega.}
\value{List, each element is an object of S3 class RefFreeCellMix, containing the last iteration of Mu and Omega.}
\references{Houseman, E. Andres, Kile, Molly L., Christiani, David C., et al. Reference-free deconvolution of DNA methylation data and mediation by cell composition effects. BMC bioinformatics, 2016, vol. 17, no 1, p. 259.
}
\author{E. Andres Houseman}
\seealso{\code{\link{RefFreeCellMix}},
\code{\link{RefFreeCellMixInitializeBySVD}}}
