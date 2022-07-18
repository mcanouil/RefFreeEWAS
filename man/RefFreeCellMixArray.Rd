#' Initialize Reference-Free Cell Mixture Projection
#'
#' Array of reference-free cell-mixture decompositions of a DNA methylation data set
#'
#' @param Y Matrix (m CpGs x n Subjects) of DNA methylation beta values
#' @param Klist List of K values (each K = assumed number of cell types)
#' @param iters Number of iterations to execute for each value of K
#' @param Yfinal Matrix (m* CpGs x n Subjects) of DNA methylation beta values on which to base final methylomes
#' @param verbose Report summary of errors after each iteration for each fit?
#' @param dist.method Method for calculating distance matrix for methylome initialization
#' @param ... Additional parameters for hclust function for methylome initialization
#' @details List of Reference-free decompositions for a range of K values.  For each value of K, the decomposition is initialized by hierarchical clutering as specified by the parameters dist.method, etc.  Note that for each K, the decomposition will be based on Y,
#'          but Yfinal (=Y by default) will be used to determine the final value of Mu based on the last iterated value of Omega.
#' @return List, each element is an object of S3 class RefFreeCellMix, containing the last iteration of Mu and Omega.
#' @references Houseman, E. Andres, Kile, Molly L., Christiani, David C., et al. Reference-free deconvolution of DNA methylation data and mediation by cell composition effects. BMC bioinformatics, 2016, vol. 17, no 1, p. 259.
#' @seealso `RefFreeCellMix`
\examples{
data(HNSCC)
Y.shortTest <- Y.HNSCC.averageBetas[1:500,]
testArray2  <- RefFreeCellMixArray(Y.shortTest,Klist=1:5,iters=5)
sapply(testArray2,deviance,Y=Y.shortTest)

\dontrun{
testBootDevs <- RefFreeCellMixArrayDevianceBoots(testArray2,Y.shortTest,R=10)

testBootDevs
apply(testBootDevs[-1,],2,mean,trim=0.25)
which.min(apply(testBootDevs[-1,],2,mean,trim=0.25))
