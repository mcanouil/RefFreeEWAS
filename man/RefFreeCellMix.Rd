#' Reference-Free Cell Mixture Projection
#'
#' Reference-free cell-mixture decomposition of DNA methylation data set
#'
#' @param Y Matrix (m CpGs x n Subjects) of DNA methylation beta values
#' @param mu0 Matrix (m CpGs x K cell types) of *initial* cell-type specific methylomes
#' @param K Number of cell types (ignored if mu0 is provided)
#' @param iters Number of iterations to execute
#' @param Yfinal Matrix (m* CpGs x n Subjects) of DNA methylation beta values on which to base final methylomes
#' @param verbose Report summary of errors after each iteration?
#' @details Reference-free decomposition of DNA methylation matrix into cell-type distributions and cell-type
#'          methylomes, Y = Mu Omega^T.  Either an initial estimate of Mu must be provided, or else the number of cell types K,
#'          in which case RefFreeCellMixInitialize will be used to initialize.  Note that the decomposition will be based on Y,
#'          but Yfinal (=Y by default) will be used to determine the final value of Mu based on the last iterated value of Omega.
#' @return Object of S3 class RefFreeCellMix, containing the last iteration of Mu and Omega.
#' @references Houseman, E. Andres, Kile, Molly L., Christiani, David C., et al. Reference-free deconvolution of DNA methylation data and mediation by cell composition effects. BMC bioinformatics, 2016, vol. 17, no 1, p. 259.
#' @seealso `RefFreeCellMixInitialize`
\examples{
data(HNSCC)

# Typical use
Y.shortTest <- Y.HNSCC.averageBetas[1:500,]
Y.shortTest.final <- Y.HNSCC.averageBetas[1:1000,]
testArray1  <- RefFreeCellMixArray(Y.shortTest,Klist=1:3,iters=5,Yfinal=Y.shortTest.final)
testArray1
lapply(testArray1,summary)
sapply(testArray1,deviance,Y=Y.shortTest.final)

# Example with explicit initialization
testKeq2  <- RefFreeCellMix(Y.shortTest, mu0=RefFreeCellMixInitialize(Y.shortTest,K=2))
testKeq2
head(testKeq2$Mu)
head(testKeq2$Omega)
