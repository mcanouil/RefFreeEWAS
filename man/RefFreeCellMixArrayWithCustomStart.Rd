#' Reference-Free Cell Mixture Projection - Custom Initialization
#'
#' Array of reference-free cell-mixture decompositions of a DNA methylation data set, with custom initialization
#'
#' @param Y Matrix (m CpGs x n Subjects) of DNA methylation beta values
#' @param mu.start matrix of starting values for Mu: number of columns must be at least the maximum in Klist
#' @param Klist List of K values (each K = assumed number of cell types)
#' @param iters Number of iterations to execute for each value of K
#' @param Yfinal Matrix (m* CpGs x n Subjects) of DNA methylation beta values on which to base final methylomes
#' @param verbose Report summary of errors after each iteration for each fit?
#' @details List of Reference-free decompositions for a range of K values.  For each value of K, the decomposition is initialized by using the first K columns of mu.start. Note that for each K, the decomposition will be based on Y, but Yfinal (=Y by default) will be used to determine the final value of Mu based on the last iterated value of Omega.
#' @return List, each element is an object of S3 class RefFreeCellMix, containing the last iteration of Mu and Omega.
#' @seealso `RefFreeCellMix`, `RefFreeCellMixInitializeBySVD`
