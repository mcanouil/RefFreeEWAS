#' One Bootstrap Sample for Reference-Free EWAS Model, Accounting for Paired Data
#'
#' Bootstrap generation procedure for reference-free method for conducting EWAS while deconvoluting DNA methylation arising as mixtures of cell types.
#' This version accounts for paired data (e.g. twin data)
#'
#' @param mod model object of class RefFreeEwasModel (generated with smallOutput=FALSE).
#' @param pairID Pair IDs (one unique value per pair).
#' @details Generates one bootstrapped data set for the reference-free method for conducting EWAS while deconvoluting DNA methylation arising as mixtures of cell types.  This version facilitates the estimation of robust standard errors to account for paired data (e.g. twin data) using a strategy similar to that employed by Generalized Estimating Equations (GEEs).  Specifically, in bootstrapping the errors, the pairs are sampled rather than individual arrays.  Typically not run by user. 
#' @return   A matrix representing a bootstrap sample of an DNA methylation assay matrix.
#' @seealso `BootRefFreeEwasModel`,`BootOneRefFreeEwasModel`
#' @keywords deconvolution,DNA methylation,EWAS,surrogate variable,cell mixture
#' @export
PairsBootOneRefFreeEwasModel <- function(mod, pairID) {
  n2 <- dim(mod$X)[1]
  iboot <- bootstrapPairs(1:n2, pairID)
  mu <- mod$Bstar %*% t(mod$X)
  return(mu + mod$dispersion * mod$E[, iboot])
}
