#' Bootstrap for Reference-Free EWAS Model, Accounting for Paired Data
#'
#' Bootstrap procedure for reference-free method for conducting EWAS while deconvoluting DNA methylation arising as mixtures of cell types.
#' This version accounts for paired data (e.g. twin data)
#'
#' @param mod model object of class RefFreeEwasModel (generated with smallOutput=FALSE).
#' @param nboot Number of bootstrap samples to generate.
#' @param pairID Pair IDs (one unique value per pair).
#' @details Generates the bootstrap samples for the reference-free method for conducting EWAS while deconvoluting DNA methylation arising as mixtures of cell types.  This paired version facilitates the estimation of robust standard errors to account for paired data (e.g. twin data) using a strategy similar to that employed by Generalized Estimating Equations (GEEs).  Specifically, in bootstrapping the errors, the pairs are sampled rather than individual arrays.  An error will be generated unless each cluster has exactly two members (i.e. exactly two observations correspond to the same unique ID given in pairID).
#' @return   An array object of class \dQuote{BootRefFreeEwasModel}. Bootstraps are generated for both Beta and Bstar.
#' @seealso `RefFreeEwasModel`,`BootRefFreeEwasModel`
#' @keywords deconvolution,DNA methylation,EWAS,surrogate variable,cell mixture,bootstrap
#' @export
#' @examples
#' 
#' data(RefFreeEWAS)
#' 
#' if (interactive()) {
#'   tmpDesign <- cbind(1, rfEwasExampleCovariate)
#'   tmpBstar <- (rfEwasExampleBetaValues %*% tmpDesign %*% solve(t(tmpDesign)%*%tmpDesign))
#' 
#'   EstDimRMT(rfEwasExampleBetaValues-tmpBstar %*% t(tmpDesign))$dim
#' }
#' 
#' test <- RefFreeEwasModel(
#'   rfEwasExampleBetaValues,
#'   cbind(1,rfEwasExampleCovariate),
#'   4
#' )
#' 
#' testBoot <- BootRefFreeEwasModel(test,10)
#' summary(testBoot)
PairsBootRefFreeEwasModel <- function(mod, nboot, pairID) {
  BetaBoot <- array(NA, dim = c(dim(mod$Beta), 2, nboot))
  dimnames(BetaBoot)[1:2] <- dimnames(mod$Beta)
  dimnames(BetaBoot)[[3]] <- c("B", "B*")
  dimnames(BetaBoot)[[4]] <- 1:nboot
  attr(BetaBoot, "nSample") <- dim(mod$X)[1]
  for (r in 1:nboot) {
    isError <- TRUE
    while (isError) {
      catchError <- try({
        Yboot <- PairsBootOneRefFreeEwasModel(mod, pairID)
        bootFit <- RefFreeEwasModel(Yboot, mod$X, dim(mod$Lambda)[2],
          smallOutput = TRUE
        )
        BetaBoot[, , 1, r] <- bootFit$Beta
        BetaBoot[, , 2, r] <- bootFit$Bstar
      })
      isError <- inherits(catchError, "try-error")
    }
    if (r %% 10 == 0) {
      cat(r, "\n")
    }
  }
  class(BetaBoot) <- "BootRefFreeEwasModel"
  BetaBoot
}
