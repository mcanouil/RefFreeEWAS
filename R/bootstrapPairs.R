#' One Bootstrap Sample for Pairs
#'
#' Bootstrap generation procedure for sampling paired data (e.g. twin data)
#'
#' @param obs Observation ids (numeric vector).
#' @param pairID Pair IDs (one unique value per pair).
#' @details Generates one bootstrapped set of ids corresponding to pairs for
#'          the method for conducting EWAS while deconvoluting DNA methylation
#'          arising as mixtures of cell types.  
#'          Typically not run by user.
#' @return   A vector of IDs corresponding to bootstrapped pairs
#' @seealso `BootRefFreeEwasModel`, `PairsBootRefFreeEwasModel`
#' @keywords deconvolution,DNA methylation,EWAS,surrogate variable,cell mixture
#' @export
bootstrapPairs <- function(obs, pairID) {
  pairIx <- split(obs, pairID)
  nobsTab <- table(sapply(pairIx, length))
  if (length(nobsTab) > 1) {
    stop("All clusters must have the same number of observations.\n")
  }

  n <- length(pairIx)
  pairBoot <- pairIx[sample(1:n, n, replace = TRUE)]

  obsBoot <- unlist(lapply(pairBoot, function(u) sample(u, replace = FALSE)))
  names(obsBoot) <- NULL
  obsBoot
}
