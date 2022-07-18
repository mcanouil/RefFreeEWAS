#' RefFreeCellMixArrayDevianceBoot
#'
#' Vector of bootstrapped deviances corresponding to an array of reference-free cell-mixture decompositions
#'
#' @param rfArray list of RefFreeCellMix objects (e.g. from RefFreeCellMixArray)
#' @param Y Methylation matrix on which x was based
#' @param EPSILON Minimum value of variance (zero variances will be reset to this value)
#' @param bootstrapIterations Number of RefFreeCellMix iterations to use in bootstrap
#' @details Vector of bootstrapped deviances corresponding to an array of reference-free cell-mixture decompositions,
#'          used to determine optimal number of cell types. This function returns one bootstrapped vector.  
#'          See `RefFreeCellMixArrayDevianceBoots` for more than one bootstrapped vector.
#'          The bootstrapped deviance is based on normal distribution applied to errors of Y after accounting for cell mixture effect, Mu   Omega^T. 
#'          See `RefFreeCellMixArray` for example.
