#' deviance.RefFreeCellMix
#'
#' Deviance method for objects of type RefFreeCellMix.
#'
#' @param object RefFreeCellMix object to summarize
#' @param Y Methylation matrix on which x was based
#' @param Y.oob Alternate ("out-of-box") methylation matrix for which to calculate deviance, based on x
#' @param EPSILON Minimum value of variance (zero variances will be reset to this value)
#' @param bootstrapIterations Number of RefFreeCellMix iterations to use in bootstrap (see details)
#' @param bootstrapIndices Bootstrap indices (see details)
#' @param ... (Unused).
#' @details Deviance based on normal distribution applied to errors of Y after accounting for cell mixture effect,
#'          Mu Omega^T. Since RefFreeCellMix does not save the original data Y in the resulting object x, Y must
#'          be supplied here.  
#'          However, deviance may be calculated for an alternative "out-of-bag" methylation matrix, Y.oob.  
#'          If bootstrapIterations=0, this is what is done.  If bootstrapIterations>0, then x$Mu is used to
#'          initialize a new value of x via RefFreeCellMix executed on a bootstrap sample of Y 
#'          with the number of indicated iterations.  If bootstrapIndices is provided, the bootstrap will be based on
#'          these indices, otherwise the indices will be sampled randomly with replacement from 1:ncol(Y).
#'          See `RefFreeCellMix` for example.
