#' Cell Mixture Projection (reference-based)
#'
#' Constrained linear projection for estimating cell mixture or related coefficients.
#'
#' @param Y Matrix (m CpGs x n Subjects) of DNA methylation beta values
#' @param Xmat Matrix (m CpGs x K cell types) of cell-type specific methylomes
#' @param nonnegative All coefficients >=0?
#' @param sumLessThanOne Coefficient rows should sum to less than one?
#' @param lessThanOne Every value should be less than one (but possibly sum to value greater than one)?
#' @details Function for projecting methylation values (Y) onto space of methyomes (Xmat), with various constraints.
#'          This is the reference-based method described in Houseman et al. (2012) and also appearing in the minfi package.
#' @return Projection coefficients resulting from constrained projection
#' @references Houseman EA, Accomando WP et al. DNA methylation arrays as surrogate measures of cell mixture distribution,
BMC Bioinformatics, 2012.
