#' One Bootstrap Sample for Pairs
#'
#' Bootstrap generation procedure for sampling paired data (e.g. twin data)
#'
#' @param obs Observation ids (numeric vector).
#' @param pairID Pair IDs (one unique value per pair).
#' @details Generates one bootstrapped set of ids corresponding to pairs for the method for conducting EWAS while deconvoluting DNA methylation arising as mixtures of cell types.  Typically not run by user. 
#' @return   A vector of IDs corresponding to bootstrapped pairs
#' @seealso `BootRefFreeEwasModel`,`PairsBootRefFreeEwasModel`
\keyword{deconvolution
\keyword{DNA methylation
\keyword{EWAS
\keyword{surrogate variable
\keyword{cell mixture
