#' Bootstrap for Reference-Free EWAS Model
#'
#' Bootstrap procedure for reference-free method for conducting EWAS while deconvoluting DNA methylation arising as mixtures of cell types.
#'
#' @param mod model object of class RefFreeEwasModel (generated with smallOutput=FALSE).
#' @param nboot Number of bootstrap samples to generate
#' @details Generates the bootstrap samples for the reference-free method for conducting EWAS while deconvoluting DNA methylation arising as mixtures of cell types.
#' @return   An array object of class \dQuote{BootRefFreeEwasModel}. Bootstraps are generated for both Beta and Bstar.
#' @seealso `RefFreeEwasModel`
\examples{

data(RefFreeEWAS)

\dontrun{
  tmpDesign <- cbind(1, rfEwasExampleCovariate)
  tmpBstar <- (rfEwasExampleBetaValues %*% tmpDesign %*% solve(t(tmpDesign)%*%tmpDesign))

  EstDimRMT(rfEwasExampleBetaValues-tmpBstar %*% t(tmpDesign))$dim


test <- RefFreeEwasModel(
  rfEwasExampleBetaValues,
  cbind(1,rfEwasExampleCovariate),
  4)

testBoot <- BootRefFreeEwasModel(test,10)
summary(testBoot)


\keyword{deconvolution
\keyword{DNA methylation
\keyword{EWAS
\keyword{surrogate variable
\keyword{cell mixture
\keyword{bootstrap
