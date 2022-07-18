#' Reference-Free EWAS Model
#'
#'  Reference-free method for conducting EWAS while deconvoluting DNA methylation arising as mixtures of cell types.
#'
#' @param Y Matrix of DNA methylation beta values (CpGs x subjects).  Missing values *are* supported.
#' @param X Design matrix (subjects x covariates).
#' @param K Latent variable dimension (d in Houseman et al., 2013, technical report)
#' @param smallOutput Smaller output?  (Should be FALSE if you intend to run bootstraps.)
#' @details Reference-free method for conducting EWAS while deconvoluting DNA methylation arising as mixtures of cell types.  This method is similar to surrogate variable analysis (SVA and ISVA), except that it makes additional use of a biological mixture assumption.  Returns mixture-adjusted Beta and unadjusted Bstar, as well as estimates of various latent quantities.
#' @return A list object of class \dQuote{RefFreeEwasModel}. The most important elements are Beta and Bstar.
#' @seealso `BootRefFreeEwasModel`
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
\keyword{svd
