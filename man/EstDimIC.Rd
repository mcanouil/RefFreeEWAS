#' Dimension estimation by AIC and BIC
#'
#' Method for estimating latent dimension by AIC and BIC.
#'
#' @param Rmat Residual matrix for which to estimate latent dimension.
#' @param Krange Vector of integers representing candidate dimensions to consider
#' @details Method for estimating latent dimension by AIC and BIC.  
#'          Inferior to the RMT method in the isva package, but it appears here because it's mentioned in our paper.
#' @return A list containing AIC and BIC for candidate dimensions, as well as the best dimension for each. 
#' @seealso `EstDimRMT`
\examples{
data(RefFreeEWAS)

\dontrun{
  tmpDesign <- cbind(1, rfEwasExampleCovariate)
  tmpBstar <- rfEwasExampleBetaValues %*% tmpDesign %*% solve(t(tmpDesign)%*%tmpDesign)

  EstDimIC(rfEwasExampleBetaValues-tmpBstar %*% t(tmpDesign))
