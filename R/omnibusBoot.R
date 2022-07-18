#' Bootstrap-based omnibus test of significance across all features
#'
#' Support for bootstrap-based omnibus test of significance accounting for correlation.
#'
#' @param est Vector of m estimates, one for each of m features.
#' @param boots Matrix (m x R) of bootstrap samples corresponding to the estimates
#' @param denDegFree Single number representing the denominator degrees-of-freedom for computing p-values
#' @details Returns one omnibus p-value based on Kolmogorov-Smirnov distance from a uniform distribution
#' @return A single number representing the p-value for the omnibus test over all features.
#' @seealso `RefFreeEwasModel`
#' @keywords bootstrap,omnibus,olmogorov-smirnov
#' @export
#' @examples
#' data(RefFreeEWAS)
#'
#' test <- RefFreeEwasModel(
#'   rfEwasExampleBetaValues,
#'   cbind(1,rfEwasExampleCovariate),
#'   4
#' )
#'
#' testBoot <- BootRefFreeEwasModel(test,10)
#' summary(testBoot)
#' omnibusBoot(test$Beta[,2], testBoot[,2,"B",], -diff(dim(test$X))) 
#' omnibusBoot(test$Bstar[,2], testBoot[,2,"B*",], -diff(dim(test$X)))
omnibusBoot <- function(est, boots, denDegFree) {
  nFeature <- length(est)
  se <- apply(boots, 1, stats::sd)
  pv <- 2 * stats::pt(-abs(est) / se, denDegFree)
  pvNull <- 2 * stats::pt(-(1 / se) * abs(sweep(boots, 1, apply(boots, 1, mean), "-")), denDegFree + 1)

  ks <- max(abs(sort(pv) - (1:nFeature - 0.5) / nFeature))
  ksNull <- apply(pvNull, 2, function(u) max(abs(sort(u) - (1:nFeature - 0.5) / nFeature)))
  mean(ks <= ksNull)
}
