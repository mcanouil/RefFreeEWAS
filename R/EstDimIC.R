#' Dimension estimation by AIC and BIC
#'
#' Method for estimating latent dimension by AIC and BIC.
#'
#' @param Rmat Residual matrix for which to estimate latent dimension.
#' @param Krange Vector of integers representing candidate dimensions to consider.
#' 
#' @details Method for estimating latent dimension by AIC and BIC.  
#'          Inferior to the RMT method in the isva package, but it appears here because it's mentioned in our paper.
#'
#' @return A list containing AIC and BIC for candidate dimensions, as well as the best dimension for each.
#'
#' @seealso `EstDimRMT`
#' @export
#' @examples
#' data(RefFreeEWAS)
#'
#' if (interactive()) {
#'   tmpDesign <- cbind(1, rfEwasExampleCovariate)
#'   tmpBstar <- rfEwasExampleBetaValues %*% tmpDesign %*% solve(t(tmpDesign)%*%tmpDesign)
#'
#'   EstDimIC(rfEwasExampleBetaValues-tmpBstar %*% t(tmpDesign))
#' }
EstDimIC <- function(
  Rmat,
  Krange = 0:25
) {
  N1 <- dim(Rmat)[1]
  N2 <- dim(Rmat)[2]
  svdRmat <- svdSafe(Rmat)
  nK <- length(Krange)
  tmpAIC <- tmpBIC <- rep(NA, nK)
  for (Ktest in Krange) {
    if (Ktest == 0) {
      tmpRminLU <- Rmat
    } else {
      tmpRminLU <- (Rmat -
        svdRmat$u[, 1:Ktest] %*% (svdRmat$d[1:Ktest] * t(svdRmat$v[, 1:Ktest])))
    }
    tmpSigSq <- apply(tmpRminLU * tmpRminLU, 1, sum) / N2
    tmpAIC[Ktest + 1] <- 2 * (N1 + Ktest * (N1 + N2)) + N1 * N2 + N2 * sum(log(tmpSigSq))
    tmpBIC[Ktest + 1] <- log(N2) * (N1 + Ktest * (N1 + N2)) + N1 * N2 + N2 * sum(log(tmpSigSq))
  }
  list(
    icTable = cbind(K = Krange, AIC = tmpAIC, BIC = tmpBIC),
    best = Krange[c(AIC = which.min(tmpAIC), BIC = which.min(tmpBIC))]
  )
}
