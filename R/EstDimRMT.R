#' Dimension estimation by Random Matrix Theory
#'
#' Method for estimating latent dimension by Random Matrix Theory.
#'
#' @param Rmat Residual matrix for which to estimate latent dimension.
#'
#' @details Method for estimating latent dimension by Random Matrix Theory.
#'          This function originated in the package isva, authored by A. Teschendorff.  
#'          Previous versions of RefFreeEWAS used the isva version of the function.  
#'          However, because of dependency issues in that package, the present
#'          version of RefFreeEWAS simply reproduces the function found in
#'          version 1.9 of isva and removes the dependency on the isva package.
#'          Documentation from isva:  Given a data matrix, it estimates
#'          the number of significant components of variation by comparing
#'          the observed distribution of spectral eigenvalues to the theoretical
#'          one under a Gaussian Orthogonal Ensemble (GOE).
#'          Specifically, a spectral decomposition of the data covariance matrix
#'          is performed and the number of eigenvalues larger than the
#'          theoretical maximum predicted by the GOE is taken as an estimate
#'          of the number of significant components.
#' @return A list with following objects:
#'         - cor Data covariance matrix.
#'         - dim Estimated intrinsic dimensionality of data.
#'         - estdens Empirical density of eigenvalues.
#'         - thdens Theoretical density of eigenvalues.
#' @export
#' @references 
#'  - Random matrix approach to cross correlations in financial data. Plerou et al. Physical Review E (2002), Vol.65.
#'  - Independent Surrogate Variable Analysis to deconvolve confounding factors in large-scale microarray profiling studies. Teschendorff AE, Zhuang JJ, Widschwendter M. Bioinformatics. 2011 Jun 1;27(11):1496-505.'
#' @examples
#' data(RefFreeEWAS)
#'
#' if (interactive()) {
#'   tmpDesign <- cbind(1, rfEwasExampleCovariate)
#'   tmpBstar <- rfEwasExampleBetaValues %*% tmpDesign %*% solve(t(tmpDesign)%*%tmpDesign)
#'   EstDimRMT(rfEwasExampleBetaValues-tmpBstar %*% t(tmpDesign))
#' }
EstDimRMT <- function(Rmat) {
  data.m <- Rmat
  M <- apply(data.m, 2, function(X) {
    (X - mean(X)) / sqrt(stats::var(X))
  })
  sigma2 <- stats::var(as.vector(M))
  Q <- nrow(data.m) / ncol(data.m)
  ns <- ncol(data.m)
  lambdaMAX <- sigma2 * (1 + 1 / Q + 2 * sqrt(1 / Q))
  lambdaMIN <- sigma2 * (1 + 1 / Q - 2 * sqrt(1 / Q))
  delta <- lambdaMAX - lambdaMIN
  roundN <- 3
  step <- round(delta / ns, roundN)
  while (step == 0) {
    roundN <- roundN + 1
    step <- round(delta / ns, roundN)
  }
  lambda.v <- seq(lambdaMIN, lambdaMAX, by = step)
  dens.v <- vector()
  ii <- 1
  for (i in lambda.v) {
    dens.v[ii] <- (Q / (2 * pi * sigma2)) * sqrt((lambdaMAX - i) * (i - lambdaMIN)) / i
    ii <- ii + 1
  }
  thdens.o <- list(
    min = lambdaMIN, max = lambdaMAX, step = step,
    lambda = lambda.v, dens = dens.v
  )
  C <- 1 / nrow(M) * t(M) %*% M
  eigen.o <- eigen(C, symmetric = TRUE)
  estdens.o <- stats::density(eigen.o$values,
    from = min(eigen.o$values),
    to = max(eigen.o$values), cut = 0
  )
  intdim <- length(which(eigen.o$values > thdens.o$max))
  evalues.v <- eigen.o$values
  return(list(
    cor = C, dim = intdim, estdens = estdens.o, thdens = thdens.o,
    evals = eigen.o$values
  ))
}
