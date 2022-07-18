#' Reference-Free EWAS Model
#'
#' Reference-free method for conducting EWAS while deconvoluting DNA methylation arising as mixtures of cell types.
#'
#' @param Y Matrix of DNA methylation beta values (CpGs x subjects).  Missing values *are* supported.
#' @param X Design matrix (subjects x covariates).
#' @param K Latent variable dimension (d in Houseman et al., 2013, technical report)
#' @param smallOutput Smaller output?  (Should be FALSE if you intend to run bootstraps.)
#' @details Reference-free method for conducting EWAS while deconvoluting DNA methylation arising as mixtures of cell types.  This method is similar to surrogate variable analysis (SVA and ISVA), except that it makes additional use of a biological mixture assumption.  Returns mixture-adjusted Beta and unadjusted Bstar, as well as estimates of various latent quantities.
#' @return A list object of class \dQuote{RefFreeEwasModel}. The most important elements are Beta and Bstar.
#' @seealso `BootRefFreeEwasModel`
#' @keywords deconvolution,DNA methylation,EWAS,surrogate variable,cell mixture,svd
#' @export
#' @examples
#' data(RefFreeEWAS)
#' 
#' if (interactive()) {
#'   tmpDesign <- cbind(1, rfEwasExampleCovariate)
#'   tmpBstar <- (rfEwasExampleBetaValues %*% tmpDesign %*% solve(t(tmpDesign)%*%tmpDesign))
#'   
#'   EstDimRMT(rfEwasExampleBetaValues-tmpBstar %*% t(tmpDesign))$dim  
#' }
#' test <- RefFreeEwasModel(
#'   rfEwasExampleBetaValues, 
#'   cbind(1,rfEwasExampleCovariate),
#'   4
#' )
#' 
#' testBoot <- BootRefFreeEwasModel(test,10)
#' summary(testBoot)
RefFreeEwasModel <- function(
  Y,
  X,
  K,
  smallOutput = FALSE
) {
  n1 <- dim(Y)[1]
  n2 <- dim(X)[1]
  pdim <- dim(X)[2]

  noMiss <- !apply(is.na(Y), 1, any)
  allObs <- all(noMiss)

  HX <- solve(t(X) %*% X)
  PX <- X %*% HX

  if (allObs) { # If no missing value
    Bstar <- Y %*% PX

    muStar <- Bstar %*% t(X)
    Estar <- Y - muStar

    sstar <- sqrt(apply(Estar * Estar, 1, sum) / (n2 - pdim))

    BetaE <- cbind(Bstar, Estar)
    svdStar <- svdSafe(BetaE)

    Lambda <- t(svdStar$d[1:K] * t(svdStar$u[, 1:K]))
    U <- svdStar$v[, 1:K]
  } else { # If missing values, do as much as possible on one fell swoop,
    # and for the rest, do on a CpG-by-CpG basis
    Bstar <- matrix(NA, n1, pdim)
    degfree <- rep(n2 - pdim, n1)
    Bstar[noMiss, ] <- Y[noMiss, ] %*% PX
    whichMiss <- which(!noMiss)
    nMiss <- length(whichMiss)

    for (j in 1:nMiss) {
      jj <- whichMiss[j]
      mflag <- !is.na(Y[jj, ])
      Xj <- X[mflag, , drop = FALSE]
      HXj <- solve(t(Xj) %*% Xj)
      PXj <- Xj %*% HXj
      Bstar[jj, ] <- Y[jj, mflag, drop = FALSE] %*% PXj
      degfree[jj] <- sum(mflag) - pdim
    }

    muStar <- Bstar %*% t(X)
    Estar <- Y - muStar

    sstar <- sqrt(apply(Estar * Estar, 1, sum, na.rm = TRUE) / degfree)
    BetaE <- cbind(Bstar, Estar)
    svdStar <- svdSafe(BetaE[noMiss, ])

    Lambda <- matrix(NA, n1, K)
    Lambda[noMiss, ] <- t(svdStar$d[1:K] * t(svdStar$u[, 1:K]))
    U <- svdStar$v[, 1:K]
    for (j in 1:nMiss) {
      jj <- whichMiss[j]
      mflag <- c(rep(TRUE, pdim), !is.na(Y[jj, ]))
      Uj <- U[mflag, , drop = FALSE]
      Lambda[jj, ] <- solve(t(Uj) %*% Uj, t(Uj) %*% BetaE[jj, mflag])
    }
  }

  LambdaProjBstar <- solve(t(Lambda) %*% Lambda, t(Lambda) %*% Bstar)
  Beta <- Bstar - Lambda %*% LambdaProjBstar

  out <- list(Bstar = Bstar, Beta = Beta, sstar = sstar, Lambda = Lambda, U = U, d = svdStar$d)

  if (!smallOutput) {
    muStar <- ifelse(muStar < 0.00001, 0.00001, muStar)
    muStar <- ifelse(muStar > 0.99999, 0.99999, muStar)
    out$dispersion <- sqrt(muStar * (1 - muStar))
    out$E <- Estar / out$dispersion
    out$X <- X
  }
  out
  class(out) <- "RefFreeEwasModel"
  out
}

#' print.RefFreeEwasModel
#'
#' Print method for objects of type RefFreeEwasModel
#'
#' @param x RefFreeEwasModel object to print
#' @param ... Unused..
#' @details See `RefFreeEwasModel` for example.
#' @export
print.RefFreeEwasModel <- function(x, ...) {
  cat("Reference Free EWAS Model\n\n")
  cat("Assay matrix: ", dim(x$Beta)[1], " features\n")
  if (!is.null(x$X)) {
    cat("Design matrix: ", dim(x$X)[1],
      " subjects x ", dim(x$X)[2], " covariates\n\n",
      sep = ""
    )
  } else {
    cat("(small output version)\n\n")
  }

  cat(dim(x$Lambda)[2], " latent variables\n\n")
}
