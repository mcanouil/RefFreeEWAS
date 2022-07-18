#################################################################################
# RefFreeEwasModel: Reference-free cell-mixture-adjusted EWAS
#################################################################################
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
