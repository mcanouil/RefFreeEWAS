#' Cell Mixture Projection (reference-based)
#'
#' Constrained linear projection for estimating cell mixture or related coefficients.
#'
#' @param Y Matrix (m CpGs x n Subjects) of DNA methylation beta values.
#' @param Xmat Matrix (m CpGs x K cell types) of cell-type specific methylomes.
#' @param nonnegative All coefficients >=0?
#' @param sumLessThanOne Coefficient rows should sum to less than one?
#' @param lessThanOne Every value should be less than one (but possibly sum to value greater than one)?
#' @details Function for projecting methylation values (Y) onto space of methyomes (Xmat), with various constraints.
#'          This is the reference-based method described in Houseman et al. (2012) and also appearing in the minfi package.
#' @return Projection coefficients resulting from constrained projection
#' @references Houseman EA, Accomando WP et al. DNA methylation arrays as surrogate measures of cell mixture distribution, BMC Bioinformatics, 2012.
#' @export
projectMix <- function(Y, Xmat, nonnegative = TRUE, sumLessThanOne = TRUE, lessThanOne = !sumLessThanOne) {
  nCol <- dim(Xmat)[2]
  nSubj <- dim(Y)[2]

  mixCoef <- matrix(0, nSubj, nCol)
  rownames(mixCoef) <- colnames(Y)
  colnames(mixCoef) <- colnames(Xmat)

  if (nonnegative) {
    if (sumLessThanOne) {
      Amat <- cbind(rep(-1, nCol), diag(nCol))
      b0vec <- c(-1, rep(0, nCol))
    } else if (lessThanOne) {
      Amat <- cbind(-diag(nCol), diag(nCol))
      b0vec <- c(rep(-1, nCol), rep(0, nCol))
    } else {
      Amat <- diag(nCol)
      b0vec <- rep(0, nCol)
    }

    for (i in 1:nSubj) {
      obs <- which(!is.na(Y[, i]))
      Dmat <- t(Xmat[obs, ]) %*% Xmat[obs, ]
      mixCoef[i, ] <- quadprog::solve.QP(Dmat, t(Xmat[obs, ]) %*% Y[obs, i], Amat, b0vec)$sol
    }
  } else {
    for (i in 1:nSubj) {
      obs <- which(!is.na(Y[, i]))
      Dmat <- t(Xmat[obs, ]) %*% Xmat[obs, ]
      mixCoef[i, ] <- solve(Dmat, t(Xmat[obs, ]) %*% Y[obs, i])
    }
  }

  return(mixCoef)
}
