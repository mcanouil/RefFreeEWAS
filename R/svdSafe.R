#################################################################################
# svdSafe: svd that traps errors and switches to QR when necessary
#################################################################################
svdSafe <- function(X) {
  sv <- try(svd(X), silent = TRUE)
  if (inherits(sv, "try-error")) {
    warning("SVD algorithm failed, using QR-decomposition instead")
    QR <- qr(X)
    sv <- list(d = rep(1, dim(X)[2]))
    sv$u <- qr.Q(QR)
    sv$v <- t(qr.R(QR))
  }
  sv
}
