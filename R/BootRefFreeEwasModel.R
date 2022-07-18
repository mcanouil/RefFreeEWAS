#################################################################################
# BootRefFreeEwasModel: Create several bootstrap samples from
#     reference-free cell-mixture-adjusted EWAS
#     and return Beta estimates
#################################################################################
# BootRefFreeEwasModel_old <- function(mod, nboot) {
#   BetaBoot <- array(NA, dim = c(dim(mod$Beta), 2, nboot))
#   dimnames(BetaBoot)[1:2] <- dimnames(mod$Beta)
#   dimnames(BetaBoot)[[3]] <- c("B", "B*")
#   dimnames(BetaBoot)[[4]] <- 1:nboot
#   attr(BetaBoot, "nSample") <- dim(mod$X)[1]

#   for (r in 1:nboot) {
#     isError <- TRUE
#     while (isError) {
#       catchError <- try({
#         Yboot <- BootOneRefFreeEwasModel(mod)
#         bootFit <- RefFreeEwasModel(Yboot, mod$X,
#           dim(mod$Lambda)[2],
#           smallOutput = TRUE
#         )
#         BetaBoot[, , 1, r] <- bootFit$Beta
#         BetaBoot[, , 2, r] <- bootFit$Bstar
#       })
#       isError <- inherits(catchError, "try-error")
#     }
#     if (r %% 10 == 0) cat(r, "\n")
#   }
#   class(BetaBoot) <- "BootRefFreeEwasModel"
#   BetaBoot
# }

BootRefFreeEwasModel <- function(mod, nboot) {
  BetaBoot <- replicate(
    n = nboot,
    expr = {
      isError <- TRUE
      niter <- 0
      while (isError & niter < 100) {
        isError <- inherits(try({
          bootFit <- RefFreeEwasModel(
            Y = BootOneRefFreeEwasModel(mod),
            X = mod$X,
            K = ncol(mod$Lambda),
            smallOutput = TRUE
          )
        }), "try-error")
        niter <- niter + 1
      }
      arr <- array(
        data = NA,
        dim = c(dim(mod$Beta), 2),
        dimnames = c(
          dimnames(mod$Beta),
          list(c("B", "B*"))
        )
      )
      arr[, , 1] <- bootFit[["Beta"]]
      arr[, , 2] <- bootFit[["Bstar"]]
      arr
    },
    simplify = "array"
  )
  dimnames(BetaBoot) <- c(
    dimnames(mod$Beta),
    list(c("B", "B*")),
    list(seq_len(nboot))
  )
  attr(BetaBoot, "nSample") <- nrow(mod$X)
  attr(BetaBoot, "class") <- "BootRefFreeEwasModel"
  BetaBoot
}

#' One Bootstrap sample for Reference-Free EWAS Model
#'
#' Bootstrap generation procedure for reference-free method for conducting EWAS
#' while deconvoluting DNA methylation arising as mixtures of cell types.
#'
#' @param mod model object of class RefFreeEwasModel
#'            (generated with `smallOutput = FALSE`).
#' @details Generates one bootstrapped data set for the reference-free method
#'          for conducting EWAS while deconvoluting DNA methylation arising as
#'          mixtures of cell types.  Typically not run by user.
#' @return A matrix representing a bootstrap sample of an DNA methylation assay matrix.
#' @seealso BootRefFreeEwasModel
#' @keywords deconvolution, DNA methylation, EWAS, surrogate variable,
#'           cell mixture
BootOneRefFreeEwasModel <- function(mod) {
  mod$Bstar %*% t(mod$X) +
    mod$dispersion * mod$E[, sample(seq_len(ncol(mod$E)), replace = TRUE)]
}

#' summary.BootRefFreeEwasModel
#'
#' Summary method for objects of type BootRefFreeEwasModel;
#' calculates bootstrap mean and standard deviation.
#'
#' @param object BootRefFreeEwasModel object to summarize.
#' @param ... Unused.
#' @details See `RefFreeEwasModel` for example.
summary.BootRefFreeEwasModel <- function(object, ...) {
  x <- object

  out <- array(NA, dim = c(dim(x)[1:3], 2))
  dimnames(out)[1:3] <- dimnames(x)[1:3]
  dimnames(out)[[4]] <- c("mean", "sd")

  out[, , , 1] <- apply(x, c(1:3), mean)
  out[, , , 2] <- apply(x, c(1:3), sd)

  class(out) <- "summaryBootRefFreeEwasModel"
  attr(out, "nBoot") <- dim(x)[4]
  attr(out, "nSample") <- attr(x, "nSample")

  out
}

print.summaryBootRefFreeEwasModel <- function(x, ...) {
  cat(attr(x, "nBoot"), "bootstrap samples, n =", attr(x, "nSample"), "subjects\n\nBeta Mean\n")
  print(x[1:6, , 1, 1], ...)
  cat("\nBeta Standard Deviation\n")
  print(x[1:6, , 1, 2], ...)
}

print.BootRefFreeEwasModel <- function(x, ...) {
  print(summary(x), ...)
}
