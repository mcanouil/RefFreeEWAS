#' Bootstrap for Reference-Free EWAS Model
#'
#' Bootstrap procedure for reference-free method for conducting EWAS while
#' deconvoluting DNA methylation arising as mixtures of cell types.
#'
#' @param mod model object of class RefFreeEwasModel (generated with smallOutput=FALSE).
#' @param nboot Number of bootstrap samples to generate.
#' @details Generates the bootstrap samples for the reference-free method
#'          for conducting EWAS while deconvoluting DNA methylation arising
#'          as mixtures of cell types.
#' @return An array object of class \dQuote{BootRefFreeEwasModel}.
#'         Bootstraps are generated for both Beta and Bstar.
#' @seealso `RefFreeEwasModel`
#' @keywords deconvolution,DNA methylation,EWAS,surrogate variable,cell mixture,bootstrap
#' @export
#' @examples
#'
#' data(RefFreeEWAS)
#'
#' if (interactive()) {
#'   tmpDesign <- cbind(1, rfEwasExampleCovariate)
#'   tmpBstar <- (rfEwasExampleBetaValues %*% tmpDesign %*% solve(t(tmpDesign)%*%tmpDesign))
#'
#'   EstDimRMT(rfEwasExampleBetaValues-tmpBstar %*% t(tmpDesign))$dim
#' }
#'
#' test <- RefFreeEwasModel(
#'   rfEwasExampleBetaValues,
#'   cbind(1,rfEwasExampleCovariate),
#'   4
#' )
#'
#' testBoot <- BootRefFreeEwasModel(test,10)
#' summary(testBoot)
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
#' @keywords deconvolution, DNA methylation, EWAS, surrogate variable, cell mixture
#' @export
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
#' @export
summary.BootRefFreeEwasModel <- function(object, ...) {
  x <- object

  out <- array(NA, dim = c(dim(x)[1:3], 2))
  dimnames(out)[1:3] <- dimnames(x)[1:3]
  dimnames(out)[[4]] <- c("mean", "sd")

  out[, , , 1] <- apply(x, c(1:3), mean)
  out[, , , 2] <- apply(x, c(1:3), stats::sd)

  class(out) <- "summaryBootRefFreeEwasModel"
  attr(out, "nBoot") <- dim(x)[4]
  attr(out, "nSample") <- attr(x, "nSample")

  out
}

#' print.summaryBootRefFreeEwasModel
#'
#' Print method for objects of type summaryBootRefFreeEwasModel
#'
#' @param x summaryBootRefFreeEwasModel object to print
#' @param ... Unused..
#' @details See `RefFreeEwasModel` for example.
#' @export
print.summaryBootRefFreeEwasModel <- function(x, ...) {
  cat(attr(x, "nBoot"), "bootstrap samples, n =", attr(x, "nSample"), "subjects\n\nBeta Mean\n")
  print(x[1:6, , 1, 1], ...)
  cat("\nBeta Standard Deviation\n")
  print(x[1:6, , 1, 2], ...)
}

#' print.BootRefFreeEwasModel
#'
#' Print method for objects of type BootRefFreeEwasModel
#'
#' @param x BootRefFreeEwasModel object to print
#' @param ... Unused..
#' @details See `RefFreeEwasModel` for example.
#' @export
print.BootRefFreeEwasModel <- function(x, ...) {
  print(summary(x), ...)
}
