#' Reference-Free Cell Mixture Projection
#'
#' Reference-free cell-mixture decomposition of DNA methylation data set
#'
#' @param Y Matrix (m CpGs x n Subjects) of DNA methylation beta values
#' @param mu0 Matrix (m CpGs x K cell types) of *initial* cell-type specific methylomes
#' @param K Number of cell types (ignored if mu0 is provided)
#' @param iters Number of iterations to execute
#' @param Yfinal Matrix (m* CpGs x n Subjects) of DNA methylation beta values on which to base final methylomes
#' @param verbose Report summary of errors after each iteration?
#' @details Reference-free decomposition of DNA methylation matrix into cell-type distributions and cell-type
#'          methylomes, Y = Mu Omega^T.  Either an initial estimate of Mu must be provided, or else the number of cell types K,
#'          in which case RefFreeCellMixInitialize will be used to initialize.  Note that the decomposition will be based on Y,
#'          but Yfinal (=Y by default) will be used to determine the final value of Mu based on the last iterated value of Omega.
#' @return Object of S3 class RefFreeCellMix, containing the last iteration of Mu and Omega.
#' @references Houseman, E. Andres, Kile, Molly L., Christiani, David C., et al. Reference-free deconvolution of DNA methylation data and mediation by cell composition effects. BMC bioinformatics, 2016, vol. 17, no 1, p. 259.
#' @seealso `RefFreeCellMixInitialize`
#' @export
#' @examples
#' data(HNSCC)
#' 
#' # Typical use
#' Y.shortTest <- Y.HNSCC.averageBetas[1:500,]
#' Y.shortTest.final <- Y.HNSCC.averageBetas[1:1000,]
#' testArray1  <- RefFreeCellMixArray(Y.shortTest,Klist=1:3,iters=5,Yfinal=Y.shortTest.final)
#' testArray1
#' lapply(testArray1,summary)
#' sapply(testArray1,stats::deviance,Y=Y.shortTest.final)
#' 
#' # Example with explicit initialization
#' testKeq2  <- RefFreeCellMix(Y.shortTest, mu0=RefFreeCellMixInitialize(Y.shortTest,K=2))
#' testKeq2
#' head(testKeq2$Mu)
#' head(testKeq2$Omega)
RefFreeCellMix <- function(Y, mu0 = NULL, K = NULL, iters = 10, Yfinal = NULL, verbose = TRUE) {
  if (is.null(mu0)) {
    if (K == 1) {
      if (!is.null(Yfinal)) Y <- Yfinal
      n <- dim(Y)[2]

      mu <- matrix(apply(Y, 1, mean, na.rm = TRUE), ncol = 1)
      omega <- matrix(1, n, 1)
      o <- list(Mu = mu, Omega = omega)
      class(o) <- "RefFreeCellMix"
      return(o)
    } else {
      mu0 <- RefFreeCellMixInitialize(Y, K = K, method = "ward")
    }
  }

  incrementalChangeSummary <- list()
  for (i in 1:iters) {
    flag <- !apply(is.na(mu0), 1, any)
    omega <- projectMix(Y[flag, ], mu0[flag, ])
    mu <- projectMix(t(Y), omega, sumLessThanOne = FALSE)
    incrementalChangeSummary[[i]] <- summary(abs(as.vector(mu - mu0)))
    if (verbose) print(incrementalChangeSummary[[i]])
    mu0 <- mu
  }
  if (!is.null(Yfinal)) {
    mu <- projectMix(t(Yfinal), omega, sumLessThanOne = FALSE)
  }

  o <- list(Mu = mu, Omega = omega, incrementalChangeSummary = incrementalChangeSummary)
  class(o) <- "RefFreeCellMix"
  o
}

#' print.RefFreeCellMix
#'
#' Print method for objects of type RefFreeCellMix
#'
#' @param x RefFreeCellMix object to print
#' @param ... Unused..
#' @details See `RefFreeCellMix` for example.
#' @export
print.RefFreeCellMix <- function(x, ...) {
  cat("Reference Free Deconvolution\n\n")
  cat("Mu: ", dim(x$Mu)[1], " cpgs x ", dim(x$Mu)[2], "cell types\n")
  cat("Omega :", dim(x$Omega)[1], " subjects x ", dim(x$Omega)[2], "cell types\n")
}

#' summary.RefFreeCellMix
#'
#' Summary method for objects of type RefFreeCellMix.
#' 
#' @param object RefFreeCellMix object to summarize.
#' @param ... Unused.
#' @details See `RefFreeCellMix()` for example.
#' @export
summary.RefFreeCellMix <- function(object, ...) {
  list(Mu = apply(object$Mu, 2, summary), Omega = apply(object$Omega, 2, summary), MuCorr = stats::cor(object$Mu))
}

#' Initialize Reference-Free Cell Mixture Projection
#'
#' Array of reference-free cell-mixture decompositions of a DNA methylation data set
#'
#' @param Y Matrix (m CpGs x n Subjects) of DNA methylation beta values
#' @param Klist List of K values (each K = assumed number of cell types)
#' @param iters Number of iterations to execute for each value of K
#' @param Yfinal Matrix (m* CpGs x n Subjects) of DNA methylation beta values on which to base final methylomes
#' @param verbose Report summary of errors after each iteration for each fit?
#' @param dist.method Method for calculating distance matrix for methylome initialization
#' @param ... Additional parameters for hclust function for methylome initialization
#' @details List of Reference-free decompositions for a range of K values.  For each value of K, the decomposition is initialized by hierarchical clutering as specified by the parameters dist.method, etc.  Note that for each K, the decomposition will be based on Y,
#'          but Yfinal (=Y by default) will be used to determine the final value of Mu based on the last iterated value of Omega.
#' @return List, each element is an object of S3 class RefFreeCellMix, containing the last iteration of Mu and Omega.
#' @references Houseman, E. Andres, Kile, Molly L., Christiani, David C., et al. Reference-free deconvolution of DNA methylation data and mediation by cell composition effects. BMC bioinformatics, 2016, vol. 17, no 1, p. 259.
#' @seealso `RefFreeCellMix`
#' @export
#' @examples
#' data(HNSCC)
#' Y.shortTest <- Y.HNSCC.averageBetas[1:500,]
#' testArray2  <- RefFreeCellMixArray(Y.shortTest,Klist=1:5,iters=5)
#' sapply(testArray2,stats::deviance,Y=Y.shortTest)
#' 
#' if (interactive()) {
#'   testBootDevs <- RefFreeCellMixArrayDevianceBoots(testArray2,Y.shortTest,R=10)
#'   testBootDevs
#'   apply(testBootDevs[-1,],2,mean,trim=0.25)
#'   which.min(apply(testBootDevs[-1,],2,mean,trim=0.25))
#' }
RefFreeCellMixArray <- function(Y, Klist = 1:5, iters = 10, Yfinal = NULL, verbose = FALSE, dist.method = "euclidean", ...) {
  D <- stats::dist(t(Y), method = dist.method)
  hc <- stats::hclust(D, ...)

  rfcmArray <- list()
  nK <- length(Klist)
  for (r in 1:nK) {
    cat("Fitting K =", Klist[r], "\n")
    if (Klist[r] == 1) {
      rfcmArray[[r]] <- RefFreeCellMix(Y, K = 1, Yfinal = Yfinal, iters = iters)
    } else {
      rfcmArray[[r]] <- RefFreeCellMix(Y,
        mu0 = RefFreeCellMixInitialize(Y, K = Klist[r], Y.Cluster = hc),
        Yfinal = Yfinal, verbose = verbose, iters = iters
      )
    }
  }
  names(rfcmArray) <- Klist
  rfcmArray
}

#' deviance.RefFreeCellMix
#'
#' Deviance method for objects of type RefFreeCellMix.
#'
#' @param object RefFreeCellMix object to summarize
#' @param Y Methylation matrix on which x was based
#' @param Y.oob Alternate ("out-of-box") methylation matrix for which to calculate deviance, based on x
#' @param EPSILON Minimum value of variance (zero variances will be reset to this value)
#' @param bootstrapIterations Number of RefFreeCellMix iterations to use in bootstrap (see details)
#' @param bootstrapIndices Bootstrap indices (see details)
#' @param ... Unused..
#'
#' @details Deviance based on normal distribution applied to errors of Y after 
#'          accounting for cell mixture effect,
#'          $Mu Omega^T$. Since RefFreeCellMix does not save the original data Y
#'          in the resulting object x, Y must be supplied here.  
#'          However, deviance may be calculated for an alternative "out-of-bag"
#'          methylation matrix, Y.oob.  
#'          If `bootstrapIterations=0`, this is what is done.  
#'          If `bootstrapIterations>0`, then `x$Mu` is used to
#'          initialize a new value of x via RefFreeCellMix executed on a bootstrap sample of Y 
#'          with the number of indicated iterations.  
#'          If bootstrapIndices is provided, the bootstrap will be based on
#'          these indices, otherwise the indices will be sampled randomly with replacement from 1:ncol(Y).
#'          See `RefFreeCellMix` for example.
#' @export
deviance.RefFreeCellMix <- function(object, Y, Y.oob = NULL, EPSILON = 1E-9,
  bootstrapIterations = 0, bootstrapIndices = NULL, ...) {
  N <- dim(Y)[2]
  if (bootstrapIterations > 0) { # Do the bootstrap and replace x (but initialize with x$Mu)
    if (is.null(bootstrapIndices)) {
      boots <- sample(1:N, N, replace = TRUE)
    } else {
      boots <- bootstrapIndices
    }
    Y.oob <- Y[, -unique(boots)]
    Y <- Y[, boots]
    if (dim(object$Mu)[2] == 1) {
      object <- RefFreeCellMix(Y, K = 1, iters = bootstrapIterations, verbose = FALSE)
    } else {
      object <- RefFreeCellMix(Y, mu0 = object$Mu, iters = bootstrapIterations, verbose = FALSE)
    }
  }

  Y.mu <- object$Mu %*% t(object$Omega)
  R <- Y - Y.mu
  Y.n <- apply(!is.na(Y), 1, sum)
  Y.SSQ <- apply(R * R, 1, sum, na.rm = TRUE)
  logSigma2 <- log(pmax(EPSILON, Y.SSQ)) - log(Y.n)

  if (!is.null(Y.oob)) {
    Omega.oob <- projectMix(Y.oob, object$Mu)
    Y.mu <- object$Mu %*% t(Omega.oob)
    R.oob <- Y.oob - Y.mu
    n.oob <- apply(!is.na(Y.oob), 1, sum)
    SSQ.oob <- apply(R.oob * R.oob, 1, sum, na.rm = TRUE)
    N <- dim(Y.oob)[2]
  } else {
    SSQ.oob <- Y.SSQ
    n.oob <- Y.n
  }

  sum(n.oob * log(2 * pi) + n.oob * logSigma2 + SSQ.oob / exp(logSigma2)) / N
}

#' Initialize Reference-Free Cell Mixture Projection
#'
#' Initializes the methylome matrix "Mu" for RefFreeCellMix
#'
#' @param Y Matrix (m CpGs x n Subjects) of DNA methylation beta values
#' @param K Number of cell types
#' @param Y.Distance Distance matrix (object of class "dist") to use for clustering.
#' @param Y.Cluster Hiearchical clustering object (from hclust function)
#' @param largeOK OK to calculate distance matrix for large number of subjects? (See details.)
#' @param dist.method Method for calculating distance matrix
#' @param \dots Additional parameters for hclust function
#' @details Initializes the methylome matrix "Mu" for RefFreeCellMix by computing the mean methylation (from Y)
#' over K clusters of Y, determined by the Y.Cluster object.  If Y.Cluster object does not exist, it will be 
#' created from Y. Distance (using additional clustering parameters if supplied).  If Y.Distance does not exist,
#' it will be created from t(Y).  As a protection against attempting to fit a very large distance matrix, the
#' program will stop if the number of columns of Y is > 2500, unless largeOK is explicitly set to TRUE.
#' @return An m x K matrix of mean methylation values.
#' @export
RefFreeCellMixInitialize <- function(Y, K = 2, Y.Distance = NULL, Y.Cluster = NULL,
  largeOK = FALSE, dist.method = "euclidean", ...) {
  if (!is.matrix(Y) | !is.numeric(Y)) {
    stop("Y is not a numeric matrix\n")
  }
  n <- dim(Y)[2]

  if (is.null(Y.Cluster)) {
    if (is.null(Y.Distance)) {
      if (n > 2500 & !largeOK) {
        stop("Y has a large number of subjects!  If this is what you really want, change 'largeOK' to TRUE\n")
      }
      Y.Distance <- stats::dist(t(Y), method = dist.method)
    }
    Y.Cluster <- stats::hclust(Y.Distance, ...)
  }

  classes <- stats::cutree(Y.Cluster, K)
  s <- split(1:n, classes)

  sapply(s, function(u) apply(Y[, u, drop = FALSE], 1, mean, na.rm = TRUE))
}

#' RefFreeCellMixArrayDevianceBoot
#'
#' Vector of bootstrapped deviances corresponding to an array of reference-free cell-mixture decompositions
#'
#' @param rfArray list of RefFreeCellMix objects (e.g. from RefFreeCellMixArray)
#' @param Y Methylation matrix on which x was based
#' @param EPSILON Minimum value of variance (zero variances will be reset to this value)
#' @param bootstrapIterations Number of RefFreeCellMix iterations to use in bootstrap
#' @details Vector of bootstrapped deviances corresponding to an array of reference-free cell-mixture decompositions,
#'          used to determine optimal number of cell types. This function returns one bootstrapped vector.  
#'          See `RefFreeCellMixArrayDevianceBoots` for more than one bootstrapped vector.
#'          The bootstrapped deviance is based on normal distribution applied to errors of Y after accounting for cell mixture effect, Mu   Omega^T. 
#'          See `RefFreeCellMixArray` for example.
#' @export
RefFreeCellMixArrayDevianceBoot <- function(rfArray, Y, EPSILON = 1E-9, bootstrapIterations = 5) {
  N <- dim(Y)[2]
  boots <- sample(1:N, N, replace = TRUE)
  sapply(rfArray, deviance.RefFreeCellMix,
    Y = Y, EPSILON = EPSILON,
    bootstrapIterations = bootstrapIterations, bootstrapIndices = boots
  )
}

#' RefFreeCellMixArrayDevianceBoots
#'
#' Matrix of bootstrapped deviances corresponding to an array of reference-free cell-mixture decompositions
#'
#' @param rfArray list of RefFreeCellMix objects (e.g. from RefFreeCellMixArray)
#' @param Y Methylation matrix on which x was based
#' @param R Number of bootstrapped vectors to return
#' @param EPSILON Minimum value of variance (zero variances will be reset to this value)
#' @param bootstrapIterations Number of RefFreeCellMix iterations to use in bootstrap
#' @details Matrix (multiple vectors) of bootstrapped deviances corresponding to an array of reference-free cell-mixture decompositions,
#'          used to determine optimal number of cell types. This function returns one bootstrapped vector.  
#'          The bootstrapped deviance is based on normal distribution applied to errors of Y after accounting for cell mixture effect, Mu   Omega^T. 
#'          See `RefFreeCellMixArray` for example.
#' @export
RefFreeCellMixArrayDevianceBoots <- function(rfArray, Y, R = 5, EPSILON = 1E-9, bootstrapIterations = 5) {
  dv <- sapply(rfArray, stats::deviance, Y = Y)
  nK <- length(dv)
  devs <- matrix(NA, R, nK)
  for (r in 1:R) {
    if (r %% 10 == 0) cat("Bootstrap", r, "\n")
    devs[r, ] <- RefFreeCellMixArrayDevianceBoot(rfArray, Y,
      EPSILON = EPSILON, bootstrapIterations = bootstrapIterations
    )
  }
  out <- rbind(dv, devs)
  rownames(out) <- 0:R
  out
}

#' Simple imputation method based on row-mean
#'
#' Simple method for imputing missing values by row-mean
#'
#' @param Y Matrix to impute.
#''
#' @return Matrix with missing values replaced by imputed values.
#' @export
ImputeByMean <- function(Y) {
  Yimpute <- Y
  whichMiss <- apply(is.na(Y), 1, which)
  nmiss <- sapply(whichMiss, sum)
  if (any(nmiss > 0)) {
    for (i in which(nmiss > 0)) {
      Yimpute[i, whichMiss[[i]]] <-
        mean(Y[i, -whichMiss[[i]]])
    }
  }
  Yimpute
}

#' SVD with missing values
#'
#' Compute singular value decomposition on a matrix with missing values,
#' using a naive/simple method for imputing missing values by row-mean.
#'
#' @param Y Matrix for which to compute SVD.
#' @details Computes singular value decomposition on a matrix with
#'          missing values, using a naive/simple method for imputing missing
#'          values by row-mean.  Not recommended for matrices with very large
#'          numbers of missing values.
#' @return singular value decomposition (as returned by svd function).
#' @export
SVDwithMissing <- function(Y) {
  svd(ImputeByMean(Y))
}

#' Reference-Free Cell Mixture Projection - Custom Initialization
#'
#' Array of reference-free cell-mixture decompositions of a DNA methylation data set, with custom initialization
#'
#' @param Y Matrix (m CpGs x n Subjects) of DNA methylation beta values
#' @param mu.start matrix of starting values for Mu: number of columns must be at least the maximum in Klist
#' @param Klist List of K values (each K = assumed number of cell types)
#' @param iters Number of iterations to execute for each value of K
#' @param Yfinal Matrix (m* CpGs x n Subjects) of DNA methylation beta values on which to base final methylomes
#' @param verbose Report summary of errors after each iteration for each fit?
#' @details List of Reference-free decompositions for a range of K values.  
#'          For each value of K, the decomposition is initialized by using the first K columns of mu.start.
#'          Note that for each K, the decomposition will be based on Y, but Yfinal (=Y by default) will be used to determine the final value of Mu based on the last iterated value of Omega.
#' @return List, each element is an object of S3 class RefFreeCellMix, containing the last iteration of Mu and Omega.
#' @seealso `RefFreeCellMix`, `RefFreeCellMixInitializeBySVD`
#' @export
RefFreeCellMixArrayWithCustomStart <- function(Y, mu.start,
  Klist = 1:5, iters = 10, Yfinal = NULL, verbose = FALSE) {
  rfcmArray <- list()
  nK <- length(Klist)
  for (r in 1:nK) {
    cat("Fitting K =", Klist[r], "\n")
    if (Klist[r] == 1) {
      rfcmArray[[r]] <- RefFreeCellMix(Y, K = 1, Yfinal = Yfinal, iters = iters)
    } else {
      rfcmArray[[r]] <- RefFreeCellMix(Y,
        mu0 = mu.start[, 1:Klist[r]],
        Yfinal = Yfinal, verbose = verbose, iters = iters
      )
    }
  }
  names(rfcmArray) <- Klist
  rfcmArray
}

#' Initialize Reference-Free Cell Mixture Projection by SVD
#'
#' Initialize Reference-Free Cell Mixture Projection by SVD
#'
#' @param Y Matrix (m CpGs x n Subjects) of DNA methylation beta values
#' @param type See details
#' @details This method initializes the reference-free cell mixture deconvolution using an ad-hoc method based on singular value decomposition.  Type=1 will attempt to discretize Mu to 0/1, Type=2 will attempt to find a continuous range using column ranks. However, neither of these strategies is guaranteed to result in stable starting values for K larger than the "true" value of K.
#' @return Matrix of starting values for Mu.
#' @seealso `RefFreeCellMix`, `RefFreeCellMixArrayWithCustomStart`
#' @examples
#' data(HNSCC)
#' Y.shortTest <- Y.HNSCC.averageBetas[1:600, ]
#' mu.start.svd <- RefFreeCellMixInitializeBySVD(Y.shortTest)
#' testArray2  <- RefFreeCellMixArrayWithCustomStart(
#'   Y.shortTest,
#'   mu.start=mu.start.svd,
#'   Klist=1:3,iters=5
#' )
#' sapply(testArray2,stats::deviance,Y=Y.shortTest)
#' 
#' if (interactive()) {
#'   testBootDevs <- RefFreeCellMixArrayBySVDDevianceBoots(testArray2,Y.shortTest,R=10)
#' 
#'   testBootDevs
#'   apply(testBootDevs[-1,],2,mean,trim=0.25)
#'   which.min(apply(testBootDevs[-1,],2,mean,trim=0.25))
#' }
#' @export
RefFreeCellMixInitializeBySVD <- function(Y, type = 1) {
  Y.svd <- SVDwithMissing(Y)

  nn <- ncol(Y.svd$u)
  Y.svd.sign <- sapply(1:nn, function(i) sign(mean(sign(Y.svd$v[, i]))))
  Y.svd.sign[Y.svd.sign == 0] <- 1
  Y.svd.u <- Y.svd.sign * t(Y.svd$u)
  if (type == 1) {
    mu.start <- apply(Y.svd.u, 1, function(x) (sign(x - stats::median(x)) + 1) / 2)
  } else {
    mu.start <- apply(Y.svd.u, 1, function(x) rank(x) - 0.5) / ncol(Y.svd.u)
  }
  mu.start
}
