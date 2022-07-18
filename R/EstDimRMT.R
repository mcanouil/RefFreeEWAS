#################################################################################
# EstDimRMT: Estimate dimension using Random Matrix Theory
#  Note:  this function was originally authored by A. Teschendorff in the
#         package isva.
#         Previous versions of RefFreeEWAS used the isva version of the function.
#         However, because of dependency issues in that package, this version
#         simply reproduces the function found in version 1.9 of isva and
#         removes the dependency on the isva package
#         Plotting functionality also removed from original source.
#################################################################################
EstDimRMT <- function(Rmat) {
  data.m <- Rmat
  M <- apply(data.m, 2, function(X) {
    (X - mean(X)) / sqrt(var(X))
  })
  sigma2 <- var(as.vector(M))
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
  estdens.o <- density(eigen.o$values,
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
