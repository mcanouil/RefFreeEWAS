bootstrapPairs <- function(obs, pairID) {
  pairIx <- split(obs, pairID)
  nobsTab <- table(sapply(pairIx, length))
  if (length(nobsTab) > 1) {
    stop("All clusters must have the same number of observations.\n")
  }

  n <- length(pairIx)
  pairBoot <- pairIx[sample(1:n, n, replace = TRUE)]

  obsBoot <- unlist(lapply(pairBoot, function(u) sample(u, replace = FALSE)))
  names(obsBoot) <- NULL
  obsBoot
}
