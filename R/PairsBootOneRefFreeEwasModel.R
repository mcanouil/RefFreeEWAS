PairsBootOneRefFreeEwasModel <- function(mod, pairID) {
  n2 <- dim(mod$X)[1]
  iboot <- bootstrapPairs(1:n2, pairID)
  mu <- mod$Bstar %*% t(mod$X)
  return(mu + mod$dispersion * mod$E[, iboot])
}
