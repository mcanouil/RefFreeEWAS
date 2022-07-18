#################################################################################
# omnibusBoot: Omnibus test of association using bootstraps
#################################################################################
omnibusBoot <- function(est, boots, denDegFree) {
  nFeature <- length(est)
  se <- apply(boots, 1, sd)
  pv <- 2 * pt(-abs(est) / se, denDegFree)
  pvNull <- 2 * pt(-(1 / se) * abs(sweep(boots, 1, apply(boots, 1, mean), "-")), denDegFree + 1)

  ks <- max(abs(sort(pv) - (1:nFeature - 0.5) / nFeature))
  ksNull <- apply(pvNull, 2, function(u) max(abs(sort(u) - (1:nFeature - 0.5) / nFeature)))
  mean(ks <= ksNull)
}
