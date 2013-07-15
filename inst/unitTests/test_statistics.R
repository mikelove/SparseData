test_calculateTStats <- function() {
  # a simple example of t-statistic calculation
  
  sds <- simulateSparseDataSet(10, c(10,10,10), nzs=1, nzg=1)
  sds <- calculateMeans(sds)
  sds <- calculateTStats(sds, offset=1)
  tStats(sds)[["c1"]]

  n_k <- sum(conditions(sds) == "c1")
  n <- ncol(sds)
  K <- nlevels(conditions(sds))
  s <- sqrt(sumSquares(sds)[["global"]]/(n - K))
  numerator <- means(sds)[["c1"]] - means(sds)[["global"]]
  denominator <- sqrt(1/n_k + 1/n) * (s + 1)

  checkEquals(as.numeric(tStats(sds)[["c1"]]), as.numeric(numerator/denominator))
}
