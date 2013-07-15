test_sparseCov <- function() {
  sds1 <- simulateSparseDataSet(100,c(10,10))
  sds2 <- simulateSparseDataSet(100, c(5,5))
  densemat1 <- as.matrix(sparseData(sds1))
  densemat2 <- as.matrix(sparseData(sds2))
  
  cormat <- sparseCov(sparseData(sds1))$cor
  dense.cormat <- cor(densemat1)
  checkEquals(cormat, dense.cormat)

  corxymat <- sparseCov(sparseData(sds1), sparseData(sds2))$cor
  dense.corxymat <- cor(densemat1, densemat2)
  checkEquals(corxymat, dense.corxymat)
}

test_sparseEuclid <- function() {
  sds1 <- simulateSparseDataSet(100,c(10,10))
  sds2 <- simulateSparseDataSet(100, c(5,5))
  densemat1 <- as.matrix(sparseData(sds1))
  densemat2 <- as.matrix(sparseData(sds2))
  
  distmat <- sparseEuclid(sparseData(sds1))
  dense.distmat <- as.matrix(dist(t(densemat1), method="euclid"))
  checkEquals(distmat, dense.distmat)

  distxymat <- sparseEuclid(sparseData(sds1), sparseData(sds2))
  dense.distxymat <- as.matrix(dist(t(cbind(densemat1,densemat2)), method="euclid"))[1:20,21:30]
  checkEquals(distxymat, dense.distxymat)
}


test_sparseCosine <- function() {
  sds1 <- simulateSparseDataSet(100,c(10,10))
  sds2 <- simulateSparseDataSet(100, c(5,5))
  densemat1 <- as.matrix(sparseData(sds1))
  densemat2 <- as.matrix(sparseData(sds2))
  
  cosinemat <- sparseCosine(sparseData(sds1))
  vec1 <- as.numeric(densemat1[,1])
  vec2 <- as.numeric(densemat1[,2])
  dense.cosine <- sum(vec1 * vec2) / (sqrt(sum(vec1^2)) * sqrt(sum(vec2^2)))
  checkEquals(cosinemat[1,2], dense.cosine)

  cosinexymat <- sparseCosine(sparseData(sds1), sparseData(sds2))
  vec1 <- as.numeric(densemat1[,1])
  vec2 <- as.numeric(densemat2[,1])
  dense.cosinexy <- sum(vec1 * vec2) / (sqrt(sum(vec1^2)) * sqrt(sum(vec2^2)))
  checkEquals(cosinexymat[1,1], dense.cosinexy)
}
