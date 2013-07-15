test_recalc <- function() {
  sds <- simulateSparseDataSet(10,c(2,2,2,2))
  sds <- calculateMeans(sds)
  
  checkEquals(sort(names(means(sds))), sort(c(levels(conditions(sds)),"global")))
  checkEquals(sort(names(sumSquares(sds))), sort(c(levels(conditions(sds)),"global")))

  sds <- calculateMeans(sds, recalc=TRUE)
  
  # recalculation removed the earlier lists
  checkEquals(sort(names(means(sds))), sort(c(levels(conditions(sds)),"global")))
  checkEquals(sort(names(sumSquares(sds))), sort(c(levels(conditions(sds)),"global")))
  
  means(sds)[["c1"]] <- NULL
  sumSquares(sds)[["c2"]] <- NULL

  # this should recalculate only the c1 and c2 means and sum of squares
  # and replace the global means and sum of squares
  sds <- calculateMeans(sds)
  checkEquals(sort(names(means(sds))), sort(c(levels(conditions(sds)),"global")))
  checkEquals(sort(names(sumSquares(sds))), sort(c(levels(conditions(sds)),"global")))
}

  

test_combine <- function() {
  x <- simulateSparseDataSet(10,c(2,2,2))
  y <- simulateSparseDataSet(10,c(2,2,2))
  sampleNames(y) <- paste("sample",(ncol(x) + 1:ncol(y)),sep="")
  pData(y)$sampleID <- sampleNames(y)
  conditions(y) <- factor(c("c1","c1","c4","c4","c5","c5"))

  x <- calculateMeans(x)
  x <- calculateTStats(x)
  y <- calculateMeans(y)
  y <- calculateTStats(y)
  
  z <- combine(x,y)
  
  # the combined statistics should be the union minus the intersection
  # as the shared conditions will need to be recalculated
  checkEquals(sort(names(means(z))), sort(setdiff(union(names(means(x)),names(means(y))), intersect(names(means(x)),names(means(y))))))
  checkEquals(sort(names(sumSquares(z))), sort(setdiff(union(names(sumSquares(x)),names(sumSquares(y))), intersect(names(sumSquares(x)),names(sumSquares(y))))))

  # t-statistics should be removed
  checkEquals(length(tStats(z)), 0)

  z <- calculateMeans(z)
  z <- calculateTStats(z)
  
  # now all means, sum of squares and t-statistics should be present
  checkEquals(sort(names(means(z))), sort(c(levels(conditions(z)),"global")))
  checkEquals(sort(names(sumSquares(z))), sort(c(levels(conditions(z)),"global")))
  checkEquals(sort(names(tStats(z))), sort(levels(conditions(z))))
}
