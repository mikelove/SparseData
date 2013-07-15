\name{calculateTStats}
\Rdversion{1.1}
\docType{methods}
\alias{calculateTStats}
\alias{calculateTStats,SparseDataSet-method}
\title{ Calculate the t-statistics for a SparseDataSet }
\description{
  Calculate the t-statistics for a SparseDataSet.  The t-statistic is a
  modified version of the typical equal-variance t-statistic, see the
  SparseData vignette for more details.  Using the notation of

  Tibshirani, R., Hastie, T., Narasimhan, B. Chu, G. "Diagnosis of
  multiple cancer types by shrunken centroids of gene
  expression". Proceedings of the National Academy of Sciences 99,
  6567-6572(2002)

  For n_k samples in condition k, n total samples, and K conditions:

  s = sqrt(sumsquare_global / (n - K))
  
  tstat = (mean_k - mean_global) / (sqrt(1/n_k + 1/n) * (s + offset))
 
}
\usage{
  \S4method{calculateTStats}{SparseDataSet}( object, offset="mean", quiet=FALSE )
}
\arguments{
  \item{object}{ a SparseDataSet, with \code{\link{calculateMeans}}
    already having been called }
  \item{offset}{ either "mean", which puts the mean of the pooled
    standard deviations as an offset to the denominator, or a numeric
    value which will be used as the offset } 
  \item{quiet}{ should the function not print out each condition name as
    it is calculated }
}
\author{ Michael Love }
\examples{
  sds <- simulateSparseDataSet(5, c(5,5))
  sds <- calculateMeans(sds)
  sds <- calculateTStats(sds)
  tStats(sds)
}