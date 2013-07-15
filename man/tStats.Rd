\name{tStats-methods}
\Rdversion{1.1}
\docType{methods}
\alias{tStats}
\alias{tStats<-}
\alias{tStats,SparseDataSet-method}
\alias{tStats<-,SparseDataSet-method}
\title{ Accessor for the tStats }
\description{
  Accessor for the tStats, a list which is created by the
  \code{\link{calculateTStats}} function, which itself requires the
  \code{\link{calculateMeans}} function to have been called previously.

}
\usage{
  \S4method{tStats}{SparseDataSet}(object)
  \S4method{tStats}{SparseDataSet}(object) <- value
}
\arguments{
  \item{object}{ a SparseDataSet}
  \item{value}{ a list }
}
\author{ Michael Love }
\examples{
  sds <- simulateSparseDataSet(5, c(5,5))
  sds <- calculateMeans(sds)
  sds <- calculateTStats(sds)
  tStats(sds)
}