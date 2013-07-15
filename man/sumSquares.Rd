\name{sumSquares-methods}
\Rdversion{1.1}
\docType{methods}
\alias{sumSquares}
\alias{sumSquares<-}
\alias{sumSquares,SparseDataSet-method}
\alias{sumSquares<-,SparseDataSet-method}
\title{ Accessor for the sum of squares }
\description{
  Accessor for the sum of squares, a list which is created by the
  \code{\link{calculateMeans}} function. 
}
\usage{
  \S4method{sumSquares}{SparseDataSet}(object)
  \S4method{sumSquares}{SparseDataSet}(object) <- value
}
\arguments{
  \item{object}{ a SparseDataSet}
  \item{value}{ a list }
}
\author{ Michael Love }
\examples{
  sds <- simulateSparseDataSet(5, c(5,5))
  sds <- calculateMeans(sds)
  sumSquares(sds)
}