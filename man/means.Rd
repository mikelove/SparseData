\name{means-methods}
\Rdversion{1.1}
\docType{methods}
\alias{means}
\alias{means<-}
\alias{means,SparseDataSet-method}
\alias{means<-,SparseDataSet-method}
\title{ Accessor for the means }
\description{
  Accessor for the means, a list which is created by the
  \code{\link{calculateMeans}} function.
  }
\usage{
  \S4method{means}{SparseDataSet}(object)
  \S4method{means}{SparseDataSet}(object) <- value
}
\arguments{
  \item{object}{ a SparseDataSet}
  \item{value}{ a list }
}
\author{ Michael Love }
\examples{
  sds <- simulateSparseDataSet(5, c(5,5))
  sds <- calculateMeans(sds)
  means(sds)
}