\name{sparseData}
\docType{methods}
\alias{sparseData}
\alias{sparseData<-}
\alias{sparseData,SparseDataSet-method}
\alias{sparseData<-,SparseDataSet,dgCMatrix-method}
\title{ Accessors for the 'sparseData' slot of a SparseDataSet object }
\description{
  The sparseData slot holds a dgCMatrix object defined by the Matrix
  package.  As an extension of the eSet, each row is a feature and each
  column is annotated with phenotypic data.
}
\usage{
  \S4method{sparseData}{SparseDataSet}(object)
  \S4method{sparseData}{SparseDataSet,dgCMatrix}(object) <- value
}
\arguments{
  \item{object}{a SparseDataSet object.}
  \item{value}{a dgCMatrix}
}
\author{
  Michael Love
}
\examples{
  sds <- simulateSparseDataSet(100, c(5,5))
  head(sparseData(sds))
}